"""
compute_incell_parameters.py
============================
Postprocessing script for a Tulip solved case.

Reads a .tulip.out.json file and, for each material with a multipolar
expansion, computes the in-cell capacitance C[i,j] (F/m) and inductance
L[i,j] (H/m) for a fixed reference conductor i and all available conductors j
in that material association.

The user must select one mode with --mode:

    * inner-region: uses getInCellCapacitanceUsingInnerRegion / 
        getInCellInductanceUsingInnerRegion formulas.
    * on-box: uses getInCellCapacitanceOnBox / getInCellInductanceOnBox via
        multipolar integration over the specified cell box.

Formulas (from InCellPotentials::getInCellCapacitanceUsingInnerRegion and
getInCellInductanceUsingInnerRegion in Results.cpp):

  C[i,j] = Q_j / (V_i|Vj=1 - <V_j>_inner) * epsilon0

  L[i,j] = (A_i|Aj=1 - <A_j>_inner) / I_j * mu0

where:
  Q_j   = electric[j].ab[0][0]              (monopole charge coefficient)
  <V_j> = electric[j].innerRegionAveragePotential
  V_i   = electric[j].conductorPotentials[i]

  I_j   = magnetic[j].ab[0][0]              (monopole current coefficient)
  <A_j> = magnetic[j].innerRegionAveragePotential
  A_i   = magnetic[j].conductorPotentials[i]

Usage:
    python scripts/compute_incell_parameters.py --mode inner-region <path/to/case.tulip.out.json> <i_ref>
    python scripts/compute_incell_parameters.py --mode on-box --cell-box xmin ymin xmax ymax <path/to/case.tulip.out.json> <i_ref>
  python scripts/compute_incell_parameters.py --run-tests
"""

import argparse
import contextlib
import io
import json
import math
import os
import sys
import unittest

# Physical constants (match tulip/src/driver/constants.h)
EPSILON0_SI = 8.8541878176e-12  # F/m
MU0_SI = 4.0e-7 * math.pi  # H/m
GRID_INTEGRATION_SAMPLING_POINTS = 100


def load_json(json_path: str) -> dict:
    """Load a Tulip output JSON file from disk."""
    abs_path = os.path.abspath(json_path)
    with open(abs_path, encoding="utf-8") as fh:
        return json.load(fh)


def capacitance(electric_solutions: list, i: int, j: int) -> float:
    """Compute C[i,j] in F/m."""
    sol_j = electric_solutions[j]
    q_j = sol_j["ab"][0][0]
    avg_v_j = sol_j["innerRegionAveragePotential"]
    v_i = sol_j["conductorPotentials"][i]
    return q_j / (v_i - avg_v_j) * EPSILON0_SI


def inductance(magnetic_solutions: list, i: int, j: int) -> float:
    """Compute L[i,j] in H/m."""
    sol_j = magnetic_solutions[j]
    i_j = sol_j["ab"][0][0]
    avg_a_j = sol_j["innerRegionAveragePotential"]
    a_i = sol_j["conductorPotentials"][i]
    return (a_i - avg_a_j) / i_j * MU0_SI


def rel_error(expected: float, computed: float) -> float:
    """Relative error used by DriverTests: |a-b| / max(|a|, |b|)."""
    denom = max(abs(expected), abs(computed))
    if denom == 0.0:
        return 0.0
    return abs(expected - computed) / denom


def make_box(xmin: float, ymin: float, xmax: float, ymax: float) -> dict:
    """Build a Box-like dict using the C++ field names."""
    return {"min": [xmin, ymin], "max": [xmax, ymax]}


def box_area(box: dict) -> float:
    """Compute area of a Box-like dict."""
    return (box["max"][0] - box["min"][0]) * (box["max"][1] - box["min"][1])


def point_in_box(point: tuple, box: dict) -> bool:
    """Check point inclusion using the same closed intervals as C++."""
    return (
        box["min"][0] <= point[0] <= box["max"][0]
        and box["min"][1] <= point[1] <= box["max"][1]
    )


def build_integration_planes_for_box(integration_box: dict, inner_region_box: dict):
    """Port of buildIntegrationPlanesForBox in Results.cpp."""
    for x in (0, 1):
        if (
            integration_box["min"][x] > inner_region_box["min"][x]
            or integration_box["max"][x] < inner_region_box["max"][x]
        ):
            raise ValueError("Integration box has to be larger than the inner region.")

    planes = []
    for x in (0, 1):
        control_points = {
            integration_box["min"][x],
            integration_box["max"][x],
            inner_region_box["min"][x],
            inner_region_box["max"][x],
        }
        sorted_points = sorted(control_points)

        axis_planes = set()
        for idx in range(1, len(sorted_points)):
            prev = sorted_points[idx - 1]
            curr = sorted_points[idx]
            step = (curr - prev) / GRID_INTEGRATION_SAMPLING_POINTS
            for k in range(GRID_INTEGRATION_SAMPLING_POINTS):
                axis_planes.add(prev + k * step)
            axis_planes.add(curr)

        planes.append(sorted(axis_planes))

    return planes


def multipolar_expansion(position: tuple, ab: list, expansion_center: list) -> float:
    """2D multipolar expansion matching multipolarExpansion.h."""
    rx = position[0] - expansion_center[0]
    ry = position[1] - expansion_center[1]
    r = math.hypot(rx, ry)
    phi = math.atan2(ry, rx)

    res = 0.0
    for n, coeff in enumerate(ab):
        an = coeff[0]
        bn = coeff[1]
        if n == 0:
            res -= an * math.log(r)
        else:
            res += (an * math.cos(n * phi) + bn * math.sin(n * phi)) / (r ** n)
    return res / (2.0 * math.pi)


def get_average_potential(solution: dict, inner_box: dict, outer_box: dict) -> float:
    """Port of getAveragePotential in Results.cpp."""
    integration_planes = build_integration_planes_for_box(outer_box, inner_box)
    outer_v = 0.0

    for m in range(1, len(integration_planes[0])):
        for n in range(1, len(integration_planes[1])):
            x_min = integration_planes[0][m - 1]
            x_max = integration_planes[0][m]
            y_min = integration_planes[1][n - 1]
            y_max = integration_planes[1][n]

            mid_point = (0.5 * (x_min + x_max), 0.5 * (y_min + y_max))
            area = (x_max - x_min) * (y_max - y_min)
            if point_in_box(mid_point, inner_box):
                continue

            outer_v += area * multipolar_expansion(
                mid_point,
                solution["ab"],
                solution["expansionCenter"],
            )

    inner_v = solution["innerRegionAveragePotential"] * box_area(inner_box)
    return (inner_v + outer_v) / box_area(outer_box)


def capacitance_on_box(
    electric_solutions: list,
    i: int,
    j: int,
    inner_region_box: dict,
    cell_box: dict,
) -> float:
    """Compute C[i,j] in F/m using getInCellCapacitanceOnBox semantics."""
    sol_j = electric_solutions[j]
    q_j = sol_j["ab"][0][0]
    avg_v_j = get_average_potential(sol_j, inner_region_box, cell_box)
    v_i = sol_j["conductorPotentials"][i]
    return q_j / (v_i - avg_v_j) * EPSILON0_SI


def inductance_on_box(
    magnetic_solutions: list,
    i: int,
    j: int,
    inner_region_box: dict,
    cell_box: dict,
) -> float:
    """Compute L[i,j] in H/m using getInCellInductanceOnBox semantics."""
    sol_j = magnetic_solutions[j]
    i_j = sol_j["ab"][0][0]
    avg_a_j = get_average_potential(sol_j, inner_region_box, cell_box)
    a_i = sol_j["conductorPotentials"][i]
    return (a_i - avg_a_j) / i_j * MU0_SI


def local_index_from_element_id(element_ids: list, element_id: int):
    """Map a physical conductor id to the local multipolar solution index."""
    try:
        return element_ids.index(element_id)
    except ValueError:
        return None


def element_ids_for_material(material_id: int, mat_assoc: list, fallback_size: int) -> list:
    """Return material-associated conductor ids, or fallback to local indices."""
    assoc = next((a for a in mat_assoc if a.get("materialId") == material_id), None)
    if assoc is not None and "elementIds" in assoc:
        return assoc["elementIds"]
    return list(range(fallback_size))


def compute_material_rows(material: dict, mat_assoc: list, i_ref: int, mode: str, cell_box: dict = None):
    """Compute (j, C[i_ref,j], L[i_ref,j]) rows for one material and mode."""
    mat_id = material.get("id")
    mat_type = material.get("type", "unknown")
    mp = material.get("multipolarExpansion", {})
    e_sols = mp.get("electric", [])
    m_sols = mp.get("magnetic", [])
    inner_region_box = mp.get("innerRegionBox")

    if not e_sols or not m_sols:
        return None, "missing electric or magnetic solutions"

    if mode == "on-box":
        if cell_box is None:
            return None, "on-box mode requires a cell box"
        if inner_region_box is None:
            return None, "on-box mode requires multipolarExpansion.innerRegionBox"
        if "expansionCenter" not in e_sols[0] or "expansionCenter" not in m_sols[0]:
            return None, "on-box mode requires expansionCenter in electric/magnetic solutions"

    element_ids = element_ids_for_material(mat_id, mat_assoc, len(e_sols))
    i_local = local_index_from_element_id(element_ids, i_ref)
    if i_local is None:
        return None, f"reference conductor id {i_ref} not in {element_ids}"

    rows = []
    for j in element_ids:
        j_local = local_index_from_element_id(element_ids, j)
        if j_local is None:
            continue
        if mode == "inner-region":
            c_val = capacitance(e_sols, i_local, j_local)
            l_val = inductance(m_sols, i_local, j_local)
        elif mode == "on-box":
            c_val = capacitance_on_box(e_sols, i_local, j_local, inner_region_box, cell_box)
            l_val = inductance_on_box(m_sols, i_local, j_local, inner_region_box, cell_box)
        else:
            return None, f"unsupported mode '{mode}'"
        rows.append((j, c_val, l_val))

    return {
        "mat_id": mat_id,
        "mat_type": mat_type,
        "element_ids": element_ids,
        "rows": rows,
        "mode": mode,
    }, None


def print_material_report(report: dict, i_ref: int):
    """Print computed rows for one material."""
    print(f"\n{'=' * 55}")
    print(f"  Material id={report['mat_id']}  type={report['mat_type']}")
    print(f"  Total conductors: {len(report['element_ids'])}")
    print(f"  Mode: {report['mode']}")
    print(f"{'=' * 55}")
    print(f"  {'j':>4}  {f'C[{i_ref},j] (F/m)':>18}  {f'L[{i_ref},j] (H/m)':>18}")
    print(f"  {'-' * 4}  {'-' * 18}  {'-' * 18}")

    for j, c_val, l_val in report["rows"]:
        print(f"  {j:>4}  {c_val:>+18.6e}  {l_val:>+18.6e}")


def parse_args(argv: list):
    """Parse command-line arguments."""
    parser = argparse.ArgumentParser(
        description="Compute in-cell C/L parameters from a .tulip.out.json file."
    )
    parser.add_argument(
        "--run-tests",
        action="store_true",
        help="Run the unit tests embedded in this file.",
    )
    parser.add_argument(
        "--mode",
        choices=["inner-region", "on-box"],
        help="Choose the computational mode. Required unless --run-tests is used.",
    )
    parser.add_argument(
        "--cell-box",
        nargs=4,
        type=float,
        metavar=("XMIN", "YMIN", "XMAX", "YMAX"),
        help="Cell box for on-box mode: xmin ymin xmax ymax.",
    )
    parser.add_argument(
        "json_path",
        nargs="?",
        help="Path to a .tulip.out.json file.",
    )
    parser.add_argument(
        "i_ref",
        nargs="?",
        type=int,
        help="Reference conductor id i.",
    )

    args = parser.parse_args(argv)
    if args.run_tests:
        return args

    if args.mode is None:
        parser.error("--mode is required unless --run-tests is used")
    if args.json_path is None or args.i_ref is None:
        parser.error("json_path and i_ref are required unless --run-tests is used")

    if args.mode == "on-box":
        if args.cell_box is None:
            parser.error("--cell-box is required when --mode on-box is selected")
        if args.cell_box[2] <= args.cell_box[0] or args.cell_box[3] <= args.cell_box[1]:
            parser.error("--cell-box must satisfy xmax > xmin and ymax > ymin")
    else:
        if args.cell_box is not None:
            parser.error("--cell-box is only valid when --mode on-box is selected")

    return args


def run_self_tests() -> int:
    """Execute same-file unit tests."""
    test_program = unittest.main(module=__name__, argv=[sys.argv[0]], exit=False)
    return 0 if test_program.result.wasSuccessful() else 1


def run(argv=None) -> int:
    """Entrypoint used by both CLI and tests."""
    args = parse_args(sys.argv[1:] if argv is None else argv)
    if args.run_tests:
        return run_self_tests()

    json_path = os.path.abspath(args.json_path)
    i_ref = args.i_ref
    mode = args.mode
    cell_box = None
    if mode == "on-box":
        cell_box = make_box(*args.cell_box)

    print(f"Output JSON path      : {json_path}")
    print(f"Reference conductor i : {i_ref}")
    print(f"Mode                  : {mode}")
    if cell_box is not None:
        print(f"Cell box              : {cell_box}")

    data = load_json(json_path)
    materials = data.get("materials", [])
    mat_assoc = data.get("materialAssociations", [])

    if not materials:
        print("No materials found in output JSON.")
        return 0

    for material in materials:
        report, warn = compute_material_rows(material, mat_assoc, i_ref, mode, cell_box)
        if warn is not None:
            print(f"[warn] Material {material.get('id')}: {warn}, skipping.")
            continue
        print_material_report(report, i_ref)
        print()

    return 0


class TestComputeInCellParameters(unittest.TestCase):
    """Unit tests for compute_incell_parameters.py."""

    def test_local_index_lookup(self):
        element_ids = [16, 30, 42]
        self.assertEqual(local_index_from_element_id(element_ids, 30), 1)
        self.assertIsNone(local_index_from_element_id(element_ids, 99))

    def test_compute_material_rows_uses_all_conductors(self):
        material = {
            "id": 7,
            "type": "dielectric",
            "multipolarExpansion": {
                "electric": [
                    {
                        "ab": [[1.0]],
                        "innerRegionAveragePotential": 0.0,
                        "conductorPotentials": [2.0, 3.0],
                    },
                    {
                        "ab": [[2.0]],
                        "innerRegionAveragePotential": 1.0,
                        "conductorPotentials": [3.0, 4.0],
                    },
                ],
                "magnetic": [
                    {
                        "ab": [[1.0]],
                        "innerRegionAveragePotential": 0.0,
                        "conductorPotentials": [2.0, 3.0],
                    },
                    {
                        "ab": [[2.0]],
                        "innerRegionAveragePotential": 1.0,
                        "conductorPotentials": [3.0, 4.0],
                    },
                ],
            },
        }
        mat_assoc = [{"materialId": 7, "elementIds": [16, 30]}]

        report, warn = compute_material_rows(material, mat_assoc, 16, mode="inner-region")
        self.assertIsNone(warn)
        self.assertIsNotNone(report)
        self.assertEqual([row[0] for row in report["rows"]], [16, 30])

    def test_reference_case_values_from_embedded_fixture(self):
        material = {
            "id": 1,
            "type": "unshieldedMultiwire",
            "multipolarExpansion": {
                "electric": [
                    {
                        "ab": [
                            [0.8419073593680602, 0.0],
                            [4.7094817324223296e-07, -7.437578541169077e-07],
                            [-2.254972536944232e-06, -1.3845721768518446e-06],
                            [-1.7502712357816652e-09, 1.6065198400637674e-09],
                        ],
                        "conductorPotentials": [1.0, 0.7525321728251249],
                        "expansionCenter": [0.0017138274803646352, -0.0027066091802577024],
                        "innerRegionAveragePotential": 0.6100537777863153,
                    },
                    {
                        "ab": [
                            [1.0308546965043754, 0.0],
                            [4.236110742075576e-07, -8.85263389452124e-07],
                            [-1.301845581534472e-07, -8.009832224471683e-08],
                            [3.712118282384746e-10, -3.4152957909780026e-10],
                        ],
                        "conductorPotentials": [0.9213984270827741, 1.0],
                        "expansionCenter": [0.0023197106176927533, -0.004847736541849076],
                        "innerRegionAveragePotential": 0.7447102973543827,
                    },
                ],
                "magnetic": [
                    {
                        "ab": [
                            [0.8419073593680602, 0.0],
                            [4.7094817324223296e-07, -7.437578541169077e-07],
                            [-2.254972536944232e-06, -1.3845721768518446e-06],
                            [-1.7502712357816652e-09, 1.6065198400637674e-09],
                        ],
                        "conductorPotentials": [1.0, 0.7525321728251249],
                        "expansionCenter": [0.0017138274803646352, -0.0027066091802577024],
                        "innerRegionAveragePotential": 0.6100537777863153,
                    },
                    {
                        "ab": [
                            [1.0308546965043754, 0.0],
                            [4.236110742075576e-07, -8.85263389452124e-07],
                            [-1.301845581534472e-07, -8.009832224471683e-08],
                            [3.712118282384746e-10, -3.4152957909780026e-10],
                        ],
                        "conductorPotentials": [0.9213984270827741, 1.0],
                        "expansionCenter": [0.0023197106176927533, -0.004847736541849076],
                        "innerRegionAveragePotential": 0.7447102973543827,
                    },
                ],
            },
        }
        mat_assoc = [{"elementIds": [16, 30], "materialId": 1}]

        report, warn = compute_material_rows(material, mat_assoc, 16, mode="inner-region")
        self.assertIsNone(warn)
        self.assertIsNotNone(report)

        rows = {j: (c_val, l_val) for j, c_val, l_val in report["rows"]}
        self.assertIn(16, rows)
        self.assertIn(30, rows)

        self.assertAlmostEqual(rows[16][0], 1.9116497250688996e-11, delta=1e-20)
        self.assertAlmostEqual(rows[16][1], 5.820365736777192e-07, delta=1e-15)
        self.assertAlmostEqual(rows[30][0], 5.165814539740492e-11, delta=1e-20)
        self.assertAlmostEqual(rows[30][1], 2.1538714707844526e-07, delta=1e-15)

    def test_on_box_mode_matches_driver_test_single_wire_values(self):
        # This fixture is equivalent to the in-cell output used by
        # DriverTest::lansink2024_single_wire_multipolar_in_cell_parameters.
        material = {
            "id": 1,
            "type": "unshieldedMultiwire",
            "multipolarExpansion": {
                "innerRegionBox": {
                    "max": [0.0040000002, 0.0040000002],
                    "min": [-0.0040000002, -0.0040000002],
                },
                "electric": [
                    {
                        "ab": [
                            [0.9766734089848975, 0.0],
                            [-5.686345063943916e-12, -3.090752789506213e-12],
                            [4.637444003229173e-11, 8.593758249532865e-12],
                            [7.014114840799095e-15, 4.218755938827563e-15],
                            [-9.241748384573654e-18, 3.9674592343454175e-18],
                            [-1.0388931187475782e-21, 3.5460599095173054e-21],
                        ],
                        "conductorPotentials": [1.0],
                        "expansionCenter": [-4.899147920455392e-08, -2.6628800768531107e-08],
                        "innerRegionAveragePotential": 0.9040784423949009,
                    }
                ],
                "magnetic": [
                    {
                        "ab": [
                            [0.9092956935239767, 0.0],
                            [4.676064105592532e-12, -1.9984633243172434e-12],
                            [1.79380465488009e-11, 3.6833862710835886e-12],
                            [-2.3149425432828747e-15, 1.8113967798566277e-15],
                            [8.898740430780706e-19, -6.653588044528165e-19],
                            [-4.9778872535021934e-23, -2.6188498558876825e-22],
                        ],
                        "conductorPotentials": [1.0],
                        "expansionCenter": [4.3634340291863276e-08, -1.8648510194811956e-08],
                        "innerRegionAveragePotential": 0.8490379205671101,
                    }
                ],
            },
        }
        mat_assoc = [{"materialId": 1, "elementIds": [0]}]
        fdtd_cell = make_box(-0.0075, -0.0075, 0.0075, 0.0075)

        report, warn = compute_material_rows(
            material,
            mat_assoc,
            0,
            mode="on-box",
            cell_box=fdtd_cell,
        )
        self.assertIsNone(warn)
        self.assertIsNotNone(report)

        rows = {j: (c_val, l_val) for j, c_val, l_val in report["rows"]}
        self.assertIn(0, rows)

        expected_c00 = 49.11e-12
        expected_l00 = 320e-9
        r_tol = 0.06

        self.assertLess(rel_error(expected_c00, rows[0][0]), r_tol)
        self.assertLess(rel_error(expected_l00, rows[0][1]), r_tol)

    def test_parse_args_requires_mode(self):
        with contextlib.redirect_stderr(io.StringIO()):
            with self.assertRaises(SystemExit):
                parse_args(["case.tulip.out.json", "16"])

    def test_parse_args_inner_region_happy_path(self):
        args = parse_args(["--mode", "inner-region", "case.tulip.out.json", "16"])
        self.assertEqual(args.json_path, "case.tulip.out.json")
        self.assertEqual(args.i_ref, 16)
        self.assertEqual(args.mode, "inner-region")

    def test_parse_args_on_box_requires_cell_box(self):
        with contextlib.redirect_stderr(io.StringIO()):
            with self.assertRaises(SystemExit):
                parse_args(["--mode", "on-box", "case.tulip.out.json", "16"])

    def test_parse_args_on_box_happy_path(self):
        args = parse_args(
            [
                "--mode",
                "on-box",
                "--cell-box",
                "-0.1",
                "-0.1",
                "0.1",
                "0.1",
                "case.tulip.out.json",
                "16",
            ]
        )
        self.assertEqual(args.mode, "on-box")
        self.assertEqual(args.cell_box, [-0.1, -0.1, 0.1, 0.1])


if __name__ == "__main__":
    sys.exit(run())
