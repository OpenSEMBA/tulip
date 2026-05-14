"""
compute_incell_parameters.py
============================
Postprocessing script for a Tulip solved case.

Reads a .tulip.out.json file and, for each material with a multipolar
expansion, computes the in-cell capacitance C[i,j] (F/m) and inductance
L[i,j] (H/m) for a fixed reference conductor i and all available conductors j
in that material association.

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
  python scripts/compute_incell_parameters.py <path/to/case.tulip.out.json> <i_ref>
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


def compute_material_rows(material: dict, mat_assoc: list, i_ref: int):
    """Compute (j, C[i_ref,j], L[i_ref,j]) rows for one material."""
    mat_id = material.get("id")
    mat_type = material.get("type", "unknown")
    mp = material.get("multipolarExpansion", {})
    e_sols = mp.get("electric", [])
    m_sols = mp.get("magnetic", [])

    if not e_sols or not m_sols:
        return None, "missing electric or magnetic solutions"

    element_ids = element_ids_for_material(mat_id, mat_assoc, len(e_sols))
    i_local = local_index_from_element_id(element_ids, i_ref)
    if i_local is None:
        return None, f"reference conductor id {i_ref} not in {element_ids}"

    rows = []
    for j in element_ids:
        j_local = local_index_from_element_id(element_ids, j)
        if j_local is None:
            continue
        c_val = capacitance(e_sols, i_local, j_local)
        l_val = inductance(m_sols, i_local, j_local)
        rows.append((j, c_val, l_val))

    return {
        "mat_id": mat_id,
        "mat_type": mat_type,
        "element_ids": element_ids,
        "rows": rows,
    }, None


def print_material_report(report: dict, i_ref: int):
    """Print computed rows for one material."""
    print(f"\n{'=' * 55}")
    print(f"  Material id={report['mat_id']}  type={report['mat_type']}")
    print(f"  Total conductors: {len(report['element_ids'])}")
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

    if args.json_path is None or args.i_ref is None:
        parser.error("json_path and i_ref are required unless --run-tests is used")
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

    print(f"Output JSON path      : {json_path}")
    print(f"Reference conductor i : {i_ref}")

    data = load_json(json_path)
    materials = data.get("materials", [])
    mat_assoc = data.get("materialAssociations", [])

    if not materials:
        print("No materials found in output JSON.")
        return 0

    for material in materials:
        report, warn = compute_material_rows(material, mat_assoc, i_ref)
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

        report, warn = compute_material_rows(material, mat_assoc, 16)
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

        report, warn = compute_material_rows(material, mat_assoc, 16)
        self.assertIsNone(warn)
        self.assertIsNotNone(report)

        rows = {j: (c_val, l_val) for j, c_val, l_val in report["rows"]}
        self.assertIn(16, rows)
        self.assertIn(30, rows)

        self.assertAlmostEqual(rows[16][0], 1.9116497250688996e-11, delta=1e-20)
        self.assertAlmostEqual(rows[16][1], 5.820365736777192e-07, delta=1e-15)
        self.assertAlmostEqual(rows[30][0], 5.165814539740492e-11, delta=1e-20)
        self.assertAlmostEqual(rows[30][1], 2.1538714707844526e-07, delta=1e-15)

    def test_parse_args_requires_positionals(self):
        with contextlib.redirect_stderr(io.StringIO()):
            with self.assertRaises(SystemExit):
                parse_args([])

    def test_parse_args_happy_path(self):
        args = parse_args(["case.tulip.out.json", "16"])
        self.assertEqual(args.json_path, "case.tulip.out.json")
        self.assertEqual(args.i_ref, 16)


if __name__ == "__main__":
    sys.exit(run())
