"""
compute_incell_parameters.py
=============================
Postprocessing script for a Tulip solved case.

Reads the .tulip.out.json file and, for each material with a multipolar
expansion, computes the in-cell capacitance C[0,j] (F/m) and inductance
L[0,j] (H/m) for a fixed reference conductor i=0, varying over the
conductors of interest listed in J_INDICES.

Formulas (from InCellPotentials::getInCellCapacitanceUsingInnerRegion and
getInCellInductanceUsingInnerRegion in Results.cpp):

  C[0,j] = Q_j / (V_0|Vj=1 - <V_j>_inner) * ε₀

  L[0,j] = (A_0|Aj=1 - <A_j>_inner) / I_j * μ₀

where:
  Q_j   = electric[j].ab[0][0]              (monopole charge coefficient)
  <V_j> = electric[j].innerRegionAveragePotential
  V_0   = electric[j].conductorPotentials[0]

  I_j   = magnetic[j].ab[0][0]              (monopole current coefficient)
  <A_j> = magnetic[j].innerRegionAveragePotential
  A_0   = magnetic[j].conductorPotentials[0]

Usage:
  python compute_incell_parameters.py [CASE_DIR]

  CASE_DIR defaults to the realistic_case_fdtd_cell_centered_in_0 directory.
"""

import json
import math
import os
import sys

# Physical constants (match tulip/src/driver/constants.h)
EPSILON0_SI = 8.8541878176e-12   # F/m
MU0_SI      = 4.0e-7 * math.pi  # H/m

# ---------------------------------------------------------------------------
# Configuration
# ---------------------------------------------------------------------------
DEFAULT_CASE_DIR = (
    "/home/luis/workspace/tulip/tmp_cases/realistic_case_just_16_and_30"
)

def load_json(case_dir: str) -> dict:
    case_name = os.path.basename(case_dir.rstrip("/"))
    json_path = os.path.join(case_dir, f"{case_name}.tulip.out.json")
    with open(json_path) as fh:
        return json.load(fh)


# ---------------------------------------------------------------------------
# Core computation (mirrors Results.cpp formulas)
# ---------------------------------------------------------------------------

def capacitance(electric_solutions: list, i: int, j: int) -> float:
    """C[i,j] in F/m."""
    sol_j   = electric_solutions[j]
    Q_j     = sol_j["ab"][0][0]
    avg_V_j = sol_j["innerRegionAveragePotential"]
    V_i     = sol_j["conductorPotentials"][i]
    return Q_j / (V_i - avg_V_j) * EPSILON0_SI


def inductance(magnetic_solutions: list, i: int, j: int) -> float:
    """L[i,j] in H/m."""
    sol_j   = magnetic_solutions[j]
    I_j     = sol_j["ab"][0][0]
    avg_A_j = sol_j["innerRegionAveragePotential"]
    A_i     = sol_j["conductorPotentials"][i]
    return (A_i - avg_A_j) / I_j * MU0_SI


def local_index_from_element_id(element_ids: list, element_id: int):
    """Map a physical conductor id to the local multipolar solution index."""
    try:
        return element_ids.index(element_id)
    except ValueError:
        return None


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    case_dir = sys.argv[1] if len(sys.argv) > 1 else DEFAULT_CASE_DIR
    case_dir = os.path.abspath(case_dir)

    # Reference conductor
    I_REF = 16

    # Conductor indices j to evaluate
    J_INDICES = [16, 30]


    print(f"Case directory : {case_dir}")
    print(f"Reference conductor i = {I_REF}")
    print(f"Conductors j          = {J_INDICES}")

    data      = load_json(case_dir)
    materials = data.get("materials", [])
    mat_assoc = data.get("materialAssociations", [])

    if not materials:
        print("No materials found in output JSON.")
        return

    for material in materials:
        mat_id   = material["id"]
        mat_type = material.get("type", "unknown")
        mp       = material.get("multipolarExpansion", {})

        e_sols = mp.get("electric", [])
        m_sols = mp.get("magnetic", [])

        if not e_sols or not m_sols:
            print(f"[warn] Material {mat_id}: missing electric or magnetic solutions, skipping.")
            continue

        assoc       = next((a for a in mat_assoc if a["materialId"] == mat_id), None)
        element_ids = assoc["elementIds"] if assoc else list(range(len(e_sols)))

        i_local = local_index_from_element_id(element_ids, I_REF)
        if i_local is None:
            print(f"[warn] Material {mat_id}: reference conductor id {I_REF} not in {element_ids}, skipping.")
            continue
        print(f"\n{'='*55}")
        print(f"  Material id={mat_id}  type={mat_type}")
        print(f"  Total conductors: {len(element_ids)}")
        print(f"{'='*55}")
        print(f"  {'j':>4}  {f'C[{I_REF},j] (F/m)':>18}  {f'L[{I_REF},j] (H/m)':>18}")
        print(f"  {'-'*4}  {'-'*18}  {'-'*18}")

        for j in J_INDICES:
            j_local = local_index_from_element_id(element_ids, j)
            if j_local is None:
                print(f"  {j:>4}  [id not in materialAssociations.elementIds]")
                continue
            C_val = capacitance(e_sols, i_local, j_local)
            L_val = inductance(m_sols, i_local, j_local)
            print(f"  {j:>4}  {C_val:>+18.6e}  {L_val:>+18.6e}")

        print()


if __name__ == "__main__":
    main()
