"""
Tulip ParaView Macro
====================
Loads all solved conductor solutions for a Tulip case and creates combined
potential views for each multipolar expansion solution in the .tulip.out.json.

Usage:
  1. Set CASE_DIR to your case directory.
  2. Run as a ParaView macro (Macros > Add New Macro or run from Python Shell).

The macro creates:
  - Layout "Individual Conductors": all PVD readers for electrostatic and
    magnetostatic solved conductors (toggle visibility per conductor).
  - For each material with a multipolarExpansion, a layout
    "Material_<id>_<type>_Modes" containing two views:
      * Electric modes: combined Phi field weighted by conductorPotentials
        for every electric solution; stored as arrays ElectricMode_000, ...
      * Magnetic modes: same for magnetic solutions (MagneticMode_000, ...).
    Switch between modes in the "Color By" selector in the Properties panel.
"""

import os
import re
import json

import paraview.simple as pvs

# ==============================================================================
# CONFIGURATION
# ==============================================================================
# Assumes that results are in this case folder.
CASE_DIR = "PATH_TO_CASE"

# ==============================================================================
# Helpers
# ==============================================================================

def _find_pvd(pv_dir: str, conductor_idx: int, kind: str) -> str | None:
    name = f"Conductor_{conductor_idx}_{kind}"
    pvd = os.path.join(pv_dir, name, f"{name}.pvd")
    return pvd if os.path.exists(pvd) else None


def _conductor_index(folder_name: str) -> int | None:
    m = re.match(r"Conductor_(\d+)_(electrostatic|magnetostatic)$", folder_name)
    return int(m.group(1)) if m else None


def _make_combined_filter(readers: list, weights_all: list[list[float]], prefix: str):
    """
    Return a ProgrammableFilter that reads Phi from each input, computes
    every weighted sum listed in weights_all, and stores results as
    <prefix>_000, <prefix>_001, ... point arrays.

    readers       – list of ParaView source objects (same mesh topology)
    weights_all   – list of weight vectors, one per solution
    prefix        – array name prefix, e.g. 'ElectricMode'
    """
    weights_repr = repr(weights_all)

    script = f"""
import numpy as np
from vtk.util import numpy_support

weights_all = {weights_repr}

output.ShallowCopy(inputs[0].VTKObject)

for sol_idx, weights in enumerate(weights_all):
    result = np.zeros(inputs[0].GetNumberOfPoints())
    for w, inp in zip(weights, inputs):
        phi_vtk = inp.PointData.GetArray('Phi')
        if phi_vtk is not None:
            result += w * numpy_support.vtk_to_numpy(phi_vtk)
    arr = numpy_support.numpy_to_vtk(result, deep=True)
    arr.SetName(f'{prefix}_{{sol_idx:03d}}')
    output.GetPointData().AddArray(arr)
"""
    pf = pvs.ProgrammableFilter(Input=readers)
    pf.Script = script
    pf.RequestInformationScript = ""
    pf.RequestUpdateExtentScript = ""
    return pf


# ==============================================================================
# Main
# ==============================================================================

pvs._DisableFirstRenderCameraReset()

case_name = os.path.basename(CASE_DIR)
json_path = os.path.join(CASE_DIR, f"{case_name}.tulip.out.json")
pv_dir = os.path.join(CASE_DIR, "ParaView")

with open(json_path) as fh:
    out_data = json.load(fh)

materials = out_data.get("materials", [])
mat_assoc = out_data.get("materialAssociations", [])

# Collect sorted conductor indices present in the ParaView folder
all_indices = sorted(
    {_conductor_index(d) for d in os.listdir(pv_dir) if _conductor_index(d) is not None}
)

# --------------------------------------------------------------------------
# Layout 1 – Individual conductor solutions
# --------------------------------------------------------------------------
layout_indiv = pvs.CreateLayout("Individual Conductors")
view_indiv = pvs.CreateRenderView()
pvs.AssignViewToLayout(view=view_indiv, layout=layout_indiv)

elec_readers: dict[int, pvs.PVDReader] = {}
mag_readers:  dict[int, pvs.PVDReader] = {}

for idx in all_indices:
    pvd_e = _find_pvd(pv_dir, idx, "electrostatic")
    if pvd_e:
        r = pvs.PVDReader(FileName=pvd_e)
        pvs.RenameSource(f"Conductor_{idx}_electrostatic", r)
        elec_readers[idx] = r
        d = pvs.Show(r, view_indiv)
        d.Visibility = 0  # hidden by default; user toggles in pipeline browser

    pvd_m = _find_pvd(pv_dir, idx, "magnetostatic")
    if pvd_m:
        r = pvs.PVDReader(FileName=pvd_m)
        pvs.RenameSource(f"Conductor_{idx}_magnetostatic", r)
        mag_readers[idx] = r
        d = pvs.Show(r, view_indiv)
        d.Visibility = 0

# Show the first electrostatic conductor as a quick default
if all_indices and all_indices[0] in elec_readers:
    first_r = elec_readers[all_indices[0]]
    disp = pvs.Show(first_r, view_indiv)
    disp.Visibility = 1
    pvs.ColorBy(disp, ("POINTS", "Phi"))
    pvs.UpdateScalarBars()

pvs.ResetCamera(view_indiv)

# --------------------------------------------------------------------------
# For each material – build combined multipolar views
# --------------------------------------------------------------------------
for material in materials:
    mat_id   = material["id"]
    mat_type = material.get("type", f"material_{mat_id}")
    mp       = material.get("multipolarExpansion", {})

    assoc = next((a for a in mat_assoc if a["materialId"] == mat_id), None)
    if assoc is None:
        print(f"  [warn] No materialAssociation found for material id={mat_id}; skipping.")
        continue

    element_ids: list[int] = assoc["elementIds"]

    e_readers = [elec_readers[i] for i in element_ids if i in elec_readers]
    m_readers = [mag_readers[i]  for i in element_ids if i in mag_readers]

    electric_solutions = mp.get("electric", [])
    magnetic_solutions = mp.get("magnetic", [])

    if not electric_solutions and not magnetic_solutions:
        continue

    layout_name = f"Material_{mat_id}_{mat_type}_Modes"
    layout_mp = pvs.CreateLayout(layout_name)

    # ---- Electric combined view ------------------------------------------
    if electric_solutions and e_readers:
        e_weights = [sol["conductorPotentials"] for sol in electric_solutions]

        e_filter = _make_combined_filter(e_readers, e_weights, "ElectricMode")
        pvs.RenameSource(f"Material_{mat_id}_{mat_type}_electric_combined", e_filter)

        view_e = pvs.CreateRenderView()
        pvs.AssignViewToLayout(view=view_e, layout=layout_mp)
        view_e.ViewSize = [800, 600]

        disp_e = pvs.Show(e_filter, view_e)
        # Default: show first mode
        pvs.ColorBy(disp_e, ("POINTS", "ElectricMode_000"))
        pvs.UpdateScalarBars()
        pvs.ResetCamera(view_e)

        print(
            f"  [material {mat_id}] Created electric combined view with "
            f"{len(electric_solutions)} modes for {len(e_readers)} conductors."
        )

    # ---- Magnetic combined view ------------------------------------------
    if magnetic_solutions and m_readers:
        m_weights = [sol["conductorPotentials"] for sol in magnetic_solutions]

        m_filter = _make_combined_filter(m_readers, m_weights, "MagneticMode")
        pvs.RenameSource(f"Material_{mat_id}_{mat_type}_magnetic_combined", m_filter)

        # Split the layout horizontally and add a second view
        if electric_solutions and e_readers:
            pvs.SplitViewHorizontal(view=view_e, layout=layout_mp)
            view_m = pvs.CreateRenderView()
            pvs.AssignViewToLayout(view=view_m, layout=layout_mp, hint=2)
        else:
            view_m = pvs.CreateRenderView()
            pvs.AssignViewToLayout(view=view_m, layout=layout_mp)

        view_m.ViewSize = [800, 600]
        disp_m = pvs.Show(m_filter, view_m)
        pvs.ColorBy(disp_m, ("POINTS", "MagneticMode_000"))
        pvs.UpdateScalarBars()
        pvs.ResetCamera(view_m)

        print(
            f"  [material {mat_id}] Created magnetic combined view with "
            f"{len(magnetic_solutions)} modes for {len(m_readers)} conductors."
        )

pvs.Render()
print("Tulip ParaView macro completed successfully.")
print("Tip: In the 'Modes' layout, use 'Color By' to switch between ElectricMode_000, "
      "ElectricMode_001, ... arrays to inspect individual multipolar solutions.")
