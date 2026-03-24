from typing import Dict
from src.mesher import ShapesClassification
from src.utils import assertListOfFloatsAlmostEqual
from src.BoundingBox import BoundingBox
import gmsh
import unittest

class testOpenCase(unittest.TestCase):
    _TEST_MODEL_NAME: str = "Test_Model"
    _DEBUG_MODE = False #Set to true if gmsh ui needed for debug

    def setUp(self):
        self.meshing_options = {
    
            "Mesh.MshFileVersion": 2.2,   # Mandatory for MFEM compatibility
            "Mesh.MeshSizeFromCurvature": 20,
            "Mesh.ElementOrder": 3,
            "Mesh.ScalingFactor": 1e-3,
            "Mesh.SurfaceFaces": 1,
            "Mesh.MeshSizeMax": 20,

            "General.DrawBoundingBoxes": 1,
            "General.Axes": 1,

            "Geometry.SurfaceType": 2,    # Diplay surfaces as solids rather than dashed lines.
            # "Geometry.OCCBoundsUseStl": 1,
            # "Geometry.OCCSewFaces": 1,
            # "Geometry.Tolerance": 1e-3,
        }

    def tearDown(self):
        if self._DEBUG_MODE:
            gmsh.fltk.run()
        gmsh.finalize()