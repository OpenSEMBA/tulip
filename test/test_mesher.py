import gmsh
import os
import unittest
import sys
import copy
from typing import Dict, List, Tuple, Optional

from pathlib import Path
from src.mesher import Mesher
from src.ShapesClassification import ShapesClassification

# Add project root to Python path
project_root = str(Path(__file__).parent.parent)
if project_root not in sys.path:
    sys.path.insert(0, project_root)


class TestMesher(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.dirPath = os.path.dirname(os.path.realpath(__file__)) + '/'
        cls.testdataPath = cls.dirPath + '/../testData/'

    @classmethod
    def tearDownClass(cls):
        del cls.dirPath
        del cls.testdataPath

    def setUp(self):
        gmsh.initialize()

    def tearDown(self):
        gmsh.finalize()

    def assertPhysicalGroup(self, expectedNames, expectedEntities):
        pGs = gmsh.model.getPhysicalGroups()
        pGNames = [gmsh.model.getPhysicalName(*pG) for pG in pGs]
        self.assertEqual(sorted(pGNames), sorted(expectedNames))

        for idx, name in enumerate(expectedNames):
            self.assertEqual(self.countEntitiesInPhysicalGroupWithName(
                name), expectedEntities[idx], name)

    def countEntitiesInPhysicalGroupWithName(self, name: str):
        return len(
            gmsh.model.getEntitiesForPhysicalGroup(
                *Mesher.getPhysicalGroupWithName(name)
            )
        )

    def inputFileFromCaseName(self, caseName):
        return self.testdataPath + caseName + '/' + caseName + ".step"

    def test_get_number_from_entity_name(self):
        self.assertEqual(
            ShapesClassification.getNumberFromName(
                'Shapes/Conductor_1',
                'Conductor_'
            ), 1
        )

        self.assertEqual(
            ShapesClassification.getNumberFromName(
                'Shapes/solid_wire_002/Conductor_002/Conductor_002',
                'Conductor_'
            ), 2
        )

    def test_mesh_from_step_with_partially_filled_coax(self):
        caseName = 'partially_filled_coax'

        meshing_options = copy.deepcopy(Mesher.DEFAULT_MESHING_OPTIONS)
        meshing_options["Mesh.ElementOrder"] = 3

        Mesher().meshFromStep(
            self.inputFileFromCaseName(caseName),
            caseName,
            meshing_options)

        # For debugging
        gmsh.write(caseName + '.msh')
        gmsh.write(caseName + '.vtk')

        pGs = gmsh.model.getPhysicalGroups()
        self.assertEqual(len(pGs), 4)

        pGNames = [gmsh.model.getPhysicalName(*pG) for pG in pGs]
        expectedNames = ['Conductor_0',
                         'Conductor_1', 'Dielectric_0', 'Vacuum_0']
        self.assertEqual(sorted(pGNames), sorted(expectedNames))

        for name in expectedNames:
            self.assertEqual(
                self.countEntitiesInPhysicalGroupWithName(name), 1)

    def test_mesh_from_step_with_empty_coax(self):
        caseName = 'empty_coax'
        Mesher().meshFromStep(self.inputFileFromCaseName(caseName), caseName)
        gmsh.write(caseName + '.msh')
        gmsh.write(caseName + '.vtk')

        pGs = gmsh.model.getPhysicalGroups()
        pGNames = [gmsh.model.getPhysicalName(*pG) for pG in pGs]
        expectedNames = ['Conductor_0', 'Conductor_1', 'Vacuum_0']
        self.assertEqual(sorted(pGNames), sorted(expectedNames))

        for name in expectedNames:
            self.assertEqual(
                self.countEntitiesInPhysicalGroupWithName(name), 1)

    def test_mesh_from_step_with_two_wires_coax(self):
        caseName = 'two_wires_coax'
        Mesher().meshFromStep(self.inputFileFromCaseName(caseName), caseName)

        pGs = gmsh.model.getPhysicalGroups()
        pGNames = [gmsh.model.getPhysicalName(*pG) for pG in pGs]
        expectedNames = ['Conductor_0',
                         'Conductor_1', 'Conductor_2', 'Vacuum_0']
        self.assertEqual(sorted(pGNames), sorted(expectedNames))

        for name in expectedNames:
            self.assertEqual(
                self.countEntitiesInPhysicalGroupWithName(name), 1)

    def test_mesh_from_step_with_two_wires_open(self):
        caseName = 'two_wires_open'

        Mesher().meshFromStep(self.inputFileFromCaseName(caseName), caseName)
       # gmsh.write(caseName + '.vtk')

        pGs = gmsh.model.getPhysicalGroups()
        pGNames = [gmsh.model.getPhysicalName(*pG) for pG in pGs]
        expectedNames = ['Conductor_0', 'Conductor_1',
                         'OpenBoundary_0', 'Vacuum_0', 'Vacuum_1']
        self.assertEqual(sorted(pGNames), sorted(expectedNames))

        expectedEntities = [1, 1, 1, 1, 1]

        for idx, name in enumerate(expectedNames):
            self.assertEqual(self.countEntitiesInPhysicalGroupWithName(
                name), expectedEntities[idx], name)

    def test_mesh_dielectric_pair_open(self) -> None:
        caseName = 'DielectricUnshieldedPair'
        expectedNames = [
            'Conductor_0', 'Conductor_1',
            'Dielectric_0', 'Dielectric_1',
            'OpenBoundary_0', 'Vacuum_0', 'Vacuum_1']
        expectedEntities = [1, 1, 1, 1, 1, 1, 1]

        Mesher().meshFromStep(self.inputFileFromCaseName(caseName), caseName)
       # gmsh.write(caseName + '.vtk')

        self.assertPhysicalGroup(expectedNames, expectedEntities)

    def test_mesh_dielectric_pair_open_defined_boundary(self) -> None:
        caseName = 'DielectricUnshieldedPairDefinedBoundary'
    
        meshing_options = copy.deepcopy(Mesher.DEFAULT_MESHING_OPTIONS)
        meshing_options["Mesh.ElementOrder"] = 1
        
        Mesher().meshFromStep(
            self.inputFileFromCaseName(caseName), 
            caseName,
            meshingOptions=meshing_options)

        gmsh.write(caseName + '.vtk')
        
        expectedNames = [
            'Conductor_0', 'Conductor_1',
            'Dielectric_0', 'Dielectric_1',
            'OpenBoundary_0', 'Vacuum_0']
        expectedEntities = [1, 1, 1, 1, 1, 1]

        self.assertPhysicalGroup(expectedNames, expectedEntities)

    def test_mesh_from_step_with_five_wires(self):
        expectedNames = [
            'Conductor_0', 'Conductor_1',
            'Conductor_2', 'Conductor_3',
            'Conductor_4', 'Conductor_5',
            'Dielectric_0', 'Dielectric_1',
            'Dielectric_2', 'Dielectric_3',
            'Dielectric_4', 'Vacuum_0'
        ]

        caseName = 'five_wires'
        Mesher().meshFromStep(self.inputFileFromCaseName(caseName), caseName)

        pGs = gmsh.model.getPhysicalGroups()
        pGNames = [gmsh.model.getPhysicalName(*pG) for pG in pGs]
        self.assertEqual(sorted(pGNames), sorted(expectedNames))

        for name in expectedNames:
            self.assertEqual(
                self.countEntitiesInPhysicalGroupWithName(name), 1)

    def test_mesh_from_step_with_three_wires_ribbon(self):
        caseName = 'three_wires_ribbon'
        Mesher().meshFromStep(self.inputFileFromCaseName(caseName), caseName)

       # gmsh.write(caseName + '.vtk')

        pGs = gmsh.model.getPhysicalGroups()
        pGNames = [gmsh.model.getPhysicalName(*pG) for pG in pGs]
        expectedNames = [
            'Conductor_0', 'Conductor_1', 'Conductor_2',
            'OpenBoundary_0',
            'Dielectric_0', 'Dielectric_1', 'Dielectric_2',
            'Vacuum_0', 'Vacuum_1'
        ]
        expectedEntities = [1, 1, 1,
                            1,
                            1, 1, 1,
                            1, 1]
        self.maxDiff = None
        self.assertEqual(sorted(pGNames), sorted(expectedNames))

        for idx, name in enumerate(expectedNames):
            self.assertEqual(self.countEntitiesInPhysicalGroupWithName(
                name), expectedEntities[idx], name)

    def test_mesh_from_step_with_nested_coax(self):
        caseName = 'nested_coax'
        Mesher().meshFromStep(self.inputFileFromCaseName(caseName), caseName)

        pGs = gmsh.model.getPhysicalGroups()
        pGNames = [gmsh.model.getPhysicalName(*pG) for pG in pGs]
        expectedNames = ['Conductor_0',
                         'Conductor_1', 'Conductor_2', 'Vacuum_0']
        self.assertEqual(sorted(pGNames), sorted(expectedNames))

        self.assertEqual(
            self.countEntitiesInPhysicalGroupWithName('Conductor_0'), 1)
        self.assertEqual(
            self.countEntitiesInPhysicalGroupWithName('Conductor_1'), 2)
        self.assertEqual(
            self.countEntitiesInPhysicalGroupWithName('Conductor_2'), 1)
        self.assertEqual(
            self.countEntitiesInPhysicalGroupWithName('Vacuum_0'), 2)

    def test_mesh_from_step_with_agrawal1981(self):
        caseName = 'agrawal1981'
        Mesher().meshFromStep(self.inputFileFromCaseName(caseName), caseName)
        # gmsh.write(caseName + '.msh')
        # gmsh.write(caseName + '.vtk')
        # gmsh.fltk.run()

        pGs = gmsh.model.getPhysicalGroups()
        pGNames = [gmsh.model.getPhysicalName(*pG) for pG in pGs]
        expectedNames = ['Conductor_0', 'Conductor_1', 'Conductor_2', 'Conductor_3',
                         'OpenBoundary_0',
                         'Dielectric_1', 'Dielectric_2', 'Dielectric_0',
                         'Vacuum_0']
        expectedEntities = [4, 1, 1, 1,
                            1,
                            1, 1, 1,
                            2]

        for idx, name in enumerate(expectedNames):
            self.assertEqual(self.countEntitiesInPhysicalGroupWithName(
                name), expectedEntities[idx], name)

    def test_unshielded_multiwire(self):
        caseName = 'unshielded_multiwire'
        Mesher().meshFromStep(self.inputFileFromCaseName(caseName), caseName)

        # gmsh.write(caseName + '.vtk')

        pGs = gmsh.model.getPhysicalGroups()
        pGNames = [gmsh.model.getPhysicalName(*pG) for pG in pGs]
        expectedNames = ['Conductor_0', 'Conductor_1',  'Dielectric_0',
                         'OpenBoundary_0',
                         'Vacuum_0', 'Vacuum_1']
        expectedEntities = [1, 1, 1,
                            1,
                            1, 1]
        self.assertEqual(sorted(pGNames), sorted(expectedNames))

        for idx, name in enumerate(expectedNames):
            self.assertEqual(self.countEntitiesInPhysicalGroupWithName(
                name), expectedEntities[idx], name)

    def test_conductor_and_outer_dielectric(self):
        caseName = 'conductor_and_outer_dielectric'
        Mesher().meshFromStep(self.inputFileFromCaseName(caseName), caseName)

        gmsh.write(caseName + '.vtk')
        gmsh.write(caseName + '.msh')

        pGs = gmsh.model.getPhysicalGroups()
        pGNames = [gmsh.model.getPhysicalName(*pG) for pG in pGs]
        expectedNames = ['Conductor_0',
                         'Dielectric_0',
                         'OpenBoundary_0',
                         'Vacuum_0', 'Vacuum_1']
        expectedEntities = [1,
                            1,
                            1,
                            2, 1]
        self.assertEqual(sorted(pGNames), sorted(expectedNames))

        for idx, name in enumerate(expectedNames):
            self.assertEqual(self.countEntitiesInPhysicalGroupWithName(
                name), expectedEntities[idx], name)

    def test_realistic_case_with_dielectrics_fdtd_cell(self):
        caseName = 'realistic_case_with_dielectrics_fdtd_cell'
        meshing_options = copy.deepcopy(Mesher.DEFAULT_MESHING_OPTIONS)
        meshing_options["Mesh.ElementOrder"] = 1

        Mesher().meshFromStep(
            self.inputFileFromCaseName(caseName),
            caseName,
            meshingOptions=meshing_options)

        # For debugging.
        gmsh.write(caseName + '.vtk')
        gmsh.write(caseName + '.msh')
        # gmsh.fltk.run()

        pGs = gmsh.model.getPhysicalGroups()
        pGNames = [gmsh.model.getPhysicalName(*pG) for pG in pGs]

        expectedNames = ['Conductor_0', 'Conductor_1', 'Conductor_2', 'Conductor_3', 'Conductor_4',
                         'Conductor_5', 'Conductor_6', 'Conductor_7', 'Conductor_8', 'Conductor_9',
                         'Conductor_10', 'Conductor_11', 'Conductor_12', 'Conductor_13', 'Conductor_14',
                         'Conductor_15', 'Conductor_16', 'Conductor_17', 'Conductor_18', 'Conductor_19',
                         'Conductor_20', 'Conductor_21', 'Conductor_22', 'Conductor_23', 'Conductor_24',
                         'Conductor_25', 'Conductor_26', 'Conductor_27', 'Conductor_28', 'Conductor_29',
                         'Conductor_30',
                         'Dielectric_0', 'Dielectric_1', 'Dielectric_2', 'Dielectric_3', 'Dielectric_4',
                         'Dielectric_5', 'Dielectric_6', 'Dielectric_7', 'Dielectric_8', 'Dielectric_9',
                         'Dielectric_10', 'Dielectric_11', 'Dielectric_12', 'Dielectric_13', 'Dielectric_14',
                         'Dielectric_15', 'Dielectric_16', 'Dielectric_17', 'Dielectric_18', 'Dielectric_19',
                         'Dielectric_20', 'Dielectric_21', 'Dielectric_22', 'Dielectric_23', 'Dielectric_24',
                         'Dielectric_25', 'Dielectric_26', 'Dielectric_27', 'Dielectric_28', 'Dielectric_29',
                         'Dielectric_30', 'Dielectric_31',
                         'OpenBoundary_0',
                         'Vacuum_0']
        self.assertEqual(sorted(pGNames), sorted(expectedNames))

    def test_lansink2024_single_wire_multipolar(self):
        caseName = 'lansink2024_single_wire_multipolar'
        Mesher().meshFromStep(self.inputFileFromCaseName(caseName), caseName)

        gmsh.write(caseName + '.msh')
        gmsh.write(caseName + '.vtk')

        pGs = gmsh.model.getPhysicalGroups()
        pGNames = [gmsh.model.getPhysicalName(*pG) for pG in pGs]
        expectedNames = ['Conductor_0',
                         'Dielectric_0',
                         'OpenBoundary_0',
                         'Vacuum_0', 'Vacuum_1']
        expectedEntities = [1, 1, 1,
                            1,
                            1, 1]
        self.assertEqual(sorted(pGNames), sorted(expectedNames))

        for idx, name in enumerate(expectedNames):
            self.assertEqual(self.countEntitiesInPhysicalGroupWithName(
                name), expectedEntities[idx], name)


if __name__ == '__main__':
    unittest.main()
