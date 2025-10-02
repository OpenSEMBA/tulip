from typing import List
import unittest
import os
import numpy as np
import gmsh
from src.mesher import Mesher
from src.AreaExporterService import AreaExporterService
class testAreaExporterService(unittest.TestCase):
    @staticmethod
    def sumAreasFromList(areas:List[float]):
        total:float = 0
        for area in areas:
            total += area
        return total

    @classmethod
    def setUpClass(cls):
        cls.dirPath = os.path.dirname(os.path.realpath(__file__)) + '/'
        cls.testdataPath = cls.dirPath + '/../testData/'

    def setUp(self):
        gmsh.initialize()

    def tearDown(self):
        gmsh.finalize()
    def inputFileFromCaseName(self, caseName) -> None:
        return self.testdataPath + caseName + '/' + caseName + ".step"

    def testAreaExporterReturnsTrueValues(self):
        caseName = 'five_wires'
        mappedElements = Mesher().meshFromStep(self.inputFileFromCaseName(caseName), caseName)
        areaExporter = AreaExporterService()
        areaExporter.addPhysicalModelOfDimension(mappedElements=mappedElements, dimension=1)
        areaExporter.addPhysicalModelOfDimension(mappedElements=mappedElements, dimension=2)
        geometries = areaExporter.computedAreas['geometries']

        internalElements = []
        for geometry in geometries:
            if geometry['geometry'] == "Conductor_0":
                totalArea = geometry['area']
            else:
                internalElements.append(geometry['area'])
        areaElements = self.sumAreasFromList(internalElements)

        self.assertAlmostEqual(totalArea, areaElements)

    def testJsonFormat(self) -> None:
        caseName = 'DielectricUnshieldedPair'
        mappedElements = Mesher().meshFromStep(self.inputFileFromCaseName(caseName), caseName)
        areaExporter = AreaExporterService()
        areaExporter.addPhysicalModelOfDimension(mappedElements=mappedElements, dimension=1)
        areaExporter.addPhysicalModelOfDimension(mappedElements=mappedElements, dimension=2)

        expectedDict = {
            'geometries': [
                {
                    'area': 201.06193,
                    'geometry': 'Conductor_1',
                    'label': 'RightConductor'
                },
                {
                    'area': 201.06193,
                    'geometry': 'Conductor_0',
                    'label': 'LeftConductor'},
                {
                    'area': 312048.117187,
                    'geometry': 'OpenBoundary_0',
                    'label': 'OpenBoundary_0'
                },
                {
                    'area': 603.185789,
                    'geometry': 'Dielectric_1',
                    'label': 'RightDielectric'
                },
                {
                    'area': 603.185789,
                    'geometry': 'Dielectric_0',
                    'label': 'LeftDielectric'
                },
                {
                    'area': 6491.504606,
                    'geometry': 'Vacuum_0',
                    'label': 'Vacuum_0'
                },
                {
                    'area': 303948.117142,
                    'geometry': 'Vacuum_1',
                    'label': 'Vacuum_1'
                }
            ]
        }
        self.maxDiff = None
        self.assertDictEqual(areaExporter.computedAreas, expectedDict)
