from typing import List
import unittest
import os
import numpy as np
import gmsh
from src.mesher import Mesher
from src.AreaExporterService import AreaExporterService


class testAreaExporterService(unittest.TestCase):
    @staticmethod
    def sumAreasFromList(areas: List[float]):
        total: float = 0
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
        mappedElements = Mesher().meshFromStep(
            self.inputFileFromCaseName(caseName), caseName)
        areaExporter = AreaExporterService()
        areaExporter.addPhysicalModelForConductors(mappedElements=mappedElements)
        geometries = areaExporter.computedAreas['geometries']

        expectedDict = {
            'geometries': [
                {
                    'area': 28.274334,
                    'geometry': 'Conductor_0',
                    'label': 'Conductor_0'
                },
                {
                    'area': 0.785398,
                    'geometry': 'Conductor_5',
                    'label': 'Conductor_1'
                },
                {
                    'area': 2.010619,
                    'geometry': 'Conductor_1',
                    'label': 'Conductor_002'
                },
                {
                    'area': 0.785398,
                    'geometry': 'Conductor_2',
                    'label': 'Conductor_003'
                },
                {
                    'area': 0.785398,
                    'geometry': 'Conductor_3',
                    'label': 'Conductor_004'
                },
                {
                    'area': 0.785398,
                    'geometry': 'Conductor_4',
                    'label': 'Conductor_005'
                },
            ]
        }

        self.maxDiff = None
        self.assertDictEqual(areaExporter.computedAreas, expectedDict)

    def testJsonFormat(self) -> None:
        caseName = 'DielectricUnshieldedPair'
        mappedElements = Mesher().meshFromStep(
            self.inputFileFromCaseName(caseName), caseName)
        areaExporter = AreaExporterService()
        areaExporter.addPhysicalModelForConductors(
            mappedElements=mappedElements)

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
                    'label': 'LeftConductor'
                }
            ]
        }
        self.maxDiff = None
        self.assertDictEqual(areaExporter.computedAreas, expectedDict)
