from copy import copy
import os
from typing import Dict, List, Tuple
import unittest
import gmsh
import json


from src.ShapesClassification import ShapesClassification


class TestShapesClassification(unittest.TestCase):

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

    def inputFileFromCaseName(self, caseName):
        return self.testdataPath + caseName + '/' + caseName + ".step"

    def initShapeClassification(self, inputFile:str) -> None:
        jsonFile = inputFile.strip('.step') + '.json'
        self.shapeClassification = ShapesClassification(
            gmsh.model.occ.importShapes(inputFile, highestDimOnly=False),
            jsonFile
        )

    def testDielectricShieldedPairClassification(self) -> None:
        case = 'DielectricShieldedPair'
        filepath = self.inputFileFromCaseName(case)
        self.initShapeClassification(filepath)
        expectedShapes = [
            (2, 1),(2, 2),(2, 3),(2, 4),(2, 5),
            (1, 1),(1, 2),(1, 3),(1, 4),(1, 5),
            (0, 1),(0, 1),(0, 2),(0, 2),(0, 3),(0, 3),(0, 4),(0, 4),(0, 5),(0, 5)
        ]
        expectedPecs = {
            'RightConductor': [(2,1)],
            'ExternalShield': [(2,2)],
            'LeftConductor': [(2,3)],
        }
        expectedDielectrics = {
            'RightDielectric': [(2,4)],
            'LeftDielectric': [(2,5)],
        }
        self.assertListEqual(self.shapeClassification.allShapes, expectedShapes)
        self.assertDictEqual(self.shapeClassification.pecs, expectedPecs)
        self.assertDictEqual(self.shapeClassification.dielectrics, expectedDielectrics)
        self.assertFalse(self.shapeClassification.isOpenCase)

    def testFusedConductors(self) -> None:
        case = 'FusedConductor'
        filepath = self.inputFileFromCaseName(case)
        self.initShapeClassification(filepath)

    def testComplexNesting(self) -> None:
        case = 'ComplexNesting'
        filepath = self.inputFileFromCaseName(case)
        self.initShapeClassification(filepath)

    def testDielectricUnshieldedPairClassification(self) -> None:
        case = 'DielectricUnshieldedPair'
        filepath = self.inputFileFromCaseName(case)
        self.initShapeClassification(filepath)
        expectedShapes = [
            (2, 1),(2, 2),(2, 3),(2, 4),
            (1, 1),(1, 2),(1, 3),(1, 4),
            (0, 1),(0, 1),(0, 2),(0, 2),(0, 3),(0, 3),(0, 4),(0, 4),
        ]
        expectedPecs = {
            'LeftConductor': [(2, 2)],
            'RightConductor': [(2, 1)],
        }
        expectedDielectrics = {
            'RightDielectric': [(2,3)],
            'LeftDielectric': [(2,4)],
        }
        self.assertListEqual(self.shapeClassification.allShapes, expectedShapes)
        self.assertDictEqual(self.shapeClassification.pecs, expectedPecs)
        self.assertDictEqual(self.shapeClassification.dielectrics, expectedDielectrics)
        self.assertTrue(self.shapeClassification.isOpenCase)