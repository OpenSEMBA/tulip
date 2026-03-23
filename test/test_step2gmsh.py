import os
from ..step2gmsh import launcher
import unittest

import sys

sys.path.insert(0, '.')
class testStep2gmsh(unittest.TestCase):
    dir_path = os.path.dirname(os.path.realpath(__file__)) + '/'
    testdata_path = dir_path + '/../testData/'

    def test_launcher(self):
        case_name = 'partially_filled_coax'
        input = self.testdata_path + case_name + '/' + case_name + '.step'
        launcher(input)
