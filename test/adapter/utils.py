import numpy as np
from typing import List, Union

def assertListOfFloatsAlmostEqual(realValues, expectedValues, tolerance = 0.0000001):
    if len(realValues) != len(expectedValues):
        raise AssertionError("List have diferent lengths: {} != {}".format(len(realValues), len(expectedValues)))
    
    if not np.allclose(realValues, expectedValues, atol=tolerance):
        raise AssertionError("List are not the same: {} != {}".format(realValues, expectedValues))

def getMinFromList(values: List[Union[float, int]]) -> Union[float, int]:
    if not values:
        raise ValueError("List cannot be empty")
    
    for value in values:
        if not isinstance(value, (int, float)):
            raise TypeError(f"List must contain only numbers. Found invalid value: {value} of type {type(value)}")
    
    return min(values)

def getMaxFromList(values: List[Union[float, int]]) -> Union[float, int]:
    if not values:
        raise ValueError("List cannot be empty")
    
    for value in values:
        if not isinstance(value, (int, float)):
            raise TypeError(f"List must contain only numbers. Found invalid value: {value} of type {type(value)}")
    
    return max(values)

