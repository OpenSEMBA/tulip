from typing import Tuple, List, Dict
import gmsh
from . import utils

class BoundingBox():
    edges: Dict[str, float]
    def __init__(self, listOfCoordinates:Tuple[float,float,float,float,float,float]):
        self.edges = {
                'XMin': listOfCoordinates[0],
                'YMin': listOfCoordinates[1],
                'ZMin': listOfCoordinates[2],
                'XMax': listOfCoordinates[3],
                'YMax': listOfCoordinates[4],
                'ZMax': listOfCoordinates[5],
            }
        
    def getOrigin(self) -> Tuple[float,float,float]:
        return (
            self.edges['XMin'],
            self.edges['YMin'],
            self.edges['ZMin']
        )

    def getCenter(self) -> Tuple[float,float,float]:
        return (
            (self.edges['XMax'] + self.edges['XMin']) / 2,
            (self.edges['YMax'] + self.edges['YMin']) / 2,
            (self.edges['ZMax'] + self.edges['ZMin']) / 2
        )
    def getDiagonal(self) -> float:
        dx = self.edges['XMax'] - self.edges['XMin']
        dy = self.edges['YMax'] - self.edges['YMin']
        dz = self.edges['ZMax'] - self.edges['ZMin']
        return (dx**2 + dy**2 + dz**2) ** 0.5
    
    def getLengths(self) -> Tuple[float, float, float]:
        return (
            self.edges['XMax'] - self.edges['XMin'],
            self.edges['YMax'] - self.edges['YMin'],
            self.edges['ZMax'] - self.edges['ZMin']
        )
    
    @staticmethod
    def _getBoundingBox(element:Tuple[int,int]) -> 'BoundingBox':
        boundingBox:BoundingBox = BoundingBox(gmsh.model.occ.get_bounding_box(*element))
        return boundingBox

    @staticmethod
    def getBoundingBoxFromGroup(elements:List[Tuple[int,int]]) -> 'BoundingBox':
        boundingBoxs:List[BoundingBox] = []
        for element in elements:
            boundingBoxs.append(BoundingBox._getBoundingBox(element))
        
        if len(boundingBoxs) != 0:
            edges: Dict[str, List[float]] = {
                'XMin': [],
                'YMin': [],
                'ZMin': [],
                'XMax': [],
                'YMax': [],
                'ZMax': [],
            }
            for boundingBox in boundingBoxs:
                for key in edges.keys():
                    edges[key].append(boundingBox.edges[key])
            
            return BoundingBox(
                (
                    utils.getMinFromList(edges['XMin']),
                    utils.getMinFromList(edges['YMin']),
                    utils.getMinFromList(edges['ZMin']),
                    utils.getMaxFromList(edges['XMax']),
                    utils.getMaxFromList(edges['YMax']),
                    utils.getMaxFromList(edges['ZMax'])
                )
            )
        else:
            return BoundingBox((0,0,0,0,0,0))