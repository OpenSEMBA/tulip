from enum import Enum
import json
import math
from time import sleep
from typing import Any, Tuple, List, Dict

import gmsh

from src.Graph import Graph
from .BoundingBox import BoundingBox
from itertools import chain
import numpy as np

class ShapesClassification:
    _ROUND_VALUE:int = 6

    isOpenCase:bool
    crossSectionData: Dict[str,List[Dict[str,any]]]
    pecs: Dict[str,List[Tuple[int,int]]]
    dielectrics: Dict[str,List[Tuple[int,int]]]
    nestedGraph: Graph

    def __init__(self, shapes, jsonFile:str):
        gmsh.model.occ.synchronize()

        self.allShapes = shapes
        with open(jsonFile) as f:
            jsonData = json.load(f)
        self.crossSectionData = jsonData['CrossSection']
        self.pecs = self.get_pecs(shapes)
        self.dielectrics = self.get_dielectrics(shapes)
        self.open = self.get_open_boundaries(shapes)
        self.vacuum = dict()
        self.nestedGraph = self.__getNestedGraph()
        self.isOpenCase = self.isOpenProblem()

    @staticmethod
    def getNumberFromName(entity_name: str, label: str):
        ini = entity_name.rindex(label) + len(label)
        num = int(entity_name[ini:])
        return num

    def get_pecs(self, entity_tags) -> Dict[str, Dict[str,any]]:
        pecNames = self.__getGeometryNamesByMaterialType('PEC')
        pecs = dict()
        for s in entity_tags:
            name = gmsh.model.get_entity_name(*s).split('/')[-1]
            if s[0] != 2 or name not in pecNames:
                continue
            pecs[name] = [s]

        return pecs
    
    def get_dielectrics(self, entity_tags) -> Dict[str, Dict[str,any]]:
        dielectricNames = self.__getGeometryNamesByMaterialType('Dielectric')
        dielectrics = dict()
        for s in entity_tags:
            name = gmsh.model.get_entity_name(*s).split('/')[-1]
            if s[0] != 2 or name not in dielectricNames:
                continue
            dielectrics[name] = [s]

        return dielectrics
    
    def get_open_boundaries(self, entity_tags) -> Dict[str, List[Tuple[int,int]]]:
        openBoundaryNames = self.__getGeometryNamesByMaterialType('OpenBoundary')
        open_boundaries = dict()
        for s in entity_tags:
            name = gmsh.model.get_entity_name(*s).split('/')[-1]
            if s[0] != 2 or name not in openBoundaryNames: 
                continue
            open_boundaries[name] = [s]

        assert len(open_boundaries) <= 1

        return open_boundaries
    
    def __getGeometryNamesByMaterialType(self, materialType:str) -> List[str]:
        names = [
            geometry['name'] 
            for geometry in self.crossSectionData 
            if geometry['material']['type'] == materialType
        ]
        return names
    
    def isOpenBoundaryDefined(self) -> bool:
        return len(self.open) > 0
    
    def isOpenProblem(self) -> bool:
        roots = self.nestedGraph.roots 
        if len(self.open) == 1:
            return True
        if len(roots) > 1: 
            return True
        if roots[0] in self.dielectrics.keys(): 
            return True
        return False
    
    def removeConductorsFromDielectrics(self):
        for num, diel in self.dielectrics.items():
            pec_surfs = []
            for num2, pec_surf in self.pecs.items():
                if (num2 in self.nestedGraph.roots) and (not self.isOpenCase):
                    continue
                pec_surfs.extend(pec_surf)
            self.dielectrics[num] = gmsh.model.occ.cut(diel, pec_surfs, removeTool=False)[0]

        gmsh.model.occ.synchronize()

    def ensureDielectricsDoNotOverlap(self):
        for currentKey in self.dielectrics.keys():

            others = list(chain(*[tag for key, tag in self.dielectrics.items() if currentKey != key]))

            if len(others) == 0:
                continue

            self.dielectrics[currentKey] = gmsh.model.occ.cut(
                self.dielectrics[currentKey], others, removeObject=True, removeTool=False
                )[0]

        gmsh.model.occ.synchronize()

    def buildVacuumDomain(self):
        if self.isOpenCase:
            self.vacuum = self._buildOpenVacuumDomain()
        else:
            self.vacuum = self._buildClosedVacuumDomain()
        return self.vacuum
    
    def _buildClosedVacuumDomain(self) -> Tuple[int, int]:
        root = self.nestedGraph.roots[0]
        dom = self.pecs[root]
        surfsToRemove = []
        for num, surf in self.pecs.items():
            if num == root:
                continue
            surfsToRemove.extend(surf)

        for _, surf in self.dielectrics.items():
            surfsToRemove.extend(surf)
        dom = gmsh.model.occ.cut(
            dom, surfsToRemove, removeObject=False, removeTool=False)[0]
        gmsh.model.occ.synchronize()
        return dict([['Vacuum_0', dom]])
    
    def _buildOpenVacuumDomain(self):
        NEAR_REGION_BOUNDING_BOX_SCALING_FACTOR = 1.15
        FAR_REGION_DISK_SCALING_FACTOR = 4.0
        nonVacuumSurfaces = []
        for _, surf in self.pecs.items():
            nonVacuumSurfaces.extend(surf)
        for _, surf in self.dielectrics.items():
            nonVacuumSurfaces.extend(surf)

        if self.isOpenBoundaryDefined():
            name, vacuum = next(iter(self.open.items()))
            
            vacuum = gmsh.model.occ.cut(vacuum, nonVacuumSurfaces,
                removeObject=True, removeTool=False)[0]
            gmsh.model.occ.synchronize()

            vacuumBoundaries = gmsh.model.getBoundary(vacuum)
            externalVacuumBoundaries = [dt for dt in vacuumBoundaries if dt[1]>0]
            self.open = dict([[name, externalVacuumBoundaries]])

            return dict([['Vacuum_0', vacuum]])
        else:            
            boundingBox = BoundingBox.getBoundingBoxFromGroup(nonVacuumSurfaces)

            bbMaxLength = np.max(boundingBox.getLengths())
            nearVacuumBoxSize = bbMaxLength*NEAR_REGION_BOUNDING_BOX_SCALING_FACTOR
            nVOrigin = tuple(
                np.subtract(boundingBox.getCenter(), 
                            (nearVacuumBoxSize/2.0, nearVacuumBoxSize/2.0, 0.0)))
            nearVacuum = [
                (2, gmsh.model.occ.addRectangle(*nVOrigin, *(nearVacuumBoxSize,)*2))
            ]

            farVacuumDiameter = FAR_REGION_DISK_SCALING_FACTOR * boundingBox.getDiagonal()
            farVacuum = [(2, gmsh.model.occ.addDisk(
                *boundingBox.getCenter(), 
                farVacuumDiameter, farVacuumDiameter))]

            gmsh.model.occ.synchronize()
            self.open = dict([['OpenBoundary_0', gmsh.model.getBoundary(farVacuum)]])

            farVacuum = gmsh.model.occ.cut(
                farVacuum, nearVacuum, removeObject=True, removeTool=False)[0]
            nearVacuum = gmsh.model.occ.cut(
                nearVacuum, nonVacuumSurfaces, removeObject=True, removeTool=False)[0]
        
            gmsh.model.occ.synchronize()
            
            # -- Set mesh size for near vacuum region
            bb = BoundingBox(
                gmsh.model.getBoundingBox(2, nearVacuum[0][1]))
            minSide = np.min(np.array([bb.getLengths()[0], bb.getLengths()[1]]))

            innerRegion = gmsh.model.getBoundary(nearVacuum, recursive=True)
            gmsh.model.mesh.setSize(innerRegion, minSide / 20)
            
            gmsh.model.occ.synchronize()

            return dict([['Vacuum_0', nearVacuum], ['Vacuum_1', farVacuum]])
        

    
    def __getNestedGraph(self):
        gmsh.model.occ.synchronize()
        graph = Graph()
        elements:Dict = {}
        elements = {**self.pecs, **self.dielectrics, **self.open}
        for key in elements:
            graph.add_node(key)
        for i, keyA in enumerate(elements):
            for j, keyB in enumerate(elements):
                if i < j:
                    inter = gmsh.model.occ.intersect(
                        elements[keyA], 
                        elements[keyB],
                        removeObject=False,
                        removeTool=False
                    )
                    if len(inter[1][0]) == 0: 
                        continue
                    else:
                        if inter[1][0] == elements[keyA]:
                            graph.add_edge(keyB, keyA)
                        elif inter[1][0] == elements[keyB]:
                            graph.add_edge(keyA, keyB)
        graph.prune_to_longest_paths()
        graph._reorderData()
        return graph
    
    def getMappedComponents(self) -> Dict[str,str]:
        
        mappedElements = []
        
        conductors = []
        sortedNodes = self.nestedGraph.getNodesByLevels()
        for node in sortedNodes:
            if node in self.pecs.keys():
                conductors.append((node, 'Conductor_{}'.format(len(conductors))))
        mappedElements.extend(conductors)
        
        dielectrics = []
        for node in self.dielectrics.keys():
            dielectrics.append((node, 'Dielectric_{}'.format(len(dielectrics))))
        mappedElements.extend(dielectrics)

        mappedComponents = {element[0]:element[1] for element in mappedElements}

        for domain in self.vacuum.keys():
            mappedComponents[domain] = domain
        for openBoundary in self.open.keys():
            mappedComponents[openBoundary] = 'OpenBoundary_0'

        return mappedComponents