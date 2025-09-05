from typing import Any, Tuple, List, Dict

import gmsh
from .BoundingBox import BoundingBox
from itertools import chain
import numpy as np

class ShapesClassification:
    isOpenCase:bool


    def __init__(self, shapes):
        gmsh.model.occ.synchronize()

        self.allShapes = shapes
        self.pecs = self.get_surfaces_with_label(shapes, "Conductor_")
        self.dielectrics = self.get_surfaces_with_label(shapes, "Dielectric_")
        self.open = self.get_surfaces_with_label(shapes, "OpenBoundary_")
        self.vacuum = dict()

        self.isOpenCase = self.isOpenProblem()

        if len(self.open) > 1:
            raise ValueError("Only one open region is allowed.")

    @staticmethod
    def getNumberFromName(entity_name: str, label: str):
        ini = entity_name.rindex(label) + len(label)
        num = int(entity_name[ini:])
        return num

    @staticmethod
    def get_surfaces_with_label(entity_tags, label: str):
        surfaces = dict()
        for s in entity_tags:
            name = gmsh.model.get_entity_name(*s)
            if s[0] != 2 or label not in name:
                continue
            num = ShapesClassification.getNumberFromName(name, label)
            surfaces[num] = [s]

        return surfaces

    def isOpenProblem(self):
        elements = list(chain(self.pecs.values()))
        for idx, element in enumerate(elements):
            for otheridx, otherElement in enumerate(elements[idx+1:]):
                if element != otherElement:
                    intersect = gmsh.model.occ.intersect(
                        element, 
                        otherElement,
                        removeObject=False,
                        tag=300+otheridx,
                        removeTool=False
                    )[0]
                    if intersect:
                        return False   
        return True
    
    def removeConductorsFromDielectrics(self):
        for num, diel in self.dielectrics.items():
            pec_surfs = []
            for num2, pec_surf in self.pecs.items():
                if num2 == 0 and not self.isOpenCase:
                    continue
                pec_surfs.extend(pec_surf)
            self.dielectrics[num] = gmsh.model.occ.cut(diel, pec_surfs, removeTool=False)[0]

        gmsh.model.occ.synchronize()

    def ensureDielectricsDoNotOverlap(self):
        for n1, diel1 in self.dielectrics.items():
            others = list(
                chain(
                    *[x[1] for x in self.dielectrics.items() if x[0] != n1]
                )
            )

            if len(others) == 0:
                continue

            self.dielectrics[n1] = gmsh.model.occ.cut(
                self.dielectrics[n1], others, removeObject=True, removeTool=False)[0]

        gmsh.model.occ.synchronize()

    def buildVacuumDomain(self):
        if self.isOpenCase and len(self.open) == 0:
            self.vacuum = self._buildDefaultVacuumDomain()
        elif self.isOpenCase and len(self.open) > 0:
            self.vacuum = self._buildVacuumDomainFromOpenBoundary()
        else:
            self.vacuum = self._buildClosedVacuumDomain()
        return self.vacuum
    
    def _buildVacuumDomainFromOpenBoundary(self) -> Dict[int, List[int]]:
        dom = self.open[0]
        
        surfsToRemove = []
        for num, surf in self.pecs.items():
            surfsToRemove.extend(surf)

        for _, surf in self.dielectrics.items():
            surfsToRemove.extend(surf)

        dom = gmsh.model.occ.cut(
            dom, surfsToRemove, removeObject=False, removeTool=False)[0]
        gmsh.model.occ.synchronize()

        return dict([[0, dom]])
    
    def _buildClosedVacuumDomain(self) -> Tuple[int, int]:
        dom = self.pecs[0]
        surfsToRemove = []
        for num, surf in self.pecs.items():
            if num == 0:
                continue
            surfsToRemove.extend(surf)

        for _, surf in self.dielectrics.items():
            surfsToRemove.extend(surf)
        dom = gmsh.model.occ.cut(
            dom, surfsToRemove, removeObject=False, removeTool=False)[0]
        gmsh.model.occ.synchronize()
        return dict([[0, dom]])
    
    def _buildDefaultVacuumDomain(self):
        NEAR_REGION_BOUNDING_BOX_SCALING_FACTOR = 1.25
        FAR_REGION_DISK_SCALING_FACTOR = 4.0
        nonVacuumSurfaces = []
        for _, surf in self.pecs.items():
            nonVacuumSurfaces.extend(surf)
        for _, surf in self.dielectrics.items():
            nonVacuumSurfaces.extend(surf)
            
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

        farVacuum = gmsh.model.occ.cut(
            farVacuum, nearVacuum, removeObject=True, removeTool=False)[0]

        nearVacuum = gmsh.model.occ.cut(
            nearVacuum, nonVacuumSurfaces, removeObject=True, removeTool=False)[0]
        
        # -- Set mesh size for near vacuum region
        bb = BoundingBox(
            gmsh.model.getBoundingBox(2, nearVacuum[0][1]))
        minSide = np.min(np.array([bb.getLengths()[0], bb.getLengths()[1]]))

        innerRegion = gmsh.model.getBoundary(nearVacuum, recursive=True)
        gmsh.model.mesh.setSize(innerRegion, minSide / 20)
        
        
        self.open = dict([[0, gmsh.model.getBoundary(farVacuum)]])
        gmsh.model.occ.synchronize()

        return dict([[0, nearVacuum], [1, farVacuum]])
    
    
