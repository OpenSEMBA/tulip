import json
import gmsh
from typing import Dict, List
import numpy as np
class AreaExporterService:
    _EMPTY_NAME_CASE = ""
    computedAreas:Dict[str,List]
    geometry: Dict
    def __init__(self):
        self.computedAreas = {
            "geometries": []
        }
    
    def addComputedArea(self, geometry:str, label:str, area:float):
        geometry:Dict ={
            "geometry": geometry,
            "label": label,
            "area": round(area,6),
        }
        self.computedAreas['geometries'].append(geometry)

    def addPhysicalModelForConductors(self, mappedElements:Dict[str,str]):
        physicalGroups = gmsh.model.getPhysicalGroups(1)
        for physicalGroup in physicalGroups:
            entityTags = gmsh.model.getEntitiesForPhysicalGroup(*physicalGroup)
            geometryName = gmsh.model.getPhysicalName(*physicalGroup)
            if not geometryName.startswith("Conductor_"):
                continue

            label = ''
            for key, geometry in mappedElements.items():
                if geometry == geometryName:
                    label = key
                    break
            
            # Find surface that has these curves as boundaries
            allSurfaces = gmsh.model.getEntities(2)
            foundSurface = None
            for surface in allSurfaces:
                boundary = gmsh.model.getBoundary([surface], oriented=False, recursive=False)
                boundaryTags = set(tag for dim, tag in boundary)
                if set(entityTags) == boundaryTags:
                    foundSurface = surface
                    break
            
            if foundSurface:
                area = gmsh.model.occ.getMass(2, foundSurface[1])
                self.addComputedArea(geometryName, label, area)
    
    def exportToJson(self, exportFileName:str):
        with open(exportFileName + ".areas.json", 'w') as f:
            json.dump(self.computedAreas, f, indent=3)
