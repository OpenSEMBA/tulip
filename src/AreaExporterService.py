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
            if geometryName.startswith("Conductor_"):
                label = ''
                for key, geometry in mappedElements.items():
                    if geometry == geometryName:
                        label = key
                        break
                
                    assert len(entityTags) == 1
                    loop = gmsh.model.geo.addCurveLoop(entityTags)
                    surface = gmsh.model.geo.addPlaneSurface([loop])

                    gmsh.model.geo.synchronize()

                    area = gmsh.model.occ.getMass(2, surface)
                    gmsh.model.occ.remove([(2, surface)])
                    if geometryName != AreaExporterService._EMPTY_NAME_CASE:
                        self.addComputedArea(geometryName, label, area)

    def exportToJson(self, exportFileName:str):
        with open(exportFileName + ".areas.json", 'w') as f:
            json.dump(self.computedAreas, f, indent=3)
