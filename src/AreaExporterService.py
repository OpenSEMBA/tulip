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

    def addPhysicalModelOfDimension(self, mappedElements:Dict[str,str],  dimension=2):
        physicalGroups = gmsh.model.getPhysicalGroups(dimension)
        for physicalGroup in physicalGroups:
            entityTags = gmsh.model.getEntitiesForPhysicalGroup(*physicalGroup)
            geometryName = gmsh.model.getPhysicalName(*physicalGroup)
            label = ''
            for key, geometry in mappedElements.items():
                if geometry == geometryName:
                    label = key
            for tag in entityTags:
                if dimension == 1:
                    rad = gmsh.model.occ.getMass(dimension, tag) / (2*np.pi)
                    area = rad*rad*np.pi
                if dimension == 2:
                    area = gmsh.model.occ.getMass(dimension, tag)
                if geometryName != AreaExporterService._EMPTY_NAME_CASE:
                    self.addComputedArea(geometryName, label, area)

    def exportToJson(self, exportFileName:str):
        with open(exportFileName + ".areas.json", 'w') as f:
            json.dump(self.computedAreas, f, indent=3)
