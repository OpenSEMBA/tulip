from typing import Tuple
import gmsh
from pathlib import Path
from typing import Dict

from src.AreaExporterService import AreaExporterService
from .ShapesClassification import ShapesClassification
from .BoundingBox import BoundingBox
import numpy as np

class Mesher():
    DEFAULT_MESHING_OPTIONS = {
    
        "Mesh.MshFileVersion": 2.2,   # Required for MFEM compatibility
        "Mesh.MeshSizeFromCurvature": 50,
        "Mesh.ElementOrder": 3,
        "Mesh.ScalingFactor": 1e-3,
        "Mesh.SurfaceFaces": 1,
        "Mesh.MeshSizeMax": 40,

        "General.DrawBoundingBoxes": 1,
        "General.Axes": 1,

        "Geometry.SurfaceType": 2,    # Display surfaces as solids rather than dashed lines.
        # "Geometry.OCCBoundsUseStl": 1,
        # "Geometry.OCCSewFaces": 1,
        # "Geometry.Tolerance": 1e-3,
    }

    def runFromInput(self, inputFile, runGui=False):
        caseName = Path(inputFile).stem

        gmsh.initialize()
        self.meshFromStep(inputFile, caseName, self.DEFAULT_MESHING_OPTIONS)
        self.exportGeometryAreas(caseName)
        gmsh.write(caseName + '.msh')
        gmsh.write(caseName + '.vtk') # vtk export is just for debugging. 
        if runGui:
            gmsh.fltk.run()

        gmsh.finalize()

    def meshFromStep(self, inputFile: str, caseName: str, meshingOptions=None):
        if meshingOptions is None:
            meshingOptions = Mesher.DEFAULT_MESHING_OPTIONS

        gmsh.model.add(caseName)
        allShapes = ShapesClassification(
            gmsh.model.occ.importShapes(inputFile, highestDimOnly=False)
        )

        # --- Geometry manipulation ---
        allShapes.ensureDielectricsDoNotOverlap()
        allShapes.removeConductorsFromDielectrics()
        vacuumDomain = allShapes.buildVacuumDomain()
        # -- Boundaries
        pecBoundaries = self.extractBoundaries(allShapes.pecs)

        self.buildPhysicalModel(
            pecBoundaries, 
            allShapes.dielectrics,
            allShapes.open,
            vacuumDomain
        )
        
        for [opt, val] in meshingOptions.items():
            gmsh.option.setNumber(opt, val)

        # --- Mesh generation ---
        
        gmsh.model.mesh.generate(2)

    def exportGeometryAreas(self, caseName:str):
        exporter = AreaExporterService()
        exporter.addPhysicalModelOfDimension(dimension=2)
        exporter.addPhysicalModelOfDimension(dimension=1)
        exporter.exportToJson(caseName)
            

    def buildPhysicalModel(self, pecBoundaries, dielectrics, openRegion, vacuumDomain):
        self._addPhysicalGroup("Conductor_", pecBoundaries, dimensionTag=1)
        self._addPhysicalGroup("OpenBoundary_", openRegion, dimensionTag=1)
        self._addPhysicalGroup("Vacuum_", vacuumDomain, dimensionTag=2)
        self._addPhysicalGroup("Dielectric_", dielectrics, dimensionTag=2)

        allEnts = gmsh.model.get_entities()
        entsInPG = []
        for pG in gmsh.model.get_physical_groups():
            ents = gmsh.model.getEntitiesForPhysicalGroup(pG[0], pG[1])
            for ent in ents:
                entsInPG.append((pG[0], ent))
        
        entsNotInPG = [x for x in allEnts if x not in entsInPG]
        gmsh.model.remove_entities(entsNotInPG, recursive=False)
        gmsh.model.occ.synchronize()


    def _addPhysicalGroup(self, physicalGroupName:str, objsDict:Dict, dimensionTag=1):
        for num, objs in objsDict.items():
            name = physicalGroupName + str(num)
            tags = [x[1] for x in objs]
            gmsh.model.addPhysicalGroup(dimensionTag, tags, name=name)
            

    @staticmethod
    def getPhysicalGroupWithName(name: str):
        pGs = gmsh.model.getPhysicalGroups()
        for pG in pGs:
            if gmsh.model.getPhysicalName(*pG) == name:
                return pG

    def extractBoundaries(self, shapes: dict):
        shapeBoundaries = dict()
        for num, surfs in shapes.items():
            bdrs = gmsh.model.getBoundary(surfs)
            shapeBoundaries[num] = bdrs
        return shapeBoundaries

    @classmethod
    def runCase(cls, folder: str, caseName: str, meshingOptions=None):
        gmsh.initialize()
        inputFile = folder + caseName + '/' + caseName + ".step"
        mesher=Mesher()
        mesher.meshFromStep(inputFile, caseName, meshingOptions)
        
        gmsh.write(caseName + '.msh')

        gmsh.finalize()

def print_entity_info(dim, tag):
    print(f"--- Entity (dim={dim}, tag={tag}) ---")
    
    # Name (if any)
    name = gmsh.model.get_entity_name(dim, tag)
    print(f"Name: {name if name else '(none)'}")
    
    # Type (e.g., 'Point', 'Line', 'Surface', 'Volume')
    entity_type = gmsh.model.get_type(dim, tag)
    print(f"Type: {entity_type}")
    
    # Bounding box
    bbox = gmsh.model.get_bounding_box(dim, tag)
    print(f"Bounding box: xmin={bbox[0]}, ymin={bbox[1]}, zmin={bbox[2]}, xmax={bbox[3]}, ymax={bbox[4]}, zmax={bbox[5]}")
    
    # Physical groups
    phys_groups = []
    for d in range(4):
        for pg in gmsh.model.getPhysicalGroups(d):
            if tag in gmsh.model.getEntitiesForPhysicalGroup(pg[0], pg[1]):
                phys_groups.append((pg[0], pg[1], gmsh.model.getPhysicalName(pg[0], pg[1])))
    print(f"Physical groups: {phys_groups if phys_groups else '(none)'}")
    
    # Parent entities
    parents = gmsh.model.get_adjacencies(dim, tag)[0]
    print(f"Parent entities: {parents if parents else '(none)'}")
    
    # Child entities
    children = gmsh.model.get_adjacencies(dim, tag)[1]
    print(f"Child entities: {children}")
    
    # Mesh nodes (if mesh exists)
    try:
        node_tags, node_coords, _ = gmsh.model.mesh.getNodes(dim, tag)
        print(f"Number of mesh nodes: {len(node_tags)}")
    except Exception as e:
        print("Mesh nodes: (not available)")
    
    print("-----------------------------\n")