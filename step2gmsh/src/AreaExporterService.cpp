#include "AreaExporterService.h"

#include <cmath>
#include <fstream>
#include <set>

#include <gmsh.h>

namespace step2gmsh {

AreaExporterService::AreaExporterService() {
    computedAreas["geometries"] = nlohmann::json::array();
}

void AreaExporterService::addComputedArea(const std::string& geometry,
                                          const std::string& label,
                                          double area)
{
    double rounded = std::round(area * 1e6) / 1e6;
    computedAreas["geometries"].push_back({
        {"geometry", geometry},
        {"label",    label},
        {"area",     rounded}
    });
}

void AreaExporterService::addPhysicalModelForConductors(
    const std::map<std::string, std::string>& mappedElements)
{
    gmsh::vectorpair physicalGroups;
    gmsh::model::getPhysicalGroups(physicalGroups, 1);

    gmsh::vectorpair allSurfaces;
    gmsh::model::getEntities(allSurfaces, 2);

    for (const auto& [pgDim, pgTag] : physicalGroups) {
        std::string geometryName;
        gmsh::model::getPhysicalName(pgDim, pgTag, geometryName);

        if (geometryName.rfind("Conductor_", 0) != 0) continue;

        std::vector<int> entityTags;
        gmsh::model::getEntitiesForPhysicalGroup(pgDim, pgTag, entityTags);
        std::set<int> entityTagSet(entityTags.begin(), entityTags.end());

        // Find the label associated with this geometry name
        std::string label;
        for (const auto& [key, geom] : mappedElements) {
            if (geom == geometryName) {
                label = key;
                break;
            }
        }

        // Find the surface whose boundary matches the curve tags
        for (const auto& surface : allSurfaces) {
            gmsh::vectorpair boundary;
            gmsh::model::getBoundary({surface}, boundary, true, false, false);
            std::set<int> boundaryTags;
            for (const auto& [dim, tag] : boundary) {
                boundaryTags.insert(tag);
            }
            if (entityTagSet == boundaryTags) {
                double area = 0.0;
                gmsh::model::occ::getMass(2, surface.second, area);
                addComputedArea(geometryName, label, area);
                break;
            }
        }
    }
}

void AreaExporterService::exportToJson(const std::string& exportFileName) const {
    std::ofstream f(exportFileName + ".areas.json");
    f << computedAreas.dump(3);
}

} // namespace step2gmsh
