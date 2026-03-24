#pragma once

#include <map>
#include <string>
#include <utility>
#include <vector>

#include <nlohmann/json.hpp>

#include "BoundingBox.h"
#include "Graph.h"

namespace step2gmsh {

using EntityTag  = std::pair<int, int>;
using EntityList = std::vector<EntityTag>;
using EntityMap  = std::map<std::string, EntityList>;

class ShapesClassification {
public:
    bool      isOpenCase;
    EntityMap pecs;
    EntityMap dielectrics;
    EntityMap open;
    EntityMap vacuum;
    EntityList allShapes;
    Graph      nestedGraph;

    ShapesClassification(const EntityList& shapes, const std::string& jsonFile);

    static int getNumberFromName(const std::string& entityName,
                                 const std::string& label);

    EntityMap getEntitiesByMaterialType(const EntityList& entityTags,
                                        const std::string& materialType,
                                        int entityDim = 2) const;

    EntityMap getPecs(const EntityList& entityTags);
    EntityMap getDielectrics(const EntityList& entityTags);
    EntityMap getOpenBoundaries(const EntityList& entityTags);

    bool isOpenBoundaryDefined() const;
    bool isOpenProblem() const;

    void removeConductorsFromDielectrics();
    void ensureDielectricsDoNotOverlap();
    EntityMap buildVacuumDomain();

    Graph getConductorOnlyGraph() const;
    std::map<std::string, std::string> getMappedComponents() const;

private:
    nlohmann::json crossSectionData_;

    std::vector<std::string> getGeometryNamesByMaterialType(
        const std::string& materialType) const;

    Graph     buildNestedGraph();
    EntityMap buildClosedVacuumDomain();
    EntityMap buildOpenVacuumDomain();
};

} // namespace step2gmsh
