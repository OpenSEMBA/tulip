#pragma once

#include <map>
#include <string>
#include <utility>
#include <vector>

#include <nlohmann/json.hpp>

#include "BoundingBox.h"
#include "Graph.h"

namespace tulip {

using EntityTag  = std::pair<int, int>;
using EntityList = std::vector<EntityTag>;
using EntityMap  = std::map<std::string, EntityList>;

class ShapesClassification {
public:
    bool isOpenCase;
    EntityMap conductors;
    EntityMap shields;  
    EntityMap dielectrics;
    EntityMap open;
    EntityMap vacuum;
    EntityList allShapes;
    Graph nestedGraph;

    ShapesClassification() = default;
    ShapesClassification(const EntityList& shapes,
                         const std::string& jsonFile,
                         double innerRegionBoxScalingFactor,
                         double farRegionDiskScalingFactor);
    ShapesClassification(const EntityList& shapes,
                         const nlohmann::json& jsonData,
                         double innerRegionBoxScalingFactor,
                         double farRegionDiskScalingFactor);

    static int getNumberFromName(const std::string& entityName,
                                 const std::string& label);

    EntityMap getEntitiesByMaterialType(const EntityList& entityTags,
                                        const std::string& materialType,
                                        int entityDim = 2) const;

    EntityMap getPECs(const EntityList& entityTags);
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
    double innerRegionBoxScalingFactor_;
    double farRegionDiskScalingFactor_;

    double getDielectricRelativePermittivity(const std::string& geometryName) const;
    bool dielectricHasPriorityOver(const std::string& lhs,
                                   const std::string& rhs) const;

    std::vector<std::string> getGeometryNamesByMaterialType(
        const std::string& materialType) const;

    Graph     buildNestedGraph();
    EntityMap buildClosedVacuumDomain();
    EntityMap buildOpenVacuumDomain();
};

} // namespace tulip
