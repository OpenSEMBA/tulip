#include "ShapesClassification.h"

#include <algorithm>
#include <cassert>
#include <cmath>
#include <fstream>
#include <stdexcept>

#include <gmsh.h>

namespace tulip {

namespace {
constexpr double innerRegionBoxScalingFactor = 1.15;
constexpr double farRegionBoxScalingFactor = 4.0;
constexpr double defaultRelativePermittivity = 1.0;
constexpr double relativePermittivityTolerance = 1e-12;

std::string toLegacyMaterialType(const std::string& type) {
    if (type == "conductor" || type == "shield") return "PEC";
    if (type == "dielectric") return "Dielectric";
    if (type == "open") return "OpenBoundary";
    return type;
}

nlohmann::json buildCrossSectionFromInputFormat(const nlohmann::json& jsonData) {
    nlohmann::json crossSection = nlohmann::json::array();
    if (!jsonData.contains("materials") || !jsonData.contains("layers") ||
        !jsonData["materials"].is_array() || !jsonData["layers"].is_array()) {
        return crossSection;
    }

    std::map<int, nlohmann::json> materialById;
    for (const auto& material : jsonData["materials"]) {
        if (material.contains("id") && material.contains("type")) {
            materialById[material["id"].get<int>()] = material;
        }
    }

    for (const auto& layer : jsonData["layers"]) {
        if (!layer.contains("name") || !layer.contains("materialId")) {
            continue;
        }
        const int materialId = layer["materialId"].get<int>();
        auto materialIt = materialById.find(materialId);
        if (materialIt == materialById.end()) {
            continue;
        }
        const auto& material = materialIt->second;
        const std::string mappedType = toLegacyMaterialType(material["type"].get<std::string>());
        const std::string geometryName = layer["name"].get<std::string>();
        nlohmann::json materialJson = {{"type", mappedType}};
        if (mappedType == "Dielectric" && material.contains("relativePermittivity")) {
            materialJson["relativePermittivity"] = material["relativePermittivity"];
        }
        crossSection.push_back({
            {"name", geometryName},
            {"material", materialJson}
        });
        // Also register the legacy "Conductor_<id>" alias for PEC layers so that
        // STEP files using that naming convention are still matched correctly.
        if (mappedType == "PEC" && layer.contains("id")) {
            const std::string alias =
                "Conductor_" + std::to_string(layer["id"].get<int>());
            if (alias != geometryName) {
                crossSection.push_back({
                    {"name", alias},
                    {"material", {{"type", mappedType}}}
                });
            }
        }
    }
    return crossSection;
}
}

ShapesClassification::ShapesClassification(const EntityList& shapes,
                                           const std::string& jsonFile)
    : ShapesClassification(shapes, [&jsonFile]() {
        std::ifstream f(jsonFile);
        if (!f.is_open()) {
            throw std::runtime_error("Cannot open JSON file: " + jsonFile);
        }
        nlohmann::json jsonData;
        f >> jsonData;
        return jsonData;
    }())
{}

ShapesClassification::ShapesClassification(const EntityList& shapes,
                                           const nlohmann::json& jsonData)
{
    gmsh::model::occ::synchronize();

    allShapes = shapes;


    crossSectionData_ = buildCrossSectionFromInputFormat(jsonData);

    conductors  = getPECs(shapes);
    dielectrics = getDielectrics(shapes);
    open        = getOpenBoundaries(shapes);
    nestedGraph = buildNestedGraph();
    isOpenCase  = isOpenProblem();
}

int ShapesClassification::getNumberFromName(const std::string& entityName,
                                            const std::string& label)
{
    auto pos = entityName.rfind(label);
    if (pos == std::string::npos) {
        throw std::runtime_error("Label '" + label + "' not found in '" + entityName + "'");
    }
    return std::stoi(entityName.substr(pos + label.size()));
}

std::vector<std::string> ShapesClassification::getGeometryNamesByMaterialType(
    const std::string& materialType) const
{
    std::vector<std::string> names;
    for (const auto& geometry : crossSectionData_) {
        if (geometry["material"]["type"] == materialType) {
            names.push_back(geometry["name"].get<std::string>());
        }
    }
    return names;
}

double ShapesClassification::getDielectricRelativePermittivity(
    const std::string& geometryName) const
{
    for (const auto& geometry : crossSectionData_) {
        if (geometry["name"].get<std::string>() != geometryName) {
            continue;
        }
        const auto& material = geometry["material"];
        if (!material.contains("type") || material["type"] != "Dielectric") {
            continue;
        }
        if (material.contains("relativePermittivity")) {
            return material["relativePermittivity"].get<double>();
        }
        return defaultRelativePermittivity;
    }
    return defaultRelativePermittivity;
}

bool ShapesClassification::dielectricHasPriorityOver(const std::string& lhs,
                                                     const std::string& rhs) const
{
    const double lhsPermittivity = getDielectricRelativePermittivity(lhs);
    const double rhsPermittivity = getDielectricRelativePermittivity(rhs);
    if (lhsPermittivity > rhsPermittivity + relativePermittivityTolerance) {
        return true;
    }
    if (std::abs(lhsPermittivity - rhsPermittivity) <= relativePermittivityTolerance) {
        return lhs > rhs;
    }
    return false;
}

EntityMap ShapesClassification::getEntitiesByMaterialType(
    const EntityList& entityTags,
    const std::string& materialType,
    int entityDim) const
{
    auto materialNames = getGeometryNamesByMaterialType(materialType);
    // Collect all matching entities with their full GMSH name
    struct EntityInfo {
        std::string fullName;
        EntityTag   tag;
    };
    std::map<std::string, std::vector<EntityInfo>> candidates;
    for (const auto& [dim, tag] : entityTags) {
        if (dim != entityDim) continue;
        std::string fullName;
        gmsh::model::getEntityName(dim, tag, fullName);
        auto pos  = fullName.rfind('/');
        std::string name = (pos != std::string::npos) ? fullName.substr(pos + 1) : fullName;
        if (std::find(materialNames.begin(), materialNames.end(), name) == materialNames.end()) {
            continue;
        }
        candidates[name].push_back({fullName, {dim, tag}});
    }

    EntityMap entities;
    for (const auto& [name, infos] : candidates) {
        // Check if any entity has a deep path (assembly instance).
        // Bare part definitions have paths like "Shapes/<name>" (1 slash).
        std::string barePath = "Shapes/" + name;
        bool hasDeepEntities = false;
        for (const auto& info : infos) {
            if (info.fullName != barePath) {
                hasDeepEntities = true;
                break;
            }
        }
        for (const auto& info : infos) {
            if (hasDeepEntities && info.fullName == barePath) {
                continue;  // Skip bare part definition
            }
            entities[name].push_back(info.tag);
        }
    }
    return entities;
}

EntityMap ShapesClassification::getPECs(const EntityList& entityTags) {
    return getEntitiesByMaterialType(entityTags, "PEC");
}

EntityMap ShapesClassification::getDielectrics(const EntityList& entityTags) {
    return getEntitiesByMaterialType(entityTags, "Dielectric");
}

EntityMap ShapesClassification::getOpenBoundaries(const EntityList& entityTags) {
    return getEntitiesByMaterialType(entityTags, "OpenBoundary");
}

bool ShapesClassification::isOpenBoundaryDefined() const {
    return !open.empty();
}

bool ShapesClassification::isOpenProblem() const {
    auto roots = nestedGraph.roots();
    if (open.size() == 1) return true;
    if (roots.size() > 1) return true;
    if (!roots.empty()) {
        const auto& root = roots[0];
        if (dielectrics.count(root)) return true;
        if (conductors.count(root)) {
            auto parentNodes = nestedGraph.getParentNodes();
            if (std::find(parentNodes.begin(), parentNodes.end(), root) == parentNodes.end()) {
                return true;
            }
        }
    }
    return false;
}

void ShapesClassification::removeConductorsFromDielectrics() {
    std::string closedCaseRootConductor;
    if (!isOpenCase) {
        const auto roots = nestedGraph.roots();
        if (!roots.empty() && conductors.count(roots.front())) {
            closedCaseRootConductor = roots.front();
        }
    }

    auto removeConductors = [this, &closedCaseRootConductor](auto& container) {
        for (auto& [name, surfs] : container) {
            EntityList pecSurfs;
            for (const auto& [condName, condSurfs] : conductors) {
                if (!isOpenCase && condName == closedCaseRootConductor) {
                    continue;
                }
                pecSurfs.insert(pecSurfs.end(), condSurfs.begin(), condSurfs.end());
            }
        
            if (!pecSurfs.empty()) {
                gmsh::vectorpair outDimTags;
                std::vector<gmsh::vectorpair> outMap;
                gmsh::model::occ::cut(surfs, pecSurfs, outDimTags, outMap, -1, true, false);
                surfs = outDimTags;
            }
        }
    };

    removeConductors(dielectrics);
    removeConductors(vacuum);

    gmsh::model::occ::synchronize();
}

void ShapesClassification::ensureDielectricsDoNotOverlap() {
    for (auto& [currentKey, currentSurfs] : dielectrics) {
        EntityList others;
        for (const auto& [key, surf] : dielectrics) {
            if (key != currentKey && dielectricHasPriorityOver(key, currentKey)) {
                others.insert(others.end(), surf.begin(), surf.end());
            }
        }
        if (others.empty()) continue;

        gmsh::vectorpair outDimTags;
        std::vector<gmsh::vectorpair> outMap;
        gmsh::model::occ::cut(currentSurfs, others, outDimTags, outMap, -1, true, false);
        currentSurfs = outDimTags;
    }

    gmsh::model::occ::synchronize();
}

EntityMap ShapesClassification::buildVacuumDomain() {
    if (isOpenCase) {
        vacuum = buildOpenVacuumDomain();
    } else {
        vacuum = buildClosedVacuumDomain();
    }
    return vacuum;
}

EntityMap ShapesClassification::buildClosedVacuumDomain() {
    const auto  roots   = nestedGraph.roots();
    const auto& root    = roots[0];
    EntityList  dom     = conductors.at(root);
    EntityList  toRemove;

    for (const auto& [name, surf] : conductors) {
        if (name == root) continue;
        toRemove.insert(toRemove.end(), surf.begin(), surf.end());
    }
    for (const auto& [name, surf] : dielectrics) {
        toRemove.insert(toRemove.end(), surf.begin(), surf.end());
    }

    gmsh::vectorpair outDimTags;
    std::vector<gmsh::vectorpair> outMap;
    gmsh::model::occ::cut(dom, toRemove, outDimTags, outMap, -1, false, false);
    gmsh::model::occ::synchronize();

    return {{"Vacuum", outDimTags}};
}

EntityMap ShapesClassification::buildOpenVacuumDomain() {
    EntityList nonVacuumSurfaces;
    for (const auto& [name, surf] : conductors) {
        nonVacuumSurfaces.insert(nonVacuumSurfaces.end(), surf.begin(), surf.end());
    }
    for (const auto& [name, surf] : dielectrics) {
        nonVacuumSurfaces.insert(nonVacuumSurfaces.end(), surf.begin(), surf.end());
    }

    if (isOpenBoundaryDefined()) {
        auto [openName, openVacuum] = *open.begin();

        gmsh::vectorpair outDimTags;
        std::vector<gmsh::vectorpair> outMap;
        gmsh::model::occ::cut(openVacuum, nonVacuumSurfaces, outDimTags, outMap, -1, true, false);
        gmsh::model::occ::synchronize();

        gmsh::vectorpair vacuumBoundaries;
        gmsh::model::getBoundary(outDimTags, vacuumBoundaries, true, true, false);
        EntityList externalBoundaries;
        for (const auto& [dim, tag] : vacuumBoundaries) {
            if (tag > 0) externalBoundaries.push_back({dim, tag});
        }
        open = {{openName, externalBoundaries}};

        return {{"Vacuum", outDimTags}};
    } else {
        auto boundingBox = BoundingBox::getBoundingBoxFromGroup(nonVacuumSurfaces);
        auto lengths     = boundingBox.getLengths();
        double bbMaxLen  = *std::max_element(lengths.begin(), lengths.end());

        double nearBoxSize = bbMaxLen * innerRegionBoxScalingFactor;
        auto   center      = boundingBox.getCenter();
        double nVX = center[0] - nearBoxSize / 2.0;
        double nVY = center[1] - nearBoxSize / 2.0;
        double nVZ = center[2];

        int nearRectTag = gmsh::model::occ::addRectangle(
            nVX, nVY, nVZ, nearBoxSize, nearBoxSize);
        EntityList nearVacuum = {{2, nearRectTag}};

        double farDiameter = farRegionBoxScalingFactor * boundingBox.getDiagonal();
        int    farDiskTag  = gmsh::model::occ::addDisk(
            center[0], center[1], center[2], farDiameter, farDiameter);
        EntityList farVacuum = {{2, farDiskTag}};

        gmsh::model::occ::synchronize();

        gmsh::vectorpair farBoundary;
        gmsh::model::getBoundary(farVacuum, farBoundary, true, true, false);
        open = {{"OpenBoundary", farBoundary}};

        {
            gmsh::vectorpair out;
            std::vector<gmsh::vectorpair> map;
            gmsh::model::occ::cut(farVacuum, nearVacuum, out, map, -1, true, false);
            farVacuum = out;
        }
        {
            gmsh::vectorpair out;
            std::vector<gmsh::vectorpair> map;
            gmsh::model::occ::cut(nearVacuum, nonVacuumSurfaces, out, map, -1, true, false);
            nearVacuum = out;
        }

        gmsh::model::occ::synchronize();

        // Set mesh size for inner region.
        BoundingBox bb = BoundingBox::getBoundingBox(2, nearVacuum[0].second);
        auto nlen = bb.getLengths();
        double minSide = std::min(nlen[0], nlen[1]);

        gmsh::vectorpair innerRegion;
        gmsh::model::getBoundary(nearVacuum, innerRegion, true, true, true);
        gmsh::model::mesh::setSize(innerRegion, minSide / 20.0);

        gmsh::model::occ::synchronize();

        return {{"InnerRegion", nearVacuum}, {"OuterRegion", farVacuum}};
    }
}

Graph ShapesClassification::buildNestedGraph() {
    gmsh::model::occ::synchronize();
    Graph graph;

    // Merge conductors, dielectrics, open into a single map
    EntityMap elements;
    for (const auto& [k, v] : conductors)    elements[k] = v;
    for (const auto& [k, v] : dielectrics) elements[k] = v;
    for (const auto& [k, v] : open)        elements[k] = v;

    for (const auto& [key, _] : elements) {
        graph.addNode(key);
    }

    std::vector<std::string> keys;
    for (const auto& [key, _] : elements) keys.push_back(key);

    for (std::size_t i = 0; i < keys.size(); ++i) {
        for (std::size_t j = i + 1; j < keys.size(); ++j) {
            const auto& keyA   = keys[i];
            const auto& keyB   = keys[j];
            const auto& elemsA = elements.at(keyA);
            const auto& elemsB = elements.at(keyB);

            BoundingBox bbA = BoundingBox::getBoundingBoxFromGroup(elemsA);
            BoundingBox bbB = BoundingBox::getBoundingBoxFromGroup(elemsB);

            bool aInB = (bbA.xMin >= bbB.xMin && bbA.xMax <= bbB.xMax &&
                         bbA.yMin >= bbB.yMin && bbA.yMax <= bbB.yMax);
            bool bInA = (bbB.xMin >= bbA.xMin && bbB.xMax <= bbA.xMax &&
                         bbB.yMin >= bbA.yMin && bbB.yMax <= bbA.yMax);

            if (aInB && !bInA) {
                graph.addEdge(keyB, keyA);  // B contains A
            } else if (bInA && !aInB) {
                graph.addEdge(keyA, keyB);  // A contains B
            }
        }
    }

    graph.pruneToLongestPaths();
    graph.reorderData();
    return graph;
}

Graph ShapesClassification::getConductorOnlyGraph() const {
    Graph conductorGraph;

    for (const auto& [name, _] : conductors) {
        auto& nodes = nestedGraph.nodes();
        if (std::find(nodes.begin(), nodes.end(), name) != nodes.end()) {
            conductorGraph.addNode(name);
        }
    }

    for (const auto& edge : nestedGraph.edges()) {
        const auto& src  = edge.first;
        const auto& dest = edge.second;

        bool srcIsPec  = conductors.count(src)  > 0;
        bool destIsPec = conductors.count(dest) > 0;
        bool destIsDiel = dielectrics.count(dest) > 0;

        if (srcIsPec && destIsPec) {
            conductorGraph.addEdge(src, dest);
        } else if (srcIsPec && destIsDiel) {
            for (const auto& childEdge : nestedGraph.edges()) {
                if (childEdge.first == dest && conductors.count(childEdge.second)) {
                    conductorGraph.addEdge(src, childEdge.second);
                }
            }
        }
    }

    conductorGraph.pruneToLongestPaths();
    conductorGraph.reorderData();
    return conductorGraph;
}

std::map<std::string, std::string> ShapesClassification::getMappedComponents() const {
    std::map<std::string, std::string> mappedComponents;

    std::vector<std::string> sortedNodes = nestedGraph.getNodesByLevels();

    int conductorIdx = 0;
    for (const auto& node : sortedNodes) {
        if (conductors.count(node)) {
            mappedComponents[node] = "Conductor_" + std::to_string(conductorIdx++);
        }
    }

    int dielectricIdx = 0;
    for (const auto& [name, _] : dielectrics) {
        mappedComponents[name] = "Dielectric_" + std::to_string(dielectricIdx++);
    }

    for (const auto& [domain, _] : vacuum) {
        mappedComponents[domain] = domain;
    }
    for (const auto& [openBoundary, _] : open) {
        mappedComponents[openBoundary] = "OpenBoundary_0";
    }

    return mappedComponents;
}


} // namespace tulip
