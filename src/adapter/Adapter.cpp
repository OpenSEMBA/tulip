#include "Adapter.h"

#include <algorithm>
#include <cassert>
#include <filesystem>
#include <map>
#include <tuple>

#include <gmsh.h>

namespace tulip {

const Adapter::MeshingOptions Adapter::DEFAULT_MESHING_OPTIONS{
    {"Mesh.MshFileVersion", 2.2},
    {"Mesh.MeshSizeFromCurvature", 50.0},
    {"Mesh.ElementOrder", 3.0},
    {"Mesh.ScalingFactor", 1e-3},
    {"Mesh.SurfaceFaces", 1.0},
    {"Mesh.MeshSizeMax", 40.0},
};

std::pair<bool,
          std::map<std::tuple<double, double, double>,
                   std::vector<std::size_t>>>
Adapter::findDuplicateNodes()
{
    std::vector<std::size_t> nodeTags;
    std::vector<double>      nodeCoords, nodeParams;
    gmsh::model::mesh::getNodes(nodeTags, nodeCoords, nodeParams);

    std::map<std::tuple<double, double, double>, std::vector<std::size_t>> groups;
    for (std::size_t idx = 0; idx < nodeTags.size(); ++idx) {
        auto key = std::make_tuple(
            nodeCoords[3 * idx],
            nodeCoords[3 * idx + 1],
            nodeCoords[3 * idx + 2]);
        groups[key].push_back(nodeTags[idx]);
    }

    // Keep only groups with duplicates
    for (auto it = groups.begin(); it != groups.end();) {
        if (it->second.size() <= 1) it = groups.erase(it);
        else ++it;
    }

    return {!groups.empty(), groups};
}

std::map<std::string, std::string> Adapter::meshFromInput(
    const std::string& inputFile,
    const std::string& caseName,
    const MeshingOptions* meshingOptions)
{
    const auto& opts = meshingOptions ? *meshingOptions : DEFAULT_MESHING_OPTIONS;
    const auto modelName = caseName.empty()
        ? std::filesystem::path(inputFile).stem().string()
        : caseName;

    gmsh::model::add(modelName);

    EntityList shapes;
    gmsh::model::occ::importShapes(inputFile, shapes, false);

    auto jsonFile = std::filesystem::path(inputFile).replace_extension(".json").string();
    ShapesClassification allShapes(shapes, jsonFile);

    // Geometry manipulation
    allShapes.ensureDielectricsDoNotOverlap();
    allShapes.removeConductorsFromDielectrics();
    allShapes.vacuum = allShapes.buildVacuumDomain();
    allShapes.pecs = extractBoundaries(allShapes.pecs);

    // Mapping
    auto mappedComponents = allShapes.getMappedComponents();
    buildPhysicalModel(allShapes, mappedComponents);

    // Meshing options
    for (const auto& [opt, val] : opts) {
        gmsh::option::setNumber(opt, val);
    }

    gmsh::model::mesh::generate(2);
    gmsh::model::mesh::removeDuplicateNodes();

    auto [hasDups, _] = findDuplicateNodes();
    assert(!hasDups);

    return mappedComponents;
}

void Adapter::buildPhysicalModel(
    ShapesClassification& shapes,
    const std::map<std::string, std::string>& labelMapping)
{
    EntityMap components;
    for (const auto& [k, v] : shapes.pecs)        components[k] = v;
    for (const auto& [k, v] : shapes.dielectrics) components[k] = v;
    for (const auto& [k, v] : shapes.open)        components[k] = v;
    for (const auto& [k, v] : shapes.vacuum)      components[k] = v;

    createPhysicalGroups(components, labelMapping);

    gmsh::vectorpair allEnts;
    gmsh::model::getEntities(allEnts);

    gmsh::vectorpair allPGs;
    gmsh::model::getPhysicalGroups(allPGs);

    gmsh::vectorpair entsInPG;
    for (const auto& [pgDim, pgTag] : allPGs) {
        std::vector<int> tags;
        gmsh::model::getEntitiesForPhysicalGroup(pgDim, pgTag, tags);
        for (int t : tags) {
            entsInPG.push_back({pgDim, t});
        }
    }

    gmsh::vectorpair entsNotInPG;
    for (const auto& ent : allEnts) {
        if (std::find(entsInPG.begin(), entsInPG.end(), ent) == entsInPG.end()) {
            entsNotInPG.push_back(ent);
        }
    }

    gmsh::model::removeEntities(entsNotInPG, false);
    gmsh::model::occ::synchronize();
}

void Adapter::createPhysicalGroups(
    const EntityMap& objsDict,
    const std::map<std::string, std::string>& labelMapping)
{
    for (const auto& [name, elements] : objsDict) {
        auto it = labelMapping.find(name);
        if (it == labelMapping.end()) continue;
        const std::string& mappedName = it->second;
        int dim = elements[0].first;
        std::vector<int> tags;
        for (const auto& [d, t] : elements) tags.push_back(t);
        gmsh::model::addPhysicalGroup(dim, tags, -1, mappedName);
    }
}

std::pair<int, int> Adapter::getPhysicalGroupWithName(const std::string& name) {
    gmsh::vectorpair pGs;
    gmsh::model::getPhysicalGroups(pGs);
    for (const auto& [dim, tag] : pGs) {
        std::string pgName;
        gmsh::model::getPhysicalName(dim, tag, pgName);
        if (pgName == name) return {dim, tag};
    }
    return {-1, -1};
}

EntityMap Adapter::extractBoundaries(const EntityMap& shapes) {
    EntityMap boundaries;
    for (const auto& [name, surfs] : shapes) {
        gmsh::vectorpair bdrs;
        gmsh::model::getBoundary(surfs, bdrs, true, true, false);
        boundaries[name] = bdrs;
    }
    return boundaries;
}


} // namespace tulip
