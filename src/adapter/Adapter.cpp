#include "Adapter.h"

#include <algorithm>
#include <cassert>
#include <filesystem>
#include <fstream>
#include <map>
#include <stdexcept>
#include <tuple>

#include <gmsh.h>
#include <nlohmann/json.hpp>

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

    for (auto it = groups.begin(); it != groups.end();) {
        if (it->second.size() <= 1) it = groups.erase(it);
        else ++it;
    }

    return {!groups.empty(), groups};
}

namespace {

bool hasSuffix(const std::string& value, const std::string& suffix)
{
    return value.size() >= suffix.size() &&
           value.compare(value.size() - suffix.size(), suffix.size(), suffix) == 0;
}

std::string getCaseNameFromInputPath(const std::filesystem::path& inputPath)
{
    const std::string fileName = inputPath.filename().string();
    constexpr const char* suffix = ".tulip.input.json";
    if (!hasSuffix(fileName, suffix)) {
        throw std::runtime_error("Input file must end with .tulip.input.json: " + fileName);
    }
    return fileName.substr(0, fileName.size() - std::string(suffix).size());
}

nlohmann::json readJsonFile(const std::filesystem::path& path)
{
    std::ifstream inputStream(path);
    if (!inputStream.is_open()) {
        throw std::runtime_error("Cannot open JSON file: " + path.string());
    }
    nlohmann::json json;
    inputStream >> json;
    return json;
}

std::map<int, std::string> buildMaterialTypeById(const nlohmann::json& inputJson)
{
    std::map<int, std::string> materialTypeById;
    if (!inputJson.contains("materials") || !inputJson["materials"].is_array()) {
        return materialTypeById;
    }
    for (const auto& material : inputJson["materials"]) {
        if (material.contains("id") && material.contains("type")) {
            materialTypeById[material["id"].get<int>()] = material["type"].get<std::string>();
        }
    }
    return materialTypeById;
}

std::map<std::string, std::string> buildLayerNameMapping(const nlohmann::json& inputJson)
{
    std::map<std::string, std::string> mapping;
    mapping["Vacuum"] = "Vacuum";
    mapping["Vacuum_0"] = "Vacuum";

    if (!inputJson.contains("layers") || !inputJson["layers"].is_array()) {
        return mapping;
    }

    const auto materialTypeById = buildMaterialTypeById(inputJson);
    for (const auto& layer : inputJson["layers"]) {
        if (!layer.contains("name") || !layer.contains("materialId")) {
            continue;
        }
        const std::string layerName = layer["name"].get<std::string>();
        const int materialId = layer["materialId"].get<int>();
        auto materialIt = materialTypeById.find(materialId);
        if (materialIt == materialTypeById.end()) {
            continue;
        }

        const std::string& materialType = materialIt->second;
        if (materialType == "conductor" || materialType == "shield" ||
            materialType == "dielectric" || materialType == "open") {
            mapping[layerName] = layerName;
            if ((materialType == "conductor" || materialType == "shield") &&
                layer.contains("id")) {
                mapping["Conductor_" + std::to_string(layer["id"].get<int>())] = layerName;
            }
        }
    }

    return mapping;
}

} // namespace

std::map<std::string, std::string> Adapter::meshFromInput(
    const std::string& inputFile)
{
    const std::filesystem::path inputPath(inputFile);
    if (!hasSuffix(inputPath.filename().string(), ".tulip.input.json")) {
        throw std::runtime_error("Unsupported input file extension: " + inputPath.string());
    }

    const nlohmann::json inputJson = readJsonFile(inputPath);
    const std::string caseName = getCaseNameFromInputPath(inputPath);

    MeshingOptions opts = DEFAULT_MESHING_OPTIONS;
    if (inputJson.contains("adapterOptions") && inputJson["adapterOptions"].contains("gmshOptions")) {
        const auto& gmshOptions = inputJson["adapterOptions"]["gmshOptions"];
        if (gmshOptions.is_object()) {
            for (auto it = gmshOptions.begin(); it != gmshOptions.end(); ++it) {
                opts[it.key()] = it.value().get<double>();
            }
        }
    }

    std::filesystem::path stepPath;
    if (inputJson.contains("adapterOptions") && inputJson["adapterOptions"].contains("stepFilename")) {
        stepPath = inputPath.parent_path() /
                   inputJson["adapterOptions"]["stepFilename"].get<std::string>();
    } else {
        stepPath = inputPath.parent_path() / (caseName + ".step");
    }

    if (!std::filesystem::exists(stepPath)) {
        throw std::runtime_error("STEP file not found: " + stepPath.string());
    }

    gmsh::model::add(caseName);

    EntityList shapes;
    gmsh::model::occ::importShapes(stepPath.string(), shapes, false);

    ShapesClassification allShapes(shapes, inputFile);

    allShapes.ensureDielectricsDoNotOverlap();
    allShapes.removeConductorsFromDielectrics();
    allShapes.vacuum = allShapes.buildVacuumDomain();
    allShapes.conductors = extractBoundaries(allShapes.conductors);

    const auto layerNameMapping = buildLayerNameMapping(inputJson);
    buildPhysicalModel(allShapes, layerNameMapping);

    for (const auto& [opt, val] : opts) {
        gmsh::option::setNumber(opt, val);
    }

    gmsh::model::mesh::generate(2);
    gmsh::model::mesh::removeDuplicateNodes();

    auto [hasDups, _] = findDuplicateNodes();
    assert(!hasDups);

    return layerNameMapping;
}

void Adapter::buildPhysicalModel(
    ShapesClassification& shapes,
    const std::map<std::string, std::string>& labelMapping)
{
    EntityMap components;
    for (const auto& [k, v] : shapes.conductors)  components[k] = v;
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
