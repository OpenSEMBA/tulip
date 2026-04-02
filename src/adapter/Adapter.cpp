#include "Adapter.h"
#include "AdapterOptions.h"

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

std::map<std::string, std::string> buildLayerTypeMapping(const nlohmann::json& inputJson)
{
    std::map<std::string, std::string> mapping;
    if (!inputJson.contains("layers") || !inputJson["layers"].is_array()) {
        return mapping;
    }

    const auto materialTypeById = buildMaterialTypeById(inputJson);
    for (const auto& layer : inputJson["layers"]) {
        if (!layer.contains("name") || !layer.contains("materialId")) {
            continue;
        }
        const int materialId = layer["materialId"].get<int>();
        auto materialIt = materialTypeById.find(materialId);
        if (materialIt == materialTypeById.end()) {
            continue;
        }
        mapping[layer["name"].get<std::string>()] = materialIt->second;
    }
    return mapping;
}

AdapterOptions parseAdapterOptions(const nlohmann::json& inputJson,
                                   const std::filesystem::path& inputPath,
                                   const std::string& caseName)
{
    AdapterOptions options;
    if (inputJson.contains("adapterOptions") && inputJson["adapterOptions"].is_object()) {
        const auto& adapterOptions = inputJson["adapterOptions"];
        if (adapterOptions.contains("innerRegionBoxScalingFactor")) {
            options.innerRegionBoxScalingFactor =
                adapterOptions["innerRegionBoxScalingFactor"].get<double>();
        }
        if (adapterOptions.contains("farRegionDiskScalingFactor")) {
            options.farRegionDiskScalingFactor =
                adapterOptions["farRegionDiskScalingFactor"].get<double>();
        }
        if (adapterOptions.contains("stepFilename")) {
            options.stepFilename = adapterOptions["stepFilename"].get<std::string>();
        }
        if (adapterOptions.contains("gmshOptions") && adapterOptions["gmshOptions"].is_object()) {
            for (auto it = adapterOptions["gmshOptions"].begin();
                 it != adapterOptions["gmshOptions"].end(); ++it) {
                options.gmshOptions[it.key()] = it.value().get<double>();
            }
        }
    }

    if (options.stepFilename.empty()) {
        options.stepFilename = (inputPath.parent_path() / (caseName + ".step")).string();
    } else {
        options.stepFilename = (inputPath.parent_path() / options.stepFilename).string();
    }

    return options;
}

nlohmann::json buildAdaptedJson(const std::string& caseName,
                                const std::map<std::string, std::string>& layerTypeByName)
{
    gmsh::vectorpair groups;
    gmsh::model::getPhysicalGroups(groups);
    std::sort(groups.begin(), groups.end(),
              [](const auto& a, const auto& b) { return a.second < b.second; });

    nlohmann::json materials = nlohmann::json::array();
    int conductorId = 0;
    for (const auto& [dim, tag] : groups) {
        std::string name;
        gmsh::model::getPhysicalName(dim, tag, name);
        auto typeIt = layerTypeByName.find(name);
        const std::string layerType =
            typeIt == layerTypeByName.end() ? "dielectric" : typeIt->second;

        if (dim == 1 || layerType == "conductor" || layerType == "shield") {
            materials.push_back({
                {"type", "conductor"},
                {"conductorId", conductorId++},
                {"attribute", tag}
            });
        } else if (layerType == "open") {
            materials.push_back({
                {"type", "openBoundary"},
                {"attribute", tag}
            });
        } else {
            materials.push_back({
                {"type", "dielectric"},
                {"attribute", tag}
            });
        }
    }

    return {
        {"DriverOptions", {{"exportFolder", "Results/" + caseName + "/"}}},
        {"model", {
            {"materials", materials},
            {"gmshFile", caseName + ".msh"}
        }}
    };
}

} // namespace

Adapter::Adapter(const std::string& inputFile)
{
    if (!gmsh::isInitialized()) {
        throw std::runtime_error("gmsh is not initialized.");
    } 

    gmsh::clear();
    
    const std::filesystem::path inputPath(inputFile);
    if (!hasSuffix(inputPath.filename().string(), ".tulip.input.json")) {
        throw std::runtime_error("Unsupported input file extension: " + inputPath.string());
    }

    const nlohmann::json inputJson = readJsonFile(inputPath);
    caseName_ = getCaseNameFromInputPath(inputPath);
    inputDir_ = inputPath.parent_path().string();

    const AdapterOptions adapterOptions = parseAdapterOptions(inputJson, inputPath, caseName_);
    if (!std::filesystem::exists(adapterOptions.stepFilename)) {
        throw std::runtime_error("STEP file not found: " + adapterOptions.stepFilename);
    }

    gmsh::model::add(caseName_);

    EntityList shapes;
    gmsh::model::occ::importShapes(adapterOptions.stepFilename, shapes, false);

    ShapesClassification allShapes(shapes, inputFile);

    allShapes.ensureDielectricsDoNotOverlap();
    allShapes.removeConductorsFromDielectrics();
    allShapes.vacuum = allShapes.buildVacuumDomain();
    allShapes.conductors = extractBoundaries(allShapes.conductors);

    const auto layerNameMapping = buildLayerNameMapping(inputJson);
    const auto layerTypeMapping = buildLayerTypeMapping(inputJson);
    buildPhysicalModel(allShapes, layerNameMapping);

    for (const auto& [opt, val] : adapterOptions.gmshOptions) {
        gmsh::option::setNumber(opt, val);
    }

    gmsh::model::mesh::generate(2);
    gmsh::model::mesh::removeDuplicateNodes();

    auto [hasDups, _] = findDuplicateNodes();
    if (hasDups) {
        throw std::runtime_error("Duplicate mesh nodes found after meshing.");
    }

    const std::filesystem::path mshPath = std::filesystem::path(inputDir_) / (caseName_ + ".msh");
    gmsh::write(mshPath.string());

    adaptedInputJSON_ = buildAdaptedJson(caseName_, layerTypeMapping);
}

nlohmann::json Adapter::getAdaptedInputJSON() const {
    return adaptedInputJSON_;
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
