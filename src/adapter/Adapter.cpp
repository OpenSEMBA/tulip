#include "Adapter.h"
#include "AdapterOptions.h"

#include <algorithm>
#include <cassert>
#include <filesystem>
#include <fstream>
#include <set>
#include <limits>
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

bool hasPrefix(const std::string& value, const std::string& prefix)
{
    return value.size() >= prefix.size() &&
           value.compare(0, prefix.size(), prefix) == 0;
}

bool contains(const std::string& value, const std::string& token)
{
    return value.find(token) != std::string::npos;
}

void enforceDielectricVacuumContinuity(EntityMap& dielectrics, EntityMap& vacuum)
{
    gmsh::vectorpair dielectricSurfaces;
    gmsh::vectorpair vacuumSurfaces;
    std::vector<std::string> dielectricOwnerBySurface;
    std::vector<std::string> vacuumOwnerBySurface;

    for (const auto& [name, entities] : dielectrics) {
        for (const auto& entity : entities) {
            if (entity.first != 2) {
                continue;
            }
            dielectricSurfaces.push_back(entity);
            dielectricOwnerBySurface.push_back(name);
        }
    }

    for (const auto& [name, entities] : vacuum) {
        for (const auto& entity : entities) {
            if (entity.first != 2) {
                continue;
            }
            vacuumSurfaces.push_back(entity);
            vacuumOwnerBySurface.push_back(name);
        }
    }

    if (dielectricSurfaces.empty() || vacuumSurfaces.empty()) {
        return;
    }

    gmsh::vectorpair outDimTags;
    std::vector<gmsh::vectorpair> outDimTagsMap;
    gmsh::model::occ::fragment(
        dielectricSurfaces, vacuumSurfaces, outDimTags, outDimTagsMap, -1, true, true);

    if (outDimTagsMap.size() != dielectricSurfaces.size() + vacuumSurfaces.size()) {
        throw std::runtime_error(
            "Unexpected gmsh::model::occ::fragment mapping size while fragmenting "
            "dielectrics and vacuum.");
    }

    EntityMap fragmentedDielectrics;
    EntityMap fragmentedVacuum;
    for (const auto& [name, _] : dielectrics) {
        fragmentedDielectrics[name] = {};
    }
    for (const auto& [name, _] : vacuum) {
        fragmentedVacuum[name] = {};
    }

    for (std::size_t idx = 0; idx < outDimTagsMap.size(); ++idx) {
        const bool isDielectricInput = idx < dielectricSurfaces.size();
        auto& targetMap = isDielectricInput ? fragmentedDielectrics : fragmentedVacuum;
        const auto& owner = isDielectricInput
            ? dielectricOwnerBySurface[idx]
            : vacuumOwnerBySurface[idx - dielectricSurfaces.size()];

        for (const auto& entity : outDimTagsMap[idx]) {
            if (entity.first != 2) {
                continue;
            }
            targetMap[owner].push_back(entity);
        }
    }

    const auto deduplicateEntities = [](EntityMap& entitiesByName) {
        for (auto& [_, entities] : entitiesByName) {
            std::set<EntityTag> uniqueEntities;
            EntityList deduplicated;
            deduplicated.reserve(entities.size());

            for (const auto& entity : entities) {
                if (uniqueEntities.insert(entity).second) {
                    deduplicated.push_back(entity);
                }
            }
            entities = std::move(deduplicated);
        }
    };

    deduplicateEntities(fragmentedDielectrics);
    deduplicateEntities(fragmentedVacuum);

    dielectrics = std::move(fragmentedDielectrics);
    vacuum = std::move(fragmentedVacuum);
    gmsh::model::occ::synchronize();
}

std::string getEntityBaseName(const std::string& fullName)
{
    const auto pos = fullName.rfind('/');
    return pos == std::string::npos ? fullName : fullName.substr(pos + 1);
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

std::map<int, nlohmann::json> buildMaterialById(const nlohmann::json& inputJson)
{
    std::map<int, nlohmann::json> materialById;
    if (!inputJson.contains("materials") || !inputJson["materials"].is_array()) {
        return materialById;
    }
    for (const auto& material : inputJson["materials"]) {
        if (material.contains("id")) {
            materialById[material["id"].get<int>()] = material;
        }
    }
    return materialById;
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
                const std::string alias =
                    "Conductor_" + std::to_string(layer["id"].get<int>());
                if (mapping.find(alias) == mapping.end()) {
                    mapping[alias] = layerName;
                }
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

std::map<std::string, nlohmann::json> buildLayerDielectricPropertiesMapping(
    const nlohmann::json& inputJson)
{
    std::map<std::string, nlohmann::json> mapping;
    if (!inputJson.contains("layers") || !inputJson["layers"].is_array()) {
        return mapping;
    }

    const auto materialById = buildMaterialById(inputJson);
    for (const auto& layer : inputJson["layers"]) {
        if (!layer.contains("name") || !layer.contains("materialId")) {
            continue;
        }
        const std::string layerName = layer["name"].get<std::string>();
        const int materialId = layer["materialId"].get<int>();
        auto materialIt = materialById.find(materialId);
        if (materialIt == materialById.end()) {
            continue;
        }

        const auto& material = materialIt->second;
        if (!material.contains("type") || material["type"] != "dielectric") {
            continue;
        }

        // Only propagate dielectric properties to explicit dielectric layers.
        if (!contains(layerName, "Dielectric")) {
            continue;
        }

        nlohmann::json properties = nlohmann::json::object();
        if (material.contains("relativePermittivity")) {
            properties["relativePermittivity"] = material["relativePermittivity"];
        }
        if (!properties.empty()) {
            mapping[layerName] = properties;
        }
    }

    return mapping;
}

double getEntityArea(const EntityList& entities)
{
    double area = 0.0;
    for (const auto& [dim, tag] : entities) {
        if (dim != 2) {
            continue;
        }
        double mass = 0.0;
        gmsh::model::occ::getMass(dim, tag, mass);
        area += mass;
    }
    return area;
}

std::map<std::string, nlohmann::json> buildLayerConductorPropertiesMapping(
    const nlohmann::json& inputJson,
    const EntityMap& conductorEntitiesByLayerName,
    double scalingFactor)
{
    std::map<std::string, nlohmann::json> mapping;
    if (!inputJson.contains("layers") || !inputJson["layers"].is_array()) {
        return mapping;
    }

    const auto materialById = buildMaterialById(inputJson);
    for (const auto& layer : inputJson["layers"]) {
        if (!layer.contains("name") || !layer.contains("materialId")) {
            continue;
        }

        const std::string layerName = layer["name"].get<std::string>();
        const int materialId = layer["materialId"].get<int>();
        auto materialIt = materialById.find(materialId);
        if (materialIt == materialById.end()) {
            continue;
        }

        const auto& material = materialIt->second;
        const std::string materialType = material.value("type", "");
        if (materialType != "conductor" && materialType != "shield") {
            continue;
        }

        const bool hasResistance = material.contains("resistancePerMeter");
        const bool hasConductivity = material.contains("conductivity");
        if (hasResistance && hasConductivity) {
            throw std::runtime_error(
                "Material '" + layerName +
                "' cannot define both resistancePerMeter and conductivity.");
        }

        nlohmann::json properties = nlohmann::json::object();
        if (materialType == "shield") {
            properties["isShield"] = true;
            nlohmann::json transferImpedance = nlohmann::json::object();
            if (material.contains("resistancePerMeter")) {
                transferImpedance["resistiveTerm"] = material["resistancePerMeter"];
            }
            if (material.contains("inductancePerMeter")) {
                transferImpedance["inductiveTerm"] = material["inductancePerMeter"];
            }
            if (material.contains("direction")) {
                transferImpedance["direction"] = material["direction"];
            }
            if (!transferImpedance.empty()) {
                properties["transferImpedancePerMeter"] = transferImpedance;
            }
        }

        if (hasResistance) {
            properties["resistancePerMeter"] = material["resistancePerMeter"];
        }
        else if (hasConductivity) {
            const double conductivity = material["conductivity"].get<double>();
            if (conductivity <= 0.0) {
                throw std::runtime_error(
                    "Material '" + layerName + "' conductivity must be positive.");
            }

            auto conductorIt = conductorEntitiesByLayerName.find(layerName);
            if (conductorIt == conductorEntitiesByLayerName.end() && layer.contains("id")) {
                const std::string alias =
                    "Conductor_" + std::to_string(layer["id"].get<int>());
                conductorIt = conductorEntitiesByLayerName.find(alias);
            }
            if (conductorIt == conductorEntitiesByLayerName.end()) {
                throw std::runtime_error(
                    "Unable to determine conductor area for layer '" + layerName + "'.");
            }

            const double areaInSquareMeters =
                getEntityArea(conductorIt->second) * scalingFactor * scalingFactor;
            if (areaInSquareMeters <= 0.0) {
                throw std::runtime_error(
                    "Conductor layer '" + layerName + "' must have positive area.");
            }

            properties["resistancePerMeter"] =
                1.0 / (conductivity * areaInSquareMeters);
        }

        if (!properties.empty()) {
            mapping[layerName] = properties;
        }
    }

    return mapping;
}

std::map<std::string, int> buildLayerConductorIdMapping(const nlohmann::json& inputJson)
{
    std::map<std::string, int> mapping;
    if (!inputJson.contains("layers") || !inputJson["layers"].is_array()) {
        return mapping;
    }

    const auto materialTypeById = buildMaterialTypeById(inputJson);
    for (const auto& layer : inputJson["layers"]) {
        if (!layer.contains("name") || !layer.contains("materialId") || !layer.contains("id")) {
            continue;
        }

        const int materialId = layer["materialId"].get<int>();
        auto materialIt = materialTypeById.find(materialId);
        if (materialIt == materialTypeById.end()) {
            continue;
        }

        const std::string& materialType = materialIt->second;
        if (materialType == "conductor" || materialType == "shield") {
            mapping[layer["name"].get<std::string>()] = layer["id"].get<int>();
        }
    }

    return mapping;
}

std::vector<std::string> buildDielectricLayerNamesById(const nlohmann::json& inputJson)
{
    std::vector<std::pair<int, std::string>> layersById;
    if (!inputJson.contains("layers") || !inputJson["layers"].is_array()) {
        return {};
    }

    const auto materialTypeById = buildMaterialTypeById(inputJson);
    for (const auto& layer : inputJson["layers"]) {
        if (!layer.contains("name") || !layer.contains("materialId") || !layer.contains("id")) {
            continue;
        }

        const int materialId = layer["materialId"].get<int>();
        auto materialIt = materialTypeById.find(materialId);
        if (materialIt == materialTypeById.end() || materialIt->second != "dielectric") {
            continue;
        }

        layersById.emplace_back(layer["id"].get<int>(), layer["name"].get<std::string>());
    }

    std::sort(layersById.begin(), layersById.end(),
              [](const auto& a, const auto& b) { return a.first < b.first; });

    std::vector<std::string> names;
    names.reserve(layersById.size());
    for (const auto& [_, name] : layersById) {
        names.push_back(name);
    }
    return names;
}

AdapterOptions parseAdapterOptions(const nlohmann::json& inputJson,
                                   const std::filesystem::path& inputDir,
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
        if (adapterOptions.contains("gmshOptions")) {
            const auto& gmshOptions = adapterOptions["gmshOptions"];
            if (gmshOptions.is_object()) {
                for (auto it = gmshOptions.begin(); it != gmshOptions.end(); ++it) {
                    options.gmshOptions[it.key()] = it.value().get<double>();
                }
            } else if (gmshOptions.is_array()) {
                for (const auto& entry : gmshOptions) {
                    if (!entry.is_object()) {
                        continue;
                    }
                    for (auto it = entry.begin(); it != entry.end(); ++it) {
                        options.gmshOptions[it.key()] = it.value().get<double>();
                    }
                }
            }
        }
    }

    if (options.stepFilename.empty()) {
        options.stepFilename = (inputDir / (caseName + ".step")).string();
    } else {
        const std::filesystem::path configuredStep(options.stepFilename);
        if (configuredStep.is_absolute()) {
            options.stepFilename = configuredStep.string();
        } else {
            options.stepFilename = (inputDir / configuredStep).string();
        }
    }

    return options;
}

void validateLayerMaterialIds(const nlohmann::json& inputJson)
{
    if (!inputJson.contains("materials") || !inputJson["materials"].is_array() ||
        !inputJson.contains("layers") || !inputJson["layers"].is_array()) {
        return;
    }

    std::set<int> materialIds;
    for (const auto& material : inputJson["materials"]) {
        if (material.contains("id")) {
            materialIds.insert(material["id"].get<int>());
        }
    }

    for (const auto& layer : inputJson["layers"]) {
        if (!layer.contains("name") || !layer.contains("materialId")) {
            continue;
        }
        const int materialId = layer["materialId"].get<int>();
        if (materialIds.find(materialId) == materialIds.end()) {
            throw std::runtime_error(
                "Layer '" + layer["name"].get<std::string>() +
                "' references unknown materialId " + std::to_string(materialId));
        }
    }
}

void validateRequiredInputSections(const nlohmann::json& inputJson)
{
    const bool hasMaterials =
        inputJson.contains("materials") && inputJson["materials"].is_array();
    const bool hasLayers = inputJson.contains("layers") && inputJson["layers"].is_array();

    if (hasMaterials && hasLayers) {
        return;
    }

    std::string missingSections;
    if (!hasMaterials) {
        missingSections += "'materials'";
    }
    if (!hasLayers) {
        if (!missingSections.empty()) {
            missingSections += " and ";
        }
        missingSections += "'layers'";
    }

    std::string message =
        "Invalid input JSON: missing required top-level array section(s): " +
        missingSections + ".";
    if (inputJson.contains("model")) {
        message += " Found 'model'; expected top-level 'materials' and 'layers'.";
    } else {
        message += " Expected top-level 'materials' and 'layers'.";
    }
    throw std::runtime_error(message);
}

std::vector<std::string> buildAcceptedStepNamesForLayer(
    const nlohmann::json& layer,
    const std::map<int, std::string>& materialTypeById)
{
    std::vector<std::string> names;
    if (!layer.contains("name") || !layer.contains("materialId")) {
        return names;
    }

    const std::string layerName = layer["name"].get<std::string>();
    const int materialId = layer["materialId"].get<int>();
    auto materialIt = materialTypeById.find(materialId);
    if (materialIt == materialTypeById.end()) {
        return names;
    }

    const std::string& materialType = materialIt->second;
    if (materialType != "conductor" && materialType != "shield" &&
        materialType != "dielectric" && materialType != "open") {
        return names;
    }

    names.push_back(layerName);
    if ((materialType == "conductor" || materialType == "shield") && layer.contains("id")) {
        const std::string alias = "Conductor_" + std::to_string(layer["id"].get<int>());
        if (alias != layerName) {
            names.push_back(alias);
        }
    }

    return names;
}

std::set<std::string> collectStepLayerNames(const EntityList& shapes)
{
    std::set<std::string> names;
    for (const auto& [dim, tag] : shapes) {
        if (dim != 2) {
            continue;
        }

        std::string fullName;
        gmsh::model::getEntityName(dim, tag, fullName);
        const std::string baseName = getEntityBaseName(fullName);
        if (!baseName.empty()) {
            names.insert(baseName);
        }
    }
    return names;
}

void validateLayerNamesMatchStep(const nlohmann::json& inputJson, const EntityList& shapes)
{
    if (!inputJson.contains("layers") || !inputJson["layers"].is_array()) {
        return;
    }

    const auto materialTypeById = buildMaterialTypeById(inputJson);
    const auto stepLayerNames = collectStepLayerNames(shapes);

    std::set<std::string> acceptedStepNames;
    for (const auto& layer : inputJson["layers"]) {
        const auto acceptedNames = buildAcceptedStepNamesForLayer(layer, materialTypeById);
        if (acceptedNames.empty()) {
            continue;
        }

        bool foundInStep = false;
        for (const auto& acceptedName : acceptedNames) {
            acceptedStepNames.insert(acceptedName);
            if (stepLayerNames.find(acceptedName) != stepLayerNames.end()) {
                foundInStep = true;
            }
        }

        if (!foundInStep) {
            throw std::runtime_error(
                "Layer '" + layer["name"].get<std::string>() +
                "' from input JSON not found in STEP file.");
        }
    }

    for (const auto& stepLayerName : stepLayerNames) {
        if (acceptedStepNames.find(stepLayerName) == acceptedStepNames.end()) {
            throw std::runtime_error(
                "Layer '" + stepLayerName + "' from STEP file not found in input JSON.");
        }
    }
}

nlohmann::json buildAdaptedJson(const std::string& caseName,
                                const nlohmann::json& inputJson,
                                const std::map<std::string, std::string>& layerTypeByName,
                                const std::map<std::string, int>& conductorIdByLayerName,
                                const std::map<std::string, nlohmann::json>&
                                    conductorPropertiesByLayerName,
                                const std::vector<std::string>& dielectricLayerNamesById,
                                const std::map<std::string, nlohmann::json>&
                                    dielectricPropertiesByLayerName)
{
    // OpenBoundary is always the 1D outer boundary of the computational domain;
    // OuterRegion is always the 2D outer vacuum domain. Only the 1D boundary
    // can carry an open (Robin/absorbing) boundary condition.
    const bool openBoundaryAsOpen = true;

    gmsh::vectorpair groups;
    gmsh::model::getPhysicalGroups(groups);
    std::sort(groups.begin(), groups.end(),
              [](const auto& a, const auto& b) { return a.second < b.second; });

    nlohmann::json materials = nlohmann::json::array();
    int fallbackConductorId = 0;
    for (const auto& [dim, tag] : groups) {
        std::string name;
        gmsh::model::getPhysicalName(dim, tag, name);

        auto typeIt = layerTypeByName.find(name);
        std::string layerType =
            typeIt == layerTypeByName.end() ? "dielectric" : typeIt->second;

        if (typeIt == layerTypeByName.end()) {
            if (hasPrefix(name, "OpenBoundary") || name == "OpenBoundary") {
                layerType = openBoundaryAsOpen ? "open" : "dielectric";
            } else if (name == "OuterRegion") {
                layerType = openBoundaryAsOpen ? "dielectric" : "open";
            }
        }

        if (layerType == "conductor" || layerType == "shield") {
            auto conductorIdIt = conductorIdByLayerName.find(name);
            const int conductorId =
                conductorIdIt == conductorIdByLayerName.end()
                    ? fallbackConductorId++
                    : conductorIdIt->second;
            materials.push_back({
                {"type", "conductor"},
                {"conductorId", conductorId},
                {"attribute", tag}
            });
            auto propertiesIt = conductorPropertiesByLayerName.find(name);
            if (propertiesIt != conductorPropertiesByLayerName.end()) {
                for (const auto& [key, value] : propertiesIt->second.items()) {
                    materials.back()[key] = value;
                }
            }
        } else if (layerType == "open") {
            materials.push_back({
                {"type", "open"},
                {"attribute", tag}
            });
        } else {
            nlohmann::json material = {
                {"type", "dielectric"},
                {"attribute", tag}
            };
            auto propertiesIt = dielectricPropertiesByLayerName.find(name);
            if (propertiesIt != dielectricPropertiesByLayerName.end()) {
                for (const auto& [key, value] : propertiesIt->second.items()) {
                    material[key] = value;
                }
            }
            materials.push_back(material);
        }
    }

    std::sort(materials.begin(), materials.end(), [](const auto& a, const auto& b) {
        const std::string typeA = a.value("type", "");
        const std::string typeB = b.value("type", "");

        const int priorityA = typeA == "conductor" ? 0 : (typeA == "dielectric" ? 1 : 2);
        const int priorityB = typeB == "conductor" ? 0 : (typeB == "dielectric" ? 1 : 2);
        if (priorityA != priorityB) {
            return priorityA < priorityB;
        }

        if (typeA == "conductor" && typeB == "conductor") {
            return a.value("conductorId", std::numeric_limits<int>::max()) <
                   b.value("conductorId", std::numeric_limits<int>::max());
        }

        return a.value("attribute", std::numeric_limits<int>::max()) <
               b.value("attribute", std::numeric_limits<int>::max());
    });

    nlohmann::json driverOptions = nlohmann::json::object();
    if (inputJson.contains("driverOptions") && inputJson["driverOptions"].is_object()) {
        driverOptions = inputJson["driverOptions"];
    }

    return {
        {"driverOptions", driverOptions},
        {"model", {
            {"materials", materials},
            {"gmshFile", caseName + ".msh"}
        }}
    };
}

} // namespace

Adapter::Adapter(const std::string& inputFile)
    : allShapes(EntityList{}, nlohmann::json::object())
{
    const std::filesystem::path inputPath(inputFile);
    if (!hasSuffix(inputPath.filename().string(), ".tulip.input.json")) {
        throw std::runtime_error("Unsupported input file extension: " + inputPath.string());
    }

    const nlohmann::json inputJson = readJsonFile(inputPath);
    initialize(inputJson, getCaseNameFromInputPath(inputPath), inputPath.parent_path().string());
}

Adapter::Adapter(const nlohmann::json& inputJson,
                 const std::string& caseName,
                 const std::string& inputDir)
    : allShapes(EntityList{}, inputJson)
{
    initialize(inputJson, caseName, inputDir);
}

void Adapter::initialize(const nlohmann::json& inputJson,
                         const std::string& caseName,
                         const std::string& inputDir)
{
    if (!gmsh::isInitialized()) {
        throw std::runtime_error("gmsh is not initialized.");
    }

    gmsh::clear();

    caseName_ = caseName;
    inputDir_ = inputDir;

    validateRequiredInputSections(inputJson);
    validateLayerMaterialIds(inputJson);

    adapterOptions_ = parseAdapterOptions(inputJson, std::filesystem::path(inputDir_), caseName_);
    if (!std::filesystem::exists(adapterOptions_.stepFilename)) {
        throw std::runtime_error("STEP file not found: " + adapterOptions_.stepFilename);
    }

    gmsh::model::add(caseName_);

    EntityList shapes;
    gmsh::model::occ::importShapes(adapterOptions_.stepFilename, shapes, false);
    gmsh::model::occ::synchronize();
    validateLayerNamesMatchStep(inputJson, shapes);

    allShapes = ShapesClassification(shapes, inputJson);

    allShapes.ensureDielectricsDoNotOverlap();
    allShapes.vacuum = allShapes.buildVacuumDomain();
    enforceDielectricVacuumContinuity(allShapes.dielectrics, allShapes.vacuum);
	const EntityMap conductorDomains = allShapes.conductors;

    allShapes.removeConductorsFromDielectrics();
    allShapes.conductors = extractBoundaries(allShapes.conductors);
    
    const auto layerNameMapping = buildLayerNameMapping(inputJson);
    const auto layerTypeMapping = buildLayerTypeMapping(inputJson);
    const auto layerDielectricPropertiesMapping =
        buildLayerDielectricPropertiesMapping(inputJson);
    const auto layerConductorIdMapping = buildLayerConductorIdMapping(inputJson);
    const double meshScalingFactor =
        adapterOptions_.gmshOptions.count("Mesh.ScalingFactor") == 0
            ? 1.0
            : adapterOptions_.gmshOptions.at("Mesh.ScalingFactor");
    const auto layerConductorPropertiesMapping =
        buildLayerConductorPropertiesMapping(
            inputJson, conductorDomains, meshScalingFactor);
    const auto dielectricLayerNamesById = buildDielectricLayerNamesById(inputJson);
    buildPhysicalModel(allShapes, layerNameMapping);

    for (const auto& [opt, val] : adapterOptions_.gmshOptions) {
        gmsh::option::setNumber(opt, val);
    }

    
    gmsh::model::mesh::generate(2);
    gmsh::model::mesh::removeDuplicateNodes();
    auto [hasDups, _] = findDuplicateNodes();
    if (hasDups) {
        throw std::runtime_error("Duplicate mesh nodes found after meshing.");
    }
    
    // For debugging.
    // gmsh::write("./debug_state.geo_unrolled");
    // gmsh::write("./debug_state.vtk");
    
    const std::filesystem::path mshPath = std::filesystem::path(inputDir_) / (caseName_ + ".msh");
    gmsh::write(mshPath.string());

    adaptedInputJSON_ = buildAdaptedJson(
        caseName_,
        inputJson,
        layerTypeMapping,
        layerConductorIdMapping,
        layerConductorPropertiesMapping,
        dielectricLayerNamesById,
        layerDielectricPropertiesMapping);
}

nlohmann::json Adapter::getAdaptedInputJSON() const {
    return adaptedInputJSON_;
}

bool Adapter::isOpenProblem() const {
    return allShapes.isOpenCase;
}

const AdapterOptions& Adapter::getAdapterOptions() const {
    return adapterOptions_;
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
        const std::string& mappedName = (it == labelMapping.end()) ? name : it->second;
        if (elements.empty()) continue;
        int dim = elements[0].first;
        std::vector<int> tags;
        for (const auto& [d, t] : elements) tags.push_back(t);
        gmsh::model::addPhysicalGroup(dim, tags, -1, mappedName);
    }
}

EntityMap Adapter::extractBoundaries(const EntityMap& shapes) {
    EntityMap boundaries;
    for (const auto& [name, surfs] : shapes) {
        gmsh::vectorpair bdrs;
        gmsh::model::getBoundary(surfs, bdrs, true, true, false);
        boundaries[name] = std::move(bdrs);
    }
    return boundaries;
}

} // namespace tulip
