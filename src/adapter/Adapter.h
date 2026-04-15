#pragma once

#include <map>
#include <string>
#include <tuple>
#include <utility>
#include <vector>

#include <nlohmann/json.hpp>

#include "ShapesClassification.h"
#include "AdapterOptions.h"

namespace tulip {

class Adapter {
public:
    using MeshingOptions = std::map<std::string, double>;

    // Baseline gmsh options used by adapter tests and CLI code.
    static const MeshingOptions DEFAULT_MESHING_OPTIONS;

    Adapter(const std::string& inputFile);
    Adapter(const nlohmann::json& inputJson,
        const std::string& caseName,
        const std::string& inputDir = ".");
    
    void buildPhysicalModel(ShapesClassification& shapes,
                            const std::map<std::string, std::string>& labelMapping);

    nlohmann::json getAdaptedInputJSON() const;
    const AdapterOptions& getAdapterOptions() const;

    bool isOpenProblem() const;
private:
    std::string caseName_;
    std::string inputDir_;
    nlohmann::json adaptedInputJSON_;
    AdapterOptions adapterOptions_;
    ShapesClassification allShapes;

    void initialize(const nlohmann::json& inputJson,
                    const std::string& caseName,
                    const std::string& inputDir);

    static std::pair<bool,
                     std::map<std::tuple<double, double, double>,
                              std::vector<std::size_t>>>
    findDuplicateNodes();

    EntityMap extractBoundaries(const EntityMap& shapes);
    void createPhysicalGroups(const EntityMap& objsDict,
                              const std::map<std::string, std::string>& labelMapping);
};

} // namespace tulip
