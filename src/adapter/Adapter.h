#pragma once

#include <map>
#include <string>
#include <tuple>
#include <utility>
#include <vector>

#include "ShapesClassification.h"
#include "AdapterOptions.h"

namespace tulip {

class Adapter {
public:
    using MeshingOptions = std::map<std::string, double>;

    // Baseline gmsh options used by adapter tests and CLI code.
    static const MeshingOptions DEFAULT_MESHING_OPTIONS;

    std::map<std::string, std::string> meshFromStep(
        const std::string& inputFile,
        const std::string& caseName = "",
        const MeshingOptions* meshingOptions = nullptr);

    void buildPhysicalModel(ShapesClassification& shapes,
                            const std::map<std::string, std::string>& labelMapping);

    static std::pair<int, int> getPhysicalGroupWithName(const std::string& name);

    EntityMap extractBoundaries(const EntityMap& shapes);

private:
    static std::pair<bool,
                     std::map<std::tuple<double, double, double>,
                              std::vector<std::size_t>>>
    findDuplicateNodes();

    void createPhysicalGroups(const EntityMap& objsDict,
                              const std::map<std::string, std::string>& labelMapping);
};

} // namespace tulip
