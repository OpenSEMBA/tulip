#pragma once

#include <map>
#include <string>
#include <tuple>
#include <utility>
#include <vector>

#include "ShapesClassification.h"
#include "AdapterOptions.h"
#include "../core/Materials.h"

namespace tulip {

struct StepLayers {
    std::string name;
    std::size_t id;
    Material material;
};

class Adapter {
public:
    Adapter(
        const std::string& inputFile,
        const AdapterOptions& options);

    
    std::map<std::string, std::string> meshFromStep();

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
