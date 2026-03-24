#pragma once

#include <map>
#include <string>
#include <tuple>
#include <utility>
#include <vector>

#include "ShapesClassification.h"

namespace step2gmsh {

class Mesher {
public:
    static const std::map<std::string, double> DEFAULT_MESHING_OPTIONS;

    static std::pair<bool,
                     std::map<std::tuple<double, double, double>,
                              std::vector<std::size_t>>>
    findDuplicateNodes();

    void runFromInput(const std::string& inputFile, bool runGui = false);

    std::map<std::string, std::string> meshFromStep(
        const std::string& inputFile,
        const std::string& caseName,
        const std::map<std::string, double>* meshingOptions = nullptr);

    void exportGeometryAreas(
        const std::string& caseName,
        const std::map<std::string, std::string>& mappedElements);

    void buildPhysicalModel(ShapesClassification& shapes,
                            const std::map<std::string, std::string>& labelMapping);

    static std::pair<int, int> getPhysicalGroupWithName(const std::string& name);

    EntityMap extractBoundaries(const EntityMap& shapes);

    static void runCase(
        const std::string& folder,
        const std::string& caseName,
        const std::map<std::string, double>* meshingOptions = nullptr);

private:
    void createPhysicalGroups(const EntityMap& objsDict,
                              const std::map<std::string, std::string>& labelMapping);
};

} // namespace step2gmsh
