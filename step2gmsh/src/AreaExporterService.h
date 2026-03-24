#pragma once

#include <map>
#include <string>
#include <vector>

#include <nlohmann/json.hpp>

namespace step2gmsh {

class AreaExporterService {
public:
    nlohmann::json computedAreas;

    AreaExporterService();

    void addComputedArea(const std::string& geometry,
                         const std::string& label,
                         double area);

    void addPhysicalModelForConductors(
        const std::map<std::string, std::string>& mappedElements);

    void exportToJson(const std::string& exportFileName) const;
};

} // namespace step2gmsh
