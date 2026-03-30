#pragma once

#include <array>
#include <utility>
#include <vector>

namespace tulip {

class BoundingBox {
public:
    double xMin, yMin, zMin;
    double xMax, yMax, zMax;

    BoundingBox(double xMin, double yMin, double zMin,
                double xMax, double yMax, double zMax);

    explicit BoundingBox(const std::array<double, 6>& coords);

    std::array<double, 3> getOrigin() const;
    std::array<double, 3> getCenter() const;
    double getDiagonal() const;
    std::array<double, 3> getLengths() const;

    static BoundingBox getBoundingBox(int dim, int tag);
    static BoundingBox getBoundingBoxFromGroup(
        const std::vector<std::pair<int, int>>& elements);
};

} // namespace tulip
