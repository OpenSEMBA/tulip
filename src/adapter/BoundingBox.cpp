#include "BoundingBox.h"

#include <algorithm>
#include <cmath>
#include <stdexcept>

#include <gmsh.h>

namespace tulip {

BoundingBox::BoundingBox(double xMin_, double yMin_, double zMin_,
                         double xMax_, double yMax_, double zMax_)
    : xMin(xMin_), yMin(yMin_), zMin(zMin_),
      xMax(xMax_), yMax(yMax_), zMax(zMax_) {}

BoundingBox::BoundingBox(const std::array<double, 6>& coords)
    : xMin(coords[0]), yMin(coords[1]), zMin(coords[2]),
      xMax(coords[3]), yMax(coords[4]), zMax(coords[5]) {}

std::array<double, 3> BoundingBox::getOrigin() const {
    return {xMin, yMin, zMin};
}

std::array<double, 3> BoundingBox::getCenter() const {
    return {
        (xMax + xMin) / 2.0,
        (yMax + yMin) / 2.0,
        (zMax + zMin) / 2.0
    };
}

double BoundingBox::getDiagonal() const {
    double dx = xMax - xMin;
    double dy = yMax - yMin;
    double dz = zMax - zMin;
    return std::sqrt(dx * dx + dy * dy + dz * dz);
}

std::array<double, 3> BoundingBox::getLengths() const {
    return {xMax - xMin, yMax - yMin, zMax - zMin};
}

BoundingBox BoundingBox::getBoundingBox(int dim, int tag) {
    double xmin, ymin, zmin, xmax, ymax, zmax;
    gmsh::model::occ::getBoundingBox(dim, tag, xmin, ymin, zmin, xmax, ymax, zmax);
    return BoundingBox(xmin, ymin, zmin, xmax, ymax, zmax);
}

BoundingBox BoundingBox::getBoundingBoxFromGroup(
    const std::vector<std::pair<int, int>>& elements)
{
    if (elements.empty()) {
        return BoundingBox(0.0, 0.0, 0.0, 0.0, 0.0, 0.0);
    }

    double xMinAll, yMinAll, zMinAll, xMaxAll, yMaxAll, zMaxAll;
    {
        auto bb = getBoundingBox(elements[0].first, elements[0].second);
        xMinAll = bb.xMin; yMinAll = bb.yMin; zMinAll = bb.zMin;
        xMaxAll = bb.xMax; yMaxAll = bb.yMax; zMaxAll = bb.zMax;
    }

    for (std::size_t i = 1; i < elements.size(); ++i) {
        auto bb = getBoundingBox(elements[i].first, elements[i].second);
        xMinAll = std::min(xMinAll, bb.xMin);
        yMinAll = std::min(yMinAll, bb.yMin);
        zMinAll = std::min(zMinAll, bb.zMin);
        xMaxAll = std::max(xMaxAll, bb.xMax);
        yMaxAll = std::max(yMaxAll, bb.yMax);
        zMaxAll = std::max(zMaxAll, bb.zMax);
    }

    return BoundingBox(xMinAll, yMinAll, zMinAll, xMaxAll, yMaxAll, zMaxAll);
}

} // namespace tulip
