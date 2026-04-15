#include <gtest/gtest.h>

#include <gmsh.h>

#include "BoundingBox.h"

using namespace tulip;

class BoundingBoxTest : public ::testing::Test {
protected:
    void TearDown() override {
        if (gmsh::isInitialized()) gmsh::finalize();
    }
};

TEST_F(BoundingBoxTest, canDefineBoundingBox) {
    BoundingBox bb(-1.0, 2.0, 3.0, 5.0, 6.0, 7.0);
    EXPECT_DOUBLE_EQ(bb.xMin, -1.0);
    EXPECT_DOUBLE_EQ(bb.yMin,  2.0);
    EXPECT_DOUBLE_EQ(bb.zMin,  3.0);
    EXPECT_DOUBLE_EQ(bb.xMax,  5.0);
    EXPECT_DOUBLE_EQ(bb.yMax,  6.0);
    EXPECT_DOUBLE_EQ(bb.zMax,  7.0);
}

TEST_F(BoundingBoxTest, constructFromArray) {
    BoundingBox bb(std::array<double, 6>{0.0, 1.0, 2.0, 3.0, 4.0, 5.0});
    EXPECT_DOUBLE_EQ(bb.xMin, 0.0);
    EXPECT_DOUBLE_EQ(bb.zMax, 5.0);
}

TEST_F(BoundingBoxTest, getOrigin) {
    BoundingBox bb(-1.0, 2.0, 3.0, 5.0, 6.0, 7.0);
    auto origin = bb.getOrigin();
    EXPECT_DOUBLE_EQ(origin[0], -1.0);
    EXPECT_DOUBLE_EQ(origin[1],  2.0);
    EXPECT_DOUBLE_EQ(origin[2],  3.0);
}

TEST_F(BoundingBoxTest, getCenter) {
    BoundingBox bb(0.0, 0.0, 0.0, 4.0, 6.0, 8.0);
    auto center = bb.getCenter();
    EXPECT_DOUBLE_EQ(center[0], 2.0);
    EXPECT_DOUBLE_EQ(center[1], 3.0);
    EXPECT_DOUBLE_EQ(center[2], 4.0);
}

TEST_F(BoundingBoxTest, getDiagonal) {
    BoundingBox bb(0.0, 0.0, 0.0, 3.0, 4.0, 0.0);
    EXPECT_DOUBLE_EQ(bb.getDiagonal(), 5.0);
}

TEST_F(BoundingBoxTest, getLengths) {
    BoundingBox bb(1.0, 2.0, 3.0, 4.0, 6.0, 9.0);
    auto len = bb.getLengths();
    EXPECT_DOUBLE_EQ(len[0], 3.0);
    EXPECT_DOUBLE_EQ(len[1], 4.0);
    EXPECT_DOUBLE_EQ(len[2], 6.0);
}

TEST_F(BoundingBoxTest, canGetBoundingBoxFromSurface) {
    gmsh::initialize();
    gmsh::model::add("Test_Model");

    int circleTag = gmsh::model::occ::addCircle(0, 0, 0, 10, 1);
    gmsh::model::occ::synchronize();

    auto bb = BoundingBox::getBoundingBox(1, circleTag);

    EXPECT_NEAR(bb.xMin, -10.0, 1e-6);
    EXPECT_NEAR(bb.yMin, -10.0, 1e-6);
    EXPECT_NEAR(bb.zMin,   0.0, 1e-6);
    EXPECT_NEAR(bb.xMax,  10.0, 1e-6);
    EXPECT_NEAR(bb.yMax,  10.0, 1e-6);
    EXPECT_NEAR(bb.zMax,   0.0, 1e-6);
}

TEST_F(BoundingBoxTest, canGetBoundingBoxFromGroupOfSurfaces) {
    gmsh::initialize();
    gmsh::model::add("Test_Model");

    int c1 = gmsh::model::occ::addCircle(0,  0,  0, 10, 1);
    int c2 = gmsh::model::occ::addCircle(0,  0,  0,  5, 2);
    int c3 = gmsh::model::occ::addCircle(25, 30, 0, 10, 3);
    gmsh::model::occ::synchronize();

    std::vector<std::pair<int,int>> group = {{1, c1}, {1, c2}, {1, c3}};
    auto bb = BoundingBox::getBoundingBoxFromGroup(group);

    EXPECT_NEAR(bb.xMin, -10.0, 1e-6);
    EXPECT_NEAR(bb.yMin, -10.0, 1e-6);
    EXPECT_NEAR(bb.zMin,   0.0, 1e-6);
    EXPECT_NEAR(bb.xMax,  35.0, 1e-6);
    EXPECT_NEAR(bb.yMax,  40.0, 1e-6);
    EXPECT_NEAR(bb.zMax,   0.0, 1e-6);
}

TEST_F(BoundingBoxTest, emptyGroupReturnsZeroBoundingBox) {
    auto bb = BoundingBox::getBoundingBoxFromGroup({});
    EXPECT_DOUBLE_EQ(bb.xMin, 0.0);
    EXPECT_DOUBLE_EQ(bb.xMax, 0.0);
}
