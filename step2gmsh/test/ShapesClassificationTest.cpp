#include <gtest/gtest.h>

#include <filesystem>

#include <gmsh.h>

#include "ShapesClassification.h"
#include "TestUtils.h"

using namespace step2gmsh;

class ShapesClassificationTest : public ::testing::Test {
protected:
    void SetUp() override { gmsh::initialize(); }
    void TearDown() override { gmsh::finalize(); }

    ShapesClassification* initShapeClassification(const std::string& inputFile) {
        EntityList shapes;
        gmsh::model::occ::importShapes(inputFile, shapes, false);
        auto jsonFile = std::filesystem::path(inputFile).replace_extension(".json").string();
        return new ShapesClassification(shapes, jsonFile);
    }
};

TEST_F(ShapesClassificationTest, getNumberFromName) {
    EXPECT_EQ(ShapesClassification::getNumberFromName("Shapes/Conductor_1", "Conductor_"), 1);
    EXPECT_EQ(ShapesClassification::getNumberFromName(
                  "Shapes/solid_wire_002/Conductor_002/Conductor_002", "Conductor_"), 2);
}

TEST_F(ShapesClassificationTest, dielectricShieldedPairClassification) {
    auto* sc = initShapeClassification(stepFileFromCaseName("DielectricShieldedPair"));

    EntityMap expectedPecs = {
        {"RightConductor", {{2, 1}}},
        {"ExternalShield",  {{2, 2}}},
        {"LeftConductor",  {{2, 3}}},
    };
    EntityMap expectedDielectrics = {
        {"RightDielectric", {{2, 4}}},
        {"LeftDielectric",  {{2, 5}}},
    };

    EXPECT_EQ(sc->pecs,        expectedPecs);
    EXPECT_EQ(sc->dielectrics, expectedDielectrics);
    EXPECT_FALSE(sc->isOpenCase);

    delete sc;
}

TEST_F(ShapesClassificationTest, dielectricUnshieldedPairClassification) {
    auto* sc = initShapeClassification(stepFileFromCaseName("DielectricUnshieldedPair"));

    EntityMap expectedPecs = {
        {"LeftConductor",  {{2, 2}}},
        {"RightConductor", {{2, 1}}},
    };
    EntityMap expectedDielectrics = {
        {"RightDielectric", {{2, 3}}},
        {"LeftDielectric",  {{2, 4}}},
    };

    EXPECT_EQ(sc->pecs,        expectedPecs);
    EXPECT_EQ(sc->dielectrics, expectedDielectrics);
    EXPECT_TRUE(sc->isOpenCase);

    delete sc;
}

TEST_F(ShapesClassificationTest, partiallyFilledCoaxHasTwoPecs) {
    auto* sc = initShapeClassification(stepFileFromCaseName("partially_filled_coax"));
    EXPECT_EQ(sc->pecs.size(), 2u);
    delete sc;
}
