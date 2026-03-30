#include <gtest/gtest.h>

#include <filesystem>
#include <set>
#include <string>

#include <gmsh.h>

#include "ShapesClassification.h"
#include "TestUtils.h"

using namespace tulip;

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

    EntityList expectedAllShapes = {
        {2, 1},{2, 2},{2, 3},{2, 4},{2, 5},
        {1, 1},{1, 2},{1, 3},{1, 4},{1, 5},
        {0, 1},{0, 1},{0, 2},{0, 2},{0, 3},{0, 3},{0, 4},{0, 4},{0, 5},{0, 5},
    };
    EntityMap expectedPecs = {
        {"RightConductor", {{2, 1}}},
        {"ExternalShield",  {{2, 2}}},
        {"LeftConductor",  {{2, 3}}},
    };
    EntityMap expectedDielectrics = {
        {"RightDielectric", {{2, 4}}},
        {"LeftDielectric",  {{2, 5}}},
    };

    EXPECT_EQ(sc->allShapes,   expectedAllShapes);
    EXPECT_EQ(sc->pecs,        expectedPecs);
    EXPECT_EQ(sc->dielectrics, expectedDielectrics);
    EXPECT_FALSE(sc->isOpenCase);

    delete sc;
}

TEST_F(ShapesClassificationTest, dielectricUnshieldedPairClassification) {
    auto* sc = initShapeClassification(stepFileFromCaseName("DielectricUnshieldedPair"));

    EntityList expectedAllShapes = {
        {2, 1},{2, 2},{2, 3},{2, 4},
        {1, 1},{1, 2},{1, 3},{1, 4},
        {0, 1},{0, 1},{0, 2},{0, 2},{0, 3},{0, 3},{0, 4},{0, 4},
    };
    EntityMap expectedPecs = {
        {"LeftConductor",  {{2, 2}}},
        {"RightConductor", {{2, 1}}},
    };
    EntityMap expectedDielectrics = {
        {"RightDielectric", {{2, 3}}},
        {"LeftDielectric",  {{2, 4}}},
    };

    EXPECT_EQ(sc->allShapes,   expectedAllShapes);
    EXPECT_EQ(sc->pecs,        expectedPecs);
    EXPECT_EQ(sc->dielectrics, expectedDielectrics);
    EXPECT_TRUE(sc->isOpenCase);

    delete sc;
}

TEST_F(ShapesClassificationTest, partiallyFilledCoaxHasTwoPecs) {
    auto* sc = initShapeClassification(stepFileFromCaseName("partially_filled_coax"));
    EXPECT_EQ(sc->pecs.size(), 2u);
    EXPECT_EQ(sc->dielectrics.size(), 1u);
    delete sc;
}

TEST_F(ShapesClassificationTest, fusedConductors) {
    auto* sc = initShapeClassification(stepFileFromCaseName("FusedConductor"));
    delete sc;
}

TEST_F(ShapesClassificationTest, complexNesting) {
    auto* sc = initShapeClassification(stepFileFromCaseName("ComplexNesting"));
    delete sc;
}

TEST_F(ShapesClassificationTest, fiveWiresStepShapes) {
    auto* sc = initShapeClassification(stepFileFromCaseName("five_wires"));
    EXPECT_EQ(sc->pecs.size(), 6u);
    EXPECT_EQ(sc->dielectrics.size(), 5u);
    delete sc;
}

TEST_F(ShapesClassificationTest, threeWiresRibbonStepShapes) {
    auto* sc = initShapeClassification(stepFileFromCaseName("three_wires_ribbon"));
    EXPECT_EQ(sc->open.size(), 0u);
    EXPECT_EQ(sc->pecs.size(), 3u);
    EXPECT_EQ(sc->dielectrics.size(), 3u);
    delete sc;
}

TEST_F(ShapesClassificationTest, conductorOnlyGraphDielectricUnshielded) {
    auto* sc = initShapeClassification(stepFileFromCaseName("DielectricUnshieldedPair"));

    const Graph& originalGraph = sc->nestedGraph;
    const Graph  conductorGraph = sc->getConductorOnlyGraph();

    std::set<std::string> conductorNames;
    for (auto& [name, _] : sc->pecs)
        conductorNames.insert(name);

    std::set<std::string> graphNodes(
        conductorGraph.nodes().begin(), conductorGraph.nodes().end());

    // All nodes in conductor graph should be conductors
    for (auto& n : graphNodes)
        EXPECT_TRUE(conductorNames.count(n) > 0);

    // The graph should contain all conductors that were in the original graph
    std::set<std::string> originalConductorNodes;
    for (auto& n : originalGraph.nodes())
        if (conductorNames.count(n))
            originalConductorNodes.insert(n);
    EXPECT_EQ(graphNodes, originalConductorNodes);

    // Verify no dielectric nodes remain
    std::set<std::string> dielectricNames;
    for (auto& [name, _] : sc->dielectrics)
        dielectricNames.insert(name);
    for (auto& n : graphNodes)
        EXPECT_TRUE(dielectricNames.count(n) == 0);

    delete sc;
}

TEST_F(ShapesClassificationTest, conductorOnlyGraphFiveWires) {
    auto* sc = initShapeClassification(stepFileFromCaseName("five_wires"));

    const Graph conductorGraph = sc->getConductorOnlyGraph();
    std::set<std::string> graphNodes(
        conductorGraph.nodes().begin(), conductorGraph.nodes().end());

    std::set<std::string> conductorNames;
    for (auto& [name, _] : sc->pecs)
        conductorNames.insert(name);

    std::set<std::string> dielectricNames;
    for (auto& [name, _] : sc->dielectrics)
        dielectricNames.insert(name);

    for (auto& n : graphNodes)
        EXPECT_TRUE(conductorNames.count(n) > 0);
    for (auto& n : graphNodes)
        EXPECT_TRUE(dielectricNames.count(n) == 0);

    delete sc;
}
