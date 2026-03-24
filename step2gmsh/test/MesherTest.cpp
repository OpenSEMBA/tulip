#include <gtest/gtest.h>

#include <gmsh.h>

#include "AreaExporterService.h"
#include "Mesher.h"
#include "TestUtils.h"

using namespace step2gmsh;

class MesherTest : public ::testing::Test {
protected:
    void SetUp() override { gmsh::initialize(); }
    void TearDown() override { gmsh::finalize(); }

    int countEntitiesInPhysicalGroupWithName(const std::string& name) {
        auto [dim, tag] = Mesher::getPhysicalGroupWithName(name);
        if (dim < 0) return 0;
        std::vector<int> tags;
        gmsh::model::getEntitiesForPhysicalGroup(dim, tag, tags);
        return (int)tags.size();
    }
};

TEST_F(MesherTest, meshFromStepWithEmptyCoax) {
    const std::string caseName = "empty_coax";
    Mesher mesher;
    mesher.meshFromStep(stepFileFromCaseName(caseName), caseName);

    gmsh::vectorpair pGs;
    gmsh::model::getPhysicalGroups(pGs);
    EXPECT_EQ(pGs.size(), 3u); // Conductor_0, Conductor_1, Vacuum_0
}

TEST_F(MesherTest, meshFromStepWithPartiallyFilledCoax) {
    const std::string caseName = "partially_filled_coax";
    Mesher mesher;
    mesher.meshFromStep(stepFileFromCaseName(caseName), caseName);

    gmsh::vectorpair pGs;
    gmsh::model::getPhysicalGroups(pGs);
    EXPECT_EQ(pGs.size(), 4u);

    std::vector<std::string> pgNames;
    for (const auto& [dim, tag] : pGs) {
        std::string name;
        gmsh::model::getPhysicalName(dim, tag, name);
        pgNames.push_back(name);
    }
    std::sort(pgNames.begin(), pgNames.end());

    std::vector<std::string> expected = {
        "Conductor_0", "Conductor_1", "Dielectric_0", "Vacuum_0"};
    std::sort(expected.begin(), expected.end());
    EXPECT_EQ(pgNames, expected);

    for (const auto& name : expected) {
        EXPECT_EQ(countEntitiesInPhysicalGroupWithName(name), 1) << name;
    }
}

TEST_F(MesherTest, areaExporterReturnsTrueValuesForDielectricUnshieldedPair) {
    const std::string caseName = "DielectricUnshieldedPair";
    auto mappedElements = Mesher().meshFromStep(
        stepFileFromCaseName(caseName), caseName);

    AreaExporterService exporter;
    exporter.addPhysicalModelForConductors(mappedElements);
    const auto& geometries = exporter.computedAreas["geometries"];

    ASSERT_EQ(geometries.size(), 2u);
    for (const auto& g : geometries) {
        EXPECT_NEAR(g["area"].get<double>(), 201.06193, 1e-4);
    }
}
