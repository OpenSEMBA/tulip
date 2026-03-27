#include <gtest/gtest.h>

#include <algorithm>
#include <string>
#include <vector>

#include <gmsh.h>

#include "AreaExporterService.h"
#include "Mesher.h"
#include "TestUtils.h"

using namespace step2gmsh;

class MesherTest : public ::testing::Test
{
protected:
    void SetUp() override
    {
        gmsh::initialize();
        gmsh::option::setNumber("General.Verbosity", 0);
    }
    void TearDown() override { gmsh::finalize(); }

    int countEntitiesInPhysicalGroupWithName(const std::string &name)
    {
        auto [dim, tag] = Mesher::getPhysicalGroupWithName(name);
        if (dim < 0) {
            return 0;
        }
        std::vector<int> tags;
        gmsh::model::getEntitiesForPhysicalGroup(dim, tag, tags);
        return (int)tags.size();
    }

    // expectedCounts[i] corresponds to expectedNames[i]. If empty, entity counts are not checked.
    void assertPhysicalGroups(const std::vector<std::string> &expectedNames,
                              const std::vector<int> &expectedCounts = {})
    {
        gmsh::vectorpair pGs;
        gmsh::model::getPhysicalGroups(pGs);

        std::vector<std::string> pgNames;
        for (const auto &[dim, tag] : pGs)
        {
            std::string n;
            gmsh::model::getPhysicalName(dim, tag, n);
            pgNames.push_back(n);
        }

        auto sortedActual = pgNames;
        auto sortedExpected = expectedNames;
        std::sort(sortedActual.begin(), sortedActual.end());
        std::sort(sortedExpected.begin(), sortedExpected.end());
        EXPECT_EQ(sortedActual, sortedExpected);

        for (std::size_t i = 0; i < expectedCounts.size(); ++i)
        {
            EXPECT_EQ(countEntitiesInPhysicalGroupWithName(expectedNames[i]),
                      expectedCounts[i])
                << expectedNames[i];
        }
    }
};

TEST_F(MesherTest, empty_coax)
{
    const std::string caseName = "empty_coax";
    Mesher().meshFromStep(stepFileFromCaseName(caseName), caseName);

    const std::vector<std::string> expectedNames = {"Conductor_0", "Conductor_1", "Vacuum_0"};
    assertPhysicalGroups(expectedNames, {1, 1, 1});
}

TEST_F(MesherTest, partially_filled_coax)
{
    const std::string caseName = "partially_filled_coax";
    Mesher mesher;
    mesher.meshFromStep(stepFileFromCaseName(caseName), caseName);

    gmsh::vectorpair pGs;
    gmsh::model::getPhysicalGroups(pGs);
    EXPECT_EQ(pGs.size(), 4u);

    std::vector<std::string> pgNames;
    for (const auto &[dim, tag] : pGs)
    {
        std::string name;
        gmsh::model::getPhysicalName(dim, tag, name);
        pgNames.push_back(name);
    }
    std::sort(pgNames.begin(), pgNames.end());

    std::vector<std::string> expected = {
        "Conductor_0", "Conductor_1", "Dielectric_0", "Vacuum_0"};
    std::sort(expected.begin(), expected.end());
    EXPECT_EQ(pgNames, expected);

    for (const auto &name : expected)
    {
        EXPECT_EQ(countEntitiesInPhysicalGroupWithName(name), 1) << name;
    }
}

TEST_F(MesherTest, two_wires_coax)
{
    const std::string caseName = "two_wires_coax";
    Mesher().meshFromStep(stepFileFromCaseName(caseName), caseName);

    const std::vector<std::string> expectedNames = {
        "Conductor_0", "Conductor_1", "Conductor_2", "Vacuum_0"};
    assertPhysicalGroups(expectedNames, {1, 1, 1, 1});
}

TEST_F(MesherTest, two_wires_open)
{
    const std::string caseName = "two_wires_open";
    Mesher().meshFromStep(stepFileFromCaseName(caseName), caseName);

    // gmsh::write(caseName + ".msh");
    // gmsh::write(caseName + ".vtk");

    const std::vector<std::string> expectedNames = {
        "Conductor_0", "Conductor_1", "OpenBoundary_0", "Vacuum_0", "Vacuum_1"};
    assertPhysicalGroups(expectedNames, {1, 1, 1, 1, 1});
}

TEST_F(MesherTest, dielectric_unshielded_pair)
{
    const std::string caseName = "DielectricUnshieldedPair";
    Mesher().meshFromStep(stepFileFromCaseName(caseName), caseName);

    const std::vector<std::string> expectedNames = {
        "Conductor_0", "Conductor_1",
        "Dielectric_0", "Dielectric_1",
        "OpenBoundary_0", "Vacuum_0", "Vacuum_1"};
    assertPhysicalGroups(expectedNames, {1, 1, 1, 1, 1, 1, 1});
}

TEST_F(MesherTest, dielectric_unshielded_pair_defined_boundary)
{
    const std::string caseName = "DielectricUnshieldedPairDefinedBoundary";
    auto opts = Mesher::DEFAULT_MESHING_OPTIONS;
    opts["Mesh.ElementOrder"] = 1;
    Mesher().meshFromStep(stepFileFromCaseName(caseName), caseName, &opts);

    const std::vector<std::string> expectedNames = {
        "Conductor_0", "Conductor_1",
        "Dielectric_0", "Dielectric_1",
        "OpenBoundary_0", "Vacuum_0"};

    assertPhysicalGroups(expectedNames, {1, 1, 1, 1, 1, 1});

    // gmsh::write(caseName + ".msh");
    // gmsh::write(caseName + ".vtk");
}

TEST_F(MesherTest, five_wires)
{
    const std::string caseName = "five_wires";
    Mesher().meshFromStep(stepFileFromCaseName(caseName), caseName);

    const std::vector<std::string> expectedNames = {
        "Conductor_0", "Conductor_1", "Conductor_2", "Conductor_3",
        "Conductor_4", "Conductor_5",
        "Dielectric_0", "Dielectric_1", "Dielectric_2", "Dielectric_3",
        "Dielectric_4", "Vacuum_0"};
    assertPhysicalGroups(expectedNames,
                         {1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1});
}

TEST_F(MesherTest, three_wires_ribbon)
{
    const std::string caseName = "three_wires_ribbon";
    Mesher().meshFromStep(stepFileFromCaseName(caseName), caseName);

    const std::vector<std::string> expectedNames = {
        "Conductor_0", "Conductor_1", "Conductor_2",
        "OpenBoundary_0",
        "Dielectric_0", "Dielectric_1", "Dielectric_2",
        "Vacuum_0", "Vacuum_1"};
    assertPhysicalGroups(expectedNames,
                         {1, 1, 1, 1, 1, 1, 1, 1, 1});

    // gmsh::write(caseName + ".msh");
    // gmsh::write(caseName + ".vtk");
}

TEST_F(MesherTest, nested_coax)
{
    const std::string caseName = "nested_coax";
    Mesher().meshFromStep(stepFileFromCaseName(caseName), caseName);

    const std::vector<std::string> expectedNames = {
        "Conductor_0", "Conductor_1", "Conductor_2", "Vacuum_0"};
    assertPhysicalGroups(expectedNames, {1, 2, 1, 2});
}

TEST_F(MesherTest, agrawal1981)
{
    const std::string caseName = "agrawal1981";
    Mesher().meshFromStep(stepFileFromCaseName(caseName), caseName);

    const std::vector<std::string> expectedNames = {
        "Conductor_0", "Conductor_1", "Conductor_2", "Conductor_3",
        "Dielectric_0", "Dielectric_1", "Dielectric_2",
        "Vacuum_0"};
    const std::vector<int> expectedCounts = {2, 1, 1, 1, 1, 1, 1, 2};
    assertPhysicalGroups(expectedNames, expectedCounts);
}

TEST_F(MesherTest, unshielded_multiwire)
{
    const std::string caseName = "unshielded_multiwire";
    Mesher().meshFromStep(stepFileFromCaseName(caseName), caseName);

    const std::vector<std::string> expectedNames = {
        "Conductor_0", "Conductor_1", "Dielectric_0",
        "OpenBoundary_0", "Vacuum_0", "Vacuum_1"};
    assertPhysicalGroups(expectedNames, {1, 1, 1, 1, 1, 1});
}

TEST_F(MesherTest, conductor_and_outer_dielectric)
{
    const std::string caseName = "conductor_and_outer_dielectric";
    Mesher().meshFromStep(stepFileFromCaseName(caseName), caseName);

    const std::vector<std::string> expectedNames = {
        "Conductor_0", "Dielectric_0",
        "OpenBoundary_0", "Vacuum_0", "Vacuum_1"};
    assertPhysicalGroups(expectedNames, {1, 1, 1, 2, 1});
}

TEST_F(MesherTest, realistic_case_with_dielectrics_fdtd_cell)
{
    const std::string caseName = "realistic_case_with_dielectrics_fdtd_cell";
    auto opts = Mesher::DEFAULT_MESHING_OPTIONS;
    opts["Mesh.ElementOrder"] = 1;
    Mesher().meshFromStep(stepFileFromCaseName(caseName), caseName, &opts);

    std::vector<std::string> expectedNames;
    for (int i = 0; i <= 30; ++i)
        expectedNames.push_back("Conductor_" + std::to_string(i));
    for (int i = 0; i <= 31; ++i)
        expectedNames.push_back("Dielectric_" + std::to_string(i));
    expectedNames.push_back("OpenBoundary_0");
    expectedNames.push_back("Vacuum_0");

    assertPhysicalGroups(expectedNames);
}

TEST_F(MesherTest, lansink2024_single_wire_multipolar)
{
    const std::string caseName = "lansink2024_single_wire_multipolar";
    Mesher().meshFromStep(stepFileFromCaseName(caseName), caseName);

    const std::vector<std::string> expectedNames = {
        "Conductor_0", "Dielectric_0",
        "OpenBoundary_0", "Vacuum_0", "Vacuum_1"};
    assertPhysicalGroups(expectedNames, {1, 1, 1, 1, 1});
}

TEST_F(MesherTest, unshielded_single_wire)
{
    const std::string caseName = "unshielded_single_wire";
    Mesher().meshFromStep(stepFileFromCaseName(caseName), caseName);

    const std::vector<std::string> expectedNames = {
        "Conductor_0", "OpenBoundary_0", "Vacuum_0", "Vacuum_1"};
    assertPhysicalGroups(expectedNames, {1, 1, 1, 1});
}

TEST_F(MesherTest, unshielded_nesting)
{
    const std::string caseName = "UnshieldedNested";
    Mesher().meshFromStep(stepFileFromCaseName(caseName), caseName);

    const std::vector<std::string> expectedNames = {
        "Conductor_0", "Conductor_1", "Conductor_2",
        "OpenBoundary_0", "Vacuum_0", "Vacuum_1"};
    assertPhysicalGroups(expectedNames, {2, 1, 1, 1, 2, 1});
}

TEST_F(MesherTest, unshielded_two_wires_with_touching_dielectric)
{
    const std::string caseName = "unshielded_two_wires_with_touching_dielectric";
    Mesher().meshFromStep(stepFileFromCaseName(caseName), caseName);

    const std::vector<std::string> expectedNames = {
        "Conductor_0", "Conductor_1",
        "Dielectric_0", "Dielectric_1",
        "OpenBoundary_0", "Vacuum_0", "Vacuum_1"};
    assertPhysicalGroups(expectedNames, {1, 1, 1, 1, 1, 1, 1});
}

TEST_F(MesherTest, area_exporter_returns_true_values_for_dielectric_unshielded_pair)
{
    const std::string caseName = "DielectricUnshieldedPair";
    auto mappedElements = Mesher().meshFromStep(
        stepFileFromCaseName(caseName), caseName);

    AreaExporterService exporter;
    exporter.addPhysicalModelForConductors(mappedElements);
    const auto &geometries = exporter.computedAreas["geometries"];

    ASSERT_EQ(geometries.size(), 2u);
    for (const auto &g : geometries)
    {
        EXPECT_NEAR(g["area"].get<double>(), 201.06193, 1e-4);
    }
}
