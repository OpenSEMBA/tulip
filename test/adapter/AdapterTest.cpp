#include <gtest/gtest.h>

#include <algorithm>
#include <cmath>
#include <fstream>
#include <stdexcept>
#include <string>
#include <vector>

#include <gmsh.h>
#include <nlohmann/json.hpp>

#include "Mesher.h"
#include "ShapesClassification.h"
#include "TestUtils.h"

using namespace tulip;

namespace {

std::string inputFolderFromCaseName(const std::string& caseName)
{
    return testDataPath() + caseName;
}

nlohmann::json readInputJsonFromCaseName(const std::string& caseName)
{
    std::ifstream inputStream(inputFileFromCaseName(caseName));
    if (!inputStream.is_open()) {
        throw std::runtime_error("Cannot open input JSON for case: " + caseName);
    }

    nlohmann::json inputJson;
    inputStream >> inputJson;
    return inputJson;
}

double getEntityMass(const EntityList& entities)
{
    double totalMass = 0.0;
    for (const auto& [dim, tag] : entities) {
        double mass = 0.0;
        gmsh::model::occ::getMass(dim, tag, mass);
        totalMass += mass;
    }
    return totalMass;
}

nlohmann::json findAdaptedConductorMaterialById(const nlohmann::json& adaptedJson,
                                                int conductorId)
{
    for (const auto& material : adaptedJson.at("model").at("materials")) {
        if (material.value("type", "") == "conductor" &&
            material.value("conductorId", -1) == conductorId) {
            return material;
        }
    }
    throw std::runtime_error(
        "Unable to find adapted conductor with id " + std::to_string(conductorId));
}

} // namespace

class AdapterTest : public ::testing::Test
{
protected:
    void SetUp() override
    {
        gmsh::initialize();
        gmsh::option::setNumber("General.Verbosity", 0);
    }

    void TearDown() override { gmsh::finalize(); }

    int countEntitiesInPhysicalGroupWithName(const std::string& name)
    {
        auto [dim, tag] = Adapter::getPhysicalGroupWithName(name);
        if (dim < 0) {
            return 0;
        }
        std::vector<int> tags;
        gmsh::model::getEntitiesForPhysicalGroup(dim, tag, tags);
        return static_cast<int>(tags.size());
    }

    // expectedCounts[i] corresponds to expectedNames[i]. If empty, entity counts are not checked.
    void assertPhysicalGroups(const std::vector<std::string>& expectedNames,
                              const std::vector<int>& expectedCounts = {})
    {
        gmsh::vectorpair pGs;
        gmsh::model::getPhysicalGroups(pGs);

        std::vector<std::string> pgNames;
        for (const auto& [dim, tag] : pGs)
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

    void assertAdaptedJsonMatchesExpected(const std::string& caseName,
                                          const Adapter& adapter)
    {
        auto adaptedJSON = adapter.getAdaptedInputJSON();
        std::ifstream expectedFile(
            testDataPath() + caseName + "/" + caseName + ".tulip.adapted.json");
        ASSERT_TRUE(expectedFile.is_open());
        nlohmann::json expectedJSON;
        expectedFile >> expectedJSON;
        EXPECT_EQ(adaptedJSON, expectedJSON);
    }
};

TEST_F(AdapterTest, empty_coax)
{
    const std::string caseName = "empty_coax";
    Adapter adapter(inputFileFromCaseName(caseName));
    assertAdaptedJsonMatchesExpected(caseName, adapter);
}

TEST_F(AdapterTest, empty_coax_propagates_conductor_resistance_per_meter)
{
    const std::string caseName = "empty_coax";
    auto inputJson = readInputJsonFromCaseName(caseName);
    inputJson["materials"][1]["resistancePerMeter"] = 12.5;

    Adapter adapter(inputJson, caseName, inputFolderFromCaseName(caseName));
    auto conductor = findAdaptedConductorMaterialById(adapter.getAdaptedInputJSON(), 1);

    ASSERT_TRUE(conductor.contains("resistancePerMeter"));
    EXPECT_DOUBLE_EQ(
        conductor.at("resistancePerMeter").get<double>(),
        12.5);
}

TEST_F(AdapterTest, empty_coax_conductivity_is_converted_to_resistance_per_meter)
{
    const std::string caseName = "empty_coax";
    auto inputJson = readInputJsonFromCaseName(caseName);
    const double conductivity = 5.8e7;
    inputJson["materials"][1]["conductivity"] = conductivity;

    Adapter adapter(inputJson, caseName, inputFolderFromCaseName(caseName));
    auto conductor = findAdaptedConductorMaterialById(adapter.getAdaptedInputJSON(), 1);

    ASSERT_TRUE(conductor.contains("resistancePerMeter"));
    const double innerRadiusMeters = 25e-3;
    const double expectedResistancePerMeter =
        1.0 / (conductivity * std::acos(-1.0) * innerRadiusMeters * innerRadiusMeters);
    EXPECT_NEAR(
        conductor.at("resistancePerMeter").get<double>(),
        expectedResistancePerMeter,
        expectedResistancePerMeter * 1e-12);
}

TEST_F(AdapterTest, partially_filled_coax)
{
    const std::string caseName = "partially_filled_coax";
    Adapter adapter(inputFileFromCaseName(caseName));
    assertAdaptedJsonMatchesExpected(caseName, adapter);
}

TEST_F(AdapterTest, two_wires_coax)
{
    const std::string caseName = "two_wires_coax";
    Adapter adapter(inputFileFromCaseName(caseName));
    assertAdaptedJsonMatchesExpected(caseName, adapter);
}

TEST_F(AdapterTest, five_wires)
{
    const std::string caseName = "five_wires";
    Adapter adapter(inputFileFromCaseName(caseName));
    assertAdaptedJsonMatchesExpected(caseName, adapter);
}

TEST_F(AdapterTest, two_wires_open)
{
    const std::string caseName = "two_wires_open";
    Adapter adapter(inputFileFromCaseName(caseName));
    assertAdaptedJsonMatchesExpected(caseName, adapter);
}

TEST_F(AdapterTest, two_wires_open_fails_if_layer_points_to_unknown_material)
{
    const std::string caseName = "two_wires_open";
    nlohmann::json inputJson = readInputJsonFromCaseName(caseName);
    ASSERT_TRUE(inputJson.contains("layers"));
    ASSERT_TRUE(inputJson["layers"].is_array());
    ASSERT_FALSE(inputJson["layers"].empty());

    inputJson["layers"][0]["materialId"] = 9999;

    EXPECT_THROW(
        Adapter adapter(inputJson, caseName, inputFolderFromCaseName(caseName));,
        std::runtime_error);
}

TEST_F(AdapterTest, dielectric_unshielded_pair_fails_if_input_layer_is_not_present_in_step)
{
    const std::string caseName = "dielectric_unshielded_pair";
    nlohmann::json inputJson = readInputJsonFromCaseName(caseName);
    ASSERT_TRUE(inputJson.contains("layers"));
    ASSERT_TRUE(inputJson["layers"].is_array());
    ASSERT_GE(inputJson["layers"].size(), 3);

    inputJson["layers"][2]["name"] = "MissingDielectric";

    EXPECT_THROW(
        Adapter adapter(inputJson, caseName, inputFolderFromCaseName(caseName));,
        std::runtime_error);
}

TEST_F(AdapterTest, dielectric_unshielded_pair_fails_if_step_layer_is_not_present_in_input)
{
    const std::string caseName = "dielectric_unshielded_pair";
    nlohmann::json inputJson = readInputJsonFromCaseName(caseName);
    ASSERT_TRUE(inputJson.contains("layers"));
    ASSERT_TRUE(inputJson["layers"].is_array());
    ASSERT_GE(inputJson["layers"].size(), 4);

    inputJson["layers"].erase(inputJson["layers"].begin() + 3);

    EXPECT_THROW(
        Adapter adapter(inputJson, caseName, inputFolderFromCaseName(caseName));,
        std::runtime_error);
}

TEST_F(AdapterTest, dielectric_unshielded_pair)
{
    const std::string caseName = "dielectric_unshielded_pair";
    Adapter adapter(inputFileFromCaseName(caseName));
    assertAdaptedJsonMatchesExpected(caseName, adapter);
}

TEST_F(AdapterTest, dielectric_unshielded_pair_defined_boundary)
{
    const std::string caseName = "dielectric_unshielded_pair_defined_boundary";
    Adapter adapter(inputFileFromCaseName(caseName));
    assertAdaptedJsonMatchesExpected(caseName, adapter);
}

TEST_F(AdapterTest, three_wires_ribbon)
{
    const std::string caseName = "three_wires_ribbon";
    Adapter adapter(inputFileFromCaseName(caseName));
    assertAdaptedJsonMatchesExpected(caseName, adapter);
}

TEST_F(AdapterTest, nested_coax)
{
    const std::string caseName = "nested_coax";
    Adapter adapter(inputFileFromCaseName(caseName));
    assertAdaptedJsonMatchesExpected(caseName, adapter);
}

TEST_F(AdapterTest, agrawal1981)
{
    const std::string caseName = "agrawal1981";
    Adapter adapter(inputFileFromCaseName(caseName));
    assertAdaptedJsonMatchesExpected(caseName, adapter);
}

TEST_F(AdapterTest, unshielded_multiwire)
{
    const std::string caseName = "unshielded_multiwire";
    Adapter adapter(inputFileFromCaseName(caseName));
    assertAdaptedJsonMatchesExpected(caseName, adapter);
}

TEST_F(AdapterTest, conductor_and_outer_dielectric)
{
    const std::string caseName = "conductor_and_outer_dielectric";
    Adapter adapter(inputFileFromCaseName(caseName));
    assertAdaptedJsonMatchesExpected(caseName, adapter);
}

TEST_F(AdapterTest, realistic_case_with_dielectrics_fdtd_cell)
{
    const std::string caseName = "realistic_case_with_dielectrics_fdtd_cell";
    Adapter adapter(testDataPath() + "realistic_case_with_dielectrics/" + caseName + ".tulip.input.json");

    auto adaptedJSON = adapter.getAdaptedInputJSON();
    std::ifstream expectedFile(
        testDataPath() + "realistic_case_with_dielectrics/" + caseName + ".tulip.adapted.json");
    ASSERT_TRUE(expectedFile.is_open());
    nlohmann::json expectedJSON;
    expectedFile >> expectedJSON;
    EXPECT_EQ(adaptedJSON, expectedJSON);
}

TEST_F(AdapterTest, realistic_case_with_dielectrics_inner_region)
{
    const std::string caseName = "realistic_case_with_dielectrics_inner_region";
    Adapter adapter(testDataPath() + "realistic_case_with_dielectrics/" + caseName + ".tulip.input.json");

    auto adaptedJSON = adapter.getAdaptedInputJSON();
    std::ifstream expectedFile(
        testDataPath() + "realistic_case_with_dielectrics/" + caseName + ".tulip.adapted.json");
    ASSERT_TRUE(expectedFile.is_open());
    nlohmann::json expectedJSON;
    expectedFile >> expectedJSON;
    EXPECT_EQ(adaptedJSON, expectedJSON);
}

TEST_F(AdapterTest, lansink2024_single_wire_multipolar)
{
    const std::string caseName = "lansink2024_single_wire_multipolar";
    Adapter adapter(inputFileFromCaseName(caseName));
    assertAdaptedJsonMatchesExpected(caseName, adapter);
}

TEST_F(AdapterTest, lansink2024_small_one_centered)
{
    const std::string caseName = "lansink2024_small_one_centered";
    Adapter adapter(inputFileFromCaseName(caseName));
    assertAdaptedJsonMatchesExpected(caseName, adapter);
}

TEST_F(AdapterTest, lansink2024_small_one_centered_fdtd_cell)
{
    const std::string caseName = "lansink2024_small_one_centered_fdtd_cell";
    Adapter adapter(inputFileFromCaseName(caseName));
    assertAdaptedJsonMatchesExpected(caseName, adapter);
}

TEST_F(AdapterTest, lansink2024_small_one_centered_fdtd_cell_parses_adapter_options)
{
    const std::string caseName = "lansink2024_small_one_centered_fdtd_cell";
    const nlohmann::json inputJson = readInputJsonFromCaseName(caseName);
    Adapter adapter(inputJson, caseName, inputFolderFromCaseName(caseName));

    const AdapterOptions& options = adapter.getAdapterOptions();
    ASSERT_TRUE(options.gmshOptions.find("Mesh.MeshSizeMax") != options.gmshOptions.end());
    EXPECT_DOUBLE_EQ(options.gmshOptions.at("Mesh.MeshSizeMax"), 10.0);
    EXPECT_DOUBLE_EQ(options.gmshOptions.at("Mesh.ElementOrder"), 3.0);
}

TEST_F(AdapterTest, single_wire)
{
    const std::string caseName = "single_wire";
    Adapter adapter(inputFileFromCaseName(caseName));
    assertAdaptedJsonMatchesExpected(caseName, adapter);
}

TEST_F(AdapterTest, coax_and_bare_wire)
{
    const std::string caseName = "coax_and_bare_wire";
    Adapter adapter(inputFileFromCaseName(caseName));
    assertAdaptedJsonMatchesExpected(caseName, adapter);
}

TEST_F(AdapterTest, two_wires_with_touching_dielectric)
{
    const std::string caseName = "two_wires_with_touching_dielectric";
    Adapter adapter(inputFileFromCaseName(caseName));
    assertAdaptedJsonMatchesExpected(caseName, adapter);
}

TEST_F(AdapterTest, overlapping_dielectrics_prioritize_higher_relative_permittivity)
{
    auto runCase = [&](double leftPermittivity, double rightPermittivity) {
        gmsh::clear();
        gmsh::model::add("dielectric_priority");

        const int leftTag = gmsh::model::occ::addRectangle(0.0, 0.0, 0.0, 4.0, 2.0);
        const int rightTag = gmsh::model::occ::addRectangle(2.0, 0.0, 0.0, 4.0, 2.0);
        gmsh::model::occ::synchronize();
        gmsh::model::setEntityName(2, leftTag, "LeftDielectric");
        gmsh::model::setEntityName(2, rightTag, "RightDielectric");

        const EntityList shapes = {{2, leftTag}, {2, rightTag}};
        const nlohmann::json inputJson = {
            {"materials", nlohmann::json::array({
                {{"id", 1}, {"type", "dielectric"}, {"relativePermittivity", leftPermittivity}},
                {{"id", 2}, {"type", "dielectric"}, {"relativePermittivity", rightPermittivity}}
            })},
            {"layers", nlohmann::json::array({
                {{"name", "LeftDielectric"}, {"id", 0}, {"materialId", 1}},
                {{"name", "RightDielectric"}, {"id", 1}, {"materialId", 2}}
            })}
        };

        ShapesClassification classification(shapes, inputJson);
        classification.ensureDielectricsDoNotOverlap();
        gmsh::model::occ::synchronize();

        return std::pair{
            getEntityMass(classification.dielectrics.at("LeftDielectric")),
            getEntityMass(classification.dielectrics.at("RightDielectric"))
        };
    };

    const auto [leftHighLeftMass, leftHighRightMass] = runCase(4.0, 2.0);
    EXPECT_NEAR(leftHighLeftMass, 8.0, 1e-9);
    EXPECT_NEAR(leftHighRightMass, 4.0, 1e-9);

    const auto [rightHighLeftMass, rightHighRightMass] = runCase(2.0, 4.0);
    EXPECT_NEAR(rightHighLeftMass, 4.0, 1e-9);
    EXPECT_NEAR(rightHighRightMass, 8.0, 1e-9);
}
