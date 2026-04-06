#include <gtest/gtest.h>

#include <algorithm>
#include <fstream>
#include <stdexcept>
#include <string>
#include <vector>

#include <gmsh.h>
#include <nlohmann/json.hpp>

#include "Mesher.h"
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

    int countEntitiesInPhysicalGroupWithName(const std::string &name)
    {
        auto [dim, tag] = Adapter::getPhysicalGroupWithName(name);
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
    assertAdaptedJsonMatchesExpected(caseName, adapter);
}

TEST_F(AdapterTest, realistic_case_with_dielectrics_inner_region)
{
    const std::string caseName = "realistic_case_with_dielectrics_inner_region";
    Adapter adapter(testDataPath() + "realistic_case_with_dielectrics/" + caseName + ".tulip.input.json");
    assertAdaptedJsonMatchesExpected(caseName, adapter);
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
    EXPECT_DOUBLE_EQ(options.gmshOptions.at("Mesh.MeshSizeMax"), 15.0);
    EXPECT_DOUBLE_EQ(options.gmshOptions.at("Mesh.ElementOrder"), 3.0);
}

TEST_F(AdapterTest, single_wire)
{
    const std::string caseName = "single_wire";
    Adapter adapter(inputFileFromCaseName(caseName));
    assertAdaptedJsonMatchesExpected(caseName, adapter);
}

TEST_F(AdapterTest, DISABLED_unshielded_nested)
{
    const std::string caseName = "unshielded_nested";
    Adapter adapter(inputFileFromCaseName(caseName));
    assertAdaptedJsonMatchesExpected(caseName, adapter);
}

TEST_F(AdapterTest, two_wires_with_touching_dielectric)
{
    const std::string caseName = "two_wires_with_touching_dielectric";
    Adapter adapter(inputFileFromCaseName(caseName));
    assertAdaptedJsonMatchesExpected(caseName, adapter);
}

