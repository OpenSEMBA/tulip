#include <gtest/gtest.h>

#include "Launcher.h"
#include "TestUtils.h"

#include <filesystem>
#include <fstream>

using namespace tulip;

namespace {

void expectFDTDOutput(
	const std::string& outputFolder,
	std::size_t expectedMaterials,
	std::size_t expectedAssociations)
{
	const auto outputFile = outputFolder + "tulip.out.json";
	ASSERT_TRUE(std::filesystem::exists(outputFile));

	const auto outJSON = readJSON(outputFile);
	ASSERT_TRUE(outJSON.contains("materials"));
	ASSERT_TRUE(outJSON.contains("materialAssociations"));
	EXPECT_EQ(expectedMaterials, outJSON["materials"].size());
	EXPECT_EQ(expectedAssociations, outJSON["materialAssociations"].size());
}

} // namespace

class LauncherTest : public ::testing::Test {};

/// Test that the Tulip launcher can run the empty_coax case successfully.
TEST_F(LauncherTest, empty_coax_from_adapted)
{
	// Create Launcher instance with empty_coax adapted JSON file
	const std::string inputFile = inputCase("empty_coax");
	const std::string outputFolder = outFolder() + "LauncherTest.empty_coax/";
	
	// Run the launcher
	Launcher tulip(inputFile, outputFolder);
	
	// This should complete without throwing an exception
	EXPECT_NO_THROW(tulip.run());
	expectFDTDOutput(outputFolder, 1, 1);
}

TEST_F(LauncherTest, empty_coax_from_input_json)
{
	const std::string caseName = "empty_coax";
	const std::string inputFile =
		casesFolder() + caseName + "/" + caseName + ".tulip.input.json";
	const std::string outputFolder =
		outFolder() + "LauncherTest.empty_coax_from_input_json/";

	Launcher tulip(inputFile, outputFolder);
	EXPECT_NO_THROW(tulip.run());
	expectFDTDOutput(outputFolder, 1, 1);
}

TEST_F(LauncherTest, coax_and_bare_wire_from_adapted)
{
	const std::string inputFile = inputCase("coax_and_bare_wire");
	const std::string outputFolder =
		outFolder() + "LauncherTest.coax_and_bare_wire/";

	Launcher tulip(inputFile, outputFolder);
	EXPECT_NO_THROW(tulip.run());

	expectFDTDOutput(outputFolder, 2, 2);
}

TEST_F(LauncherTest, nested_shield_transfer_impedance_written_to_output_json)
{
	const std::string caseName = "coax_and_bare_wire";
	auto inputJson = readJSON(
		casesFolder() + caseName + "/" + caseName + ".tulip.input.json");
	for (auto& material : inputJson["materials"]) {
		if (material.value("type", "") == "shield") {
			material["resistancePerMeter"] = 1.5e-3;
			material["inductancePerMeter"] = 3.0e-9;
			material["direction"] = "outwards";
		}
	}

	const std::string tempInputFile =
		outFolder() + "LauncherTest.coax_and_bare_wire_with_shield_transfer.tulip.input.json";
	std::ofstream out(tempInputFile);
	ASSERT_TRUE(out.is_open());
	out << inputJson;
	out.close();

	const std::string outputFolder =
		outFolder() + "LauncherTest.coax_and_bare_wire_with_shield_transfer/";

	Launcher tulip(tempInputFile, outputFolder);
	EXPECT_NO_THROW(tulip.run());

	const auto outJson = readJSON(outputFolder + "tulip.out.json");
	ASSERT_EQ(2, outJson["materials"].size());
	ASSERT_EQ(2, outJson["materialAssociations"].size());

	int containedMaterialId = -1;
	for (const auto& association : outJson["materialAssociations"]) {
		if (association.contains("containedWithinElementId")) {
			containedMaterialId = association["materialId"].get<int>();
			break;
		}
	}
	ASSERT_NE(-1, containedMaterialId);

	const nlohmann::json* containedMaterial = nullptr;
	for (const auto& material : outJson["materials"]) {
		if (material["id"] == containedMaterialId) {
			containedMaterial = &material;
			break;
		}
	}
	ASSERT_NE(nullptr, containedMaterial);
	ASSERT_TRUE(containedMaterial->contains("transferImpedancePerMeter"));
	EXPECT_DOUBLE_EQ(
		1.5e-3,
		(*containedMaterial)["transferImpedancePerMeter"]["resistiveTerm"]);
	EXPECT_DOUBLE_EQ(
		3.0e-9,
		(*containedMaterial)["transferImpedancePerMeter"]["inductiveTerm"]);
	EXPECT_EQ(
		"outwards",
		(*containedMaterial)["transferImpedancePerMeter"]["direction"]);
}
