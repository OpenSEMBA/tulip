#include <gtest/gtest.h>

#include "Launcher.h"
#include "TestUtils.h"

#include <filesystem>
#include <fstream>
#include <sstream>

using namespace tulip;

namespace {

std::string readTextFile(const std::string& path)
{
	std::ifstream in(path);
	std::stringstream buffer;
	buffer << in.rdbuf();
	return buffer.str();
}

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

void expectParaViewCollections(
	const std::string& outputFolder,
	bool expectFloatingCollection)
{
	const auto allConductors = outputFolder + "ParaView/all_conductors.pvd";
	ASSERT_TRUE(std::filesystem::exists(allConductors));
	const auto allContent = readTextFile(allConductors);
	EXPECT_NE(std::string::npos, allContent.find("_electrostatic"));
	EXPECT_NE(std::string::npos, allContent.find("_magnetostatic"));

	const auto floating = outputFolder + "ParaView/floating_conductors.pvd";
	if (expectFloatingCollection) {
		ASSERT_TRUE(std::filesystem::exists(floating));
		const auto floatingContent = readTextFile(floating);
		EXPECT_NE(std::string::npos, floatingContent.find("_floating"));
	}
	else {
		EXPECT_FALSE(std::filesystem::exists(floating));
	}
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
	expectParaViewCollections(outputFolder, false);
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
	expectParaViewCollections(outputFolder, false);
}

TEST_F(LauncherTest, coax_and_bare_wire_from_adapted)
{
	const std::string inputFile = inputCase("coax_and_bare_wire");
	const std::string outputFolder =
		outFolder() + "LauncherTest.coax_and_bare_wire/";

	Launcher tulip(inputFile, outputFolder);
	EXPECT_NO_THROW(tulip.run());

	expectFDTDOutput(outputFolder, 2, 2);
	expectParaViewCollections(outputFolder, true);
}
