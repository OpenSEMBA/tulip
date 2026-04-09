#include <gtest/gtest.h>

#include "Launcher.h"
#include "TestUtils.h"

using namespace tulip;

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
}
