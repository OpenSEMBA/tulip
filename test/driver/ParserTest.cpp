#include <gtest/gtest.h>

#include "Parser.h"
#include "TestUtils.h"
#include "constants.h"

using namespace tulip;

class ParserTest : public ::testing::Test {};

TEST_F(ParserTest, empty_coax)
{
	const std::string CASE{ "empty_coax" };
	Parser parser{ 
		casesFolder() + CASE + "/" + CASE + ".tulip.in.json" 
	};

	auto opts{ parser.readDriverOptions() };
	EXPECT_EQ(3, opts.solverOptions.order);
	EXPECT_EQ(true, opts.exportParaViewSolution);

	auto model{ parser.readModel() };
	const auto& pecs{ model.getMaterials().pecs };
	ASSERT_EQ(2, pecs.size());
	EXPECT_EQ(pecs[0].name, "Conductor_0");
	EXPECT_EQ(pecs[0].attribute, 1);
	EXPECT_EQ(pecs[1].name, "Conductor_1");
	EXPECT_EQ(pecs[1].attribute, 2);

	const auto& diels{ model.getMaterials().dielectrics };
	ASSERT_EQ(1, diels.size());
	EXPECT_EQ(diels[0].name, "Vacuum_0");
	EXPECT_EQ(diels[0].attribute, 3);
	EXPECT_EQ(diels[0].relativePermittivity, VACUUM_RELATIVE_PERMITTIVITY);

	EXPECT_NE(0, model.getMesh()->GetNE());
	EXPECT_EQ(2, model.getMesh()->bdr_attributes.Size());
}

TEST_F(ParserTest, partially_filled_coax)
{
	const std::string CASE{ "partially_filled_coax" };
	Parser parser{
		casesFolder() + CASE + "/" + CASE + ".tulip.in.json"
	};

	auto opts{ parser.readDriverOptions() };
	EXPECT_EQ(3, opts.solverOptions.order);
	EXPECT_EQ(true, opts.exportParaViewSolution);

	auto model{ parser.readModel() };
	const auto& pecs{ model.getMaterials().pecs };
	ASSERT_EQ(2, pecs.size());
	EXPECT_EQ(pecs[0].name, "Conductor_0");
	EXPECT_EQ(pecs[0].attribute, 1);
	EXPECT_EQ(pecs[1].name, "Conductor_1");
	EXPECT_EQ(pecs[1].attribute, 2);

	const auto& diels{ model.getMaterials().dielectrics };
	ASSERT_EQ(2, diels.size());
	EXPECT_EQ(diels[0].name, "Dielectric_1");
	EXPECT_EQ(diels[0].attribute, 4);
	EXPECT_EQ(diels[0].relativePermittivity, 4.0);
	EXPECT_EQ(diels[1].name, "Vacuum_0");
	EXPECT_EQ(diels[1].attribute, 3);
	EXPECT_EQ(diels[1].relativePermittivity, VACUUM_RELATIVE_PERMITTIVITY);

	EXPECT_NE(0, model.getMesh()->GetNE());
	EXPECT_EQ(2, model.getMesh()->bdr_attributes.Size());
}

TEST_F(ParserTest, hasDielectrics_returns_false_for_vacuum_only_case)
{
	const std::string CASE{ "empty_coax" };
	Parser parser{ casesFolder() + CASE + "/" + CASE + ".tulip.in.json" };
	auto model{ parser.readModel() };
	EXPECT_FALSE(model.getMaterials().hasDielectrics());
}

TEST_F(ParserTest, hasDielectrics_returns_true_when_dielectric_present)
{
	const std::string CASE{ "partially_filled_coax" };
	Parser parser{ casesFolder() + CASE + "/" + CASE + ".tulip.in.json" };
	auto model{ parser.readModel() };
	EXPECT_TRUE(model.getMaterials().hasDielectrics());
}

