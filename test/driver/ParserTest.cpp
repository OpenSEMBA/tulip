#include <gtest/gtest.h>

#include "AdaptedInputParser.h"
#include "TestUtils.h"
#include "constants.h"

using namespace tulip;

class ParserTest : public ::testing::Test {};

TEST_F(ParserTest, empty_coax)
{
	const std::string CASE{ "empty_coax" };
	Parser parser{ 
		casesFolder() + CASE + "/" + CASE + ".tulip.adapted.json" 
	};

	auto opts{ parser.readDriverOptions() };
	EXPECT_EQ(3, opts.solverOptions.order);
	EXPECT_EQ(true, opts.exportParaViewSolution);

	auto model{ parser.readModel() };
	auto pecs = model.getMaterials().getConductors();
	ASSERT_EQ(2, pecs.size());

	auto it = pecs.begin();
	ASSERT_NE(it, pecs.end());
	EXPECT_EQ((*it)->getConductorId(), 0);
	EXPECT_EQ((*it)->getAttribute(), 1);

	++it;
	ASSERT_NE(it, pecs.end());
	EXPECT_EQ((*it)->getConductorId(), 1);
	EXPECT_EQ((*it)->getAttribute(), 2);

	auto diels = model.getMaterials().getDielectrics();
	ASSERT_EQ(1, diels.size());
	auto dIt = diels.begin();
	ASSERT_NE(dIt, diels.end());
	EXPECT_EQ((*dIt)->getAttribute(), 3);
	EXPECT_EQ((*dIt)->getRelativePermittivity(), VACUUM_RELATIVE_PERMITTIVITY);

	EXPECT_NE(0, model.getMesh()->GetNE());
	EXPECT_EQ(2, model.getMesh()->bdr_attributes.Size());
}

TEST_F(ParserTest, partially_filled_coax)
{
	const std::string CASE{ "partially_filled_coax" };
	Parser parser{
		casesFolder() + CASE + "/" + CASE + ".tulip.adapted.json"
	};

	auto opts{ parser.readDriverOptions() };
	EXPECT_EQ(3, opts.solverOptions.order);
	EXPECT_EQ(true, opts.exportParaViewSolution);

	auto model{ parser.readModel() };
	auto pecs = model.getMaterials().getConductors();
	ASSERT_EQ(2, pecs.size());

	auto it = pecs.begin();
	ASSERT_NE(it, pecs.end());
	EXPECT_EQ((*it)->getConductorId(), 0);
	EXPECT_EQ((*it)->getAttribute(), 1);

	++it;
	ASSERT_NE(it, pecs.end());
	EXPECT_EQ((*it)->getConductorId(), 1);
	EXPECT_EQ((*it)->getAttribute(), 2);

	auto diels = model.getMaterials().getDielectrics();
	ASSERT_EQ(2, diels.size());

	bool foundDielectric = false;
	bool foundVacuum = false;
	for (const auto* d : diels) {
		if (d->getAttribute() == 3 && d->getRelativePermittivity() == 4.0) {
			foundDielectric = true;
		}
		if (d->getAttribute() == 4 &&
			d->getRelativePermittivity() == VACUUM_RELATIVE_PERMITTIVITY) {
			foundVacuum = true;
		}
	}
	EXPECT_TRUE(foundDielectric);
	EXPECT_TRUE(foundVacuum);

	EXPECT_NE(0, model.getMesh()->GetNE());
	EXPECT_EQ(2, model.getMesh()->bdr_attributes.Size());
}

TEST_F(ParserTest, hasDielectrics_returns_false_for_vacuum_only_case)
{
	const std::string CASE{ "empty_coax" };
	Parser parser{ casesFolder() + CASE + "/" + CASE + ".tulip.adapted.json" };
	auto model{ parser.readModel() };
	EXPECT_FALSE(model.getMaterials().hasDielectrics());
}

TEST_F(ParserTest, hasDielectrics_returns_true_when_dielectric_present)
{
	const std::string CASE{ "partially_filled_coax" };
	Parser parser{ casesFolder() + CASE + "/" + CASE + ".tulip.adapted.json" };
	auto model{ parser.readModel() };
	EXPECT_TRUE(model.getMaterials().hasDielectrics());
}

