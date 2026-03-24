#include "gtest/gtest.h"

#include "TestUtils.h"

#include "Domain.h"
#include "Parser.h"

using namespace pulmtln;

class DomainTest : public ::testing::Test {};

TEST_F(DomainTest, build_domains_for_empty_coax)
{
	auto model{ Parser{inputCase("empty_coax") }.readModel() };
	
	auto domains{ Domain::buildDomains(model) };

	ASSERT_EQ(1, domains.size());
	EXPECT_EQ(1, domains.count(0));

	EXPECT_EQ(0,                         domains.at(0).ground);
	EXPECT_EQ(IdSet({0,1}),              domains.at(0).conductorIds);
	EXPECT_EQ(model.getMesh()->GetNE(),  domains.at(0).elems.size());
	EXPECT_EQ(model.getMesh()->GetNBE(), domains.at(0).bdrElems.size());
}

TEST_F(DomainTest, nested_coax_domains)
{
	auto model{ Parser{ inputCase("nested_coax")}.readModel() };

	auto domains{ Domain::buildDomains(model) };

	ASSERT_EQ(2, domains.size());
	ASSERT_EQ(1, domains.count(0));
	ASSERT_EQ(1, domains.count(1));

	EXPECT_EQ(0,                         domains.at(0).ground);
	EXPECT_EQ(IdSet({0,1}),              domains.at(0).conductorIds);
	EXPECT_GT(model.getMesh()->GetNE(),  domains.at(0).elems.size());
	EXPECT_GT(model.getMesh()->GetNBE(), domains.at(0).bdrElems.size());

	EXPECT_EQ(1,                         domains.at(1).ground);
	EXPECT_EQ(IdSet({1,2}),              domains.at(1).conductorIds);
	EXPECT_GT(model.getMesh()->GetNE(),  domains.at(1).elems.size());
	EXPECT_GT(model.getMesh()->GetNBE(), domains.at(1).bdrElems.size());
}

TEST_F(DomainTest, nested_coax_build_model_for_domain)
{
	auto model{ Parser{inputCase("nested_coax")}.readModel() };

	auto domains{ Domain::buildDomains(model) };
	ASSERT_EQ(2, domains.size());
	ASSERT_EQ(1, domains.count(0));
	ASSERT_EQ(1, domains.count(1));

	{
		const auto& dom{ domains.at(0) };
		auto sm{ 
			Domain::buildModelForDomain(
				*model.getMesh(), 
				model.getMaterials(), 
				dom) 
		};

		EXPECT_EQ(dom.elems.size(),    sm.getMesh()->GetNE());
		EXPECT_EQ(dom.bdrElems.size(), sm.getMesh()->GetNBE());
	}

	{
		const auto& dom{ domains.at(1) };
		auto sm{
			Domain::buildModelForDomain(
				*model.getMesh(),
				model.getMaterials(),
				dom)
		};

		EXPECT_EQ(dom.elems.size(), sm.getMesh()->GetNE());
		EXPECT_EQ(dom.bdrElems.size(), sm.getMesh()->GetNBE());
	}
}

