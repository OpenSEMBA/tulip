#include <gtest/gtest.h>

#include "Model.h"
#include "AdaptedInputParser.h"
#include "TestUtils.h"

using namespace tulip;

class ModelTest : public ::testing::Test {};

TEST_F(ModelTest, empty_coax_is_closed)
{
	EXPECT_EQ(
		Parser{ inputCase("empty_coax") }.readModel().determineOpenness(),
		Model::Openness::closed
	);
}

TEST_F(ModelTest, three_wires_ribbon_is_open)
{
	EXPECT_EQ(
		Parser{ inputCase("three_wires_ribbon") }.readModel().determineOpenness(),
		Model::Openness::open
	);
}

TEST_F(ModelTest, agrawal1981_conductors_in_mesh)
{
	auto m{ Parser{ inputCase("agrawal1981") }.readModel() };
	
	ASSERT_EQ(4, m.getMaterials().getConductors().size());
}

TEST_F(ModelTest, lansink2024_fdtd_cell_material_areas)
{
	const std::string CASE{ "lansink2024_fdtd_cell" };
	auto model{ Parser{ inputCase(CASE) }.readModel() };
	const double rTol{ 1e-8 };

	// Conductor areas (enclosed area computed from boundary integrals)
	Conductor conductor0{ 1, 0 };
	Conductor conductor1{ 2, 1 };
	EXPECT_NEAR(0.0, relError(3.14159253264265e-06, model.getAreaOfMaterial(&conductor0)), rTol);
	EXPECT_NEAR(0.0, relError(3.14159253264264e-04, model.getAreaOfMaterial(&conductor1)), rTol);

	// Dielectric domain areas
	Dielectric dielectric4{ 4 };
	Dielectric dielectric5{ 5 };
	EXPECT_NEAR(0.0, relError(0.0617876225248815, model.getAreaOfMaterial(&dielectric4)), rTol);
	EXPECT_NEAR(0.0, relError(0.0396826990779316, model.getAreaOfMaterial(&dielectric5)), rTol);
}

TEST_F(ModelTest, bounding_box_realistic_case_with_dielectrics_fdtd_cell)
{
	const std::string CASE{ "realistic_case_with_dielectrics_fdtd_cell" };
	auto fn =  casesFolder() + "realistic_case_with_dielectrics/" + CASE + ".tulip.adapted.json";
	auto model{ Parser{ fn }.readModel() };
	
	auto bb = model.getInnerRegionBoundingBox();
	EXPECT_NEAR(30e-3*30e-3, bb.area(), 1e-12);
}

TEST_F(ModelTest, bounding_box_realistic_case_with_dielectrics_inner_region)
{
	const std::string CASE{ "realistic_case_with_dielectrics_inner_region" };
	auto fn =  casesFolder() + "realistic_case_with_dielectrics/" + CASE + ".tulip.adapted.json";
	auto model{ Parser{ fn }.readModel() };
	
	auto bb = model.getInnerRegionBoundingBox();
	EXPECT_NEAR(16.5e-3*16.5e-3, bb.area(), 1e-12);
}

TEST_F(ModelTest, domains_for_empty_coax)
{
	auto model{ Parser{inputCase("empty_coax") }.readModel() };
	
	auto domains{ model.getDomains() };
	
	ASSERT_EQ(1, domains.size());
	EXPECT_EQ(1, domains.count(0));
	
	EXPECT_EQ(Domain::UNDEFINED_GROUND,  domains.at(0).ground);
	EXPECT_EQ(IdSet({0,1}),              domains.at(0).conductorIds);
	
	DomainTree dT(domains);
	EXPECT_EQ(IdSet({0,1}), dT.getConductorsInDomain(0));
	EXPECT_EQ(IdSet({0}), dT.getConductorsInsideConductor(0));
	EXPECT_EQ(IdSet({1}), dT.getConductorsInsideConductor(1));
}

TEST_F(ModelTest, domains_for_nested_coax)
{
	auto model{ Parser{ inputCase("nested_coax")}.readModel() };

	auto domains{ model.getDomains() };

	ASSERT_EQ(2, domains.size());
	ASSERT_EQ(1, domains.count(0));
	ASSERT_EQ(1, domains.count(1));

	EXPECT_EQ(Domain::UNDEFINED_GROUND,  domains.at(0).ground);
	EXPECT_EQ(IdSet({0,1}),              domains.at(0).conductorIds);

	EXPECT_EQ(1,                         domains.at(1).ground);
	EXPECT_EQ(IdSet({1,2}),              domains.at(1).conductorIds);

	DomainTree dT(domains);
	EXPECT_EQ(IdSet({0,1}), dT.getConductorsInDomain(0));
	EXPECT_EQ(IdSet({0,1,2}), dT.getConductorsInsideConductor(0));
	EXPECT_EQ(IdSet({1,2}), dT.getConductorsInsideConductor(1));
	EXPECT_EQ(IdSet({2}), dT.getConductorsInsideConductor(2));
}

TEST_F(ModelTest, domains_for_realistic_case_with_dielectrics_fdtd_cell)
{
	const std::string CASE{ "realistic_case_with_dielectrics_fdtd_cell" };
	auto fn =  casesFolder() + "realistic_case_with_dielectrics/" + CASE + ".tulip.adapted.json";
	auto model{ Parser{ fn }.readModel() };

	auto domains = model.getDomains();
	ASSERT_EQ(1, domains.size());
	ASSERT_EQ(1, domains.count(0));

	EXPECT_EQ(Domain::UNDEFINED_GROUND, domains.at(0).ground);
	EXPECT_EQ(31,                       domains.at(0).conductorIds.size());
}

TEST_F(ModelTest, domains_for_coax_and_bare_wire)
{
	auto model{ Parser{ inputCase("coax_and_bare_wire")}.readModel() };

	auto domains{ model.getDomains() };

	ASSERT_EQ(2, domains.size());
	ASSERT_EQ(1, domains.count(0));
	ASSERT_EQ(1, domains.count(1));

	DomainTree dT(domains);
	EXPECT_EQ(IdSet({1, 2}), dT.getConductorsInDomain(0));
	EXPECT_EQ(IdSet({0, 1}), dT.getConductorsInDomain(1));
	
	EXPECT_EQ(IdSet({0}), dT.getConductorsInsideConductor(0));
	EXPECT_EQ(IdSet({0,1}), dT.getConductorsInsideConductor(1));
	EXPECT_EQ(IdSet({2}), dT.getConductorsInsideConductor(2));
}
