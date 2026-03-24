#include <gtest/gtest.h>

#include "multipolarExpansion.h"

using namespace mfem;
using namespace pulmtln;

class multipolarExpansionTest : public ::testing::Test {
};


TEST_F(multipolarExpansionTest, dipole)
{
	Vector expansionCenter({ 0.0, 0.0 });
	double d = 0.1;
	multipolarCoefficients ab{ 
		{0.0, 0.0},
		{d, 0.0}
	};

	{
		double r = 1.0;
		Vector pos({ r, 0.0 });
		double vComputed = multipolarExpansion(pos, ab, expansionCenter);
		double vExpected = 1.0 / (2.0 * M_PI) * std::log((r + d / 2.0) / (r - d / 2.0));

		EXPECT_NEAR(vExpected, vComputed, 1e-4);
	}

	{
		double r = 1.0;
		Vector pos({ 0.0, r });
		double vComputed = multipolarExpansion(pos, ab, expansionCenter);
		double vExpected = 0.0;

		EXPECT_NEAR(vExpected, vComputed, 1e-8);
	}
}