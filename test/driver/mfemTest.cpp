#include <gtest/gtest.h>

#include "FES.h"

using namespace mfem;
using namespace pulmtln;

class mfemTest : public ::testing::Test {
};


TEST_F(mfemTest, surface_integral)
{
	auto mesh{
		Mesh::MakeCartesian2D(10, 10, Element::QUADRILATERAL, 1.0, 1.0)
	};
	RT_FESpace hdiv(&mesh, 4, mesh.Dimension());
	
	GridFunction d(&hdiv);
	VectorConstantCoefficient xU(Vector({1.0, 0.0}));
	d.ProjectCoefficient(xU);
	
	{
		LinearForm surf_int(&hdiv);
		surf_int.AddBoundaryIntegrator(new VectorFEBoundaryFluxLFIntegrator);
		surf_int.Assemble();
		
		auto res{ surf_int(d) };

		EXPECT_NEAR(0.0, res, 1e-12);
	}

	{
		LinearForm left_bdr_int(&hdiv);
		
		Array<int> attr(1);
		attr[0] = 4;

		Array<Coefficient*> coefs(1);
		ConstantCoefficient one{ 1.0 };
		coefs[0] = &one;

		PWCoefficient pwcoeff{ attr, coefs };

		left_bdr_int.AddBoundaryIntegrator(
			new VectorFEBoundaryFluxLFIntegrator(pwcoeff));
		left_bdr_int.Assemble();
		
		auto res{ left_bdr_int(d) };

		EXPECT_NEAR(-1.0, res, 1e-12);
	}

	{
		LinearForm right_bdr_int(&hdiv);

		Array<int> attr(1);
		attr[0] = 2;

		Array<Coefficient*> coefs(1);
		ConstantCoefficient one{ 1.0 };
		coefs[0] = &one;

		PWCoefficient pwcoeff{ attr, coefs };

		right_bdr_int.AddBoundaryIntegrator(
			new VectorFEBoundaryFluxLFIntegrator(pwcoeff));
		right_bdr_int.Assemble();

		auto res{ right_bdr_int(d) };

		EXPECT_NEAR(1.0, res, 1e-12);
	}
}