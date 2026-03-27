#include "FES.h"

namespace pulmtln {

H1_FESpace::H1_FESpace(Mesh *m,
                             const int p, const int space_dim, const int type,
                             int vdim, int order)
   : FiniteElementSpace(m, new H1_FECollection(p,space_dim,type),vdim,order)
{
   FEC_ = this->FiniteElementSpace::fec;
}

H1_FESpace::~H1_FESpace()
{
    delete FEC_;
}

ND_FESpace::ND_FESpace(Mesh *m, const int p, const int space_dim,
                             int vdim, int order)
   : FiniteElementSpace(m, new ND_FECollection(p,space_dim),vdim,order)
{
   FEC_ = this->FiniteElementSpace::fec;
}

ND_FESpace::~ND_FESpace()
{
   delete FEC_;
}

RT_FESpace::RT_FESpace(Mesh *m, const int p, const int space_dim,
                             int vdim, int order)
   : FiniteElementSpace(m, new RT_FECollection(p-1,space_dim),vdim,order)
{
   FEC_ = this->FiniteElementSpace::fec;
}

RT_FESpace::~RT_FESpace()
{
   delete FEC_;
}

L2_FESpace::L2_FESpace(Mesh *m, const int p, const int space_dim,
                             int vdim, int order)
   : FiniteElementSpace(m, new L2_FECollection(p,space_dim),vdim,order)
{
   FEC_ = this->FiniteElementSpace::fec;
}

L2_FESpace::~L2_FESpace()
{
   delete FEC_;
}


DiscreteInterpolationOperator::~DiscreteInterpolationOperator()
{}

DiscreteGradOperator::DiscreteGradOperator(FiniteElementSpace* dfes,
    FiniteElementSpace* rfes)
    : DiscreteInterpolationOperator(dfes, rfes)
{
    this->AddDomainInterpolator(new GradientInterpolator);
}

DiscreteDivOperator::DiscreteDivOperator(FiniteElementSpace* dfes,
    FiniteElementSpace* rfes)
    : DiscreteInterpolationOperator(dfes, rfes)
{
    this->AddDomainInterpolator(new DivergenceInterpolator);
}

}