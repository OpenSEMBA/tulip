#pragma once

#define MFEM_USE_DOUBLE
#include "mfem.hpp"

#include <cstddef>

namespace pulmtln {

using namespace mfem;

class H1_FESpace : public FiniteElementSpace
{
public:
   H1_FESpace(Mesh *m,
                 const int p, const int space_dim = 3,
                 const int type = BasisType::GaussLobatto,
                 int vdim = 1, int order = Ordering::byNODES);
   ~H1_FESpace();
private:
   const FiniteElementCollection *FEC_;
};

class ND_FESpace : public FiniteElementSpace
{
public:
   ND_FESpace(Mesh *m, const int p, const int space_dim,
                 int vdim = 1, int order = Ordering::byNODES);
   ~ND_FESpace();
private:
   const FiniteElementCollection *FEC_;
};

class RT_FESpace : public FiniteElementSpace
{
public:
   RT_FESpace(Mesh *m, const int p, const int space_dim,
                 int vdim = 1, int order = Ordering::byNODES);
   ~RT_FESpace();
private:
   const FiniteElementCollection *FEC_;
};

class L2_FESpace: public FiniteElementSpace
{
public:
   L2_FESpace(Mesh *m, const int p, const int space_dim,
                 int vdim = 1, int order = Ordering::byNODES);
   ~L2_FESpace();
private:
   const FiniteElementCollection *FEC_;
};

class DiscreteInterpolationOperator : public DiscreteLinearOperator {
public:
    DiscreteInterpolationOperator(
        FiniteElementSpace* dfes,
        FiniteElementSpace* rfes)
        : DiscreteLinearOperator(dfes, rfes) {}
    virtual ~DiscreteInterpolationOperator();
};


class DiscreteGradOperator : public DiscreteInterpolationOperator {
public:
    DiscreteGradOperator(FiniteElementSpace* dfes, FiniteElementSpace* rfes);
};

class DiscreteDivOperator : public DiscreteInterpolationOperator {
public:
    DiscreteDivOperator(FiniteElementSpace* dfes, FiniteElementSpace* rfes);
};


}