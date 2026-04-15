#include "Solver.h"

#include "constants.h"
#include "AttrToValueMap.h"

namespace tulip {

Array<int> AttrToMarker(const Mesh& mesh, const Array<int>& attrs)
{
    if (attrs.Size() == 0) {
        return Array<int>();
    }

    int maxAttrInConditions{*std::max_element(attrs.begin(), attrs.end())};
    int maxAttrInMesh{ mesh.bdr_attributes.Max() };
    int maxAttr{ std::max({maxAttrInConditions, maxAttrInMesh}) };

    Array<int> marker(maxAttr);
    marker = 0;

    
    MFEM_ASSERT(attrs.Max() <= maxAttr, "Invalid attribute number present.");

    if (attrs.Size() == 1 && attrs[0] == -1)
    {
        marker = 1;
    }
    else
    {
        marker = 0;
        for (int j = 0; j < attrs.Size(); j++)
        {
            int attr = attrs[j];
            MFEM_VERIFY(attr > 0, "Attribute number less than one!");
            marker[attr - 1] = 1;
        }
    }

    return marker;
}

Array<int> toArray(const std::vector<int>& v)
{
    Array<int> r{ (int)v.size() };
    for (int i{ 0 }; i < r.Size(); i++) {
        r[i] = v[i];
    }
    return r;
}

double firstOrderABC(const Vector& rVec)
{
    double r{ rVec.Norml2() };

    return - EPSILON0_NATURAL / (r * std::log(r));
}


Solver::Solver(
    Mesh& mesh,
    const SolverInputs& parameters,
    const SolverOptions opts) : 
    opts_(opts),
    mesh_(&mesh),
    parameters_(parameters),
    H1FESpace_(NULL),
    HCurlFESpace_(NULL),
    HDivFESpace_(NULL),
    L2FESpace_(NULL),
    divEpsGrad_(NULL),
    hDivMass_(NULL),
    hCurlHDivEps_(NULL),
    h1SurfMass_(NULL),
    rhod_(NULL),
    grad_(NULL),
    phi_(NULL),
    sigma_src_(NULL),
    e_(NULL),
    d_(NULL),
    oneCoef_(1.0)
{
    
    // Define compatible parallel finite element spaces on the parallel
    // mesh. Here we use arbitrary order H1, Nedelec, and Raviart-Thomas finite
    // elements.
    auto order{ opts_.order };
    H1FESpace_ = new H1_FESpace(mesh_, order, mesh_->Dimension());
    HCurlFESpace_ = new ND_FESpace(mesh_, order, mesh_->Dimension());
    HDivFESpace_ = new RT_FESpace(mesh_, order, mesh_->Dimension());
    L2FESpace_ = new L2_FESpace(mesh_, order - 1, mesh_->Dimension());

    // Select surface attributes for Dirichlet BCs
    ess_bdr_ = AttrToMarker(
        *mesh_,
        toArray(parameters.dirichletBoundaries.getAttributesAsArray())
    );

    // Setup domain permittivity coefficients.
    if (parameters_.domainPermittivities.empty()) {
        epsCoef_ = new ConstantCoefficient(EPSILON0_NATURAL);
    } else {
        mfem::Vector eps(mesh_->attributes.Max());
        eps = EPSILON0_NATURAL;
        for (const auto& [attr, epsr] : parameters_.domainPermittivities) {
            const auto epsSize{ eps.Size() };
            assert(attr <= epsSize);
            eps[attr-1] *= epsr;
        }
        epsCoef_ = new PWConstCoefficient(eps);
    }
    
    // Bilinear Forms
    divEpsGrad_ = new BilinearForm(H1FESpace_);
    divEpsGrad_->AddDomainIntegrator(new DiffusionIntegrator(*epsCoef_));

    // Setup open regions.
    std::unique_ptr<Coefficient> openBoundaryCoeff;
    if (!parameters.openBoundaries.empty()) {
        open_bdr_ = AttrToMarker(*mesh_, toArray(parameters.openBoundaries));
        openBoundaryCoeff.reset(new FunctionCoefficient(firstOrderABC));
        divEpsGrad_->AddBoundaryIntegrator(
            new BoundaryMassIntegrator(*openBoundaryCoeff),
            open_bdr_
        );
    }

    hDivMass_ = new BilinearForm(HDivFESpace_);
    hDivMass_->AddDomainIntegrator(new VectorFEMassIntegrator);
        
    h1SurfMass_ = new BilinearForm(H1FESpace_);
    h1SurfMass_->AddBoundaryIntegrator(new MassIntegrator);
    
    hCurlHDivEps_ = new MixedBilinearForm(HCurlFESpace_, HDivFESpace_);
    hCurlHDivEps_->AddDomainIntegrator(new VectorFEMassIntegrator(*epsCoef_));

    rhod_ = new LinearForm(H1FESpace_);

    // Discrete derivative operator
    grad_ = new DiscreteGradOperator(H1FESpace_, HCurlFESpace_);
    
    // Build grid functions
    phi_ = new GridFunction(H1FESpace_);
    d_ = new GridFunction(HDivFESpace_);
    e_ = new GridFunction(HCurlFESpace_);
    sigma_src_ = new GridFunction(H1FESpace_);

    Assemble();

    // Apply Neumann conditions on sigma_src.
    *sigma_src_ = 0.0;
    applyBoundaryConstantValuesToGridFunction(parameters_.neumannBoundaries, *sigma_src_);
}

Solver::~Solver()
{   
    delete phi_;
    delete rhod_;
    delete d_;
    delete e_;
    delete sigma_src_;
    
    delete grad_;

    delete divEpsGrad_;
    delete hDivMass_;
    delete hCurlHDivEps_;
    delete h1SurfMass_;
    
    delete H1FESpace_;
    delete HCurlFESpace_;
    delete HDivFESpace_;
    delete L2FESpace_;

    delete epsCoef_;
}

void Solver::Assemble()
{
    divEpsGrad_->Assemble();
    divEpsGrad_->Finalize();

    hDivMass_->Assemble();
    hDivMass_->Finalize();

    hCurlHDivEps_->Assemble();
    hCurlHDivEps_->Finalize();

    h1SurfMass_->Assemble();
    h1SurfMass_->Finalize();

    rhod_->Assemble();

    grad_->Assemble();
    grad_->Finalize();
}

void Solver::applyBoundaryConstantValuesToGridFunction(
    const AttrToValueMap& bdrValues,
    GridFunction& gf) const
{
    auto attributes{ bdrValues.getAttributesAsArray() };
    auto values{ bdrValues.getValuesAsArray() };
    if (attributes.empty()) {
        return;
    }

    // Apply piecewise constant boundary condition
    Array<int> bdr_attr(mesh_->bdr_attributes.Max());
    for (size_t i = 0; i < attributes.size(); i++)
    {
        ConstantCoefficient val(values[i]);
        bdr_attr = 0;
        if (attributes[i] <= bdr_attr.Size())
        {
            bdr_attr[attributes[i] - 1] = 1;
        }
        gf.ProjectBdrCoefficient(val, bdr_attr);
    }
}

void Solver::Solve()
{
    // Initialize the surface charge density (From Neumann boundaries).
    h1SurfMass_->Mult(*sigma_src_, *rhod_);
    
    // Solves phi (electrostatic potential).
    {
        *phi_ = 0.0;
        auto dbcs{ parameters_.dirichletBoundaries.getAttributesAsArray() };
        applyBoundaryConstantValuesToGridFunction(parameters_.dirichletBoundaries, *phi_);

        // Determine the essential BC degrees of freedom
        if (dbcs.size() > 0) {
            H1FESpace_->GetEssentialTrueDofs(ess_bdr_, ess_bdr_tdofs_);
        }
        else {
            ess_bdr_tdofs_.SetSize(1);
            ess_bdr_tdofs_[0] = 0;
        }


        // Apply essential BC and form linear system
        SparseMatrix DivEpsGrad;
        Vector Phi;
        Vector RHS;

        divEpsGrad_->FormLinearSystem(
            ess_bdr_tdofs_, *phi_, *rhod_, 
            DivEpsGrad, Phi, RHS);

        GSSmoother M(DivEpsGrad);
        PCG(DivEpsGrad, M, RHS, Phi, opts_.printIterations, 500, 1e-12, 0.0);

        divEpsGrad_->RecoverFEMSolution(Phi, *rhod_, *phi_);
    }

    // Computes E.
    grad_->Mult(*phi_, *e_); *e_ *= -1.0;
    
    // Computes D.
    {
        GridFunction ed(HDivFESpace_);
        hCurlHDivEps_->Mult(*e_, ed);
    
        SparseMatrix MassHDiv;
        Vector ED, D;

        Array<int> dbc_dofs_d;
        hDivMass_->FormLinearSystem(dbc_dofs_d, *d_, ed, MassHDiv, D, ED);

        GSSmoother M(MassHDiv);
        PCG(MassHDiv, M, ED, D, opts_.printIterations, 500, 1e-12);
        
        hDivMass_->RecoverFEMSolution(D, ed, *d_);
    }
        
}

void Solver::setDirichletConditions(const AttrToValueMap& dbcs)
{
    parameters_.dirichletBoundaries = dbcs;
}

void Solver::setNeumannCondition(
    const int bdrAttribute, 
    Coefficient& chargeDensity)
{
    Array<int> bdr_attr(mesh_->bdr_attributes.Max());
    bdr_attr = 0;
    bdr_attr[bdrAttribute - 1] = 1;
    
    sigma_src_->ProjectBdrCoefficient(chargeDensity, bdr_attr);
}

double Solver::getTotalCharge() const
{
    LinearForm rt_surf_int{HDivFESpace_};
    rt_surf_int.AddBoundaryIntegrator(new VectorFEBoundaryFluxLFIntegrator);
    rt_surf_int.Assemble();
    return rt_surf_int(*d_);
}

double Solver::getChargeMomentComponent(
    int n, int d, const Vector& center) const
{
    // Computes the 2^n-polar term component for direction d around an expansion center.
    // d == 0, corresponds with $a_n$
    // d == 1, corresponds with $b_n$
    
    Array<int> chargedBoundaries;
    for (const auto& [attr, val] : parameters_.dirichletBoundaries) {
        chargedBoundaries.Append(attr);
    }
    for (const auto& [attr, val] : parameters_.neumannBoundaries) {
        chargedBoundaries.Append(attr);
    }

    if (n == 0) {
        if (d == 0) {
            double Qt = 0.0;
            for (auto b : chargedBoundaries) {
                Qt += getChargeInBoundary(b);
            }
            return Qt; // a_0
        }
        else {
            return 0.0;  // b_0
        }
    }

    std::function<double(Vector)> xComponent =
        std::bind(momentComponent, std::placeholders::_1, n, d, center);
    FunctionCoefficient xComponentFunctionCoeff(xComponent);
    ProductCoefficient weight{ -1.0, xComponentFunctionCoeff };

    Array<Coefficient*> coefsArray(chargedBoundaries.Size());
    for (int i = 0; i < coefsArray.Size(); i++) {
        coefsArray[i] = &weight;
    }
    PWCoefficient pwcoeff{ chargedBoundaries, coefsArray };

    LinearForm surf_int(HDivFESpace_);
    surf_int.AddBoundaryIntegrator(new VectorFEBoundaryFluxLFIntegrator(pwcoeff));
    surf_int.Assemble();

    return surf_int(*d_);
}

Vector Solver::getCenterOfCharge() const
{
    // "Center of Charge" is the place where the dipole moment is zero.
    // It can only be defined for open and non-neutral systems.

    if (parameters_.openBoundaries.size() !=  1) {
        throw std::runtime_error("Not implemented for closed problems.");
    }
    const int& openBoundaryAttr = parameters_.openBoundaries.front();
    auto Q{ this->getTotalCharge() - getChargeInBoundary(openBoundaryAttr)};

    Vector origin(2);
    origin = 0.0;

    Vector res(2);
    for (int x = 0; x < 2; ++x) {
        res(x) = getChargeMomentComponent(1,x, origin) / Q;
    }
    
    return res;
}

multipolarCoefficients Solver::getMultipolarCoefficients(
    std::size_t order) const
{
    auto centerOfCharge{ getCenterOfCharge() };

    multipolarCoefficients ab(order + 1);
    
    for (int n = 0; n < order + 1; n++) {
        ab[n] = {
            getChargeMomentComponent(n, 0, centerOfCharge),
            getChargeMomentComponent(n, 1, centerOfCharge)
        };
    }
    
    return ab;
}

std::unique_ptr<LinearForm> buildHDivBoundaryIntegrator(
    RT_FESpace* fes, 
    int bdrAttribute,
    Coefficient& coeff)
{
    mfem::Array<int> attr(1);
    attr[0] = bdrAttribute;

    Array<Coefficient*> coefsArray(attr.Size());
    for (int i = 0; i < coefsArray.Size(); i++) {
        coefsArray[i] = &coeff;
    }
    PWCoefficient pwcoeff{ attr, coefsArray };

    auto surf_int =  std::make_unique<LinearForm>(fes);
    surf_int->AddBoundaryIntegrator(new VectorFEBoundaryFluxLFIntegrator(pwcoeff));
    surf_int->Assemble();

    return surf_int;
}

std::unique_ptr<LinearForm> buildH1BoundaryIntegrator(
    H1_FESpace* fes, 
    int bdrAttribute,
    Coefficient& coeff)
{
    mfem::Array<int> attr(1);
    attr[0] = bdrAttribute;

    Array<Coefficient*> coeffsArray(attr.Size());
    for (int i = 0; i < coeffsArray.Size(); i++) {
        coeffsArray[i] = &coeff;
    }
    PWCoefficient pwcoeff{ attr, coeffsArray };

    auto surf_int = std::make_unique<LinearForm>(fes);
    surf_int->AddBoundaryIntegrator(new BoundaryLFIntegrator(pwcoeff));
    surf_int->Assemble();

    return surf_int;
}

std::unique_ptr<GridFunction> cloneGridFunction(GridFunction* gf)
{
    auto res = std::make_unique<GridFunction>(gf->FESpace());
    for (int i = 0; i < gf->Size(); ++i) {
        (*res)[i] = (*gf)[i];
    }
    return std::move(res);
}

void copyGridFunctionValues(GridFunction* tgt, const GridFunction* src)
{
    assert(tgt->Size() == src->Size());
    for (int i = 0; i < src->Size(); ++i) {
        (*tgt)[i] = (*src)[i];
    }
}

SolverSolution Solver::getSolution() const
{
    SolverSolution res;
    res.phi = cloneGridFunction(phi_);
    res.e = cloneGridFunction(e_);
    res.d = cloneGridFunction(d_);
    return res;
}

void Solver::setSolution(const SolverSolution& s)
{
    copyGridFunctionValues(phi_, s.phi.get());
    copyGridFunctionValues(e_, s.e.get());
    copyGridFunctionValues(d_, s.d.get());
}

double Solver::getChargeInBoundary(int bdrAttribute) const
{
    ConstantCoefficient minusOne{ -1.0 };
    auto surf_int{ buildHDivBoundaryIntegrator(HDivFESpace_, bdrAttribute, minusOne) };
    return (*surf_int)(*d_);
}

double Solver::getAveragePotentialInDomain(int domainAttribute) const
{
    mfem::Array<int> attr(1);
    attr[0] = domainAttribute;

    Array<Coefficient*> coeffsArray(attr.Size());
    ConstantCoefficient one;
    for (int i = 0; i < coeffsArray.Size(); i++) {
        coeffsArray[i] = &one;
    }
    PWCoefficient pwcoeff{ attr, coeffsArray };

    LinearForm domain_int(H1FESpace_);
    domain_int.AddDomainIntegrator(new DomainLFIntegrator(pwcoeff));
    domain_int.Assemble();

    GridFunction ones(H1FESpace_);
    ones = 1.0;
    
    double totalPotential = domain_int(*phi_);
    double area = domain_int(ones);

    return totalPotential / area;
}

double Solver::getAveragePotentialInBoundary(int bdrAttribute) const
{
    ConstantCoefficient one{ 1.0 };
    auto surf_int{ buildH1BoundaryIntegrator(H1FESpace_, bdrAttribute, one) };

    GridFunction ones(H1FESpace_);
    ones = 1.0;
    
    auto totalPotential = (*surf_int)(*phi_);
    auto totalLength = (*surf_int)(ones);

    return totalPotential / totalLength;
}

double Solver::getTotalEnergy() const
{
    Array<int> domainAttributes = mesh_->attributes;
    int maxAttr = domainAttributes.Max();

    Vector domainPermittivities(maxAttr);
    domainPermittivities = EPSILON0_NATURAL;
    for (const auto& [domainAttr, relPerm] : parameters_.domainPermittivities) {
        domainPermittivities[domainAttr - 1] *= relPerm;
    }
    PWConstCoefficient permittivities(domainPermittivities);
    
    BilinearForm mass{ HCurlFESpace_ };
    mass.AddDomainIntegrator(new VectorFEMassIntegrator(permittivities));
    mass.Assemble();
    
    GridFunction aux(e_->FESpace());
    mass.Mult(*e_, aux);

    double energy{ 0.0 };
    energy = 0.5*InnerProduct(*e_, aux);

    return energy;
}

void Solver::writeParaViewFields(ParaViewDataCollection& pv) const
{
    pv.SetHighOrderOutput(true);
    pv.SetLevelsOfDetail(3);
    pv.RegisterField("Phi", phi_);
    pv.RegisterField("D", d_);
    pv.RegisterField("E", e_);

    pv.Save();
}

}