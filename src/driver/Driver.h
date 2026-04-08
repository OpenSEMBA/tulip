#pragma once

#include "DriverOptions.h"
#include "Model.h"
#include "Results.h"
#include "Solver.h"

#include <vector>

namespace tulip {

struct SolvedProblem {
    std::unique_ptr<Solver> solver;
    std::map<ConductorId,SolverSolution> solutions;
};

class Driver {
public:
    Driver(Model&& model, const DriverOptions& opts);

    PULParameters getPULMTL();
    InCellPotentials getInCellPotentials();
    
    MultiwireParametersByDomain getMultiwireParametersByDomains();
    std::map<ConductorId, double> getFloatingPotentials(
	    ConductorId prescribedId, 
        bool ignoreDielectrics);

    void setExportFolder(const std::string folder) { opts_.exportFolder = folder; }

    static SolverInputs buildSolverInputsFromModel(
        const Model& model,
        bool ignoreDielectrics);
    DenseMatrix getCFromGeneralizedC(
        const mfem::DenseMatrix& gC,
        const Model::Openness&);

    void run();

    DenseMatrix getGeneralizedCMatrix(bool ignoreDielectrics = false);
    DenseMatrix getCMatrix();
    DenseMatrix getLMatrix();

    const Model& getModel() const { return model_; }

    static Driver loadFromAdaptedFile(const std::string& filename);
    static Driver adaptFromFile(const std::string& filename);

    SolvedProblem*
        getElectrostaticSolvedProblem() { return &electric_; }
    SolvedProblem* 
        getMagnetostaticSolvedProblem() { return &magnetic_; }
    
    void loadFloatingPotentials(
        SolvedProblem* sP,
        const std::map<ConductorId, double>& fp) const;
private:
    Model model_;
    DriverOptions opts_;
    SolvedProblem electric_, magnetic_;

    SolvedProblem solveForAllConductors(
        bool ignoreDielectrics);
    PULParameters buildPULParametersForModel();
    PULParameters buildGeneralizedLCMatrices();
    double getInnerRegionAveragePotential(
        const Solver& s,
        bool includeConductors);
    std::map<ConductorId, FieldReconstruction> getFieldParameters(
        bool ignoreDielectrics);
    std::map<ConductorId, double> computeFloatingPotentialsFromC(
        ConductorId prescribedId,
        const mfem::DenseMatrix& C) const;

};

}