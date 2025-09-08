#pragma once

#include "SolverOptions.h"
#include "Model.h"
#include "Results.h"
#include "ElectrostaticSolver.h"

#include <vector>

namespace pulmtln {

struct SolvedProblem {
    std::unique_ptr<ElectrostaticSolver> solver;
    std::vector<SolverSolution> solutions;
};

class Driver {
public:
    Driver(Model&& model, const DriverOptions& opts);

    PULParameters getPULMTL();
    PULParametersByDomain getPULMTLByDomains();
    InCellPotentials getInCellPotentials();
    FloatingPotentials getFloatingPotentials();

    void setExportFolder(const std::string folder) { opts_.exportFolder = folder; }

    static SolverInputs buildSolverInputsFromModel(
        const Model& model,
        bool ignoreDielectrics);
    static DenseMatrix getCFromGeneralizedC(
        const mfem::DenseMatrix& gC,
        const Model::Openness&);

    void run();

    DenseMatrix getGeneralizedCMatrix(bool ignoreDielectrics = false);
    DenseMatrix getCMatrix();
    DenseMatrix getLMatrix();
    DenseMatrix getFloatingPotentialsMatrix(const bool ignoreDielectrics);

    static Driver loadFromFile(const std::string& filename);

private:
    Model model_;
    DriverOptions opts_;
    SolvedProblem electric_, magnetic_;

    SolvedProblem solveForAllConductors(
        bool ignoreDielectrics);
    PULParameters buildPULParametersForModel();
    double getInnerRegionAveragePotential(
        const ElectrostaticSolver& s,
        bool includeConductors);
    std::map<MaterialId, FieldReconstruction> getFieldParameters(
        bool ignoreDielectrics);
};

}