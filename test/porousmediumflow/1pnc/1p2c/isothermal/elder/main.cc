// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup OnePNCTests
 * \brief Elder problem with the Boussinesq approximation: solute mass and fluxes compared with reference values
 */
#include <config.h>

#include <algorithm>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <memory>
#include <string>
#include <vector>

#include <dune/common/parallel/mpihelper.hh>

#include <dumux/common/dumuxmessage.hh>
#include <dumux/common/initialize.hh>
#include <dumux/common/parameters.hh>
#include <dumux/common/properties.hh>
#include <dumux/experimental/assembly/multistagefvassembler.hh>
#include <dumux/experimental/timestepping/multistagemethods.hh>
#include <dumux/experimental/timestepping/multistagetimestepper.hh>
#include <dumux/io/grid/gridmanager_yasp.hh>
#include <dumux/io/vtkoutputmodule.hh>
#include <dumux/linear/istlsolvers.hh>
#include <dumux/linear/linearalgebratraits.hh>
#include <dumux/linear/linearsolvertraits.hh>
#include <dumux/nonlinear/newtonsolver.hh>

#include "properties.hh"

namespace Dumux {

/*!
 * \brief Solute mass fraction without buoyancy (Ra = 0) in the half layer [0, length] x [0, 1] in units of the height
 *
 * Each cosine mode \f$ f_m \cos(k_m x) \f$, \f$ k_m = m \pi / L \f$, of the top boundary data diffuses into the
 * layer, initially free of solute, as
 * \f$ \sinh(k_m z)/\sinh(k_m) - \sum_n 2 (-1)^{n+1} n \pi / (k_m^2 + n^2 \pi^2) \sin(n \pi z) e^{-(k_m^2 + n^2 \pi^2) t} \f$.
 */
template<class Scalar>
class ElderSolutionWithoutBuoyancy
{
public:
    template<class BoundaryData>
    ElderSolutionWithoutBuoyancy(const BoundaryData& topMassFraction, Scalar length, int numModes, int numQuadraturePoints = 40001)
    : length_(length), coefficients_(numModes, 0.0)
    {
        const Scalar h = length/(numQuadraturePoints - 1);
        for (int m = 0; m < numModes; ++m)
        {
            const Scalar k = wavenumber_(m);
            Scalar integral = 0.0;
            for (int i = 0; i < numQuadraturePoints; ++i)
            {
                const Scalar x = i*h;
                const Scalar weight = (i == 0 || i == numQuadraturePoints - 1) ? 0.5*h : h;
                integral += weight*topMassFraction(x)*std::cos(k*x);
            }
            coefficients_[m] = (m == 0 ? 1.0 : 2.0)/length*integral;
        }
    }

    Scalar operator()(Scalar x, Scalar z, Scalar t) const
    {
        const int numTerms = static_cast<int>(std::ceil(std::sqrt(40.0/(M_PI*M_PI*t)))) + 2;
        Scalar value = 0.0;
        for (std::size_t m = 0; m < coefficients_.size(); ++m)
        {
            const Scalar k = wavenumber_(m);
            Scalar profile = (m == 0) ? z : std::exp(k*(z - 1.0))*(1.0 - std::exp(-2.0*k*z))/(1.0 - std::exp(-2.0*k));
            for (int n = 1; n <= numTerms; ++n)
            {
                const Scalar rate = k*k + n*n*M_PI*M_PI;
                profile -= 2.0*(n % 2 == 1 ? 1.0 : -1.0)*n*M_PI/rate*std::sin(n*M_PI*z)*std::exp(-rate*t);
            }
            value += coefficients_[m]*std::cos(k*x)*profile;
        }
        return value;
    }

private:
    Scalar wavenumber_(int m) const
    { return m*M_PI/length_; }

    Scalar length_;
    std::vector<Scalar> coefficients_;
};

template<class Scalar>
std::shared_ptr<const Experimental::MultiStageMethod<Scalar>> makeTimeSteppingMethod(const std::string& name)
{
    if (name == "ImplicitEuler")
        return std::make_shared<Experimental::MultiStage::ImplicitEuler<Scalar>>();
    if (name == "DIRK2")
        return std::make_shared<Experimental::MultiStage::DIRKSecondOrderAlexander<Scalar>>();
    if (name == "DIRK3")
        return std::make_shared<Experimental::MultiStage::DIRKThirdOrderAlexander<Scalar>>();
    DUNE_THROW(Dune::InvalidStateException, "Unknown time stepping method " << name);
}

} // end namespace Dumux

int main(int argc, char** argv)
{
    using namespace Dumux;
    using TypeTag = Properties::TTag::TYPETAG;

    Dumux::initialize(argc, argv);
    const auto& mpiHelper = Dune::MPIHelper::instance();
    if (mpiHelper.rank() == 0)
        DumuxMessage::print(/*firstCall=*/true);

    Parameters::init(argc, argv);

    GridManager<GetPropType<TypeTag, Properties::Grid>> gridManager;
    gridManager.init();
    const auto& leafGridView = gridManager.grid().leafGridView();

    using GridGeometry = GetPropType<TypeTag, Properties::GridGeometry>;
    auto gridGeometry = std::make_shared<GridGeometry>(leafGridView);

    using Problem = GetPropType<TypeTag, Properties::Problem>;
    auto problem = std::make_shared<Problem>(gridGeometry);

    using SolutionVector = GetPropType<TypeTag, Properties::SolutionVector>;
    SolutionVector x;
    problem->applyInitialSolution(x);
    auto xOld = x;

    using GridVariables = GetPropType<TypeTag, Properties::GridVariables>;
    auto gridVariables = std::make_shared<GridVariables>(problem, gridGeometry);
    gridVariables->init(x);

    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using FluidSystem = GetPropType<TypeTag, Properties::FluidSystem>;
    using Indices = typename GetPropType<TypeTag, Properties::ModelTraits>::Indices;
    constexpr int soluteIdx = FluidSystem::comp1Idx;
    constexpr int soluteEqIdx = Indices::conti0EqIdx + FluidSystem::comp1Idx;

    // Rayleigh number Ra = K drho g H / (mu D_eff) and time scale phi H^2 / D_eff of the homogeneous layer
    const Scalar height = gridGeometry->bBoxMax()[1] - gridGeometry->bBoxMin()[1];
    const Scalar length = gridGeometry->bBoxMax()[0] - gridGeometry->bBoxMin()[0];
    auto fvGeometry = localView(*gridGeometry);
    auto elemVolVars = localView(gridVariables->curGridVolVars());
    auto elemFluxVarsCache = localView(gridVariables->gridFluxVarsCache());
    const auto element = *elements(leafGridView).begin();
    fvGeometry.bind(element);
    elemVolVars.bind(element, fvGeometry, x);
    const auto& volVars = elemVolVars[*scvs(fvGeometry).begin()];
    const Scalar density = volVars.density(0);
    const Scalar porosity = volVars.porosity();
    const Scalar effectiveDiffusivity = volVars.effectiveDiffusionCoefficient(0, FluidSystem::comp0Idx, soluteIdx);
    const Scalar g = -problem->spatialParams().gravity(gridGeometry->bBoxMin())[1];
    const Scalar rayleigh = volVars.permeability()*problem->densityDifference()*g*height
                            /(volVars.viscosity(0)*effectiveDiffusivity);
    const Scalar timeScale = porosity*height*height/effectiveDiffusivity;

    // dimensionless solute mass int X dA / H^2 and mean vertical gradients of X at the top and the bottom times H;
    // the solute mass of the Dirichlet control volumes is constant, so their balance gives the boundary fluxes
    using LocalResidual = GetPropType<TypeTag, Properties::LocalResidual>;
    const LocalResidual localResidual(problem.get());
    using ElementBoundaryTypes = GetPropType<TypeTag, Properties::ElementBoundaryTypes>;
    ElementBoundaryTypes elemBcTypes;
    struct Diagnostics { Scalar mass, sherwood, bottom; };
    const auto diagnostics = [&]()
    {
        Diagnostics d{0.0, 0.0, 0.0};
        for (const auto& element : elements(leafGridView))
        {
            fvGeometry.bind(element);
            for (const auto& scv : scvs(fvGeometry))
                d.mass += x[scv.dofIndex()][soluteIdx]*scv.volume();

            if (!fvGeometry.hasBoundaryScvf())
                continue;

            elemVolVars.bind(element, fvGeometry, x);
            elemFluxVarsCache.bind(element, fvGeometry, elemVolVars);
            elemBcTypes.update(*problem, element, fvGeometry);
            const auto residual = localResidual.evalFluxAndSource(element, fvGeometry, elemVolVars, elemFluxVarsCache, elemBcTypes);
            for (const auto& scv : scvs(fvGeometry))
            {
                if (problem->onTop(scv.dofPosition()))
                    d.sherwood += residual[scv.localDofIndex()][soluteEqIdx];
                else if (problem->onBottom(scv.dofPosition()))
                    d.bottom -= residual[scv.localDofIndex()][soluteEqIdx];
            }
        }
        d.mass /= height*height;
        const Scalar fluxScale = density*effectiveDiffusivity*length/height;
        d.sherwood /= fluxScale;
        d.bottom /= fluxScale;
        return d;
    };

    const Scalar dt = getParam<Scalar>("TimeLoop.DimensionlessTimeStepSize")*timeScale;
    auto outputTimes = getParam<std::vector<Scalar>>("TimeLoop.DimensionlessOutputTimes");
    std::sort(outputTimes.begin(), outputTimes.end());

    const auto timeSteppingMethod = makeTimeSteppingMethod<Scalar>(getParam<std::string>("TimeLoop.Method", "ImplicitEuler"));
    using Assembler = Experimental::MultiStageFVAssembler<TypeTag, DiffMethod::numeric>;
    auto assembler = std::make_shared<Assembler>(problem, gridGeometry, gridVariables, timeSteppingMethod, xOld);

    using LinearSolver = UMFPackIstlSolver<SeqLinearSolverTraits, LinearAlgebraTraitsFromAssembler<Assembler>>;
    auto linearSolver = std::make_shared<LinearSolver>();
    using NonLinearSolver = NewtonSolver<Assembler, LinearSolver>;
    auto nonLinearSolver = std::make_shared<NonLinearSolver>(assembler, linearSolver);
    Experimental::MultiStageTimeStepper<NonLinearSolver> timeStepper(nonLinearSolver, timeSteppingMethod);

    VtkOutputModule vtkWriter(*gridVariables, x, problem->name());
    using VelocityOutput = GetPropType<TypeTag, Properties::VelocityOutput>;
    vtkWriter.addVelocityOutput(std::make_shared<VelocityOutput>(*gridVariables));
    using IOFields = GetPropType<TypeTag, Properties::IOFields>;
    IOFields::initOutputModule(vtkWriter);
    vtkWriter.write(0.0);

    std::ofstream logFile(problem->name() + "_diagnostics.csv");
    logFile << std::setprecision(15) << "time,mass,sherwood,bottom\n";
    const auto initial = diagnostics();
    logFile << 0.0 << "," << initial.mass << "," << initial.sherwood << "," << initial.bottom << "\n";

    // Ra = 0: L2 error of the mass fraction against the exact solution relative to the L2 norm of the exact solution,
    // with the sub-control volumes as quadrature
    std::unique_ptr<ElderSolutionWithoutBuoyancy<Scalar>> exactSolution;
    if (problem->densityDifference() == 0.0)
        exactSolution = std::make_unique<ElderSolutionWithoutBuoyancy<Scalar>>(
            [&](Scalar xi){ return problem->sourceMassFraction(xi*height + gridGeometry->bBoxMin()[0]); },
            length/height, getParam<int>("Problem.NumExactSolutionModes", 400));
    const auto errorWithoutBuoyancy = [&](Scalar t)
    {
        Scalar error = 0.0, norm = 0.0;
        for (const auto& element : elements(leafGridView))
        {
            fvGeometry.bindElement(element);
            for (const auto& scv : scvs(fvGeometry))
            {
                const auto pos = scv.dofPosition() - gridGeometry->bBoxMin();
                const Scalar exact = (*exactSolution)(pos[0]/height, pos[1]/height, t);
                const Scalar difference = x[scv.dofIndex()][soluteIdx] - exact;
                error += difference*difference*scv.volume();
                norm += exact*exact*scv.volume();
            }
        }
        return std::sqrt(error/norm);
    };

    // each interval between output times is split into equal steps of at most dt; the time is computed from the
    // step index, so the last step of an interval ends exactly at the output time
    std::vector<Diagnostics> results;
    std::vector<Scalar> exactErrors;
    Scalar intervalStart = 0.0;
    for (const auto outputTime : outputTimes)
    {
        const Scalar intervalEnd = outputTime*timeScale;
        const int numSteps = std::max(1, static_cast<int>(std::ceil((intervalEnd - intervalStart)/dt - 1e-6)));
        const Scalar stepSize = (intervalEnd - intervalStart)/numSteps;
        for (int step = 0; step < numSteps; ++step)
        {
            timeStepper.step(x, intervalStart + step*stepSize, stepSize);
            xOld = x;
            gridVariables->advanceTimeStep();

            const bool atOutputTime = step + 1 == numSteps;
            const Scalar time = atOutputTime ? intervalEnd : intervalStart + (step + 1)*stepSize;
            const auto d = diagnostics();
            logFile << time/timeScale << "," << d.mass << "," << d.sherwood << "," << d.bottom << "\n";
            if (atOutputTime)
            {
                results.push_back(d);
                vtkWriter.write(time);
                if (exactSolution)
                    exactErrors.push_back(errorWithoutBuoyancy(outputTime));
            }
        }
        std::cout << "Reached t = " << outputTime << " with " << numSteps << " steps of size "
                  << stepSize/timeScale << std::endl;
        intervalStart = intervalEnd;
    }

    std::cout << std::setprecision(12) << "Ra = " << rayleigh << ", time scale " << timeScale << " s" << std::endl;
    for (std::size_t i = 0; i < results.size(); ++i)
    {
        std::cout << "t = " << outputTimes[i] << ": mass " << results[i].mass << ", Sherwood " << results[i].sherwood
                  << ", bottom " << results[i].bottom;
        if (exactSolution)
            std::cout << ", relative L2 error without buoyancy " << exactErrors[i];
        std::cout << std::endl;
    }

    if (mpiHelper.rank() == 0)
        DumuxMessage::print(/*firstCall=*/false);

    // relative errors of the solute mass and the Sherwood number at the output times against reference values
    if (hasParam("Problem.ReferenceMass"))
    {
        const auto referenceMass = getParam<std::vector<Scalar>>("Problem.ReferenceMass");
        const auto referenceSherwood = getParam<std::vector<Scalar>>("Problem.ReferenceSherwood");
        const auto tolerance = getParam<Scalar>("Problem.RelativeTolerance");
        if (referenceMass.size() != results.size() || referenceSherwood.size() != results.size())
            DUNE_THROW(Dune::InvalidStateException, "One reference value per output time is required");
        for (std::size_t i = 0; i < results.size(); ++i)
        {
            const Scalar massError = std::abs(results[i].mass - referenceMass[i])/referenceMass[i];
            const Scalar sherwoodError = std::abs(results[i].sherwood - referenceSherwood[i])/referenceSherwood[i];
            std::cout << "t = " << outputTimes[i] << ": relative error of the mass " << massError
                      << ", of the Sherwood number " << sherwoodError << std::endl;
            if (!(massError < tolerance) || !(sherwoodError < tolerance))
                DUNE_THROW(Dune::Exception, "Error against the reference exceeds the tolerance " << tolerance);
        }
    }

    if (hasParam("Problem.ExactSolutionTolerance"))
    {
        const auto tolerance = getParam<Scalar>("Problem.ExactSolutionTolerance");
        for (const auto error : exactErrors)
            if (!(error < tolerance))
                DUNE_THROW(Dune::Exception, "Relative L2 error " << error << " against the solution without buoyancy exceeds the tolerance");
    }

    return 0;
}
