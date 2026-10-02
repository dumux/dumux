// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup TwoPTests
 * \brief Buckley-Leverett test solved with an IMPES scheme
 *        (implicit pressure, explicit saturation).
 *
 * In each time step \f$ t^n \to t^{n+1} \f$:
 *  1. Solve the pressure equation \f$ -\nabla\cdot(K\lambda_t(S_w^n)\nabla p^{n+1}) = 0 \f$
 *     with the single-phase model (see spatialparams_impes.hh).
 *  2. Compute the total volume flux \f$ q_{t,\sigma} \f$ over each face \f$ \sigma \f$.
 *  3. Update the saturation explicitly with first-order upwinding
 *     \f[ S_{w,i}^{n+1} = S_{w,i}^n - \frac{\Delta t}{\phi |V_i|}
 *         \sum_\sigma q_{t,\sigma} f_w(S^n_{w,\text{up}(\sigma)}). \f]
 *  The time step size is restricted by the CFL condition.
 *
 * The two-phase problem, spatial parameters, analytical solution and output
 * of the fully implicit test (main.cc) are reused.
 */
#include <config.h>

#include <algorithm>
#include <cmath>
#include <iostream>
#include <limits>
#include <sstream>
#include <type_traits>
#include <vector>

#include <dune/common/exceptions.hh>
#include <dune/common/fvector.hh>

#include <dumux/assembly/fvassembler.hh>
#include <dumux/common/initialize.hh>
#include <dumux/common/integrate.hh>
#include <dumux/common/parameters.hh>
#include <dumux/common/properties.hh>
#include <dumux/common/timeloop.hh>
#include <dumux/io/grid/gridmanager_yasp.hh>
#include <dumux/io/vtkoutputmodule.hh>
#include <dumux/linear/istlsolvers.hh>
#include <dumux/linear/linearalgebratraits.hh>
#include <dumux/linear/linearsolvertraits.hh>

#include "analyticsolution.hh"
#include "properties.hh"
#include "properties_impes.hh"

int main(int argc, char** argv)
{
    using namespace Dumux;
    using TwoPTypeTag = Properties::TTag::TwoPBuckleyLeverettTpfa;
    using PressureTypeTag = Properties::TTag::TwoPBuckleyLeverettImpesPressureTpfa;

    Dumux::initialize(argc, argv);
    Parameters::init(argc, argv);

    if (Dune::MPIHelper::getCommunication().size() > 1)
        DUNE_THROW(Dune::NotImplemented, "The explicit saturation update is implemented for sequential runs only");

    GridManager<GetPropType<TwoPTypeTag, Properties::Grid>> gridManager;
    gridManager.init();

    const auto& leafGridView = gridManager.grid().leafGridView();

    // both models share the grid geometry
    using GridGeometry = GetPropType<TwoPTypeTag, Properties::GridGeometry>;
    static_assert(std::is_same_v<GridGeometry, GetPropType<PressureTypeTag, Properties::GridGeometry>>);
    auto gridGeometry = std::make_shared<GridGeometry>(leafGridView);

    // the two-phase problem provides initial/boundary saturations, the material laws and the output
    using TwoPProblem = GetPropType<TwoPTypeTag, Properties::Problem>;
    auto problem = std::make_shared<TwoPProblem>(gridGeometry);

    using PressureProblem = GetPropType<PressureTypeTag, Properties::Problem>;
    auto pressureProblem = std::make_shared<PressureProblem>(gridGeometry);

    using Scalar = GetPropType<TwoPTypeTag, Properties::Scalar>;
    const auto tEnd = getParam<Scalar>("TimeLoop.TEnd");
    const auto maxDt = getParam<Scalar>("TimeLoop.MaxTimeStepSize");
    const auto dt = getParam<Scalar>("TimeLoop.DtInitial");
    const auto cflFactor = getParam<Scalar>("Impes.CFLFactor");
    const auto pressureUpdateInterval = getParam<int>("Impes.PressureUpdateInterval");

    // the two-phase solution vector (pw, Sn) is only used for output and the error evaluation
    using ModelTraits = GetPropType<TwoPTypeTag, Properties::ModelTraits>;
    constexpr auto pressureIdx = ModelTraits::Indices::pressureIdx;
    constexpr auto saturationIdx = ModelTraits::Indices::saturationIdx;
    static_assert(ModelTraits::priVarFormulation() == TwoPFormulation::p0s1, "Expects the (pw, Sn) formulation");

    using SolutionVector = GetPropType<TwoPTypeTag, Properties::SolutionVector>;
    SolutionVector sol(gridGeometry->numDofs());
    problem->applyInitialSolution(sol);

    // the wetting-phase saturation is the transported quantity
    std::vector<Scalar> sw(gridGeometry->numDofs());
    for (std::size_t i = 0; i < sw.size(); ++i)
        sw[i] = 1.0 - sol[i][saturationIdx];

    using GridVariables = GetPropType<TwoPTypeTag, Properties::GridVariables>;
    auto gridVariables = std::make_shared<GridVariables>(problem, gridGeometry);
    gridVariables->init(sol);

    // Pressure model
    //
    // Adding the two incompressible phase mass balances (divided by the constant densities)
    // and using S_w + S_n = 1 eliminates the time derivative:
    //     phi d(S_w + S_n)/dt + div(v_w + v_n) = div(v_t) = 0,   v_t = -K lambda_t(S_w) grad p.
    // Hence, the pressure equation is elliptic and stationary: at each time step, p is fully
    // determined by the current saturation (through lambda_t) and the boundary conditions.
    // It has no initial condition; the vector p is only the starting point of the linear solve.
    //
    // The pressure equation has the same form as the single-phase equation
    //     div(rho/mu_ref K_eff grad p) = 0,
    // so we solve it with the DuMux single-phase model (OneP) by setting the effective
    // permeability K_eff = K lambda_t mu_ref in the spatial params (see spatialparams_impes.hh).
    using PressureSolutionVector = GetPropType<PressureTypeTag, Properties::SolutionVector>;
    PressureSolutionVector p(gridGeometry->numDofs());
    p = 0.0;

    // The grid variables hold the volume variables (density, viscosity, permeability per cell)
    // and the flux variables cache, which stores the TPFA transmissibility of each face.
    using PressureGridVariables = GetPropType<PressureTypeTag, Properties::GridVariables>;
    auto pressureGridVariables = std::make_shared<PressureGridVariables>(pressureProblem, gridGeometry);
    pressureGridVariables->init(p);

    // The assembler builds the global system from the element-local residuals.
    // No time loop is passed, so the stationary (storage-free) residual is assembled.
    // DiffMethod::analytic: the incompressible single-phase local residual
    // (OnePIncompressibleLocalResidual) provides the exact derivatives d(flux)/dp,
    // which are simply the face transmissibilities.
    using PressureAssembler = FVAssembler<PressureTypeTag, DiffMethod::analytic>;
    auto pressureAssembler = std::make_shared<PressureAssembler>(pressureProblem, gridGeometry, pressureGridVariables);

    // A: Jacobian matrix dR/dp (one row/column per cell), r: residual vector R(p)
    using JacobianMatrix = GetPropType<PressureTypeTag, Properties::JacobianMatrix>;
    auto A = std::make_shared<JacobianMatrix>();
    auto r = std::make_shared<PressureSolutionVector>();
    pressureAssembler->setLinearSystem(A, r);

    // The TPFA pressure matrix is symmetric and, thanks to the Dirichlet boundary, positive
    // definite. We can therefore use the conjugate gradient (CG) method, preconditioned with
    // algebraic multigrid (AMG), which makes the number of iterations nearly independent of
    // the number of cells. Stopping criterion: residual reduction by the factor
    // LinearSolver.ResidualReduction (default 1e-13).
    using PressureLinearSolver = AMGCGIstlSolver<LinearSolverTraits<GridGeometry>, LinearAlgebraTraitsFromAssembler<PressureAssembler>>;
    auto pressureLinearSolver = std::make_shared<PressureLinearSolver>(gridGeometry->gridView(), gridGeometry->dofMapper());

    // fluid and material properties for the fractional flow function
    using FluidSystem = GetPropType<TwoPTypeTag, Properties::FluidSystem>;
    using FluidState = GetPropType<TwoPTypeTag, Properties::FluidState>;
    using GridView = typename GridGeometry::GridView;
    using Element = typename GridView::template Codim<0>::Entity;
    using GlobalPosition = typename Element::Geometry::GlobalCoordinate;

    const Scalar referencePressure = getParam<Scalar>("Problem.ReferencePressure");
    FluidState fluidState;
    fluidState.setTemperature(problem->spatialParams().temperatureAtPos(GlobalPosition{}));
    fluidState.setPressure(FluidSystem::phase0Idx, referencePressure);
    fluidState.setPressure(FluidSystem::phase1Idx, referencePressure);
    const Scalar viscosityW = FluidSystem::viscosity(fluidState, FluidSystem::phase0Idx);
    const Scalar viscosityN = FluidSystem::viscosity(fluidState, FluidSystem::phase1Idx);

    // viscosity of the fluid used by the pressure model
    using PressureFluidSystem = GetPropType<PressureTypeTag, Properties::FluidSystem>;
    using PressureFluidState = GetPropType<PressureTypeTag, Properties::FluidState>;
    PressureFluidState pressureFluidState;
    pressureFluidState.setTemperature(pressureProblem->spatialParams().temperatureAtPos(GlobalPosition{}));
    pressureFluidState.setPressure(0, referencePressure);
    const Scalar referenceViscosity = PressureFluidSystem::viscosity(pressureFluidState, 0);

    // Phase mobilities lambda_alpha = kr_alpha/mu_alpha and the fractional flow function
    //     f_w(S_w) = lambda_w/(lambda_w + lambda_n),
    // which gives the part of the total volume flux carried by the wetting phase: q_w = f_w q_t.
    // (Valid without gravity and capillary pressure, where both phases are driven by the same
    // pressure gradient.) The material laws are homogeneous, so they are evaluated once.
    const auto fluidMatrixInteraction = problem->spatialParams().fluidMatrixInteractionAtPos(GlobalPosition{});
    auto mobilityW = [&](Scalar s) { return fluidMatrixInteraction.krw(s)/viscosityW; };
    auto mobilityN = [&](Scalar s) { return fluidMatrixInteraction.krn(s)/viscosityN; };
    auto fractionalFlowW = [&](Scalar s) { return mobilityW(s)/(mobilityW(s) + mobilityN(s)); };

    // The saturation equation phi dS_w/dt + div(f_w(S_w) v_t) = 0 is hyperbolic: saturation
    // values travel with the characteristic speed f_w'(S_w) v_t/phi. The fastest possible speed,
    // determined by max f_w' over [S_wr, 1 - S_nr], enters the CFL condition below.
    // It is estimated here by sampling forward difference quotients (for the parameters of
    // params_impes.input: max f_w' = 5.29 at S_w = 0.53).
    const auto& effToAbs = fluidMatrixInteraction.pcSwCurve().effToAbsParams();
    const Scalar swMin = effToAbs.swr();
    const Scalar swMax = 1.0 - effToAbs.snr();
    Scalar maxDFractionalFlowW = 0.0;
    {
        constexpr int numSamples = 1000;
        const Scalar ds = (swMax - swMin)/numSamples;
        for (int k = 0; k < numSamples; ++k)
        {
            const Scalar s = swMin + k*ds;
            maxDFractionalFlowW = std::max(maxDFractionalFlowW, (fractionalFlowW(s + ds) - fractionalFlowW(s))/ds);
        }
    }

    // output
    VtkOutputModule<GridVariables, SolutionVector> vtkWriter(*gridVariables, sol, problem->name());
    using IOFields = GetPropType<TwoPTypeTag, Properties::IOFields>;
    using VelocityOutput = GetPropType<TwoPTypeTag, Properties::VelocityOutput>;
    vtkWriter.addVelocityOutput(std::make_shared<VelocityOutput>(*gridVariables));
    IOFields::initOutputModule(vtkWriter);

    BuckleyLeverettAnalyticSolution<TwoPTypeTag> analyticSolution(problem);
    vtkWriter.addField(analyticSolution.values(), "Sw_exact");
    vtkWriter.addField(pressureProblem->spatialParams().relativeTotalMobility(), "lambda_t*mu_ref");
    vtkWriter.write(0.0);

    // total volume flux over each scvf (positive if leaving the element)
    std::vector<Scalar> volumeFlux(gridGeometry->numScvf(), 0.0);

    // step 1 + 2: implicit pressure solve and computation of the total volume fluxes
    auto solvePressure = [&]()
    {
        // (a) Freeze the coefficients: evaluate the total mobility with the OLD saturation S_w^n.
        //     This is what makes the scheme "IMplicit Pressure, Explicit Saturation": the pressure
        //     equation becomes linear in p (no Newton iterations over p and S_w as in main.cc).
        //     Per cell i, the effective permeability is K_eff,i = K lambda_t(S_w,i^n) mu_ref.
        std::vector<Scalar> relativeTotalMobility(sw.size());
        for (std::size_t i = 0; i < sw.size(); ++i)
            relativeTotalMobility[i] = (mobilityW(sw[i]) + mobilityN(sw[i]))*referenceViscosity;
        pressureProblem->spatialParams().setRelativeTotalMobility(relativeTotalMobility);

        // (b) Recompute the TPFA transmissibilities with the new K_eff. For a face sigma between
        //     cells i and j with cell-center-to-face distances d_i, d_j and face area |sigma|:
        //         T_sigma = |sigma| / (d_i/K_eff,i + d_j/K_eff,j)     (harmonic averaging)
        //     and on a Dirichlet boundary face: T_sigma = |sigma| K_eff,i / d_i.
        //     Note that lambda_t is thus averaged harmonically, not upwinded.
        //     The transmissibilities are cached and do not depend on p, so GridVariables::update
        //     would not recompute them; init forces the recomputation.
        pressureGridVariables->init(p);

        // (c) Assemble and solve the discrete system. The residual of cell i is
        //         R_i(p) = sum_{sigma in interior faces} rho/mu_ref T_sigma (p_i - p_j)
        //                + sum_{sigma in Dirichlet faces} rho/mu_ref T_sigma (p_i - p_inj)
        //                + sum_{sigma in Neumann faces}   |sigma| (prescribed outward mass flux),
        //     i.e. the net mass outflow of cell i, which has to vanish. R is affine in p:
        //         R(p) = A p - b,   with   A = dR/dp,
        //     where A contains the scaled transmissibilities and b the boundary data.
        //     DuMux solves nonlinear problems with Newton's method: A deltaP = R(p), p <- p - deltaP.
        //     For an affine residual, a single Newton step yields the exact solution from any
        //     starting value (here, p = 0 in the first step, then the previous pressure):
        //         p - A^{-1}(A p - b) = A^{-1} b.
        //     (up to the tolerance of the iterative linear solver)
        pressureAssembler->assembleJacobianAndResidual(p);
        PressureSolutionVector deltaP(p.size());
        deltaP = 0.0;
        pressureLinearSolver->solve(*A, deltaP, *r);
        p -= deltaP;
        pressureGridVariables->update(p);

        // (d) Post-process the pressure into total volume fluxes per face (in m^3/s, positive if
        //     leaving the element). These fluxes are conservative: the flux computed from cell i
        //     equals minus the flux computed from cell j. The explicit saturation update below uses
        //     exactly these fluxes, which guarantees mass conservation of the transport step.
        //     Interior and Dirichlet faces: Darcy flux q = T_sigma (p_i - p_j) * (1/mu_ref), where
        //     the "upwind term" 1/mu_ref turns the transmissibility into a volume flux
        //     (K_eff/mu_ref = K lambda_t). Neumann faces: the prescribed mass flux divided by rho.
        using FluxVariables = GetPropType<PressureTypeTag, Properties::FluxVariables>;
        auto upwindTerm = [](const auto& volVars) { return volVars.mobility(0); };
        auto fvGeometry = localView(*gridGeometry);
        auto elemVolVars = localView(pressureGridVariables->curGridVolVars());
        auto elemFluxVarsCache = localView(pressureGridVariables->gridFluxVarsCache());
        for (const auto& element : elements(leafGridView))
        {
            fvGeometry.bind(element);
            elemVolVars.bind(element, fvGeometry, p);
            elemFluxVarsCache.bind(element, fvGeometry, elemVolVars);

            for (const auto& scvf : scvfs(fvGeometry))
            {
                if (scvf.boundary() && pressureProblem->boundaryTypes(element, scvf).hasNeumann())
                {
                    // prescribed mass flux -> volume flux
                    const auto& insideVolVars = elemVolVars[scvf.insideScvIdx()];
                    const auto massFlux = pressureProblem->neumannAtPos(scvf.ipGlobal())[0];
                    volumeFlux[scvf.index()] = massFlux/insideVolVars.density(0)*scvf.area();
                }
                else
                {
                    FluxVariables fluxVars;
                    fluxVars.init(*pressureProblem, element, fvGeometry, elemVolVars, scvf, elemFluxVarsCache);
                    volumeFlux[scvf.index()] = fluxVars.advectiveFlux(0, upwindTerm);
                }
            }
        }
    };

    // CFL condition (Courant-Friedrichs-Lewy)
    //
    // Consider the 1D update of cell i with inflow from cell i-1 and constant flux q > 0
    // (see updateSaturation below):
    //     S_i^{n+1} = S_i^n - dt q/(phi |V_i|) (f_w(S_i^n) - f_w(S_{i-1}^n)).
    // With the mean value theorem, f_w(S_i) - f_w(S_{i-1}) = f_w'(xi) (S_i - S_{i-1}), so
    //     S_i^{n+1} = (1 - c) S_i^n + c S_{i-1}^n,   c = dt q f_w'(xi)/(phi |V_i|) (Courant number).
    // For 0 <= c <= 1, the new value is a convex combination of old values. Then the scheme is
    // monotone: S_w stays within [S_wr, 1 - S_nr] and no oscillations occur at the front.
    // For c > 1, errors are amplified and the solution becomes unstable.
    // Since xi is unknown, the bound uses max f_w' and the total outflow of the cell:
    //     dt <= phi |V_i| / (max f_w' * sum_{outflow faces} q_sigma)   for all cells i.
    // Physically: in one time step, no saturation value may travel further than one cell.
    //
    // Remarks for experiments (parameter Impes.CFLFactor; L1 error = integral of |S_w - S_w,exact|
    // at t_end, 100 cells as in params_impes.input):
    //  - Within the guaranteed range CFLFactor <= 1, explicit upwinding adds numerical diffusion.
    //    For linear advection with speed v it is D_num = v dx/2 (1 - c), i.e. it DEcreases for
    //    larger Courant numbers c: L1 error 0.42 m (CFLFactor 0.5) vs. 0.31 m (0.9).
    //    For implicit Euler with upwinding, D_num = v dx/2 (1 + c) INcreases with c, which is why
    //    the fully implicit scheme with its large time steps (main.cc) smears the front more (1.18 m).
    //  - The bound is not sharp: the exact solution jumps over the saturations around S_w = 0.53
    //    where f_w' is maximal (the shock speed is (f_w(S_shock) - f_w(S_wr))/(S_shock - S_wr)
    //    = 2.11 < 5.29); only the few smeared cells at the front take such values. Therefore
    //    CFLFactor slightly above 1 does not blow up immediately.
    //  - Beyond 1, however, monotonicity is lost and the results become unreliable: L1 error
    //    0.44 m (1.5), 0.73 m (1.6), 2.06 m (2.0), where a wrong saturation plateau forms behind
    //    the front. The error does not grow monotonically (e.g. 0.21 m at 2.5 by chance), so a
    //    single good result above 1 proves nothing. With 3.0 the error check fails.
    auto cflTimeStepSize = [&]()
    {
        Scalar dtCfl = std::numeric_limits<Scalar>::max();
        auto fvGeometry = localView(*gridGeometry);
        for (const auto& element : elements(leafGridView))
        {
            fvGeometry.bindElement(element);
            Scalar outflow = 0.0;
            for (const auto& scvf : scvfs(fvGeometry))
                outflow += std::max(volumeFlux[scvf.index()], 0.0);

            if (outflow > 0.0)
            {
                const auto& scv = fvGeometry.scv(fvGeometry.gridGeometry().elementMapper().index(element));
                const auto porosity = problem->spatialParams().porosityAtPos(scv.center());
                dtCfl = std::min(dtCfl, porosity*scv.volume()/(maxDFractionalFlowW*outflow));
            }
        }
        return cflFactor*dtCfl;
    };

    // Step 3: explicit upwind saturation update
    //
    // Finite-volume discretization of phi dS_w/dt + div(f_w(S_w) v_t) = 0 on cell i, integrated
    // over [t^n, t^{n+1}] with the explicit (forward) Euler method:
    //     phi |V_i| (S_i^{n+1} - S_i^n)/dt + sum_sigma q_sigma f_w(S_up(sigma)^n) = 0.
    //  - Explicit: all fluxes use the OLD saturation S^n. Each cell is updated independently;
    //    no system of equations has to be solved (in contrast to the pressure step).
    //    Hence the new values are written to a copy (swNew); updating sw in place would mix
    //    old and new values and make the result depend on the order of the cells.
    //  - Upwind: f_w is evaluated in the cell the total flux comes from. This is the correct
    //    choice here because f_w' >= 0 and both phases flow in the direction of v_t (no gravity,
    //    no capillarity). With gravity, the phases may flow in opposite directions
    //    (counter-current flow) and each phase flux has to be upwinded separately.
    //    Central differences instead of upwinding would produce oscillations at the front.
    //  - Boundary faces: outflow faces use the cell's own saturation; inflow faces use the
    //    Dirichlet saturation of the two-phase problem (S_w = 1 - S_nr, i.e. f_w = 1, pure
    //    wetting phase is injected).
    //  - Mass conservation: q_sigma f_w(S_up) is the same value seen from both cells adjacent to
    //    sigma (with opposite signs), so the sum over all cells cancels in the interior and the
    //    total wetting-phase mass only changes through the boundary fluxes. This holds because
    //    q_sigma itself is conservative (see step (d) of the pressure solve).
    // The scheme is first-order accurate in space and time.
    auto updateSaturation = [&](Scalar dt)
    {
        std::vector<Scalar> swNew(sw);
        auto fvGeometry = localView(*gridGeometry);
        for (const auto& element : elements(leafGridView))
        {
            fvGeometry.bindElement(element);
            const auto eIdx = gridGeometry->elementMapper().index(element);

            Scalar fluxW = 0.0;
            for (const auto& scvf : scvfs(fvGeometry))
            {
                const auto q = volumeFlux[scvf.index()];
                Scalar swUpwind = sw[eIdx];
                if (q < 0.0)
                {
                    if (!scvf.boundary())
                        swUpwind = sw[scvf.outsideScvIdx()];
                    else if (problem->boundaryTypes(element, scvf).isDirichlet(saturationIdx))
                        swUpwind = 1.0 - problem->dirichletAtPos(scvf.ipGlobal())[saturationIdx];
                    else
                        DUNE_THROW(Dune::NotImplemented, "Inflow over Neumann boundary at " << scvf.ipGlobal());
                }
                fluxW += q*fractionalFlowW(swUpwind);
            }

            const auto& scv = fvGeometry.scv(eIdx);
            const auto porosity = problem->spatialParams().porosityAtPos(scv.center());
            swNew[eIdx] = sw[eIdx] - dt/(porosity*scv.volume())*fluxW;
        }
        sw = std::move(swNew);
    };

    auto timeLoop = std::make_shared<TimeLoop<Scalar>>(0.0, dt, tEnd);
    timeLoop->setMaxTimeStepSize(maxDt);

    // IMPES time loop: pressure (implicit) -> fluxes -> CFL time step -> saturation (explicit)
    timeLoop->start();
    do {
        // The pressure (and with it q_t) is updated only every n-th step (Impes.PressureUpdateInterval).
        // Re-solving the elliptic pressure equation is the expensive part; the total velocity
        // usually changes much more slowly than the saturation, so it may be reused for several
        // cheap transport steps. In this 1D setup, q_t is fixed by the Neumann boundary condition
        // and does not change at all, so the interval does not affect the result.
        if (timeLoop->timeStepIndex() % pressureUpdateInterval == 0)
            solvePressure();

        // The time step size is not chosen by accuracy or Newton convergence (as in main.cc) but
        // by stability. The time loop additionally limits it to MaxTimeStepSize and the remaining time.
        timeLoop->setTimeStepSize(cflTimeStepSize());
        updateSaturation(timeLoop->timeStepSize());

        // copy the IMPES solution into the two-phase solution vector for output
        for (std::size_t i = 0; i < sol.size(); ++i)
        {
            sol[i][pressureIdx] = p[i][0];
            sol[i][saturationIdx] = 1.0 - sw[i];
        }
        gridVariables->update(sol);

        timeLoop->advanceTimeStep();

        analyticSolution.update(timeLoop->time());
        vtkWriter.write(timeLoop->time());

        timeLoop->reportTimeStep();
    } while (!timeLoop->finished());

    timeLoop->finalize(leafGridView.comm());

    // compute relative error in wetting-phase center of mass and total mass
    // assumptions: solutions are (pseudo)1D, densities and porosity are constant, TPFA discretization
    const auto densityW = FluidSystem::density(fluidState, FluidSystem::phase0Idx);
    const auto porosity = problem->spatialParams().porosityAtPos(GlobalPosition{});
    const auto xMin = gridGeometry->bBoxMin()[0];
    const auto xMax = gridGeometry->bBoxMax()[0];
    const auto domainWidth = gridGeometry->bBoxMax()[1] - gridGeometry->bBoxMin()[1];

    // analytic computation of center of mass and total mass
    auto swFunc = [&analyticSolution, &tEnd](Scalar x)
    {
        return analyticSolution.computeSaturation(x, tEnd);
    };

    auto firstMomentOfMassWIntegrand = [&densityW, &porosity, &domainWidth, &swFunc](Scalar x)
    {
        return x*densityW*swFunc(x)*porosity*domainWidth;
    };
    Scalar firstMomentOfMassWAnalytic = integrateScalarFunction(firstMomentOfMassWIntegrand, xMin, xMax);

    auto totalMassWIntegrand = [&densityW, &porosity, &domainWidth, &swFunc](Scalar x)
    {
        return densityW*swFunc(x)*porosity*domainWidth;
    };
    Scalar totalMassWAnalytic = integrateScalarFunction(totalMassWIntegrand, xMin, xMax);

    // numeric computation of center of mass and total mass
    Scalar firstMomentOfMassWNumeric = 0.0;
    Scalar totalMassWNumeric = 0.0;
    for (const auto& element : elements(gridGeometry->gridView(), Dune::Partitions::interior))
    {
        const auto globalPos = element.geometry().center();
        const auto eIdx = gridGeometry->elementMapper().index(element);
        const auto volume = element.geometry().volume();

        Scalar localMassWNumeric = densityW * porosity * volume * sw[eIdx];

        firstMomentOfMassWNumeric += globalPos[0] * localMassWNumeric;
        totalMassWNumeric += localMassWNumeric;
    }

    Scalar centerOfMassNumeric = firstMomentOfMassWNumeric/totalMassWNumeric;
    Scalar centerOfMassAnalytic = firstMomentOfMassWAnalytic/totalMassWAnalytic;

    // compute relative errors
    Scalar distanceCenterOfMass = std::abs(centerOfMassAnalytic - centerOfMassNumeric);
    Scalar relErrorCenterOfMass = distanceCenterOfMass/centerOfMassAnalytic;
    Scalar differenceTotalMass = std::abs(totalMassWAnalytic - totalMassWNumeric);
    Scalar relErrorTotalMass = differenceTotalMass/totalMassWAnalytic;
    const auto maxError = getParam<Scalar>("Problem.MaxRelError");
    const bool centerOfMassFailed = relErrorCenterOfMass > maxError;
    const bool totalMassFailed = relErrorTotalMass > maxError;
    if (centerOfMassFailed || totalMassFailed)
    {
        std::ostringstream msg;
        msg << "Analytical Buckley-Leverett check failed.";

        if (centerOfMassFailed)
            msg << " Relative wetting-phase center of mass error "
                << relErrorCenterOfMass << " exceeds the threshold " << maxError << ".";

        if (totalMassFailed)
            msg << " Relative total wetting-phase mass error "
                << relErrorTotalMass << " exceeds the threshold " << maxError << ".";

        DUNE_THROW(Dune::InvalidStateException, msg.str());
    }

    std::cout << "max dfw/dSw: " << maxDFractionalFlowW << std::endl;
    std::cout << "numeric center of mass: " << centerOfMassNumeric << std::endl;
    std::cout << "analytic center of mass: " << centerOfMassAnalytic << std::endl;
    std::cout << "integrated numeric mass: " << totalMassWNumeric << std::endl;
    std::cout << "integrated analytics mass: " << totalMassWAnalytic << std::endl;
    Parameters::print();

    return 0;
}
