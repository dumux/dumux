// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup NavierStokesTests
 * \brief P2-P1-P2-P2 Taylor-Hood Cahn-Hilliard / Navier-Stokes rising-bubble test (static mesh).
 */
#include <config.h>

#include <iostream>
#include <vector>
#include <array>
#include <cmath>
#include <optional>

#include <dune/common/parallel/mpihelper.hh>
#include <dune/common/timer.hh>
#include <dune/common/fvector.hh>

#include <dumux/common/initialize.hh>
#include <dumux/common/dumuxmessage.hh>
#include <dumux/common/parameters.hh>
#include <dumux/common/properties.hh>

#include <dumux/io/vtkoutputmodule.hh>
#include <dumux/io/grid/gridmanager_alu.hh>
#include <dumux/io/gridwriter.hh>

#include <dumux/linear/istlsolvers.hh>
#include <dumux/linear/linearsolvertraits.hh>
#include <dumux/linear/linearalgebratraits.hh>

#include <dumux/multidomain/fvassembler.hh>
#include <dumux/multidomain/assembler.hh>
#include <dumux/multidomain/traits.hh>
#include <dumux/multidomain/newtonsolver.hh>
#include <dumux/discretization/cvfe/interpolate.hh>

#include <dumux/freeflow/navierstokes/momentum/velocityoutput.hh>
#include <dumux/discretization/elementsolution.hh>
#include <dumux/discretization/evalsolution.hh>
#include <dune/geometry/quadraturerules.hh>
#include <dune/geometry/referenceelements.hh>

#include <dumux/adaptive/adapt.hh>
#include <dumux/adaptive/markelements.hh>

#include "properties.hh"
// Reuses the P1 baseline's AMR machinery (indicator + data transfer classes), which is
// basis-agnostic / templated on TypeTag and PrimaryVariables size: the phase-field jump
// indicator only calls elementSolution/evalSolution, and the (non-vertex-dof) element-local
// data transfer already handles an arbitrary number of non-CV local dofs per element (written
// for the PQ1Bubble centroid dof, it generalizes unmodified to the PQ2 edge dofs here). The
// pressure-only P1 mass transfer degenerates to plain interpolation automatically because its
// "mass-projected component range" [1, PrimaryVariables::size()) is empty for a 1-component
// (pressure-only) primary variable vector.
#include "../adapt.hh"

namespace {

//! Dune-functions-compatible wrapper extracting a single scalar component of a
//! CVFE solution, so GridFormat writes it as a 1-component (scalar) point field
//! evaluated at the higher-order Lagrange nodes (rather than a padded 3-vector).
template<class GridGeometry, class SolutionVector>
class CVFEScalarComponent
{
    using GridView = typename GridGeometry::GridView;
    using Element = typename GridView::template Codim<0>::Entity;
public:
    class LocalFunction
    {
    public:
        LocalFunction(const GridGeometry& gg, const SolutionVector& sol, int comp)
        : gg_(&gg), sol_(&sol), comp_(comp) {}
        void bind(const Element& element) { element_.emplace(element); }
        void unbind() { element_.reset(); }
        auto operator()(const typename Element::Geometry::LocalCoordinate& localPos) const
        {
            const auto elemSol = Dumux::elementSolution(*element_, *sol_, *gg_);
            return Dumux::evalSolutionAtLocalPos(
                *element_, element_->geometry(), *gg_, elemSol, localPos
            )[comp_];
        }
    private:
        const GridGeometry* gg_; const SolutionVector* sol_; int comp_;
        std::optional<Element> element_;
    };

    CVFEScalarComponent(const GridGeometry& gg, const SolutionVector& sol, int comp)
    : gg_(&gg), sol_(&sol), comp_(comp) {}

    friend LocalFunction localFunction(const CVFEScalarComponent& f)
    { return LocalFunction(*f.gg_, *f.sol_, f.comp_); }
private:
    const GridGeometry* gg_; const SolutionVector* sol_; int comp_;
};

} // end anonymous namespace

int main(int argc, char** argv)
{
    using namespace Dumux;

    using MomentumTypeTag = Properties::TTag::TYPETAG_MOMENTUM;
    using MassTypeTag = Properties::TTag::TYPETAG_MASS;

    initialize(argc, argv);
    const auto& mpiHelper = Dune::MPIHelper::instance();
    if (mpiHelper.rank() == 0)
        DumuxMessage::print(/*firstCall=*/true);

    Parameters::init(argc, argv);
    Dune::Timer timer;

    // grid (static; global refinement via Grid.Refinement)
    using Grid = GetPropType<MassTypeTag, Properties::Grid>;
    GridManager<Grid> gridManager;
    gridManager.init();
    const auto& leafGridView = gridManager.grid().leafGridView();

    using MomentumGridGeometry = GetPropType<MomentumTypeTag, Properties::GridGeometry>;
    auto momentumGridGeometry = std::make_shared<MomentumGridGeometry>(leafGridView);
    using MassGridGeometry = GetPropType<MassTypeTag, Properties::GridGeometry>;
    auto massGridGeometry = std::make_shared<MassGridGeometry>(leafGridView);

    using CouplingManager = GetPropType<MomentumTypeTag, Properties::CouplingManager>;
    auto couplingManager = std::make_shared<CouplingManager>();

    using MomentumProblem = GetPropType<MomentumTypeTag, Properties::Problem>;
    auto momentumProblem = std::make_shared<MomentumProblem>(momentumGridGeometry, couplingManager);
    using MassProblem = GetPropType<MassTypeTag, Properties::Problem>;
    auto massProblem = std::make_shared<MassProblem>(massGridGeometry, couplingManager);

    using Traits = MultiDomainTraits<MomentumTypeTag, MassTypeTag>;
    using Scalar = typename Traits::Scalar;
    static constexpr int dim = MomentumGridGeometry::GridView::dimension;

    const auto tEnd = getParam<Scalar>("TimeLoop.TEnd");
    const auto maxDt = getParam<Scalar>("TimeLoop.MaxTimeStepSize");
    auto dt = getParam<Scalar>("TimeLoop.DtInitial");
    auto timeLoop = std::make_shared<TimeLoop<Scalar>>(0.0, dt, tEnd);
    timeLoop->setMaxTimeStepSize(maxDt);

    constexpr auto momentumIdx = CouplingManager::freeFlowMomentumIndex;
    constexpr auto massIdx = CouplingManager::freeFlowMassIndex;
    using SolutionVector = typename Traits::SolutionVector;
    SolutionVector x;
    x[momentumIdx].resize(momentumGridGeometry->numDofs());
    x[massIdx].resize(massGridGeometry->numDofs());
    std::cout << "Total number of dofs: momentum(P2)=" << momentumGridGeometry->numDofs()
              << " x " << (dim+2) << " + mass(P1)=" << massGridGeometry->numDofs() << std::endl;

    // initial condition: nodal interpolation at ALL dof positions (vertices + P2 edge dofs).
    // Nodal-basis (DofPositionEvaluation) mode - no L2 projection needed for the PQ2/P1 Lagrange bases.
    interpolate(*momentumGridGeometry, x[momentumIdx],
                [&](const auto& pos){ return momentumProblem->initialAtPos(pos); });
    interpolate(*massGridGeometry, x[massIdx],
                [&](const auto& pos){ return massProblem->initialAtPos(pos); });
    auto xOld = x;

    using MomentumGridVariables = GetPropType<MomentumTypeTag, Properties::GridVariables>;
    auto momentumGridVariables = std::make_shared<MomentumGridVariables>(momentumProblem, momentumGridGeometry);
    using MassGridVariables = GetPropType<MassTypeTag, Properties::GridVariables>;
    auto massGridVariables = std::make_shared<MassGridVariables>(massProblem, massGridGeometry);

    couplingManager->init(momentumProblem, massProblem, std::make_tuple(momentumGridVariables, massGridVariables), x);
    massGridVariables->init(x[massIdx]);
    momentumGridVariables->init(x[momentumIdx]);

    /////////////////////////////////////////////////////////
    // adaptive grid refinement (P2 momentum+CH domain drives the indicator: phi lives there,
    // not on the pressure-only P1 mass domain). The mass-side transfer is the existing box-P1
    // one (degenerates to plain interpolation for a pressure-only, 1-component primary variable
    // vector); the momentum-side transfer is the existing "velocity" element-local transfer,
    // which already handles an arbitrary number of non-CV (here: P2 edge) local dofs per element
    // and, with PrimaryVariables::size()==dim+2, carries phi/mu along with velocity generically.
    /////////////////////////////////////////////////////////
    const Scalar refineTol = getParam<Scalar>("Adaptive.RefineTolerance");
    const Scalar coarsenTol = getParam<Scalar>("Adaptive.CoarsenTolerance");
    TwoPhaseCahnHilliardGridAdaptIndicator<MomentumTypeTag> indicator(momentumGridGeometry);
    using GridDataTransfer = TwoDomainOneGridDataTransfer<Grid, TwoPhaseCahnHilliardMassGridDataTransfer<MassTypeTag>, TwoPhaseCahnHilliardVelocityGridDataTransfer<MomentumTypeTag>>;
    GridDataTransfer dataTransfer(
        std::make_shared<TwoPhaseCahnHilliardMassGridDataTransfer<MassTypeTag>>(massProblem, massGridGeometry, massGridVariables, x[massIdx]),
        std::make_shared<TwoPhaseCahnHilliardVelocityGridDataTransfer<MomentumTypeTag>>(momentumProblem, momentumGridGeometry, momentumGridVariables, x[momentumIdx])
    );

    // Initial adaptation: grow the mesh incrementally from the (coarse) base grid, re-sampling
    // the analytic IC via `interpolate` (not the transfer) after every pass, so the indicator
    // always refines against the exact sharp interface and no projection smoothing accumulates.
    const std::size_t initMaxLevel = getParam<std::size_t>("Adaptive.InitMaxLevel", 10);
    for (std::size_t i = 0; i < initMaxLevel; ++i)
    {
        indicator.calculate(x[momentumIdx], refineTol, coarsenTol);
        if (!markElements(gridManager.grid(), indicator))
            break;
        if (!adapt(gridManager.grid(), dataTransfer))
            break;

        interpolate(*momentumGridGeometry, x[momentumIdx],
                    [&](const auto& pos){ return momentumProblem->initialAtPos(pos); });
        interpolate(*massGridGeometry, x[massIdx],
                    [&](const auto& pos){ return massProblem->initialAtPos(pos); });
        xOld = x;
        couplingManager->updateSolution(x);
        massGridVariables->updateAfterGridAdaption(x[massIdx]);
        momentumGridVariables->updateAfterGridAdaption(x[momentumIdx]);
    }
    std::cout << "\033[1;34m[adapt] after init adaptation: leafCells=" << gridManager.grid().leafGridView().size(0)
              << "  maxLevel=" << gridManager.grid().maxLevel()
              << "  massDofs=" << massGridGeometry->numDofs()
              << "  momentumDofs=" << momentumGridGeometry->numDofs() << "\033[0m" << std::endl;

    // reinitialize coupling manager & grid variables after grid adaption
    couplingManager->init(momentumProblem, massProblem, std::make_tuple(momentumGridVariables, massGridVariables), x);
    massGridVariables->updateAfterGridAdaption(x[massIdx]);
    momentumGridVariables->updateAfterGridAdaption(x[momentumIdx]);

    /////////////////////////////////////////////////////////

    // vtk: pressure (mass) + velocity + phase field (interpolated at the shared vertices)
    using IOFields = GetPropType<MassTypeTag, Properties::IOFields>;
    VtkOutputModule vtkWriter(*massGridVariables, x[massIdx], massProblem->name());
    IOFields::initOutputModule(vtkWriter);
    vtkWriter.addVelocityOutput(std::make_shared<NavierStokesVelocityOutput<MassGridVariables>>());
    std::vector<Scalar> phiField(massGridGeometry->numDofs(), 0.0);
    std::vector<Scalar> muField(massGridGeometry->numDofs(), 0.0);
    auto updateScalarFields = [&]{
        // the first numVertices momentum (P2) dofs are the grid vertices == mass (P1) dofs.
        // Resize every call: AMR changes the dof count between remeshes.
        phiField.resize(massGridGeometry->numDofs());
        muField.resize(massGridGeometry->numDofs());
        for (std::size_t i = 0; i < phiField.size(); ++i)
        {
            phiField[i] = x[momentumIdx][i][dim];
            muField[i]  = x[momentumIdx][i][dim+1];
        }
    };
    updateScalarFields();
    vtkWriter.addField(phiField, "phi");
    vtkWriter.addField(muField, "mu");
    vtkWriter.write(0.0);

    // Higher-order (P2) output via GridFormat. The phase field, chemical potential
    // and velocity live on the P2 momentum basis; writing them only at the P1
    // vertices (above) drops all edge DOFs and aliases the interface -> garbage.
    // Here we write true order-2 Lagrange VTK so ParaView shows the full P2 solution.
#if DUMUX_HAVE_GRIDFORMAT
    std::vector<Dune::FieldVector<Scalar, dim>> velHO(momentumGridGeometry->numDofs());
    std::vector<Scalar> phiHO(momentumGridGeometry->numDofs(), 0.0), muHO(momentumGridGeometry->numDofs(), 0.0);
    auto updateHOFields = [&]{
        // P2 momentum solution is DOF-indexed identically to the order-2 Lagrange grid.
        // Resize every call: AMR changes the dof count between remeshes.
        const auto numMomentumDofs = momentumGridGeometry->numDofs();
        velHO.resize(numMomentumDofs);
        phiHO.resize(numMomentumDofs);
        muHO.resize(numMomentumDofs);
        for (std::size_t i = 0; i < numMomentumDofs; ++i)
        {
            for (int d = 0; d < dim; ++d) velHO[i][d] = x[momentumIdx][i][d];
            phiHO[i] = x[momentumIdx][i][dim];     // phaseFieldIdx
            muHO[i]  = x[momentumIdx][i][dim+1];   // chemicalPotentialIdx
        }
    };
    updateHOFields();
    // Use base64-inlined, uncompressed VTU: the widely-compatible encoding that every
    // ParaView version reads reliably (the gridformat default is LZ4-compressed appended
    // data, which older/stricter readers can mis-decode into garbage-looking fields).
    IO::GridWriter hoWriter{
        IO::Format::pvd_with(IO::Format::vtu.with({
            .encoder = IO::Encoding::base64,
            .compressor = IO::Compression::none,
            .data_format = IO::VTK::DataFormat::inlined
        })), leafGridView,
        massProblem->name() + "_ho", IO::order<2>
    };
    // P2 momentum fields: index the DOF-indexed containers directly (no FE eval)
    hoWriter.setPointField("velocity", velHO);
    hoWriter.setPointField("phi", phiHO);
    hoWriter.setPointField("mu", muHO);
    // P1 pressure sampled at the order-2 Lagrange nodes via FE interpolation
    hoWriter.setPointField("p", CVFEScalarComponent(*massGridGeometry, x[massIdx], 0));
    hoWriter.write(0.0);
#endif

    using Assembler = Experimental::MultiDomainAssembler<Traits, CouplingManager, DiffMethod::numeric>;
    auto assembler = std::make_shared<Assembler>(std::make_tuple(momentumProblem, massProblem),
                                                 std::make_tuple(momentumGridGeometry, massGridGeometry),
                                                 std::make_tuple(momentumGridVariables, massGridVariables),
                                                 couplingManager, timeLoop, xOld);

    using LinearSolver = UMFPackIstlSolver<SeqLinearSolverTraits, LinearAlgebraTraitsFromAssembler<Assembler>>;
    auto linearSolver = std::make_shared<LinearSolver>();
    using NewtonSolver = MultiDomainNewtonSolver<Assembler, LinearSolver, CouplingManager>;
    NewtonSolver nonLinearSolver(assembler, linearSolver, couplingManager);

    // write visualization output every OutputInterval of simulated time (default: only the
    // final state). The per-step Hysing QoIs below are printed to stdout on every step regardless.
    const auto outputInterval = getParam<Scalar>("TimeLoop.OutputInterval", tEnd);
    Scalar nextOutputTime = outputInterval;

    timeLoop->start(); do
    {
        nonLinearSolver.solve(x, *timeLoop);
        xOld = x;
        momentumGridVariables->advanceTimeStep();
        massGridVariables->advanceTimeStep();
        timeLoop->advanceTimeStep();

        updateScalarFields();
        if (timeLoop->time() >= nextOutputTime - 1e-10 || timeLoop->finished())
        {
            vtkWriter.write(timeLoop->time());
#if DUMUX_HAVE_GRIDFORMAT
            updateHOFields();
            hoWriter.write(timeLoop->time());
#endif
            while (timeLoop->time() >= nextOutputTime - 1e-10)
                nextOutputTime += outputInterval;
        }
        timeLoop->reportTimeStep();
        timeLoop->setTimeStepSize(nonLinearSolver.suggestTimeStepSize(timeLoop->timeStepSize()));

        // adaptive mesh refinement: keep the fine mesh on the moving interface (see the P1
        // baseline's adapt.hh / main.cc for the rationale on AdaptEveryNSteps > 1).
        static const int adaptEvery = getParam<int>("Adaptive.AdaptEveryNSteps", 1);
        if (getParam<bool>("Adaptive.EnableInTimeStep", true)
            && (timeLoop->timeStepIndex() % adaptEvery == 0))
        {
            indicator.calculate(x[momentumIdx], refineTol, coarsenTol);
            if (markElements(gridManager.grid(), indicator)
                && adapt(gridManager.grid(), dataTransfer))
            {
                xOld = x;
                couplingManager->updateSolution(x);
                massGridVariables->updateAfterGridAdaption(x[massIdx]);
                momentumGridVariables->updateAfterGridAdaption(x[momentumIdx]);
                couplingManager->init(momentumProblem, massProblem,
                                      std::make_tuple(momentumGridVariables, massGridVariables), x);
                assembler->updateAfterGridAdaption();
#if DUMUX_HAVE_GRIDFORMAT
                // The GridFormat order-2 Lagrange grid wrapper caches its node/element topology at
                // construction and does not track subsequent grid adaptation; update() rebuilds
                // that topology from the (already-adapted) leafGridView. The registered point
                // fields stay valid: they close over gridGeometry_/x[.]/the raw containers by
                // reference, which are all updated in place, not by re-registering.
                hoWriter.update();
#endif
            }
        }

        // ---- Hysing benchmark QoIs on the co-located P2 momentum solution (phi = component dim) ----
        // Ported from the pq1bubble parent (problem.hh printMassBalanceSummary): marching-triangle
        // circularity/perimeter from the phi=0 contour, SHARP (phi>0) and diffuse ((1+phi)/2) rise
        // velocity and centroid, and the effective interface thickness eps_eff = sqrt(2) int|grad phi|
        // / (3 int|grad phi|^2). Printed in the parent's exact line format so postprocessing/
        // plot_hysing_benchmark.py parses this log unchanged. Mesh is the LEFT-HALF domain
        // x in [0,0.5] (symmetry at x=0.5), so full-bubble area & perimeter are 2x the half values.
        using GP = Dune::FieldVector<Scalar, dim>;
        Scalar riseNum = 0.0, riseDen = 0.0, riseNumSharp = 0.0, riseDenSharp = 0.0;
        Scalar bubbleVol = 0.0, bubbleYmoment = 0.0;      // diffuse (1+phi)/2 volume & y-moment
        Scalar areaHalfPos = 0.0, yMomentHalfPos = 0.0;   // sharp {phi>0} area & y-moment (marching tri)
        Scalar contourLenHalf = 0.0, tvHalf = 0.0, gradEnergyHalf = 0.0;
        {
            const auto refTri = Dune::referenceElement<Scalar, dim>(Dune::GeometryTypes::simplex(dim));
            for (const auto& element : elements(leafGridView))
            {
                const auto elemSol = elementSolution(element, x[momentumIdx], *momentumGridGeometry);
                const auto& geo = element.geometry();

                // Integrate phi and u_y over the element by P1 (vertex) interpolation. NOTE:
                // evalSolution() with the standard P2 Lagrange basis MIS-READS this hybrid CVFE grid's
                // edge dofs at non-nodal points (spurious O(100) phi/u values); the VERTEX dofs are
                // nodal and correct, so interpolate linearly from them. The benchmark QoI is P1-accurate
                // (the pq1bubble reference is P1), so dropping the P2 edge correction is acceptable.
                const auto& quad = Dune::QuadratureRules<Scalar, dim>::rule(geo.type(), 3);
                for (const auto& qp : quad)
                {
                    const auto& lp = qp.position();
                    const Scalar l0 = 1.0 - lp[0] - lp[1], l1 = lp[0], l2 = lp[1]; // barycentric (vertex 0,1,2)
                    const Scalar phi = l0*elemSol[0][dim] + l1*elemSol[1][dim] + l2*elemSol[2][dim];
                    const Scalar uy  = l0*elemSol[0][1]   + l1*elemSol[1][1]   + l2*elemSol[2][1];
                    const Scalar y = geo.global(lp)[1];
                    const Scalar w = qp.weight()*geo.integrationElement(lp);
                    const Scalar wb = 0.5*(1.0 + phi)*w;
                    riseNum += wb*uy; riseDen += wb; bubbleVol += wb; bubbleYmoment += wb*y;
                    if (phi > 0.0) { riseNumSharp += w*uy; riseDenSharp += w; }
                }

                // marching triangles on the 3 corners (phi at the P1 vertex dofs)
                if (geo.corners() == 3)
                {
                    std::array<Scalar, 3> phiC;
                    std::array<GP, 3> posC;
                    for (int i = 0; i < 3; ++i)
                    {
                        phiC[i] = elemSol[i][dim];   // P2 vertex dofs are local dofs 0..2
                        posC[i] = geo.corner(i);
                    }
                    auto triArea = [](const GP& a, const GP& b, const GP& c)
                    { return 0.5*std::abs((b[0]-a[0])*(c[1]-a[1]) - (c[0]-a[0])*(b[1]-a[1])); };
                    auto edgeCross = [&](int a, int b)
                    { const Scalar t = phiC[a]/(phiC[a]-phiC[b]); GP p = posC[a]; p.axpy(t, posC[b]-posC[a]); return p; };
                    const int npos = int(phiC[0]>0.0) + int(phiC[1]>0.0) + int(phiC[2]>0.0);
                    const Scalar fullA = triArea(posC[0], posC[1], posC[2]);
                    const Scalar ycFull = (posC[0][1]+posC[1][1]+posC[2][1])/3.0;
                    if (npos == 3) { areaHalfPos += fullA; yMomentHalfPos += fullA*ycFull; }
                    else if (npos == 1)
                    {
                        const int i = phiC[0]>0.0 ? 0 : (phiC[1]>0.0 ? 1 : 2);
                        const auto Xij = edgeCross(i, (i+1)%3), Xik = edgeCross(i, (i+2)%3);
                        const Scalar aT = triArea(posC[i], Xij, Xik);
                        areaHalfPos += aT; yMomentHalfPos += aT*(posC[i][1]+Xij[1]+Xik[1])/3.0;
                        contourLenHalf += (Xik - Xij).two_norm();
                    }
                    else if (npos == 2)
                    {
                        const int k = !(phiC[0]>0.0) ? 0 : (!(phiC[1]>0.0) ? 1 : 2);
                        const auto Xki = edgeCross(k, (k+1)%3), Xkj = edgeCross(k, (k+2)%3);
                        const Scalar aS = triArea(posC[k], Xki, Xkj);
                        areaHalfPos += fullA - aS;
                        yMomentHalfPos += fullA*ycFull - aS*(posC[k][1]+Xki[1]+Xkj[1])/3.0;
                        contourLenHalf += (Xkj - Xki).two_norm();
                    }
                    const Scalar twoA = (posC[1][0]-posC[0][0])*(posC[2][1]-posC[0][1])
                                      - (posC[2][0]-posC[0][0])*(posC[1][1]-posC[0][1]);
                    if (std::abs(twoA) > 1e-30)
                    {
                        const Scalar dphidx = ((phiC[1]-phiC[0])*(posC[2][1]-posC[0][1])
                                             - (phiC[2]-phiC[0])*(posC[1][1]-posC[0][1]))/twoA;
                        const Scalar dphidy = ((phiC[2]-phiC[0])*(posC[1][0]-posC[0][0])
                                             - (phiC[1]-phiC[0])*(posC[2][0]-posC[0][0]))/twoA;
                        tvHalf += std::hypot(dphidx, dphidy)*fullA;
                        gradEnergyHalf += (dphidx*dphidx + dphidy*dphidy)*fullA;
                    }
                }
            }
        }
        const Scalar areaFull = 2.0*areaHalfPos, perimeterFull = 2.0*contourLenHalf;
        const Scalar epsEff = gradEnergyHalf > 0.0 ? std::sqrt(2.0)*tvHalf/(3.0*gradEnergyHalf) : 0.0;
        const Scalar circularity = perimeterFull > 0.0 ? 2.0*std::sqrt(M_PI*areaFull)/perimeterFull : 0.0;

        Scalar maxV = 0.0, maxVface = 0.0, maxVbulk = 0.0;
        for (const auto& dof : x[momentumIdx]) // dof = [vx, vy, phi, mu]; classify by co-located phi
        {
            const Scalar sp = std::hypot(dof[0], dof[1]);
            maxV = std::max(maxV, sp);
            if (std::abs(dof[dim]) < 0.9) maxVface = std::max(maxVface, sp); else maxVbulk = std::max(maxVbulk, sp);
        }
        Scalar maxPhi = 0.0; for (const auto& p : phiField) maxPhi = std::max(maxPhi, std::abs(p));
        Scalar muMin = 1e30, muMax = -1e30; // mu = momentum component dim+1
        for (const auto& d : x[momentumIdx]) { muMin = std::min(muMin, d[dim+1]); muMax = std::max(muMax, d[dim+1]); }
        Scalar pMin = 1e30, pMax = -1e30;    // pressure = mass component 0
        for (const auto& d : x[massIdx]) { pMin = std::min(pMin, d[0]); pMax = std::max(pMax, d[0]); }
        std::cout << "\033[1;36m[fields] mu=[" << muMin << "," << muMax << "]  p=[" << pMin << "," << pMax
                  << "]  (Laplace sigma/R=" << getParam<Scalar>("Problem.SurfaceTension",24.5)/getParam<Scalar>("Problem.BubbleRadius",0.25)
                  << ")\033[0m" << std::endl;

        // parent-format QoI lines (parsed by postprocessing/plot_hysing_benchmark.py). NOTE: the
        // parser starts a record on the "time:" line and FLUSHES it when rise_velocity matches, so
        // rise_velocity MUST be printed LAST (after max|v|, centroid_y, circularity, eps_eff).
        std::cout << "\033[1;33m[vel] max|v| = " << maxV << " m/s  (face=" << maxVface
                  << " bulk=" << maxVbulk << " max|phi|=" << maxPhi << ")\033[0m" << std::endl;
        std::cout << "\033[1;35m[bubble] centroid_y = " << (bubbleVol>0.0 ? bubbleYmoment/bubbleVol : 0.0)
                  << " m  |  centroid_y_sharp = " << (areaHalfPos>0.0 ? yMomentHalfPos/areaHalfPos : 0.0)
                  << " m (benchmark def., phi>0)  |  volume = " << (2.0*bubbleVol) << " m^2/m\033[0m" << std::endl;
        std::cout << "\033[1;32m[bubble] circularity = " << circularity
                  << "  |  perimeter = " << perimeterFull << " m  |  area = " << areaFull << " m^2\033[0m" << std::endl;
        std::cout << "\033[1;33m[bubble] eps_eff = " << epsEff << " m  (IC eps = "
                  << getParam<Scalar>("Problem.InterfaceThickness", 0.005) << " m)\033[0m" << std::endl;
        std::cout << "\033[1;32m[bubble] rise_velocity = " << (riseDen>0.0 ? riseNum/riseDen : 0.0)
                  << " m/s  |  rise_velocity_sharp = " << (riseDenSharp>0.0 ? riseNumSharp/riseDenSharp : 0.0)
                  << " m/s (benchmark def., phi>0)\033[0m" << std::endl;
    } while (!timeLoop->finished());

    timeLoop->finalize(leafGridView.comm());
    if (mpiHelper.rank() == 0)
    {
        Parameters::print();
        DumuxMessage::print(/*firstCall=*/false);
    }
    return 0;
}
