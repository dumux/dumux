// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup EmbeddedTests
 * \brief Convergence test for the 1d-3d embedded coupling with a vessel that ends inside the tissue.
 */
#include <config.h>

#include <cmath>
#include <iomanip>
#include <iostream>

#include <dune/common/parallel/mpihelper.hh>

#include <dumux/common/initialize.hh>
#include <dumux/common/properties.hh>
#include <dumux/common/parameters.hh>
#include <dumux/common/dumuxmessage.hh>
#include <dumux/geometry/diameter.hh>
#include <dumux/linear/seqsolverbackend.hh>
#include <dumux/io/vtkoutputmodule.hh>
#include <dumux/io/grid/gridmanager_yasp.hh>
#include <dumux/io/grid/gridmanager_foam.hh>

#include <dumux/multidomain/traits.hh>
#include <dumux/multidomain/fvassembler.hh>

#include "properties.hh"

/*!
 * \brief Relative discrete \f$ L^2 \f$ error with the exact solution evaluated at the degrees of freedom
 *        and the control volumes as weights
 */
template<class Problem, class SolutionVector>
double relativeDiscreteL2Error(const Problem& problem, const SolutionVector& sol)
{
    double errorSquared = 0.0, normSquared = 0.0;
    const auto& exactSolution = problem.exactSolutionAtDofs();
    auto fvGeometry = localView(problem.gridGeometry());
    for (const auto& element : elements(problem.gridGeometry().gridView()))
    {
        fvGeometry.bindElement(element);
        for (const auto& scv : scvs(fvGeometry))
        {
            const double exact = exactSolution[scv.dofIndex()];
            const double error = sol[scv.dofIndex()][0] - exact;
            errorSquared += error*error*scv.volume();
            normSquared += exact*exact*scv.volume();
        }
    }
    return std::sqrt(errorSquared/normSquared);
}

template<class GridView>
double maxElementDiameter(const GridView& gridView)
{
    double hMax = 0.0;
    for (const auto& element : elements(gridView))
        hMax = std::max(hMax, Dumux::diameter(element.geometry()));
    return hMax;
}

int main(int argc, char** argv)
{
    using namespace Dumux;

    Dumux::initialize(argc, argv);
    const auto& mpiHelper = Dune::MPIHelper::instance();
    if (mpiHelper.rank() == 0)
        DumuxMessage::print(/*firstCall=*/true);

    Parameters::init(argc, argv);

    using BulkTypeTag = Properties::TTag::BULKTYPETAG;
    using LowDimTypeTag = Properties::TTag::LOWDIMTYPETAG;

    GridManager<GetPropType<BulkTypeTag, Properties::Grid>> bulkGridManager;
    bulkGridManager.init("Tissue");
    GridManager<GetPropType<LowDimTypeTag, Properties::Grid>> lowDimGridManager;
    lowDimGridManager.init("Vessel");

    const auto& bulkGridView = bulkGridManager.grid().leafGridView();
    const auto& lowDimGridView = lowDimGridManager.grid().leafGridView();

    using BulkGridGeometry = GetPropType<BulkTypeTag, Properties::GridGeometry>;
    auto bulkGridGeometry = std::make_shared<BulkGridGeometry>(bulkGridView);
    using LowDimGridGeometry = GetPropType<LowDimTypeTag, Properties::GridGeometry>;
    auto lowDimGridGeometry = std::make_shared<LowDimGridGeometry>(lowDimGridView);

    // the vessel starts on the symmetry plane z = 0 and ends at the tip z = L
    const auto& vesselMin = lowDimGridGeometry->bBoxMin();
    const auto& vesselMax = lowDimGridGeometry->bBoxMax();
    if (std::hypot(vesselMin[0], vesselMin[1]) > 1e-12 || std::hypot(vesselMax[0], vesselMax[1]) > 1e-12
        || std::abs(vesselMin[2]) > 1e-12 || std::abs(bulkGridGeometry->bBoxMin()[2]) > 1e-12)
        DUNE_THROW(Dune::InvalidStateException, "The vessel has to lie on the z-axis and start at z = 0 on the tissue boundary");

    using Scalar = GetPropType<BulkTypeTag, Properties::Scalar>;
    const Scalar radius = getParam<Scalar>("SpatialParams.Radius");
    auto exactSolution = std::make_shared<VesselTipExactSolution<Scalar>>(
        radius, getParam<Scalar>("MixedDimension.KernelWidthFactor")*radius, vesselMax[2],
        getParam<Scalar>("Problem.ExchangeRate"), getParam<Scalar>("SpatialParams.VesselConductivity"),
        getParam<Scalar>("Problem.TipPressure")
    );

    using Traits = MultiDomainTraits<BulkTypeTag, LowDimTypeTag>;
    constexpr auto bulkIdx = Traits::template SubDomain<0>::Index();
    constexpr auto lowDimIdx = Traits::template SubDomain<1>::Index();

    using CouplingManager = GetPropType<BulkTypeTag, Properties::CouplingManager>;
    auto couplingManager = std::make_shared<CouplingManager>(bulkGridGeometry, lowDimGridGeometry);

    using BulkProblem = GetPropType<BulkTypeTag, Properties::Problem>;
    auto bulkProblem = std::make_shared<BulkProblem>(bulkGridGeometry, couplingManager, exactSolution);
    using LowDimProblem = GetPropType<LowDimTypeTag, Properties::Problem>;
    auto lowDimProblem = std::make_shared<LowDimProblem>(lowDimGridGeometry, couplingManager, exactSolution);

    Traits::SolutionVector sol;
    sol[bulkIdx].resize(bulkGridGeometry->numDofs());
    sol[lowDimIdx].resize(lowDimGridGeometry->numDofs());
    bulkProblem->applyInitialSolution(sol[bulkIdx]);
    lowDimProblem->applyInitialSolution(sol[lowDimIdx]);

    couplingManager->init(bulkProblem, lowDimProblem, sol);
    bulkProblem->computeCouplingCoefficients();
    lowDimProblem->computeCouplingCoefficients();
    bulkProblem->computePointSourceMap();
    lowDimProblem->computePointSourceMap();

    using BulkGridVariables = GetPropType<BulkTypeTag, Properties::GridVariables>;
    auto bulkGridVariables = std::make_shared<BulkGridVariables>(bulkProblem, bulkGridGeometry);
    bulkGridVariables->init(sol[bulkIdx]);
    using LowDimGridVariables = GetPropType<LowDimTypeTag, Properties::GridVariables>;
    auto lowDimGridVariables = std::make_shared<LowDimGridVariables>(lowDimProblem, lowDimGridGeometry);
    lowDimGridVariables->init(sol[lowDimIdx]);

    using Assembler = MultiDomainFVAssembler<Traits, CouplingManager, DiffMethod::numeric>;
    auto assembler = std::make_shared<Assembler>(std::make_tuple(bulkProblem, lowDimProblem),
                                                 std::make_tuple(bulkGridGeometry, lowDimGridGeometry),
                                                 std::make_tuple(bulkGridVariables, lowDimGridVariables),
                                                 couplingManager);

    // the problem is linear: one Newton step from the initial solution solves it
    couplingManager->updateSolution(sol);
    assembler->assembleJacobianAndResidual(sol);
    auto deltaSol = sol;
    BlockDiagILU0BiCGSTABSolver linearSolver;
    if (!linearSolver.solve(assembler->jacobian(), deltaSol, assembler->residual()))
        DUNE_THROW(Dune::MathError, "Linear solver did not converge");
    sol -= deltaSol;
    couplingManager->updateSolution(sol);
    assembler->updateGridVariables(sol);

    // the exchange leaving the vessel has to enter the tissue
    const Scalar tissueSource = bulkProblem->computeSourceIntegral(sol[bulkIdx], *bulkGridVariables);
    const auto [vesselSource, errorExchange] = lowDimProblem->computeSourceIntegral(sol[lowDimIdx], *lowDimGridVariables);
    const Scalar exactSource = exactSolution->exchangeRate()*exactSolution->length();
    std::cout << std::setprecision(10)
              << "Total exchange: vessel " << vesselSource << ", tissue " << tissueSource
              << ", exact " << exactSource << '\n';
    if (std::abs(vesselSource - tissueSource) > 1e-8*exactSource)
        DUNE_THROW(Dune::InvalidStateException, "The exchange is not conservative: vessel " << vesselSource << ", tissue " << tissueSource);

    if (getParam<bool>("Vtk.EnableVtkOutput", true))
    {
        using BulkSolutionVector = std::decay_t<decltype(sol[bulkIdx])>;
        VtkOutputModule<BulkGridVariables, BulkSolutionVector> bulkVtkWriter(*bulkGridVariables, sol[bulkIdx], bulkProblem->name());
        GetPropType<BulkTypeTag, Properties::IOFields>::initOutputModule(bulkVtkWriter);
        bulkProblem->addVtkOutputFields(bulkVtkWriter);
        bulkVtkWriter.write(1.0);

        using LowDimSolutionVector = std::decay_t<decltype(sol[lowDimIdx])>;
        VtkOutputModule<LowDimGridVariables, LowDimSolutionVector> lowDimVtkWriter(*lowDimGridVariables, sol[lowDimIdx], lowDimProblem->name());
        GetPropType<LowDimTypeTag, Properties::IOFields>::initOutputModule(lowDimVtkWriter);
        lowDimProblem->addVtkOutputFields(lowDimVtkWriter);
        lowDimVtkWriter.write(1.0);
    }

    const double errorTissue = relativeDiscreteL2Error(*bulkProblem, sol[bulkIdx]);
    const double errorVessel = relativeDiscreteL2Error(*lowDimProblem, sol[lowDimIdx]);
    if (!std::isfinite(errorTissue) || !std::isfinite(errorVessel) || !std::isfinite(errorExchange))
        DUNE_THROW(Dune::InvalidStateException, "Error norms are not finite");

    std::cout << std::setprecision(6) << std::scientific
              << "Errors: hmax_3d " << maxElementDiameter(bulkGridView)
              << " hmax_1d " << maxElementDiameter(lowDimGridView)
              << " p_t " << errorTissue << " p_v " << errorVessel << " q " << errorExchange << std::endl;

    if (mpiHelper.rank() == 0)
    {
        Parameters::print();
        DumuxMessage::print(/*firstCall=*/false);
    }

    return 0;
}
