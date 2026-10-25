// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup OnePTests
 * \brief Convergence test for the hybrid control-volume finite element schemes on a line
 *        and on a junction of lines embedded in three dimensions
 */
#include <config.h>

#include <algorithm>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <vector>

#include <dune/common/parallel/mpihelper.hh>
#include <dune/geometry/quadraturerules.hh>

#include <dumux/common/initialize.hh>
#include <dumux/common/properties.hh>
#include <dumux/common/parameters.hh>
#include <dumux/common/metadata.hh>

#include <dumux/nonlinear/newtonsolver.hh>
#include <dumux/linear/istlsolvers.hh>
#include <dumux/linear/linearsolvertraits.hh>
#include <dumux/linear/linearalgebratraits.hh>
#include <dumux/assembly/assembler.hh>

#include <dumux/io/grid/gridmanager_oned.hh>
#if HAVE_DUNE_FOAMGRID
#include <dumux/io/grid/gridmanager_foam.hh>
#endif

#include "properties.hh"

/*!
 * \brief L2 error of the pressure and H1 seminorm error of its gradient, integrated with
 *        the local basis of each element, and the largest element length
 */
template<class Problem, class SolutionVector>
void writeErrors(const Problem& problem, const SolutionVector& x)
{
    using namespace Dumux;
    const auto& gridGeometry = problem.gridGeometry();
    using GridView = std::decay_t<decltype(gridGeometry.gridView())>;
    using GlobalPosition = typename GridView::template Codim<0>::Geometry::GlobalCoordinate;
    static constexpr int dim = GridView::dimension;

    double l2 = 0.0, h1 = 0.0, hMax = 0.0;
    std::vector<Dune::FieldVector<double, 1>> values;
    std::vector<Dune::FieldMatrix<double, 1, dim>> jacobians;
    auto fvGeometry = localView(gridGeometry);
    for (const auto& element : elements(gridGeometry.gridView()))
    {
        fvGeometry.bind(element);
        const auto geometry = element.geometry();
        hMax = std::max(hMax, geometry.volume());
        const auto& localBasis = fvGeometry.feLocalBasis();
        for (const auto& qp : Dune::QuadratureRules<double, dim>::rule(geometry.type(), 12))
        {
            const auto& local = qp.position();
            localBasis.evaluateFunction(local, values);
            localBasis.evaluateJacobian(local, jacobians);
            const auto jacInvT = geometry.jacobianInverseTransposed(local);

            double p = 0.0;
            GlobalPosition gradP(0.0);
            for (const auto& localDof : localDofs(fvGeometry))
            {
                const double coefficient = x[localDof.dofIndex()][0];
                p += coefficient*values[localDof.index()][0];
                GlobalPosition gradPhi;
                jacInvT.mv(jacobians[localDof.index()][0], gradPhi);
                gradP.axpy(coefficient, gradPhi);
            }

            const auto globalPos = geometry.global(local);
            const double weight = qp.weight()*geometry.integrationElement(local);
            const double dp = p - problem.analyticalPressure(globalPos);
            l2 += weight*dp*dp;
            gradP -= problem.analyticalPressureGradient(globalPos);
            h1 += weight*gradP.two_norm2();
        }
    }

    l2 = std::sqrt(l2);
    h1 = std::sqrt(h1);
    std::cout << std::setprecision(10) << std::scientific
              << "** errors for " << gridGeometry.numDofs() << " dofs, h = " << hMax
              << ": L2 = " << l2 << ", H1 = " << h1 << std::endl;

    MetaData::Collector collector;
    if (MetaData::jsonFileExists(problem.name()))
        MetaData::readJsonFile(collector, problem.name());
    collector["numDofs"].push_back(gridGeometry.numDofs());
    collector["h"].push_back(hMax);
    collector["L2errors"].push_back(l2);
    collector["H1errors"].push_back(h1);
    MetaData::writeJsonFile(collector, problem.name());
}

int main(int argc, char** argv)
{
    using namespace Dumux;
    using TypeTag = Properties::TTag::TYPETAG;

    Dumux::initialize(argc, argv);
    Parameters::init(argc, argv);

    using Grid = GetPropType<TypeTag, Properties::Grid>;
    GridManager<Grid> gridManager;
    gridManager.init();
    const auto& leafGridView = gridManager.grid().leafGridView();

    using GridGeometry = GetPropType<TypeTag, Properties::GridGeometry>;
    auto gridGeometry = std::make_shared<GridGeometry>(leafGridView);

    using Problem = GetPropType<TypeTag, Properties::Problem>;
    auto problem = std::make_shared<Problem>(gridGeometry);

    using SolutionVector = GetPropType<TypeTag, Properties::SolutionVector>;
    SolutionVector x(gridGeometry->numDofs());

    using GridVariables = GetPropType<TypeTag, Properties::GridVariables>;
    auto gridVariables = std::make_shared<GridVariables>(problem, gridGeometry);
    gridVariables->init(x);

    using Assembler = Experimental::Assembler<TypeTag, DiffMethod::numeric>;
    auto assembler = std::make_shared<Assembler>(problem, gridGeometry, gridVariables);

    using LinearSolver = UMFPackIstlSolver<SeqLinearSolverTraits, LinearAlgebraTraitsFromAssembler<Assembler>>;
    auto linearSolver = std::make_shared<LinearSolver>();

    using NewtonSolver = Dumux::NewtonSolver<Assembler, LinearSolver>;
    NewtonSolver nonLinearSolver(assembler, linearSolver);
    nonLinearSolver.solve(x);

    writeErrors(*problem, x);

    if (Dune::MPIHelper::instance().rank() == 0)
        Parameters::print();

    return 0;
}
