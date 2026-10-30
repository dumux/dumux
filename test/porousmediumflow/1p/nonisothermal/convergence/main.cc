// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup OnePTests
 * \brief Convergence test of the non-isothermal single-phase model with a manufactured solution
 */
#include <config.h>

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
#include <dumux/io/grid/gridmanager_yasp.hh>
#include <dumux/io/grid/gridmanager_alu.hh>
#include <dumux/assembly/assembler.hh>
#include <dumux/linear/istlsolvers.hh>
#include <dumux/linear/linearsolvertraits.hh>
#include <dumux/linear/linearalgebratraits.hh>
#include <dumux/nonlinear/newtonsolver.hh>

#include "properties.hh"

/*!
 * \brief The L2 norm of the difference between the discrete and the exact solution
 * \note The discrete solution is evaluated with the basis functions of the scheme at
 *       quadrature points exact for polynomials of order 10.
 */
template<class GridDiscretization, class SolutionVector, class Exact>
double l2Error(const GridDiscretization& gridDiscretization, const SolutionVector& x, int pvIdx, const Exact& exact)
{
    constexpr int dim = GridDiscretization::GridView::dimension;
    double error = 0.0;
    auto elemDisc = localView(gridDiscretization);
    std::vector<Dune::FieldVector<double, 1>> shapeValues;
    for (const auto& element : elements(gridDiscretization.gridView()))
    {
        elemDisc.bindElement(element);
        const auto geometry = element.geometry();
        for (const auto& qp : Dune::QuadratureRules<double, dim>::rule(geometry.type(), 10))
        {
            elemDisc.feLocalBasis().evaluateFunction(qp.position(), shapeValues);
            double value = 0.0;
            for (const auto& localDof : localDofs(elemDisc))
                value += shapeValues[localDof.index()][0]*x[localDof.dofIndex()][pvIdx];

            const double delta = value - exact(geometry.global(qp.position()));
            error += qp.weight()*geometry.integrationElement(qp.position())*delta*delta;
        }
    }
    using std::sqrt;
    return sqrt(error);
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

    using GridGeometry = GetPropType<TypeTag, Properties::GridGeometry>;
    auto gridGeometry = std::make_shared<GridGeometry>(gridManager.grid().leafGridView());

    using Problem = GetPropType<TypeTag, Properties::Problem>;
    auto problem = std::make_shared<Problem>(gridGeometry);

    using Indices = typename GetPropType<TypeTag, Properties::ModelTraits>::Indices;
    using SolutionVector = GetPropType<TypeTag, Properties::SolutionVector>;
    SolutionVector x(gridGeometry->numDofs());
    for (auto& priVars : x)
    {
        priVars[Indices::pressureIdx] = getParam<double>("Problem.ReferencePressure");
        priVars[Indices::temperatureIdx] = getParam<double>("Problem.ReferenceTemperature");
    }

    using GridVariables = GetPropType<TypeTag, Properties::GridVariables>;
    auto gridVariables = std::make_shared<GridVariables>(problem, gridGeometry);
    gridVariables->init(x);

    using Assembler = Experimental::Assembler<TypeTag, DiffMethod::numeric>;
    auto assembler = std::make_shared<Assembler>(problem, gridGeometry, gridVariables);

    using LinearSolver = UMFPackIstlSolver<SeqLinearSolverTraits, LinearAlgebraTraitsFromAssembler<Assembler>>;
    auto linearSolver = std::make_shared<LinearSolver>();

    NewtonSolver<Assembler, LinearSolver> nonLinearSolver(assembler, linearSolver);
    nonLinearSolver.solve(x);

    const double errorT = l2Error(*gridGeometry, x, Indices::temperatureIdx,
                                  [&](const auto& pos) { return problem->exactTemperature(pos); });
    const double errorP = l2Error(*gridGeometry, x, Indices::pressureIdx,
                                  [&](const auto& pos) { return problem->exactPressure(pos); });

    std::cout << std::setprecision(8) << std::scientific
              << "** L2 errors for " << gridGeometry->numDofs() << " dofs: temperature " << errorT
              << ", pressure " << errorP << std::endl;

    MetaData::Collector collector;
    if (MetaData::jsonFileExists(problem->name()))
        MetaData::readJsonFile(collector, problem->name());
    collector["numDofs"].push_back(gridGeometry->numDofs());
    collector["numElements"].push_back(gridGeometry->gridView().size(0));
    collector["L2errorsTemperature"].push_back(errorT);
    collector["L2errorsPressure"].push_back(errorP);
    MetaData::writeJsonFile(collector, problem->name());

    Parameters::print();

    return 0;
}
