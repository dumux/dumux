// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*****************************************************************************
 *   See the file COPYING for full copying permissions.                      *
 *                                                                           *
 *   This program is free software: you can redistribute it and/or modify    *
 *   it under the terms of the GNU General Public License as published by    *
 *   the Free Software Foundation, either version 2 of the License, or       *
 *   (at your option) any later version.                                     *
 *                                                                           *
 *   This program is distributed in the hope that it will be useful,         *
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of          *
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the            *
 *   GNU General Public License for more details.                            *
 *                                                                           *
 *   You should have received a copy of the GNU General Public License       *
 *   along with this program.  If not, see <http://www.gnu.org/licenses/>.   *
 *****************************************************************************/
/*!
 * \file
 * \brief Convergence test for single-phase flow
 */
#include <config.h>

#include <iostream>
#include <fstream>
#include <utility>
#include <algorithm>
#include <type_traits>

#include <dune/functions/gridfunctions/discreteglobalbasisfunction.hh>
#include <dune/functions/gridfunctions/analyticgridviewfunction.hh>
#include <dune/common/parallel/mpihelper.hh>

#include <dumux/common/parameters.hh>
#include <dumux/common/properties.hh>
#include <dumux/common/dumuxmessage.hh>

#include "solver.hh"

#include <dumux/common/integrate.hh>
#include <dumux/multidomain/glue.hh>

#include <dumux/discretization/functionspacebasis.hh>
#include <dumux/discretization/projection/projector.hh>

// main program
int main(int argc, char** argv)
{
    using namespace Dumux;
    using TypeTag = Properties::TTag::TYPETAG;
    static constexpr auto dm = GetPropType<TypeTag, Properties::GridGeometry>::discMethod;
    static constexpr bool isBox = dm == DiscretizationMethod::box;

    // initialize MPI, finalize is done automatically on exit
    const auto& mpiHelper = Dune::MPIHelper::instance(argc, argv);

    // print dumux start message
    if (mpiHelper.rank() == 0)
        DumuxMessage::print(/*firstCall=*/true);

    // initialize parameter tree
    Parameters::init(argc, argv);

    // obtain the sequence of discretizations
    auto cellsVector = getParam< std::vector<int> >("Grid.CellsSequence");

    std::sort(cellsVector.begin(), cellsVector.end());
    cellsVector.erase(std::unique(cellsVector.begin(), cellsVector.end()), cellsVector.end());

    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    std::vector<Scalar> dxVector(cellsVector.size());
    std::transform(cellsVector.begin(),
                   cellsVector.end(),
                   dxVector.begin(),
                   [] (auto cells) { return 1.0/cells; });

    // run the finest-resolved simulation first
    std::cout << "\n --- Solving finest solution (dx = " << 1.0/cellsVector.back() << ") --- \n\n";
    const auto finestStorage = solveRefinementLevel<TypeTag>(cellsVector.back());
    const auto& finestGridGeometry = *finestStorage.gridGeometry;
    const auto& finestSolution = *finestStorage.solution;
    const auto& finestBasis = getFunctionSpaceBasis(finestGridGeometry);

    // grid function with the discrete solution
    using namespace Dune::Functions;
    using BlockType = typename std::decay_t<decltype(finestSolution)>::block_type;
    const auto gfFine = makeDiscreteGlobalBasisFunction<BlockType>(finestBasis, finestSolution);

    // grid function with the exact solution
    using Problem = GetPropType<TypeTag, Properties::Problem>;
    auto analyticSol = [] (const auto x) { return Problem::exact(x); };
    const auto gfAnalyticFine = makeAnalyticGridViewFunction(analyticSol, finestBasis.gridView());

    // run simulations with coarser grids and compute errors
    std::vector<Scalar> discreteErrors; discreteErrors.reserve(cellsVector.size()-1);
    std::vector<Scalar> analyticErrors; analyticErrors.reserve(cellsVector.size());

    cellsVector.pop_back();
    const auto intOrder = getParam<int>("L2Norm.IntegrationOrder");
    for (auto dx : cellsVector)
    {
        std::cout << "\n --- Solving with dx = " << dx << " --- \n\n";
        const auto storage = solveRefinementLevel<TypeTag>(dx);
        const auto& gridGeometry = *storage.gridGeometry;
        const auto& solution = *storage.solution;
        const auto& basis = getFunctionSpaceBasis(gridGeometry);

        const auto gf = makeDiscreteGlobalBasisFunction<BlockType>(basis, solution);
        const auto gfAnalytic = makeAnalyticGridViewFunction(analyticSol, basis.gridView());

        std::cout << "Projecting solution into reference space" << std::endl;
        const auto glue = makeGlue(gridGeometry, finestGridGeometry);
        const auto projector = makeProjector(basis, finestBasis, glue);

        auto params = projector.defaultParams();
        params.residualReduction = 1e-16;

        const auto projectedSolution = projector.project(solution, params);
        const auto gfProjected = makeDiscreteGlobalBasisFunction<BlockType>(finestBasis, projectedSolution);

        std::cout << "Computing error norms" << std::endl;
        discreteErrors.push_back( integrateL2Error(finestBasis.gridView(), gfFine, gfProjected, intOrder) );
        analyticErrors.push_back( integrateL2Error(basis.gridView(), gf, gfAnalytic, intOrder) );
    }

    // compute analytical error for finest solution
    std::cout << "Computing error norm for finest grid w.r.t. analytical solution" << std::endl;
    analyticErrors.push_back( integrateL2Error(finestBasis.gridView(), gfFine, gfAnalyticFine, intOrder) );

    using std::log;
    std::cout << "\nComputed the following errors/rates w.r.t. the analytical solution:\n";
    for (unsigned int i = 0; i < analyticErrors.size(); ++i)
    {
        std::cout << analyticErrors[i];
        if (i == 0)
            std::cout << "\n";
        else
        {
            const auto rate = (log(analyticErrors[i]) - log(analyticErrors[i-1]))
                              /(log(dxVector[i]) - log(dxVector[i-1]));

            std::cout << " \t--->\t " << rate << std::endl;
            if (i == analyticErrors.size()-1)
                if ( (isBox && rate < 1.95) || (!isBox && rate < 0.95) )
                    DUNE_THROW(Dune::InvalidStateException, "Computed rate is below expected value: " << rate);
        }
    }

    std::cout << "\nComputed the following errors/rates w.r.t. the discrete solution:\n";
    for (unsigned int i = 0; i < discreteErrors.size(); ++i)
    {
        std::cout << discreteErrors[i];
        if (i == 0)
            std::cout << "\n";
        else
        {
            const auto rate = (log(discreteErrors[i]) - log(discreteErrors[i-1]))
                              /(log(dxVector[i]) - log(dxVector[i-1]));

            std::cout << " \t--->\t " << rate << std::endl;
            if (i == discreteErrors.size()-1)
                if ( (isBox && rate < 1.95) || (!isBox && rate < 0.95) )
                    DUNE_THROW(Dune::InvalidStateException, "Computed rate is below expected value: " << rate);
        }
    }

    // maybe write the errors to disk
    if (getParam<bool>("L2Norm.WriteLogFiles", false))
    {
        std::ofstream analyticFile(getParam<std::string>("L2Norm.AnalyticErrorsFile"));
        std::ofstream discreteFile(getParam<std::string>("L2Norm.DiscreteErrorsFile"));

        for (unsigned int i = 0; i < analyticErrors.size(); ++i)
            analyticFile << dxVector[i] << "," << analyticErrors[i] << "\n";
        for (unsigned int i = 0; i < discreteErrors.size(); ++i)
            discreteFile << dxVector[i] << "," << discreteErrors[i] << "\n";
    }

    // show parameters used
    Parameters::print();

    // print dumux good bye message
    if (mpiHelper.rank() == 0)
        DumuxMessage::print(/*firstCall=*/false);

    return 0;
}
