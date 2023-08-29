// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*****************************************************************************
 *   See the file COPYING for full copying permissions.                      *
 *                                                                           *
 *   This program is free software: you can redistribute it and/or modify    *
 *   it under the terms of the GNU General Public License as published by    *
 *   the Free Software Foundation, either version 3 of the License, or       *
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
 *
 * \brief test for the two-phase pore-network model
 */
#include <config.h>

#include <ctime>
#include <iostream>

#include <dune/common/parallel/mpihelper.hh>
#include <dune/common/timer.hh>
#include <dune/grid/io/file/dgfparser/dgfexception.hh>
#include <dune/grid/io/file/vtk.hh>

#include <dumux/assembly/fvassembler.hh>
#include <dumux/common/properties.hh>
#include <dumux/common/parameters.hh>
#include <dumux/common/dumuxmessage.hh>
#include <dumux/io/grid/porenetwork/gridmanager.hh>
#include <dumux/linear/seqsolverbackend.hh>
#include <dumux/porenetwork/common/pnmvtkoutputmodule.hh>
#include <dumux/porenetwork/2p/newtonsolver.hh>
#include <dumux/porenetwork/2p/invasionstate.hh>

#include <vector>

#include "properties.hh"

int main(int argc, char** argv)
{
    using namespace Dumux;

    using TypeTag = Properties::TTag::DrainageProblem;
    const auto& mpiHelper = Dune::MPIHelper::instance(argc, argv);
    if (mpiHelper.rank() == 0)
        DumuxMessage::print(/*firstCall=*/true);
    Parameters::init(argc, argv);
    using GridManager = Dumux::PoreNetwork::GridManager<1>;
    GridManager gridManager;
    gridManager.init();
    const auto& leafGridView = gridManager.grid().leafGridView();
    auto gridData = gridManager.getGridData();

    // create the finite volume grid geometry
    using GridGeometry = GetPropType<TypeTag, Properties::GridGeometry>;
    auto gridGeometry = std::make_shared<GridGeometry>(leafGridView, *gridData);

    // the spatial parameters
    using SpatialParams = GetPropType<TypeTag, Properties::SpatialParams>;
    auto spatialParams = std::make_shared<SpatialParams>(gridGeometry);

    // the problem (boundary conditions)
    using Problem = GetPropType<TypeTag, Properties::Problem>;
    auto problem = std::make_shared<Problem>(gridGeometry, spatialParams);

    // the solution vector
    using GridView = typename GridGeometry::GridView;
    using SolutionVector = GetPropType<TypeTag, Properties::SolutionVector>;
    SolutionVector x(leafGridView.size(GridView::dimension));

    // the grid variables
    using GridVariables = GetPropType<TypeTag, Properties::GridVariables>;
    auto gridVariables = std::make_shared<GridVariables>(problem, gridGeometry);
    gridVariables->init(x);

    // need an throat and a pore
    auto fvGeometry = localView(*gridGeometry);
    auto dummyElement = gridGeometry->element(0);
    fvGeometry.bind(dummyElement);
    auto dummyElemVolVars = localView(gridVariables->curGridVolVars());
    dummyElemVolVars.bindElement(dummyElement, fvGeometry, x);
    auto dummyScv = fvGeometry.scv(0);
    auto dummyScvf = fvGeometry.scvf(0);
    auto dummyElemSol = elementSolution(dummyElement, x, fvGeometry.gridGeometry());
    auto fluidMatrixInteraction = problem->spatialParams().fluidMatrixInteraction(dummyElement, dummyScv, dummyElemSol);

    // get the plot interval
    double SwMin = getParam<double>("PlotInterval.SwMin", 0.0);
    double SwMax = getParam<double>("PlotInterval.SwMax", 1.0);
    double PcMin, PcMax;
    PcMax = fluidMatrixInteraction.pc(SwMin);
    PcMin = fluidMatrixInteraction.pc(SwMax);

    // pc-Sw
    unsigned int numIntervals = getParam<int>("PlotInterval.SamplePoints", 50000);
    std::vector<double> Sw(numIntervals + 1);
    std::vector<double> Pc(numIntervals + 1);
    std::vector<double> entryKw(numIntervals + 1);
    std::vector<double> entryKn(numIntervals + 1);
    std::vector<double> snapoffKw(numIntervals + 1);
    std::vector<double> snapoffKn(numIntervals + 1);

    double SwRange = SwMax - SwMin;
    std::cout << "The conductivity data for Sw between " << SwMin << " and " << SwMax << " will be ploted" << std::endl;
    std::cout << "pcMin is: " << PcMin << " pcMax is: " << PcMax << "\n";
    std::cout << "pc (Sw = 1) is : " << fluidMatrixInteraction.pc(1.0) << "\n";

    std::ofstream outputFile("Pc_Sw.log");
    for (int i = 0; i <= numIntervals; i++)
    {
        Sw[i] = SwMin + SwRange * double(i) /double(numIntervals);
        Pc[i] = fluidMatrixInteraction.pc(Sw[i]);
        using std::max;
        using std::min;
        PcMin = min(PcMin, Pc[i]);
        PcMax = max(PcMax, Pc[i]);
        outputFile << Sw[i] << " " <<  Pc[i] << std::endl;
    }

    // plot kn - sw
    auto elemFluxVarsCache = localView(gridVariables->gridFluxVarsCache());
    // ! Enable to make this work, we need to comment the caching in properties.hh
    auto fluxVarsCache = elemFluxVarsCache[dummyScvf];

    std::ofstream outputFile2("K-Sw.log");

    for (int i = 0; i <= numIntervals; i++)
    {
        fluxVarsCache.update(*problem, dummyElement, fvGeometry,
                              dummyElemVolVars, dummyScvf, false, Pc[i]);
        entryKw[i] = fluxVarsCache.transmissibility(0);
        entryKn[i] = fluxVarsCache.transmissibility(1);

        fluxVarsCache.update(*problem, dummyElement, fvGeometry,
                              dummyElemVolVars, dummyScvf, true, Pc[i]);
        snapoffKw[i] = fluxVarsCache.transmissibility(0);
        snapoffKn[i] = fluxVarsCache.transmissibility(1);

        outputFile2 << entryKw[i] << " "
                    << entryKn[i] << " "
                    << snapoffKw[i] << " "
                    << snapoffKn[i] << " "
                    << std::endl;
    }
}
