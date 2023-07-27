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
 * \ingroup Pore Scale Simulations
 * \ingroup Basic Pore Structure
 * \brief A test for the evaluation of gradients
 */

#include <config.h>
#include <cstdlib>
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <cmath>

#include <dune/common/parallel/mpihelper.hh>
#include <dune/common/timer.hh>

#include <dune/grid/io/file/dgfparser/dgfexception.hh>
#include <dune/grid/io/file/gmshwriter.hh>
#include <dune/grid/io/file/vtk.hh>
#include <dune/grid/yaspgrid.hh>
#include <dune/istl/io.hh>

#include <dumux/common/initialize.hh>
#include <dumux/common/parameters.hh>
#include <dumux/common/properties.hh>

#include <dumux/io/grid/gridmanager_yasp.hh>
#include <dumux/io/grid/gridmanager_sub.hh>
#include <dumux/io/vtk/function.hh>
#include <dumux/io/vtk/intersectionwriter.hh>
#include <dumux/io/vtkoutputmodule.hh>

#include <dumux/discretization/method.hh>
#include <dumux/discretization/cellcentered/tpfa/fvgridgeometry.hh>

#include <dumux/geometry/volumeaveraging.hh>

int main(int argc, char** argv)
{
    using namespace Dumux;

    // Initialize MPI, and print dumux start message
    Dumux::initialize(argc, argv);

    // parse command line arguments and input file
    Parameters::init(argc, argv);

    // Define dimension and type containers
    constexpr int dim = 2;
    using Scalar = double;
    using GlobalPosition = Dune::FieldVector<Scalar, dim>;

    ////// Set up the host-grids
    using Grid = Dune::YaspGrid<dim,Dune::EquidistantOffsetCoordinates<Scalar,dim>>;
    using GridManager = GridManager<Grid>;
    using GridGeometry = CCTpfaFVGridGeometry<typename Grid::LeafGridView>;

    // Create a coarse base grid, gridGeometry, and VTKOutput Module
    GridManager baseGridManager;
    baseGridManager.init("Base");
    auto& baseGrid = baseGridManager.grid();
    const auto& baseGridView = baseGrid.leafGridView();
    auto baseGridGeometry = std::make_shared<GridGeometry>(baseGridView);
    std::string baseGridOutputName = "BaseGrid";
    VtkOutputModuleBase<GridGeometry> baseGridVtkWriter(*baseGridGeometry, baseGridOutputName , "Double");

    // Create a finer grid, gridGeometry, and VTKOutput Module
    GridManager extendedGridManager;
    extendedGridManager.init("Extended");
    auto& extendedGrid = extendedGridManager.grid();
    const auto& extendedGridView = extendedGrid.leafGridView();
    auto extendedGridGeometry = std::make_shared<GridGeometry>(extendedGridView);
    std::string extendedGridOutputName = "ExtendedGrid";
    VtkOutputModuleBase<GridGeometry> extendedGridVtkWriter(*extendedGridGeometry, extendedGridOutputName , "Double");

    // Set up the Unit Averaging Manager Coarse -> Fine
    using UnitAveragingManager = VolumeAveraging::UnitAveragingManager<GridGeometry, GridGeometry>;
    std::shared_ptr<UnitAveragingManager> coarseFineUnitAveragingManager = nullptr;
    coarseFineUnitAveragingManager = std::make_shared<UnitAveragingManager>(baseGridGeometry, extendedGridGeometry);

    // Set up the Unit Averaging Manager Fine -> Coarse
    auto fineCoarseUnitAveragingManager = std::make_shared<UnitAveragingManager>(extendedGridGeometry, baseGridGeometry);


    // Set up the convolutional Averaging Manager
    Dune::Timer convSetUpTimer;
    const Scalar filterSize = getParam<Scalar>("Averaging.FilterSize");
    const bool useCircleFilter = getParam<bool>("Averaging.UseCircleFilter");
    using ConvolutionalAveragingManager = VolumeAveraging::ConvolutionalAveragingManager<GridGeometry, GridGeometry>;

    // Open the file for writing
    const bool readIndexData = getParam<bool>("IO.ReadIndexData");
    const bool writeIndexData = getParam<bool>("IO.WriteIndexData");
    const std::string dataFileName = getParam<std::string>("IO.FilePath");
    std::shared_ptr<ConvolutionalAveragingManager> convolutionalAveragingManager = nullptr;
    if (readIndexData)
        convolutionalAveragingManager = std::make_shared<ConvolutionalAveragingManager>(dataFileName, extendedGridGeometry, extendedGridGeometry, filterSize, useCircleFilter);
    else
        convolutionalAveragingManager = std::make_shared<ConvolutionalAveragingManager>(extendedGridGeometry, extendedGridGeometry, filterSize, useCircleFilter);

    if (writeIndexData)
        convolutionalAveragingManager->writeIndexDataToFile(dataFileName);

    convSetUpTimer.stop();
    std::cout << "Establishing the Conv Averaging Manager took: " << convSetUpTimer.elapsed() << "\n";

    // Create some coarse solution
    std::vector<Scalar> scalarField(baseGridGeometry->numDofs(), 0.0);
    std::vector<GlobalPosition> vectorField(baseGridGeometry->numDofs(), GlobalPosition(0.0));
    for (const auto& element : elements(baseGridGeometry->gridView()))
    {
        int eIdx = baseGridGeometry->elementMapper().index(element);
        GlobalPosition globalPos = element.geometry().center();
        scalarField[eIdx] = (globalPos - GlobalPosition(0.0)).two_norm();
        vectorField[eIdx][0] = (globalPos - GlobalPosition(0.0))[0];
        vectorField[eIdx][1] = (globalPos - GlobalPosition(0.0))[1];
    }
    baseGridVtkWriter.addField(scalarField, "scalarField", Vtk::Precision::float64);
    baseGridVtkWriter.addField(vectorField, "vectorField", Vtk::Precision::float64);

    // Extend the base solution the fine grid
    auto extendedScalarField = coarseFineUnitAveragingManager->unitAverageSolution(scalarField, false);
    auto extendedVectorField = coarseFineUnitAveragingManager->unitAverageSolution(vectorField, false);
    extendedGridVtkWriter.addField(extendedScalarField, "extendedScalarField", Vtk::Precision::float64);
    extendedGridVtkWriter.addField(extendedVectorField, "extendedVectorField", Vtk::Precision::float64);

    auto averageScalarField = convolutionalAveragingManager->averageSolution(extendedScalarField);
    auto averageVectorField = convolutionalAveragingManager->averageSolution(extendedVectorField);
    extendedGridVtkWriter.addField(averageScalarField, "averageScalarField", Vtk::Precision::float64);
    extendedGridVtkWriter.addField(averageVectorField, "averageVectorField", Vtk::Precision::float64);

    auto indexSetSize = convolutionalAveragingManager->indexCountPerCell();
    extendedGridVtkWriter.addField(indexSetSize, "IndexCount", Vtk::Precision::float64);

    // Reduce the extended Fields to the corase grid
    auto averagedExtendedScalarField = fineCoarseUnitAveragingManager->unitAverageSolution(extendedScalarField);
    auto averagedExtendedVectorField = fineCoarseUnitAveragingManager->unitAverageSolution(extendedVectorField);
    baseGridVtkWriter.addField(averagedExtendedScalarField, "avg_ext_scalarField", Vtk::Precision::float64);
    baseGridVtkWriter.addField(averagedExtendedScalarField, "avg_ext_vectorField", Vtk::Precision::float64);

    // Reduce the convolutional Fields to the corase grid
    auto averagedConvolutionalScalarField = fineCoarseUnitAveragingManager->unitAverageSolution(averageScalarField);
    auto averagedConvolutionalVectorField = fineCoarseUnitAveragingManager->unitAverageSolution(averageVectorField);
    baseGridVtkWriter.addField(averagedConvolutionalScalarField, "avg_conv_scalarField", Vtk::Precision::float64);
    baseGridVtkWriter.addField(averagedConvolutionalVectorField, "avg_conv_vectorField", Vtk::Precision::float64);

    baseGridVtkWriter.write(0.0);
    extendedGridVtkWriter.write(0.0);

    return 0;

} // end main
