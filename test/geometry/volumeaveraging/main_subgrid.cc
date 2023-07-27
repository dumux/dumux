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
    using IntArray = std::array<int, dim>;
    using GlobalPosition = Dune::FieldVector<Scalar, dim>;

    ////// Set up the host-grids
    using HostGrid = Dune::YaspGrid<dim,Dune::EquidistantOffsetCoordinates<Scalar,dim>>;
    using HostGridManager = GridManager<HostGrid>;
    HostGridManager sourceHostGridManager;
    sourceHostGridManager.init("Source");
    auto& sourceHostGrid = sourceHostGridManager.grid();

    HostGridManager targetHostGridManager;
    targetHostGridManager.init("Target");
    auto& targetHostGrid = targetHostGridManager.grid();

    // Collect the image type, number and set up the stamped images.
    const std::string imgFileName = getParam<std::string>("Grid.ImagePath");
    const bool marked = getParam<bool>("Grid.Marker", 0);
    const auto numRepeats = getParam<IntArray>("Grid.Repeat");
    const auto img = NetPBMReader::readPBM(imgFileName);

    // rewrite img file according to stamp request
    std::vector<bool> stampedImg;
    const int numRows = img.header().nRows;
    const int numCols = img.header().nCols;

    // Reorder the single img file to fill the stamped vector
    for (int stampY = 0; stampY < numRepeats[1]; stampY++)
        for (int j = 0; j < numRows; j++)
            for (int stampX = 0; stampX < numRepeats[0]; stampX++)
                for (int i = 0; i < numCols; i++)
                    stampedImg.push_back(img[ (j * numCols + i) ]);

    GlobalPosition sourceBBoxMin = getParam<GlobalPosition>("Source.Grid.CutMin");
    GlobalPosition sourceBBoxMax = getParam<GlobalPosition>("Source.Grid.CutMax");

    // Set up how to cut up the Source grid.
    auto sourceGridSelector = [&sourceHostGrid, &sourceBBoxMin, &sourceBBoxMax, &stampedImg, &marked](const auto& element)
    {
        Scalar eps = 1e-6;
        GlobalPosition globalPos = element.geometry().center();
        if ((globalPos[0] < (sourceBBoxMin[0] - eps)) || (globalPos[1] < (sourceBBoxMin[1] - eps)) ||
            (globalPos[0] > (sourceBBoxMax[0] + eps)) || (globalPos[1] > (sourceBBoxMax[1] + eps)) )
            return false;
        else
            return stampedImg[sourceHostGrid.leafGridView().indexSet().index(element) ] == marked;
    };

    //! Build the Source grid.
    using SubGrid = Dune::SubGrid<dim, HostGrid>;
    using SubGridManager = Dumux::GridManager<SubGrid>;
    SubGridManager sourceSubGridManager;
    sourceSubGridManager.init(sourceHostGrid, sourceGridSelector, "SubGridSource");
    const auto& sourceSubGridView = sourceSubGridManager.grid().leafGridView();

    //! Set up how to cut up the Tracer grid.
    auto targetGridSelector = [&targetHostGrid, &img, &marked](const auto& element)
    { return img[targetHostGrid.leafGridView().indexSet().index(element)] == marked; };

    //! Build the Tracer grid.
    SubGridManager targetSubGridManager;
    targetSubGridManager.init(targetHostGrid, targetGridSelector, "SubGridTarget");
    const auto& targetSubGridView = targetSubGridManager.grid().leafGridView();

    //! Create the Grid Geometries
    using PoreScaleGridGeometry = CCTpfaFVGridGeometry<typename SubGrid::LeafGridView>;
    auto sourceSubGridGeometry = std::make_shared<PoreScaleGridGeometry>(sourceSubGridView);
    auto targetSubGridGeometry = std::make_shared<PoreScaleGridGeometry>(targetSubGridView);

    //! Set up I/O
    std::string sourceGridOutputName = "SourceGrid";
    std::string targetGridOutputName = "TargetGrid";
    VtkOutputModuleBase<PoreScaleGridGeometry> sourceGridVtkWriter(*sourceSubGridGeometry, sourceGridOutputName , "Double");
    VtkOutputModuleBase<PoreScaleGridGeometry> targetGridVtkWriter(*targetSubGridGeometry, targetGridOutputName , "Double");

    // Open the file for writing
    const bool readIndexData = getParam<bool>("IO.ReadIndexData");
    const bool writeIndexData = getParam<bool>("IO.WriteIndexData");
    const std::string dataFileName = getParam<std::string>("Grid.ImagePath");

    Dune::Timer convSetUpTimer;
    const Scalar filterSize = 1.0;
    using ConvolutionalAveragingManager = VolumeAveraging::ConvolutionalAveragingManager<PoreScaleGridGeometry,PoreScaleGridGeometry>;
    std::shared_ptr<ConvolutionalAveragingManager> convolutionalAveragingManager = nullptr;
    if (readIndexData)
        convolutionalAveragingManager = std::make_shared<ConvolutionalAveragingManager>(dataFileName + ".bin", targetSubGridGeometry, sourceSubGridGeometry, filterSize);
    else
        convolutionalAveragingManager = std::make_shared<ConvolutionalAveragingManager>(targetSubGridGeometry, sourceSubGridGeometry, filterSize);

    convSetUpTimer.stop();
    std::cout << "Establishing the Conv Averaging Manager took: " << convSetUpTimer.elapsed() << "\n";

    if (writeIndexData)
        convolutionalAveragingManager->writeIndexDataToFile(dataFileName + ".bin");

    auto indexSetSize = convolutionalAveragingManager->indexCountPerCell();
    targetGridVtkWriter.addField(indexSetSize, "IndexCount", Vtk::Precision::float64);

    using std::sin;
    using std::cos;
    std::vector<Scalar> tracer(sourceSubGridGeometry->numDofs(), 0.0);
    std::vector<GlobalPosition> velocity(sourceSubGridGeometry->numDofs(), GlobalPosition(0.0));
    for (const auto& element : elements(sourceSubGridGeometry->gridView()))
    {
        int eIdx = sourceSubGridGeometry->elementMapper().index(element);
        GlobalPosition globalPos = element.geometry().center();
        tracer[eIdx] = sin(2.0*M_PI*globalPos[0]) * cos(2.0*M_PI*globalPos[1]);
        velocity[eIdx][0] = 2.0*M_PI*cos(2.0*M_PI*globalPos[0]) * cos(2.0*M_PI*globalPos[1]);
        velocity[eIdx][1] = -2.0*M_PI*sin(2.0*M_PI*globalPos[0]) * sin(2.0*M_PI*globalPos[1]);
    }
    sourceGridVtkWriter.addField(tracer, "tracer", Vtk::Precision::float64);
    sourceGridVtkWriter.addField(velocity, "velocity", Vtk::Precision::float64);


    auto averageTracer = convolutionalAveragingManager->averageSolution(tracer);
    auto averageVelocity = convolutionalAveragingManager->averageSolution(velocity);
    targetGridVtkWriter.addField(averageTracer, "averageTracer", Vtk::Precision::float64);
    targetGridVtkWriter.addField(averageVelocity, "averageVelocity", Vtk::Precision::float64);

    targetGridVtkWriter.write(0.0);
    sourceGridVtkWriter.write(0.0);

    return 0;

} // end main
