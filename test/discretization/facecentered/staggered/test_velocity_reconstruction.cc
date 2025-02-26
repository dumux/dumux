// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \brief Test for finite volume element geometry, sub control volume, and sub
          control volume faces
 */
#include <config.h>

#include <dune/grid/spgrid.hh>
#include <dune/grid/io/file/vtk/vtkwriter.hh>

#include <dumux/common/initialize.hh>
#include <dumux/common/parameters.hh>

#include <dumux/io/grid/gridmanager_sp.hh>
#include <dumux/io/vtkoutputmodule.hh>

#include <dumux/discretization/cellcentered/tpfa/fvgridgeometry.hh>
#include <dumux/discretization/facecentered/staggered/fvgridgeometry.hh>

#include <dumux/freeflow/navierstokes/momentum/velocityreconstruction.hh>

namespace Dumux::Test {

template<class GlobalPosition>
double xSolution(const GlobalPosition& globalPos)
{ return std::sin((2.0 * M_PI) * globalPos[0]); }  // One period across [0, 1]

template<class GlobalPosition>
double ySolution(const GlobalPosition& globalPos)
{ return std::cos((2.0 * M_PI) * globalPos[1]); }  // One period across [0, 1]

template<class GlobalPosition>
GlobalPosition solution(const GlobalPosition& globalPos)
{
    GlobalPosition velocityVector(0.0);
    velocityVector[0] = xSolution(globalPos);
    velocityVector[1] = ySolution(globalPos);
    return velocityVector;
}

} // end namespace Dumux::Test

int main (int argc, char *argv[])
{
    using namespace Dumux;

    // maybe initialize MPI and/or multithreading backend
    initialize(argc, argv);

    // parse command line arguments and input file
    Parameters::init(argc, argv);
    using Scalar = double;
    constexpr int dim = 2;
    using IntArray = std::array<int, dim>;
    std::string groupName = "VelocityReconstruction";
    std::string outputName = getParamFromGroup<std::string>(groupName, "VTK.OutputName");
    bool verbose = getParamFromGroup<bool>(groupName, "Problem.Verbose");

    // Theoretical sinusoidal linear reconstruction error: 2.0/(cells per period)
    Scalar thresholdX = 2.0 / getParamFromGroup<IntArray>(groupName, "Grid.Cells")[0];
    Scalar thresholdY = 2.0 / getParamFromGroup<IntArray>(groupName, "Grid.Cells")[1];
    Scalar threshold = std::sqrt(thresholdX*thresholdX + thresholdY*thresholdY);

    // Create a periodic grid for evaluation
    using Grid = Dune::SPGrid<Scalar, dim>;
    using GlobalPosition = typename Dune::FieldVector<Grid::ctype, dim>;
    using GridManager = GridManager<Grid>;
    GridManager gridManager;
    gridManager.init("VelocityReconstruction");
    const auto& leafGridView = gridManager.grid().leafGridView();
    using MomentumGridGeometry = FaceCenteredStaggeredFVGridGeometry<typename Grid::LeafGridView, /*enable caching=*/ true>;
    auto momentumGridGeometry = std::make_shared<MomentumGridGeometry>(leafGridView);
    using MassGridGeometry = CCTpfaFVGridGeometry<typename Grid::LeafGridView, /*enable caching=*/ true>;
    auto massGridGeometry = std::make_shared<MassGridGeometry>(leafGridView);
    std::cout << "\nMomentum Grid Geometry is periodic: " << std::boolalpha << momentumGridGeometry->isPeriodic() << "\n";
    std::cout << "Mass Grid Geometry is periodic: " << std::boolalpha << massGridGeometry->isPeriodic() << "\n\n";


    std::vector<GlobalPosition> velocityVectorSolutionFace(momentumGridGeometry->numScv(), GlobalPosition(0.0));
    std::vector<GlobalPosition> velocityVectorSolutionCell(massGridGeometry->numScv(), GlobalPosition(0.0));
    std::cout << "Num momentum scvs: " << momentumGridGeometry->numScv() << "\n";
    std::cout << "Num mass scvs: " << massGridGeometry->numScv() << "\n";
    std::cout << "Num mass scvfs: " << massGridGeometry->numScvf() << "\n";
    std::vector<Scalar> faceOnlySolution(momentumGridGeometry->numScv(), 0.0);

    // Store a full vector solution at each face, and the face only solution at each face
    // Iterate over elements to access each element, and also each SCV
    auto momentumFVGeometry = localView(*momentumGridGeometry);
    for (const auto& element : elements(leafGridView))
    {
        momentumFVGeometry.bind(element);
        auto eIdx = momentumGridGeometry->elementMapper().index(element);
        velocityVectorSolutionCell[eIdx] = Test::solution(element.geometry().center());
        if (verbose)
            std::cout << "\nElement #" << eIdx << ":\n";
        for (auto&& scv : scvs(momentumFVGeometry))
        {
            if (verbose)
            {
                std::cout << "    SCV Index: " << scv.index()
                        << ", Location: (" << scv.dofPosition()[0] << ", " << scv.dofPosition()[1] << ")"
                        << ", Dof Index: " << scv.dofIndex() <<", Dof Axis: " << std::to_string(scv.dofAxis()) << "\n";
            }

            velocityVectorSolutionFace[scv.index()] = Test::solution(scv.dofPosition());
            faceOnlySolution[scv.index()] = velocityVectorSolutionFace[scv.index()][scv.dofAxis()];
        }
    }

    std::cout << "Checking the reconstruction of the full velocity vectors at the cell centers" << std::endl;
    // Test the face reconstruction methods
    std::vector<GlobalPosition> cellReconstructedVelocityVector(massGridGeometry->numScv(), GlobalPosition(0.0));
    int cellErrorCount = 0;
    auto massFVGeometry = localView(*massGridGeometry);
    for (const auto& element : elements(leafGridView))
    {
        auto eIdx = massGridGeometry->elementMapper().index(element);
        massFVGeometry.bind(element);
        for (auto&& scv : scvs(massFVGeometry))
        {
            // ConsistantlyOrientedgrids massSCVF == momentumSCV
            const auto getFaceVelocity = [&](const auto& fvGeometry, const auto& scvf){ return velocityVectorSolutionFace[scvf.index()]; };
            cellReconstructedVelocityVector[eIdx] = StaggeredVelocityReconstruction::cellCenterVelocity(getFaceVelocity, massFVGeometry);

            if ((cellReconstructedVelocityVector[eIdx] - velocityVectorSolutionCell[eIdx]).two_norm() > threshold)
            {
                std::cout << "error in element: " << eIdx << " at scv(dof): " << scv.dofIndex() << ". "
                          << "True sol: " <<  velocityVectorSolutionCell[scv.dofIndex()] << " vs. "
                          << "Reconstructed sol: " << cellReconstructedVelocityVector[scv.dofIndex()] << "\n";
                cellErrorCount++;
            }
        }
    }

    if (cellErrorCount > 0)
        DUNE_THROW(Dune::Exception, "The reconstructed velocities at the cell centers are incorrect.");
    else
        std::cout << "Cell Centered reconstruction successful. \n\n";

    // Test the face reconstruction methods
    std::cout << "Checking the reconstruction of the full velocity vectors at the face centers" << std::endl;
    std::vector<GlobalPosition> faceReconstructedVelocityVector(momentumGridGeometry->numScv(), GlobalPosition(0.0));
    int faceErrorCount = 0;
    // Iterate through all faces and reconstruct a velocity vector at each face.
    for (const auto& element : elements(leafGridView))
    {
        auto eIdx = momentumGridGeometry->elementMapper().index(element);
        momentumFVGeometry.bind(element);
        for (auto&& scv : scvs(momentumFVGeometry))
        {
            const auto getVelocitySCV = [&](const auto& scv){ return faceOnlySolution[scv.index()]; };
            faceReconstructedVelocityVector[scv.index()] = StaggeredVelocityReconstruction::faceVelocity(scv, momentumFVGeometry, getVelocitySCV);
            if ((faceReconstructedVelocityVector[scv.index()] - velocityVectorSolutionFace[scv.index()]).two_norm() > threshold)
            {
                std::cout << "error in element: " << eIdx << " at scv(dof, idx): " << scv.dofIndex() << ", " << scv.index() << ". "
                          << "True sol: " <<  velocityVectorSolutionFace[scv.index()] << " vs. "
                          << "Reconstructed sol: " << faceReconstructedVelocityVector[scv.index()] << "\n";
                faceErrorCount++;
            }
        }
    }

    if (faceErrorCount > 0)
        DUNE_THROW(Dune::Exception, "The reconstructed velocities at the faces are incorrect.");
    else
        std::cout << "Face Centered reconstruction successful. \n\n";

    VtkOutputModuleBase<MassGridGeometry> vtkWriter(*massGridGeometry, outputName , "Double");
    vtkWriter.addField(velocityVectorSolutionCell, "CCSolution");
    vtkWriter.addField(cellReconstructedVelocityVector, "CCReconstructed");
    vtkWriter.write(0.0);

    return 0;
}
