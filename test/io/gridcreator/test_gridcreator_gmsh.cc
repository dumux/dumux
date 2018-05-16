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
 *
 * \brief Test for gmsh interface of the grid creator
 */
#include <config.h>
#include <iostream>

#include <dune/geometry/referenceelements.hh>
#include <dune/grid/io/file/vtk.hh>
#include <dune/grid/common/mcmgmapper.hh>
#include <dune/common/parallel/mpihelper.hh>

#include <dumux/common/properties.hh>
#include <dumux/common/parameters.hh>
#include <dumux/io/gridcreator.hh>
#include <dumux/discretization/methods.hh>

namespace Dumux {

namespace Properties
{
    NEW_TYPE_TAG(GridCreatorGmshTest);
    SET_TYPE_PROP(GridCreatorGmshTest, Grid, Dune::UGGrid<3>);
    SET_TYPE_PROP(GridCreatorGmshTest, Scalar, double);
    SET_PROP(GridCreatorGmshTest, FVGridGeometry) {
        struct GG {
            static constexpr DiscretizationMethod discMethod = DiscretizationMethod::cctpfa;
        };
        using type = GG;
    };
}

template<class TypeTag>
class GridCreatorGmshTest
{
    using Scalar = typename GET_PROP_TYPE(TypeTag, Scalar);
    using Grid = typename GET_PROP_TYPE(TypeTag, Grid);
    static const int dim = Grid::dimension;
    using GridCreator = typename Dumux::GridCreator<TypeTag>;
    using ReferenceElements = typename Dune::ReferenceElements<Scalar, dim>;
    using VertexMapper = typename Dune::MultipleCodimMultipleGeomTypeMapper<typename Grid::LeafGridView>;

public:

    static void getBoundaryDomainMarkers(std::vector<int>& boundaryMarker)
    {
        const auto& gridView = GridCreator::grid().leafGridView();
        VertexMapper vertexMapper(gridView, Dune::mcmgVertexLayout());

        boundaryMarker.clear();
        boundaryMarker.resize(gridView.size(dim));
        for(auto eIt = gridView.template begin<0>(); eIt != gridView.template end<0>(); ++eIt)
        {
            for(auto isIt = gridView.ibegin(*eIt); isIt != gridView.iend(*eIt); ++isIt)
            {
                if(!isIt->boundary())
                    continue;

                const auto refElement = ReferenceElements::general(eIt->geometry().type());

                // loop over vertices of the intersection facet
                for(int vIdx = 0; vIdx < refElement.size(isIt->indexInInside(), 1, dim); vIdx++)
                {
                    // get local vertex index with respect to the element
                    int vIdxLocal = refElement.subEntity(isIt->indexInInside(), 1, vIdx, dim);
                    int vIdxGlobal = vertexMapper.subIndex(*eIt, vIdxLocal, dim);

                    // make sure we always take the lowest non-zero marker (problem dependent!)
                    if (boundaryMarker[vIdxGlobal] == 0)
                        boundaryMarker[vIdxGlobal] = GridCreator::getBoundaryDomainMarker(isIt->boundarySegmentIndex());
                    else
                    {
                        if (boundaryMarker[vIdxGlobal] > GridCreator::getBoundaryDomainMarker(isIt->boundarySegmentIndex()))
                            boundaryMarker[vIdxGlobal] = GridCreator::getBoundaryDomainMarker(isIt->boundarySegmentIndex());
                    }
                }
            }
        }
    }
};

}

int main(int argc, char** argv) try
{

    // initialize MPI, finalize is done automatically on exit
    Dune::MPIHelper::instance(argc, argv);

    // using declarations
    using TypeTag = TTAG(GridCreatorGmshTest);
    using Grid = typename GET_PROP_TYPE(TypeTag, Grid);
    using GridCreator = typename Dumux::GridCreator<TypeTag>;

    // Read the parameters from the input file
    Dumux::Parameters::init(argc, argv, "test_gridcreator_gmsh.input");

    // Make the grid
    GridCreator::makeGrid("Bifurcation");

    // Read the boundary markers and convert them to vertex flags (e.g. for use in a box method)
    // Write a map from vertex position to boundaryMarker
    std::vector<int> boundaryMarker;
    Dumux::GridCreatorGmshTest<TypeTag>::getBoundaryDomainMarkers(boundaryMarker);

    // construct a vtk output writer and attach the boundaryMakers
    Dune::VTKSequenceWriter<Grid::LeafGridView> vtkWriter(GridCreator::grid().leafGridView(), "bifurcation", ".", "");
    vtkWriter.addVertexData(boundaryMarker, "boundaryMarker");
    vtkWriter.write(0);

    // refine grid once. Due to parametrized boundaries this will result in a grid closer to the orginal geometry.
    GridCreator::grid().globalRefine(1);
    Dumux::GridCreatorGmshTest<TypeTag>::getBoundaryDomainMarkers(boundaryMarker);
    vtkWriter.write(1);
}
catch (Dumux::ParameterException &e) {
    std::cerr << e << ". Abort!\n";
    return 1;
}
catch (Dune::Exception &e) {
    std::cerr << "Dune reported error: " << e << std::endl;
    return 3;
}
catch (...) {
    std::cerr << "Unknown exception thrown!\n";
    return 4;
}
