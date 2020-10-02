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
 * \brief Test for finite volume element geometry, sub control volume, and sub
          control volume faces
 */
#include <config.h>

#include <iostream>
#include <utility>
#include <iomanip>

#include <dune/common/test/iteratortest.hh>
#include <dune/grid/utility/structuredgridfactory.hh>
#include <dune/grid/yaspgrid.hh>
#include <dumux/common/parameters.hh>

#include <dumux/common/intersectionmapper.hh>
#include <dumux/common/defaultmappertraits.hh>
#include <dumux/discretization/cellcentered/subcontrolvolume.hh>
#include <dumux/discretization/staggered/fvelementgeometry.hh>
#include <dumux/discretization/staggered/freeflow/fvgridgeometrytraits.hh>
#include <dumux/discretization/staggered/fvgridgeometry.hh>
#include <dumux/discretization/staggered/freeflow/subcontrolvolumeface.hh>
#include <dumux/discretization/staggered/freeflow/staggeredgeometryhelper.hh>
#include <dumux/discretization/staggered/freeflow/connectivitymap.hh>

#ifndef DOXYGEN
namespace Dumux {
namespace Detail {
template<class T>
class NoopFunctor {
public:
  NoopFunctor() {}
  void operator()(const T& t){}
};
} // end namespace Detail

} // end namespace Dumux
#endif

int main (int argc, char *argv[])
{
    using namespace Dumux;

    // maybe initialize mpi
    Dune::MPIHelper::instance(argc, argv);
    std::cout << "Checking the FVGeometries, SCVs and SCV faces" << std::endl;

    // parse command line arguments and input file
    Parameters::init(argc, argv);

    using Grid = Dune::YaspGrid<3>;

    constexpr int dim = Grid::dimension;

    static constexpr int upwindSchemeOrder = 2;

    using GridGeometry = StaggeredFVGridGeometry<typename Grid::LeafGridView, /*enable caching=*/ true,
                                                 StaggeredFreeFlowDefaultFVGridGeometryTraits<typename Grid::LeafGridView, upwindSchemeOrder> >;
    using FVElementGeometry = typename GridGeometry::LocalView;
    using SubControlVolume = typename FVElementGeometry::SubControlVolume;
    using SubControlVolumeFace = typename FVElementGeometry::SubControlVolumeFace;
    using GlobalPosition = typename SubControlVolume::GlobalPosition;

    // make a grid
    GlobalPosition lower(0.0);
    GlobalPosition upper(10.0);
    std::array<unsigned int, dim> els{{5, 5, 5}};
    std::shared_ptr<Grid> grid = Dune::StructuredGridFactory<Grid>::createCubeGrid(lower, upper, els);

    auto leafGridView = grid->leafGridView();
    GridGeometry gridGeometry(leafGridView);
    gridGeometry.update();

    std::cout << "Abbreviatons:\n"
              << "pos - postition of face center\n"
              << "fIdx - face index\n"
              << "dofIdx (self/oppo.) - dofIdx on face (self/opposite)\n"
              << "norm in/out - dofIdx on face normal to own face (within own element / in adjacent element)\n"
              << "par - dofIdx on face parallel to own face\n"
              << "normDist - distance bewteen the dofs on the faces normal the own face\n"
              << "parDist - distance bewteen the dof on the parallel face and the one on the own face\n"
              << "norm in/out - dofIdx on face normal to own face (within own element / in adjacent element)" << std::endl;

    // iterate over elements. For every element get fv geometry and loop over scvs and scvfaces
    for (const auto& element : elements(leafGridView))
    {
        auto eIdx = gridGeometry.elementMapper().index(element);
        if(eIdx == 12 || eIdx == 0)
        {
            std::cout << std::endl << "Checking fvGeometry of element " << eIdx << std::endl;
            auto fvGeometry = localView(gridGeometry);
            fvGeometry.bind(element);

            auto range = scvs(fvGeometry);
            Detail::NoopFunctor<SubControlVolume> op;
            if(0 != testForwardIterator(range.begin(), range.end(), op))
                DUNE_THROW(Dune::Exception, "Iterator does not fulfill the forward iterator concept");

            for (auto&& scv : scvs(fvGeometry))
            {
                std::cout << "-- scv " << scv.dofIndex() << " center at: " << scv.center() << std::endl;
            }

            auto range2 = scvfs(fvGeometry);
            Detail::NoopFunctor<SubControlVolumeFace> op2;
            if(0 != testForwardIterator(range2.begin(), range2.end(), op2))
                DUNE_THROW(Dune::Exception, "Iterator does not fulfill the forward iterator concept");


            for (auto&& scvf : scvfs(fvGeometry))
            {
                std::cout << "\n";
                std::cout << " Scvf Index " << scvf.index()  << " , local Index " << scvf.localFaceIdx() << "\n";
                std::cout << " Center pos "<< scvf.ipGlobal() << "\n";

                if (scvf.boundary())
                {std::cout << " (on boundary)" << "\n";}
                else
                {std::cout << "\n";}

                std::cout <<  std::fixed << std::left << std::setprecision(2);

                std::cout << " On Axis Dof Index: \n";
                if(gridGeometry.upwindStencilOrder() > 1)
                {std::cout << " | Forward dofIdx : " << std::setw(3) << scvf.axisData().inAxisForwardDofs[0] << "\n";}
                std::cout << " | Self dofIdx : " << std::setw(3) << scvf.dofIndex() << "\n";
                std::cout << " | Opposite dofIdx : " << std::setw(3) << scvf.dofIndexOpposingFace() << "\n";
                if(gridGeometry.upwindStencilOrder() > 1)
                {std::cout << " | Backward dofIdx : " << std::setw(3) << scvf.axisData().inAxisBackwardDofs[0] << "\n";}

                std::cout << " Normal Dof Index: \n";
                for(int i = 0; i < scvf.pairData().size(); i++)
                {
                    std::cout << " | normal inner dofIdx "<< i <<": " << std::setw(3) << scvf.pairData(i).lateralPair.first << "\n";
                    std::cout << " | normal outer dofIdx "<< i <<": " << std::setw(3) << scvf.pairData(i).lateralPair.second << "\n";
                }

                std::cout << " Parallel Dof Index: \n";
                for(int i = 0; i < scvf.pairData().size(); i++)
                {
                    for(int j = 0; j < gridGeometry.upwindStencilOrder(); j++)
                    {
                        std::cout << " | Parallel Dof "<< j << " on axis " << i << ": "<<  std::setw(3) << scvf.pairData(i).parallelDofs[j] << "\n";
                    }
                }

                std::cout << " Distances: \n";
                if(gridGeometry.upwindStencilOrder() > 1)
                {std::cout << " | Opposite To Backwards Face Dist : " << std::setw(3) << scvf.axisData().inAxisBackwardDistances[0] << "\n";}
                std::cout << " | self To Opposite Dist : " << std::setw(3) << scvf.selfToOppositeDistance() << "\n";
                if(gridGeometry.upwindStencilOrder() > 1)
                {std::cout << " | self To Forwards Face Dist : " << std::setw(3) << scvf.axisData().inAxisForwardDistances[0] << "\n";}

                for(int i = 0; i < scvf.pairData().size(); i++)
                {
                    for(int j = 0; j < gridGeometry.upwindStencilOrder(); j++)
                    {
                        std::cout << " | Parallel Cell Widths "<< j << " on axis " << i << ": "<< std::setw(3) << scvf.pairData(i).parallelCellWidths[j] << "\n";
                    }
                }

                for(int i = 0; i < scvf.pairData().size(); i++)
                {
                    for(int j = 0; j < gridGeometry.upwindStencilOrder(); j++)
                    {
                        std::cout << " | Cell Centered Parallel Distance "<< j << " on axis " << i << ": "<< std::setw(3) << scvf.parallelDofsDistance(i,j) << "\n";
                    }
                }
                std::cout << std::endl;
                std::cout << std::endl;
            }
        }
    }
}
