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

#include <dune/common/test/iteratortest.hh>
#include <dune/grid/utility/structuredgridfactory.hh>
#include <dune/grid/yaspgrid.hh>

#include <dumux/common/basicproperties.hh>
#include <dumux/implicit/tpfa/fvelementgeometryvector.hh>
#include <dumux/implicit/fvelementgeometry.hh>
#include <dumux/implicit/subcontrolvolume.hh>
#include <dumux/implicit/subcontrolvolumeface.hh>

namespace Dumux
{
namespace Properties
{
NEW_PROP_TAG(SubControlVolume);
NEW_PROP_TAG(SubControlVolumeFace);
NEW_PROP_TAG(FVElementGeometry);
NEW_PROP_TAG(FVElementGeometryVector);

NEW_TYPE_TAG(TestFVGeometry, INHERITS_FROM(NumericModel));

SET_TYPE_PROP(TestFVGeometry, Grid, Dune::YaspGrid<2>);

SET_TYPE_PROP(TestFVGeometry, GridView, typename GET_PROP_TYPE(TypeTag, Grid)::LeafGridView);

SET_PROP(TestFVGeometry, SubControlVolume)
{
private:
    using Grid = typename GET_PROP_TYPE(TypeTag, Grid);
    using ScvGeometry = typename Grid::template Codim<0>::Geometry;
    using IndexType = typename Grid::LeafGridView::IndexSet::IndexType;
public:
    typedef Dumux::SubControlVolume<ScvGeometry, IndexType> type;
};

SET_PROP(TestFVGeometry, SubControlVolumeFace)
{
private:
    using Grid = typename GET_PROP_TYPE(TypeTag, Grid);
    using ScvfGeometry = typename Grid::template Codim<1>::Geometry;
    using IndexType = typename Grid::LeafGridView::IndexSet::IndexType;
public:
    typedef Dumux::SubControlVolumeFace<ScvfGeometry, IndexType> type;

};

SET_TYPE_PROP(TestFVGeometry, FVElementGeometry, FVElementGeometry<TypeTag>);

SET_TYPE_PROP(TestFVGeometry, FVElementGeometryVector, TpfaFVElementGeometryVector<TypeTag>);

}
}

template<class T>
class NoopFunctor {
public:
  NoopFunctor() {}
  void operator()(const T& t){}
};

int main (int argc, char *argv[]) try
{
    // maybe initialize mpi
    Dune::MPIHelper::instance(argc, argv);

    std::cout << "Checking the FVGeometries, SCVs and SCV faces" << std::endl;

    // aliases
    using TypeTag = TTAG(TestFVGeometry);
    using Grid = typename GET_PROP_TYPE(TypeTag, Grid);
    using GridView = typename Grid::LeafGridView;

    constexpr int dim = GridView::dimension;
    constexpr int dimworld = GridView::dimensionworld;

    using Scalar = typename GET_PROP_TYPE(TypeTag, Scalar);
    using GlobalPosition = Dune::FieldVector<Scalar, dimworld>;
    using SubControlVolume = typename GET_PROP_TYPE(TypeTag, SubControlVolume);
    using SubControlVolumeFace = typename GET_PROP_TYPE(TypeTag, SubControlVolumeFace);
    using FVElementGeometryVector = typename GET_PROP_TYPE(TypeTag, FVElementGeometryVector);

    // make a grid
    GlobalPosition lower(0.0);
    GlobalPosition upper(1.0);
    std::array<unsigned int, dim> els{{2, 2}};
    std::shared_ptr<Grid> grid = Dune::StructuredGridFactory<Grid>::createCubeGrid(lower, upper, els);
    auto leafGridView = grid->leafGridView();

    FVElementGeometryVector fvGeometries(leafGridView);
    fvGeometries.update();

    // iterate over elements. For every element get fv geometry and loop over scvs and scvfaces
    for (const auto& element : elements(leafGridView))
    {
        std::cout << std::endl << "Checking fvGeometry of element " << leafGridView.indexSet().index(element) << std::endl;
        auto fvGeometry = fvGeometries.fvGeometry(element);

        auto range = fvGeometry.scvs();
        NoopFunctor<SubControlVolume> op;
        if(0 != testForwardIterator(range.begin(), range.end(), op))
            DUNE_THROW(Dune::Exception, "Iterator does not fulfill the forward iterator concept");

        for (auto&& scv : fvGeometry.scvs())
        {
            std::cout << "-- scv center at: " << scv.center() << std::endl;
        }

        auto range2 = fvGeometry.scvfs();
        NoopFunctor<SubControlVolumeFace> op2;
        if(0 != testForwardIterator(range2.begin(), range2.end(), op2))
            DUNE_THROW(Dune::Exception, "Iterator does not fulfill the forward iterator concept");

        for (auto&& scvf : fvGeometry.scvfs())
        {
            std::cout << "-- scvf center at: " << scvf.center() << std::endl;
        }
    }
}
// //////////////////////////////////
//   Error handler
// /////////////////////////////////
catch (Dune::Exception e) {

    std::cout << e << std::endl;
    return 1;
}
