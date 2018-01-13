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
 * \ingroup CCTpfaDiscretization
 * \brief Properties for all models using cell-centered finite volume scheme with TPFA
 * \note Inherit from these properties to use a cell-centered finite volume scheme with TPFA
 */

#ifndef DUMUX_CC_TPFA_PROPERTIES_HH
#define DUMUX_CC_TPFA_PROPERTIES_HH

#include <dune/common/fvector.hh>
#include <dune/common/reservedvector.hh>
#include <dune/geometry/multilineargeometry.hh>

#include <dumux/common/properties.hh>
#include <dumux/common/boundaryflag.hh>

#include <dumux/assembly/cclocalresidual.hh>

#include <dumux/discretization/methods.hh>
#include <dumux/discretization/fvproperties.hh>

#include <dumux/discretization/cellcentered/gridvolumevariables.hh>
#include <dumux/discretization/cellcentered/subcontrolvolume.hh>

#include <dumux/discretization/cellcentered/elementboundarytypes.hh>
#include <dumux/discretization/cellcentered/connectivitymap.hh>
#include <dumux/discretization/cellcentered/elementsolution.hh>
#include <dumux/discretization/cellcentered/tpfa/fvgridgeometry.hh>
#include <dumux/discretization/cellcentered/tpfa/gridfluxvariablescache.hh>
#include <dumux/discretization/cellcentered/tpfa/fvelementgeometry.hh>
#include <dumux/discretization/cellcentered/tpfa/elementvolumevariables.hh>
#include <dumux/discretization/cellcentered/tpfa/elementfluxvariablescache.hh>
#include <dumux/discretization/cellcentered/tpfa/subcontrolvolumeface.hh>

namespace Dumux
{
namespace Properties
{
//! Type tag for the cell-centered tpfa scheme.
NEW_TYPE_TAG(CCTpfaModel, INHERITS_FROM(FiniteVolumeModel));

//! Set the corresponding discretization method property
SET_PROP(CCTpfaModel, DiscretizationMethod)
{
    static const DiscretizationMethods value = DiscretizationMethods::CCTpfa;
};

//! Set the default for the global finite volume geometry
SET_TYPE_PROP(CCTpfaModel, FVGridGeometry, CCTpfaFVGridGeometry<TypeTag, GET_PROP_VALUE(TypeTag, EnableFVGridGeometryCache)>);

//! The global flux variables cache vector class
SET_TYPE_PROP(CCTpfaModel, GridFluxVariablesCache, CCTpfaGridFluxVariablesCache<TypeTag, GET_PROP_VALUE(TypeTag, EnableGridFluxVariablesCache)>);

//! Set the default for the local finite volume geometry
SET_TYPE_PROP(CCTpfaModel, FVElementGeometry, CCTpfaFVElementGeometry<TypeTag, GET_PROP_VALUE(TypeTag, EnableFVGridGeometryCache)>);

//! The global previous volume variables vector class
SET_TYPE_PROP(CCTpfaModel, ElementVolumeVariables, CCTpfaElementVolumeVariables<TypeTag, GET_PROP_VALUE(TypeTag, EnableGridVolumeVariablesCache)>);

//! The local flux variables cache vector class
SET_TYPE_PROP(CCTpfaModel, ElementFluxVariablesCache, CCTpfaElementFluxVariablesCache<TypeTag, GET_PROP_VALUE(TypeTag, EnableGridFluxVariablesCache)>);

//! The global current volume variables vector class
SET_TYPE_PROP(CCTpfaModel, GridVolumeVariables, CCGridVolumeVariables<TypeTag, GET_PROP_VALUE(TypeTag, EnableGridVolumeVariablesCache)>);

//! Set the maximum admissible number of branches per scvf
SET_PROP(CCTpfaModel, MaxNumNeighborsPerScvf)
{
private:
    using GridView = typename GET_PROP_TYPE(TypeTag, GridView);
    static constexpr int dim = GridView::dimension;
    static constexpr int dimWorld = GridView::dimensionworld;

public:
    // Per default, we allow for 8 neighbors on network/surface grids
    static constexpr std::size_t value = dim < dimWorld ? 9 : 2;
};

//! The sub control volume
SET_PROP(CCTpfaModel, SubControlVolume)
{
private:
    using Grid = typename GET_PROP_TYPE(TypeTag, Grid);
    struct ScvGeometryTraits
    {
        using Geometry = typename Grid::template Codim<0>::Geometry;
        using GridIndexType = typename Grid::LeafGridView::IndexSet::IndexType;
        using LocalIndexType = unsigned int;
        using Scalar = typename Grid::ctype;
        using GlobalPosition = Dune::FieldVector<Scalar, Grid::dimensionworld>;
    };
public:
    using type = CCSubControlVolume<ScvGeometryTraits>;
};

//! The sub control volume face
SET_PROP(CCTpfaModel, SubControlVolumeFace)
{
private:
    using Grid = typename GET_PROP_TYPE(TypeTag, Grid);
    static constexpr int dim = Grid::dimension;
    static constexpr int dimWorld = Grid::dimensionworld;

    // we use geometry traits that use static corner vectors to and a fixed geometry type
    template <class ct>
    struct ScvfMLGTraits : public Dune::MultiLinearGeometryTraits<ct>
    {
        // we use static vectors to store the corners as we know
        // the number of corners in advance (2^(dim-1) corners (1<<(dim-1))
        template< int mydim, int cdim >
        struct CornerStorage
        {
            using Type = Dune::ReservedVector< Dune::FieldVector< ct, cdim >, (1<<(dim-1)) >;
        };
    };

    struct ScvfGeometryTraits
    {
        using GridIndexType = typename Grid::LeafGridView::IndexSet::IndexType;
        using LocalIndexType = unsigned int;
        using Scalar = typename Grid::ctype;
        using Geometry = Dune::MultiLinearGeometry<Scalar, dim-1, dimWorld, ScvfMLGTraits<Scalar> >;
        using CornerStorage = typename ScvfMLGTraits<Scalar>::template CornerStorage<dim-1, dimWorld>::Type;
        using GlobalPosition = typename CornerStorage::value_type;
        using BoundaryFlag = Dumux::BoundaryFlag<Grid>;
    };
public:
    using type = Dumux::CCTpfaSubControlVolumeFace<ScvfGeometryTraits>;
};

//! Set the solution vector type for an element
SET_TYPE_PROP(CCTpfaModel, ElementSolutionVector, CCElementSolution<TypeTag>);

//! Set the default for the ElementBoundaryTypes
SET_TYPE_PROP(CCTpfaModel, ElementBoundaryTypes, CCElementBoundaryTypes<TypeTag>);

//! Set the BaseLocalResidual to CCLocalResidual
SET_TYPE_PROP(CCTpfaModel, BaseLocalResidual, CCLocalResidual<TypeTag>);

//! Set the AssemblyMap to the SimpleAssemblyMap per default
SET_TYPE_PROP(CCTpfaModel, AssemblyMap, CCSimpleConnectivityMap<TypeTag>);

} // namespace Properties
} // namespace Dumux

#endif
