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
 * \ingroup Properties
 * \file
 *
 * \brief Defines a type tag and some properties for models using the box scheme.
 */

#ifndef DUMUX_BOX_PROPERTIES_HH
#define DUMUX_BOX_PROPERTIES_HH

#include <dune/common/fvector.hh>
#include <dune/geometry/multilineargeometry.hh>

#include <dumux/common/properties.hh>
#include <dumux/common/boundaryflag.hh>

#include <dumux/assembly/boxlocalresidual.hh>

#include <dumux/discretization/methods.hh>
#include <dumux/discretization/fvproperties.hh>

#include <dumux/discretization/box/subcontrolvolume.hh>
#include <dumux/discretization/box/subcontrolvolumeface.hh>
#include <dumux/discretization/box/elementsolution.hh>
#include <dumux/discretization/box/elementboundarytypes.hh>
#include <dumux/discretization/box/gridfluxvariablescache.hh>
#include <dumux/discretization/box/elementfluxvariablescache.hh>
#include <dumux/discretization/box/gridvolumevariables.hh>
#include <dumux/discretization/box/elementvolumevariables.hh>
#include <dumux/discretization/box/fvgridgeometry.hh>
#include <dumux/discretization/box/fvelementgeometry.hh>

namespace Dumux
{
namespace Properties
{
//! Type tag for the box scheme.
NEW_TYPE_TAG(BoxModel, INHERITS_FROM(FiniteVolumeModel));

//! Set the corresponding discretization method property
SET_PROP(BoxModel, DiscretizationMethod)
{
    static const DiscretizationMethods value = DiscretizationMethods::Box;
};

//! Set the default for the FVElementGeometry vector
SET_TYPE_PROP(BoxModel, FVGridGeometry, BoxFVGridGeometry<TypeTag,
                            GET_PROP_VALUE(TypeTag, EnableFVGridGeometryCache)>);

//! Set the default for the FVElementGeometry vector
SET_TYPE_PROP(BoxModel, FVElementGeometry, BoxFVElementGeometry<TypeTag,
                            GET_PROP_VALUE(TypeTag, EnableFVGridGeometryCache)>);

//! The sub control volume
SET_PROP(BoxModel, SubControlVolume)
{
private:
    using Grid = typename GET_PROP_TYPE(TypeTag, Grid);
    static const int dim = Grid::dimension;
    static const int dimWorld = Grid::dimensionworld;

    // we use geometry traits that use static corner vectors to and a fixed geometry type
    template <class ct>
    struct ScvfMLGTraits : public Dune::MultiLinearGeometryTraits<ct>
    {
        // we use static vectors to store the corners as we know
        // the number of corners in advance (2^(dim) corners (1<<(dim))
        template< int mydim, int cdim >
        struct CornerStorage
        {
            using Type = std::array< Dune::FieldVector< ct, cdim >, (1<<(dim)) >;
        };

        // we know all scvfs will have the same geometry type
        template< int mydim >
        struct hasSingleGeometryType
        {
            static const bool v = true;
            static const unsigned int topologyId = Dune::Impl::CubeTopology< mydim >::type::id;
        };
    };

    struct ScvGeometryTraits
    {
        using GridIndexType = typename Grid::LeafGridView::IndexSet::IndexType;
        using LocalIndexType = unsigned int;
        using Scalar = typename Grid::ctype;
        using Geometry = Dune::MultiLinearGeometry<Scalar, dim, dimWorld, ScvfMLGTraits<Scalar>>;
        using CornerStorage = typename ScvfMLGTraits<Scalar>::template CornerStorage<dim, dimWorld>::Type;
        using GlobalPosition = typename CornerStorage::value_type;
    };
public:
    using type = BoxSubControlVolume<ScvGeometryTraits>;
};

SET_PROP(BoxModel, SubControlVolumeFace)
{
private:
    using Grid = typename GET_PROP_TYPE(TypeTag, Grid);
    static const int dim = Grid::dimension;
    static const int dimWorld = Grid::dimensionworld;

    // we use geometry traits that use static corner vectors to and a fixed geometry type
    template <class ct>
    struct ScvfMLGTraits : public Dune::MultiLinearGeometryTraits<ct>
    {
        // we use static vectors to store the corners as we know
        // the number of corners in advance (2^(dim-1) corners (1<<(dim-1))
        template< int mydim, int cdim >
        struct CornerStorage
        {
            using Type = std::array< Dune::FieldVector< ct, cdim >, (1<<(dim-1)) >;
        };

        // we know all scvfs will have the same geometry type
        template< int mydim >
        struct hasSingleGeometryType
        {
            static const bool v = true;
            static const unsigned int topologyId = Dune::Impl::CubeTopology< mydim >::type::id;
        };
    };

    struct ScvfGeometryTraits
    {
        using GridIndexType = typename Grid::LeafGridView::IndexSet::IndexType;
        using LocalIndexType = unsigned int;
        using Scalar = typename Grid::ctype;
        using Geometry = Dune::MultiLinearGeometry<Scalar, dim-1, dimWorld, ScvfMLGTraits<Scalar>>;
        using CornerStorage = typename ScvfMLGTraits<Scalar>::template CornerStorage<dim-1, dimWorld>::Type;
        using GlobalPosition = typename CornerStorage::value_type;
        using BoundaryFlag = Dumux::BoundaryFlag<Grid>;
    };
public:
    using type = BoxSubControlVolumeFace<ScvfGeometryTraits>;
};

//! Set the solution vector type for an element
SET_TYPE_PROP(BoxModel, ElementSolutionVector, BoxElementSolution<TypeTag>);

//! Set the default for the ElementBoundaryTypes
SET_TYPE_PROP(BoxModel, ElementBoundaryTypes, BoxElementBoundaryTypes<TypeTag>);

//! The global volume variables vector class
SET_TYPE_PROP(BoxModel, GridVolumeVariables, BoxGridVolumeVariables<TypeTag,
                            GET_PROP_VALUE(TypeTag, EnableGridVolumeVariablesCache)>);

//! The element volume variables vector class
SET_TYPE_PROP(BoxModel, ElementVolumeVariables, BoxElementVolumeVariables<TypeTag,
                            GET_PROP_VALUE(TypeTag, EnableGridVolumeVariablesCache)>);

//! The global flux variables cache vector class
SET_TYPE_PROP(BoxModel, GridFluxVariablesCache, BoxGridFluxVariablesCache<TypeTag,
                            GET_PROP_VALUE(TypeTag, EnableGridFluxVariablesCache)>);

//! The local flux variables cache vector class
SET_TYPE_PROP(BoxModel, ElementFluxVariablesCache, BoxElementFluxVariablesCache<TypeTag,
                            GET_PROP_VALUE(TypeTag, EnableGridFluxVariablesCache)>);

//! Set the BaseLocalResidual to BoxLocalResidual
SET_TYPE_PROP(BoxModel, BaseLocalResidual, BoxLocalResidual<TypeTag>);

} // namespace Properties
} // namespace Dumux

#endif
