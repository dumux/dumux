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
 * \ingroup CCMpfaDiscretization
 * \brief Properties for all models using cell-centered finite volume scheme with mpfa
 * \note Inherit from these properties to use a cell-centered finite volume scheme with mpfa
 */
#ifndef DUMUX_CC_MPFA_PROPERTIES_HH
#define DUMUX_CC_MPFA_PROPERTIES_HH

#include <dune/common/fvector.hh>
#include <dune/geometry/multilineargeometry.hh>

#include <dumux/common/properties.hh>

#include <dumux/assembly/cclocalresidual.hh>

#include <dumux/discretization/methods.hh>
#include <dumux/discretization/fvproperties.hh>

#include <dumux/discretization/cellcentered/gridvolumevariables.hh>
#include <dumux/discretization/cellcentered/subcontrolvolume.hh>
#include <dumux/discretization/cellcentered/elementsolution.hh>
#include <dumux/discretization/cellcentered/elementboundarytypes.hh>

#include <dumux/discretization/cellcentered/mpfa/connectivitymap.hh>
#include <dumux/discretization/cellcentered/mpfa/fvgridgeometry.hh>
#include <dumux/discretization/cellcentered/mpfa/gridfluxvariablescache.hh>
#include <dumux/discretization/cellcentered/mpfa/fvelementgeometry.hh>
#include <dumux/discretization/cellcentered/mpfa/elementvolumevariables.hh>
#include <dumux/discretization/cellcentered/mpfa/elementfluxvariablescache.hh>
#include <dumux/discretization/cellcentered/mpfa/subcontrolvolumeface.hh>
#include <dumux/discretization/cellcentered/mpfa/dualgridindexset.hh>
#include <dumux/discretization/cellcentered/mpfa/helper.hh>

#include <dumux/discretization/cellcentered/mpfa/omethod/interactionvolume.hh>

namespace Dumux
{
namespace Properties
{
//! Type tag for the cell-centered mpfa scheme.
NEW_TYPE_TAG(CCMpfaModel, INHERITS_FROM(FiniteVolumeModel));

//! Set the corresponding discretization method property
SET_PROP(CCMpfaModel, DiscretizationMethod)
{
    static const DiscretizationMethods value = DiscretizationMethods::CCMpfa;
};

//! Extract the used mpfa method from the primary interaction volume
SET_PROP(CCMpfaModel, MpfaMethod)
{
    static const MpfaMethods value = GET_PROP_TYPE(TypeTag, PrimaryInteractionVolume)::MpfaMethod;
};

//! Set the maximum admissible number of branches per scvf
SET_PROP(CCMpfaModel, MaxNumBranchesPerScvf)
{
private:
    using GridView = typename GET_PROP_TYPE(TypeTag, GridView);
    static constexpr int dim = GridView::dimension;
    static constexpr int dimWorld = GridView::dimensionworld;

public:
    // Per default, we allow for 8 neighbors on network/surface grids
    static const std::size_t value = dim < dimWorld ? 9 : 2;
};

//! Set the index set type used on the dual grid nodes
SET_PROP(CCMpfaModel, DualGridNodalIndexSet)
{
private:
    using GridView = typename GET_PROP_TYPE(TypeTag, GridView);

    // Use the grid view's index type for grid indices
    using GI = typename GridView::IndexSet::IndexType;

    // per default, use uint8_t as iv-local index type
    using LI = std::uint8_t;

    // the specified maximum admissible number of branches per scvf
    static constexpr int maxB = GET_PROP_VALUE(TypeTag, MaxNumBranchesPerScvf);

    // maximum admissible number of elements around a node
    // if for a given grid this number is still not high enough,
    // overwrite this property in your problem with a higher number
    static constexpr int dim = GridView::dimension;
    static constexpr int maxE = dim == 3 ? 45 : 15;

public:
    using type = CCMpfaDualGridNodalIndexSet<GI, LI, dim, maxE, maxB>;
};

//! The mpfa helper class
SET_TYPE_PROP(CCMpfaModel, MpfaHelper, CCMpfaHelper<TypeTag>);

//! Per default, we use the dynamic mpfa-o interaction volume
SET_PROP(CCMpfaModel, PrimaryInteractionVolume)
{
private:
    //! use the default traits
    using Traits = CCMpfaODefaultInteractionVolumeTraits< TypeTag >;
public:
    using type = CCMpfaOInteractionVolume< Traits >;
};

//! Per default, we use the dynamic mpfa-o interaction volume as secondary type
SET_PROP(CCMpfaModel, SecondaryInteractionVolume)
{
private:
    //! use the default traits
    using Traits = CCMpfaODefaultInteractionVolumeTraits< TypeTag >;
public:
    using type = CCMpfaOInteractionVolume< Traits >;
};

//! Set the default for the global finite volume geometry
SET_TYPE_PROP(CCMpfaModel,
              FVGridGeometry,
              CCMpfaFVGridGeometry<TypeTag, GET_PROP_VALUE(TypeTag, EnableFVGridGeometryCache)>);

//! Set the default for the local finite volume geometry
SET_TYPE_PROP(CCMpfaModel,
              FVElementGeometry,
              CCMpfaFVElementGeometry<TypeTag, GET_PROP_VALUE(TypeTag, EnableFVGridGeometryCache)>);

//! The global flux variables cache vector class
SET_TYPE_PROP(CCMpfaModel,
              GridFluxVariablesCache,
              CCMpfaGridFluxVariablesCache<TypeTag, GET_PROP_VALUE(TypeTag, EnableGridFluxVariablesCache)>);

//! The local flux variables cache vector class
SET_TYPE_PROP(CCMpfaModel,
              ElementFluxVariablesCache,
              CCMpfaElementFluxVariablesCache<TypeTag, GET_PROP_VALUE(TypeTag, EnableGridFluxVariablesCache)>);


//! The global previous volume variables vector class
SET_TYPE_PROP(CCMpfaModel,
              ElementVolumeVariables,
              CCMpfaElementVolumeVariables<TypeTag, GET_PROP_VALUE(TypeTag, EnableGridVolumeVariablesCache)>);

//! The global current volume variables vector class
SET_TYPE_PROP(CCMpfaModel, GridVolumeVariables, CCGridVolumeVariables<TypeTag, GET_PROP_VALUE(TypeTag, EnableGridVolumeVariablesCache)>);

//! The sub control volume
SET_PROP(CCMpfaModel, SubControlVolume)
{
private:
    using Grid = typename GET_PROP_TYPE(TypeTag, Grid);
    struct ScvGeometryTraits
    {
        using Geometry = typename Grid::template Codim<0>::Geometry;
        using GridIndexType = typename Grid::LeafGridView::IndexSet::IndexType;
        using LocalIndexType = unsigned int;
        using Scalar = typename Geometry::ctype;
        using GlobalPosition = Dune::FieldVector<Scalar, Geometry::coorddimension>;
    };
public:
    using type = CCSubControlVolume<ScvGeometryTraits>;
};

//! The sub-control volume face class
SET_PROP(CCMpfaModel, SubControlVolumeFace)
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
        template< int dim >
        struct hasSingleGeometryType
        {
            static const bool v = true;
            static const unsigned int topologyId = Dune::Impl::CubeTopology< dim >::type::id;
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
    };

public:
    using type = Dumux::CCMpfaSubControlVolumeFace< ScvfGeometryTraits>;
};

//! Set the solution vector type for an element
SET_TYPE_PROP(CCMpfaModel, ElementSolutionVector, CCElementSolution<TypeTag>);

//! Set the default for the ElementBoundaryTypes
SET_TYPE_PROP(CCMpfaModel, ElementBoundaryTypes, CCElementBoundaryTypes<TypeTag>);

//! Set the BaseLocalResidual to CCLocalResidual
SET_TYPE_PROP(CCMpfaModel, BaseLocalResidual, CCLocalResidual<TypeTag>);

//! Set the AssemblyMap property
SET_TYPE_PROP(CCMpfaModel, AssemblyMap, Dumux::CCMpfaConnectivityMap<TypeTag>);

} // namespace Properties
} // namespace Dumux

#endif
