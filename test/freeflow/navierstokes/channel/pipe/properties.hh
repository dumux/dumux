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

#ifndef DUMUX_TEST_FREEFLOW_PIPE_PROPERTIES_HH
#define DUMUX_TEST_FREEFLOW_PIPE_PROPERTIES_HH

#include <dune/grid/yaspgrid.hh>

#include <dumux/common/properties.hh>

#include <dumux/discretization/staggered/freeflow/properties.hh>
#include <dumux/discretization/rotationsymmetricgridgeometrytraits.hh>
#include <dumux/freeflow/navierstokes/model.hh>
#include <dumux/material/fluidsystems/1pgas.hh>
#include <dumux/material/components/air.hh>

#include "problem.hh"

namespace Dumux::Properties {

// Create new type tags
namespace TTag {
struct PipeFlow { using InheritsFrom = std::tuple<NavierStokes, StaggeredFreeFlowModel>; };
} // end namespace TTag

// the fluid system
template<class TypeTag>
struct FluidSystem<TypeTag, TTag::PipeFlow>
{
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using type = FluidSystems::OnePGas<Scalar, Dumux::Components::Air<Scalar> > ;
};

// Set the grid type
template<class TypeTag>
struct Grid<TypeTag, TTag::PipeFlow>
{ using type = Dune::YaspGrid<2, Dune::TensorProductCoordinates<GetPropType<TypeTag, Properties::Scalar>, 2> >; };

// Set the problem property
template<class TypeTag>
struct Problem<TypeTag, TTag::PipeFlow>
{ using type = FreeFlowPipeProblem<TypeTag> ; };

template<class TypeTag>
struct EnableGridGeometryCache<TypeTag, TTag::PipeFlow> { static constexpr bool value = true; };
template<class TypeTag>
struct EnableGridFluxVariablesCache<TypeTag, TTag::PipeFlow> { static constexpr bool value = true; };
template<class TypeTag>
struct EnableGridVolumeVariablesCache<TypeTag, TTag::PipeFlow> { static constexpr bool value = true; };

// rotation-symmetric grid geometry forming a cylinder channel
template<class TypeTag>
struct GridGeometry<TypeTag, TTag::PipeFlow>
{
    static constexpr auto upwindSchemeOrder = getPropValue<TypeTag, Properties::UpwindSchemeOrder>();
    static constexpr bool enableCache = getPropValue<TypeTag, Properties::EnableGridGeometryCache>();
    using GridView = typename GetPropType<TypeTag, Properties::Grid>::LeafGridView;

    using DefaultTraits = StaggeredFreeFlowDefaultFVGridGeometryTraits<GridView, upwindSchemeOrder>;

    struct GGTraits : public DefaultTraits
    {
        using SubControlVolume = RotationSymmetricSubControlVolume<typename DefaultTraits::SubControlVolume, RotationPolicy::toroid>;
        using SubControlVolumeFace = RotationSymmetricSubControlVolumeFace<typename DefaultTraits::SubControlVolumeFace, RotationPolicy::toroid>;

        struct PublicTraits
        {
            using CellSubControlVolume = SubControlVolume;
            using CellSubControlVolumeFace = SubControlVolumeFace;
            using FaceSubControlVolume = RotationSymmetricSubControlVolume<typename DefaultTraits::PublicTraits::FaceSubControlVolume, RotationPolicy::toroid>;
            using FaceLateralSubControlVolumeFace = RotationSymmetricSubControlVolumeFace<typename DefaultTraits::PublicTraits::FaceLateralSubControlVolumeFace, RotationPolicy::toroid>;
            using FaceFrontalSubControlVolumeFace = RotationSymmetricSubControlVolumeFace<typename DefaultTraits::PublicTraits::FaceFrontalSubControlVolumeFace, RotationPolicy::toroid>;
        };
    };

    using type = StaggeredFVGridGeometry<GridView, enableCache, GGTraits>;
};

} // end namespace Dumux::Properties

#endif
