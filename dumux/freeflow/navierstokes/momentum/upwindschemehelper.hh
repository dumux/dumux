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
 * \ingroup NavierStokesModel
 * \copydoc Dumux::MomentumUpwindSchemeHelper
 */
#ifndef DUMUX_MOMENTUM_UPWIND_SCHEME_HELPER_HH
#define DUMUX_MOMENTUM_UPWIND_SCHEME_HELPER_HH

#include <array>
#include <optional>

#include <dumux/common/math.hh>
#include <dumux/common/exceptions.hh>
#include <dumux/common/parameters.hh>
#include <dumux/common/properties.hh>
#include <dumux/common/typetraits/problem.hh>

#include <dumux/discretization/method.hh>
#include <dumux/freeflow/staggeredupwindmethods.hh>
#include "velocitygradients.hh"

namespace Dumux {

/*!
 * \ingroup NavierStokesModel
 * \brief The upwinding variables class for the Navier-Stokes model using the staggered grid discretization.
 */
template<class TypeTag, int upwindSchemeOrder>
class FaceCenteredStaggeredUpwindHelper
{
    using Problem = GetPropType<TypeTag, Properties::Problem>;
    using GridVariables = GetPropType<TypeTag, Properties::GridVariables>;
    using GridGeometry = typename GridVariables::GridGeometry;
    using FVElementGeometry = typename GridGeometry::LocalView;
    using GridVolumeVariables = typename GridVariables::GridVolumeVariables;
    using ElementVolumeVariables = typename GridVolumeVariables::LocalView;
    using VolumeVariables = typename GridVolumeVariables::VolumeVariables;
    using GridView = typename GridGeometry::GridView;
    using Element = typename GridView::template Codim<0>::Entity;
    using SubControlVolume = typename FVElementGeometry::SubControlVolume;
    using SubControlVolumeFace = typename FVElementGeometry::SubControlVolumeFace;
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using GlobalPosition = typename Element::Geometry::GlobalCoordinate;
    using ElementBoundaryTypes = GetPropType<TypeTag, Properties::ElementBoundaryTypes>;

    using UpwindScheme = StaggeredUpwindMethods<Scalar, upwindSchemeOrder>;
    static constexpr bool useHigherOrder = upwindSchemeOrder > 1;
    static_assert(upwindSchemeOrder <= 2, "Not implemented: Order higher than 2!");

public:
    FaceCenteredStaggeredUpwindHelper(const Element& element,
                                      const FVElementGeometry& fvGeometry,
                                      const Problem& problem,
                                      const SubControlVolumeFace& scvf,
                                      const ElementVolumeVariables& elemVolVars,
                                      const ElementBoundaryTypes& elemBcTypes,
                                      const UpwindScheme& upwindScheme)
    : element_(element)
    , fvGeometry_(fvGeometry)
    , problem_(problem)
    , scvf_(scvf)
    , elemVolVars_(elemVolVars)
    , elemBcTypes_(elemBcTypes)
    , upwindScheme_(upwindScheme)
    {}

private:
    const Element& element_;
    const FVElementGeometry& fvGeometry_;
    const Problem& problem_;
    const SubControlVolumeFace& scvf_;
    const ElementVolumeVariables& elemVolVars_;
    const ElementBoundaryTypes& elemBcTypes_;
    const UpwindScheme& upwindScheme_;
};

} // end namespace Dumux

#endif
