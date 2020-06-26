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
 * \copydoc Dumux::NavierStokesFluxVariablesImpl
 */
#ifndef DUMUX_NAVIERSTOKES_MASS_AND_ENERGY_FLUXVARIABLES_HH
#define DUMUX_NAVIERSTOKES_MASS_AND_ENERGY_FLUXVARIABLES_HH

#include <array>
#include <optional>

#include <dumux/common/math.hh>
#include <dumux/common/exceptions.hh>
#include <dumux/common/parameters.hh>
#include <dumux/common/properties.hh>

#include <dumux/flux/fluxvariablesbase.hh>
#include <dumux/discretization/method.hh>


namespace Dumux {

/*!
 * \ingroup NavierStokesModel
 * \brief The flux variables class for the Navier-Stokes model using the staggered grid discretization.
 */
template<class TypeTag>
class NavierStokesMassAndEnergyFluxVariables
: public FluxVariablesBase<GetPropType<TypeTag, Properties::Problem>,
                           typename GetPropType<TypeTag, Properties::GridGeometry>::LocalView,
                           typename GetPropType<TypeTag, Properties::GridVolumeVariables>::LocalView,
                           typename GetPropType<TypeTag, Properties::GridFluxVariablesCache>::LocalView>
{
    using GridVariables = GetPropType<TypeTag, Properties::GridVariables>;

    using GridVolumeVariables = typename GridVariables::GridVolumeVariables;
    using ElementVolumeVariables = typename GridVolumeVariables::LocalView;
    using VolumeVariables = typename GridVolumeVariables::VolumeVariables;

    using GridFluxVariablesCache = typename GridVariables::GridFluxVariablesCache;
    using FluxVariablesCache = typename GridFluxVariablesCache::FluxVariablesCache;

    using GridGeometry = typename GridVariables::GridGeometry;
    using FVElementGeometry = typename GridGeometry::LocalView;
    using SubControlVolume = typename FVElementGeometry::SubControlVolume;
    using SubControlVolumeFace = typename FVElementGeometry::SubControlVolumeFace;


    using Problem = GetPropType<TypeTag, Properties::Problem>;
    using GridView = typename GridGeometry::GridView;
    using Element = typename GridView::template Codim<0>::Entity;
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using ModelTraits = GetPropType<TypeTag, Properties::ModelTraits>;
    using Indices = typename ModelTraits::Indices;

    using GlobalPosition = typename Element::Geometry::GlobalCoordinate;

    using NumEqVector = GetPropType<TypeTag, Properties::NumEqVector>;


public:


};

} // end namespace Dumux

#endif // DUMUX_NAVIERSTOKES_STAGGERED_FLUXVARIABLES_HH
