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
 * \ingroup ShallowWaterModel
 * \brief A helper class to have a dummy filler for the flux variable cache
 */
#ifndef DUMUX_FREEFLOW_SHALLOW_WATER_FLUXVARSCACHE_FILLER_HH
#define DUMUX_FREEFLOW_SHALLOW_WATER_FLUXVARSCACHE_FILLER_HH

#include <dumux/common/properties.hh>
#include <dumux/discretization/method.hh>

namespace Dumux {

/*!
* \ingroup ShallowWaterModel
* \brief A helper class to have a dummy filler for the flux variables cache
*/
template<class TypeTag>
class ShallowWaterFluxVariablesCacheFiller
{
    using ModelTraits = GetPropType<TypeTag, Properties::ModelTraits>;
    using Problem = GetPropType<TypeTag, Properties::Problem>;
    using GridView = GetPropType<TypeTag, Properties::GridView>;
    using FVElementGeometry = typename GetPropType<TypeTag, Properties::FVGridGeometry>::LocalView;
    using SubControlVolume = typename FVElementGeometry::SubControlVolume;
    using SubControlVolumeFace = typename FVElementGeometry::SubControlVolumeFace;
    using ElementVolumeVariables = typename GetPropType<TypeTag, Properties::GridVolumeVariables>::LocalView;
    using FluxVariablesCache = GetPropType<TypeTag, Properties::FluxVariablesCache>;

    using Element = typename GridView::template Codim<0>::Entity;

public:
    static constexpr bool isSolDependent = (true);


    //! The constructor. Sets the problem pointer
    ShallowWaterFluxVariablesCacheFiller(const Problem& problem) : problemPtr_(&problem) {}

    /*!
     * \brief function to fill the flux variables caches
     *
     * \param fluxVarsCacheContainer Either the element or global flux variables cache
     * \param scvfFluxVarsCache The flux var cache to be updated corresponding to the given scvf
     * \param element The finite element
     * \param fvGeometry The finite volume geometry
     * \param elemVolVars The element volume variables
     * \param scvf The corresponding sub-control volume face
     * \param forceUpdateAll if true, forces all caches to be updated (even the solution-independent ones)
     */
    template<class FluxVariablesCacheContainer>
    void fill(FluxVariablesCacheContainer& fluxVarsCacheContainer,
              FluxVariablesCache& scvfFluxVarsCache,
              const Element& element,
              const FVElementGeometry& fvGeometry,
              const ElementVolumeVariables& elemVolVars,
              const SubControlVolumeFace& scvf,
              bool forceUpdateAll = false)
    {
    }

private:

    const Problem& problem() const
    { return *problemPtr_; }

    const Problem* problemPtr_;
};

} // end namespace Dumux

#endif
