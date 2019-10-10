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
 * \ingroup ShallowWaterModel
 * \copydoc Dumux::ShallowWaterFluxVariables
 */
#ifndef DUMUX_FREEFLOW_SHALLOW_WATER_FLUXVARIABLES_HH
#define DUMUX_FREEFLOW_SHALLOW_WATER_FLUXVARIABLES_HH

#include <dumux/common/properties.hh>
#include <dumux/flux/fluxvariablesbase.hh>

namespace Dumux {

/*!
 * \ingroup ShallowWaterModel
 * \brief The flux variables class for the shallow water model.
 *
 */
template<class TypeTag>
class ShallowWaterFluxVariables
: public FluxVariablesBase<GetPropType<TypeTag, Properties::Problem>,
                           typename GetPropType<TypeTag, Properties::GridGeometry>::LocalView,
                           typename GetPropType<TypeTag, Properties::GridVolumeVariables>::LocalView,
                           typename GetPropType<TypeTag, Properties::GridFluxVariablesCache>::LocalView>
{
    using Problem = GetPropType<TypeTag, Properties::Problem>;
    using ModelTraits = GetPropType<TypeTag, Properties::ModelTraits>;
    using NumEqVector = GetPropType<TypeTag, Properties::NumEqVector>;
    using AdvectionType = GetPropType<TypeTag, Properties::AdvectionType>;
    //using DiffusionType = GetPropType<TypeTag, Properties::DiffusionType>;

    using GridVariables = GetPropType<TypeTag, Properties::GridVariables>;
    using GridVolumeVariables = typename GridVariables::GridVolumeVariables;
    using GridGeometry = GetPropType<TypeTag, Properties::GridGeometry>;
    using FVElementGeometry = typename GridGeometry::LocalView;
    using GridView = typename GridGeometry::GridView;

    using Element = typename GridView::template Codim<0>::Entity;
    using ElementVolumeVariables = typename GridVolumeVariables::LocalView;
    using SubControlVolumeFace = typename FVElementGeometry::SubControlVolumeFace;

    static constexpr bool enableAdvection = ModelTraits::enableAdvection();
    static constexpr bool enableDiffusion = ModelTraits::enableDiffusion();

public:

    /*!
     * \brief Returns the advective flux computed by the Riemann solver
     *
     */
    NumEqVector advectiveFlux(const Problem& problem,
                              const Element& element,
                              const FVElementGeometry& fvGeometry,
                              const ElementVolumeVariables& elemVolVars,
                              const SubControlVolumeFace& scvf) const
    {
        if (enableAdvection)
            return AdvectionType::flux(problem, element, fvGeometry, elemVolVars, scvf);

        return NumEqVector(0.0);
    }

    /*!
     * \brief Returns the diffusive flux (e.g. diffusion of tracer)
     *
     */
    NumEqVector diffusiveFlux(const Problem& problem,
                              const Element& element,
                              const FVElementGeometry& fvGeometry,
                              const ElementVolumeVariables& elemVolVars,
                              const SubControlVolumeFace& scvf) const
    {
        // TODO: add diffusive flux (e.g. tracer and viscosity)
        if (enableDiffusion)
            return NumEqVector(0.0);

        return NumEqVector(0.0);
    }
};

} // end namespace Dumux

#endif
