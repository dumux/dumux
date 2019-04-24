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

namespace Dumux{

/*!
 * \ingroup ShallowWaterModel
 * \brief The flux variables class for the shallow water model.
 *
 */
template<class TypeTag>
class ShallowWaterFluxVariables
: public FluxVariablesBase<GetPropType<TypeTag,
                           Properties::Problem>,
                           typename GetPropType<TypeTag, Properties::FVGridGeometry>::LocalView,
                           typename GetPropType<TypeTag, Properties::GridVolumeVariables>::LocalView,
                           typename GetPropType<TypeTag, Properties::GridFluxVariablesCache>::LocalView>
{
    //using ParentType = FluxVariablesBase<TypeTag>;
    using Problem = GetPropType<TypeTag, Properties::Problem>;
    using ModelTraits = GetPropType<TypeTag, Properties::ModelTraits>;
    using NumEqVector = GetPropType<TypeTag, Properties::NumEqVector>;
    using AdvectionType = GetPropType<TypeTag, Properties::AdvectionType>;
    //using DiffusionType = GetPropType<TypeTag, Properties::DiffusionType>;

    static constexpr bool enableAdvection = ModelTraits::enableAdvection();
    static constexpr bool enableDiffusion = ModelTraits::enableDiffusion();


public:

    /*!
     * \brief Returns the advective flux computed by the Riemann solver
     *
     */
    NumEqVector advectiveFlux() const
    {

        NumEqVector fluxVector(0.0);
        if (enableAdvection)
        {

             return AdvectionType::flux(this->problem(),
                                        this->element(),
                                        this->fvGeometry(),
                                        this->elemVolVars(),
                                        this->scvFace(),
                                        this->elemFluxVarsCache());
        }
        else
        {
            return fluxVector;
        }
    }

    /*!
     * \brief Returns the diffusive flux (e.g. diffusion of tracer)
     *
     */
    NumEqVector diffusiveFlux() const
    {
        NumEqVector fluxVector(0.0);
        if (enableDiffusion)
        {
             /*return DiffusionType::flux(this->problem(),
                                          this->element(),
                                          this->fvGeometry(),
                                          this->elemVolVars(),
                                          this->scvFace(),
                                          this->elemFluxVarsCache());
            */
            return fluxVector;
        }

        else
        {
            return fluxVector;
        }
    }
};

} // end namespace Dumux

#endif
