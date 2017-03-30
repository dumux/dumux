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
 *
 * \brief Contains the supplemental quantities, which are constant within a
 *        finite volume in the non-isothermal Stokes staggered grid model.
 */
#ifndef DUMUX_STAGGEREDNI_VOLUME_VARIABLES_HH
#define DUMUX_STAGGEREDNI_VOLUME_VARIABLES_HH

#include <dumux/freeflow/staggered/volumevariables.hh>
#include "properties.hh"

namespace Dumux
{

/*!
 * \ingroup NavierStokesModel
 * \ingroup ImplicitVolumeVariables
 * \brief Contains the quantities which are are constant within a
 *        finite volume in the non-isothermal Stokes staggered grid model.
 */
template <class TypeTag>
class NavierStokesNIVolumeVariables : public NavierStokesVolumeVariables<TypeTag>
{
    using ParentType = ImplicitVolumeVariables<TypeTag>;
    using Implementation = typename GET_PROP_TYPE(TypeTag, VolumeVariables);
    using Scalar = typename GET_PROP_TYPE(TypeTag, Scalar);
    using Problem = typename GET_PROP_TYPE(TypeTag, Problem);
    using FVElementGeometry = typename GET_PROP_TYPE(TypeTag, FVElementGeometry);
    using SubControlVolume = typename GET_PROP_TYPE(TypeTag, SubControlVolume);
    using PrimaryVariables = typename GET_PROP_TYPE(TypeTag, PrimaryVariables);
    using ElementSolutionVector = typename GET_PROP_TYPE(TypeTag, ElementSolutionVector);
    using Indices = typename GET_PROP_TYPE(TypeTag, Indices);
    using FluidSystem = typename GET_PROP_TYPE(TypeTag, FluidSystem);
    using FluidState = typename GET_PROP_TYPE(TypeTag, FluidState);
    using GridView = typename GET_PROP_TYPE(TypeTag, GridView);
    using Element = typename GridView::template Codim<0>::Entity;
    using CellCenterPrimaryVariables = typename GET_PROP_TYPE(TypeTag, CellCenterPrimaryVariables);

    enum { numPhases = FluidSystem::numPhases,
           phaseIdx = Indices::phaseIdx,
           temperatureIdx = Indices::temperatureIdx
    };

    static constexpr bool useMoles = GET_PROP_VALUE(TypeTag, UseMoles);

public:
    /*!
     * \copydoc ImplicitVolumeVariables::update()
     */
    void update(const ElementSolutionVector &elemSol,
                const Problem &problem,
                const Element &element,
                const SubControlVolume& scv)
    {
        ParentType::update(elemSol, problem, element, scv);

        completeFluidState(elemSol, problem, element, scv, fluidState_);
    }

    /*!
     * \copydoc ImplicitModel::completeFluidState
     */
    static void completeFluidState(const ElementSolutionVector& elemSol,
                                   const Problem& problem,
                                   const Element& element,
                                   const SubControlVolume& scv,
                                   FluidState& fluidState)
    {
        fluidState.setTemperature(elemSol[0][Indices::temperatureIdx]);
        fluidState.setPressure(/*phaseIdx=*/0, elemSol[0][Indices::pressureIdx]);

        // saturation in a single phase is always 1 and thus redundant
        // to set. But since we use the fluid state shared by the
        // immiscible multi-phase models, so we have to set it here...
        fluidState.setSaturation(/*phaseIdx=*/0, 1.0);

        typename FluidSystem::ParameterCache paramCache;
        paramCache.updatePhase(fluidState, /*phaseIdx=*/0);

        Scalar value = FluidSystem::density(fluidState, paramCache, /*phaseIdx=*/0);
        fluidState.setDensity(/*phaseIdx=*/0, value);

        value = FluidSystem::viscosity(fluidState, paramCache, /*phaseIdx=*/0);
        fluidState.setViscosity(/*phaseIdx=*/0, value);

        // compute and set the enthalpy
        value = FluidSystem::enthalpy(fluidState, paramCache, /*phaseIdx=*/0);
        fluidState.setEnthalpy(/*phaseIdx=*/0, value);
    }

    /*!
     * \brief Return the mass density \f$\mathrm{[kg/m^3]}\f$ of a given phase within the
     *        control volume.
     */
    Scalar density(int phaseIdx = 0) const
    { return this->fluidState_.density(phaseIdx); }

    /*!
    * \brief Returns the mass density of a given phase within the
    *        control volume.
    *
    * \param phaseIdx The phase index
    */
   Scalar molarDensity() const
   {
       return this->fluidState_.molarDensity(phaseIdx);
   }

    /*!
     * \brief Returns the total internal energy of the fluid phase in the
     *        sub-control volume.
     */
    Scalar internalEnergy() const
    { return this->fluidState_.internalEnergy(phaseIdx); }

    /*!
     * \brief Returns the total enthalpy of the fluid phase in the sub-control
     *        volume.
     */
    Scalar enthalpy() const
    { return this->fluidState_.enthalpy(phaseIdx); }

    /*!
     * \brief Return the specific isobaric heat capacity \f$\mathrm{[J/(kg*K)]}\f$
     *        in the sub-control volume.
     */
    Scalar heatCapacity() const
    { return FluidSystem::heatCapacity(this->fluidState_, phaseIdx); }

    /*!
     * \brief Returns the thermal conductivity \f$\mathrm{[W/(m*K)]}\f$
     *        of the fluid phase in the sub-control volume.
     */
    Scalar thermalConductivity() const
    { return FluidSystem::thermalConductivity(this->fluidState_, phaseIdx); }


protected:
    FluidState fluidState_;

    // this method gets called by the parent class. since this method
    // is protected, we are friends with our parent...
    friend class NavierStokesVolumeVariables<TypeTag>;

    static Scalar temperature_(const CellCenterPrimaryVariables &priVars,
                            const Problem& problem,
                            const Element &element,
                            const FVElementGeometry &fvGeometry,
                            const int scvIdx)
    {
        return priVars[temperatureIdx];
    }

    template<class ParameterCache>
    static Scalar enthalpy_(const FluidState& fluidState,
                            const ParameterCache& paramCache,
                            const int phaseIdx)
    {
        return FluidSystem::enthalpy(fluidState, paramCache, phaseIdx);
    }
};

} // end namespace

#endif
