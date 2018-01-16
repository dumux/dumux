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
 * \ingroup NavierStokesModel
 *
 * \copydoc Dumux::NavierStokesVolumeVariables
 */
#ifndef DUMUX_NAVIERSTOKES_VOLUME_VARIABLES_HH
#define DUMUX_NAVIERSTOKES_VOLUME_VARIABLES_HH

#include <dumux/common/properties.hh>
#include <dumux/material/fluidstates/immiscible.hh>

namespace Dumux
{

// forward declaration
template <class TypeTag, bool enableEnergyBalance>
class NavierStokesVolumeVariablesImplementation;

/*!
 * \ingroup NavierStokesModel
 * \brief Volume variables for the single-phase Navier-Stokes model.
 *        The class is specialized for isothermal and non-isothermal models.
 */
template <class TypeTag>
using NavierStokesVolumeVariables = NavierStokesVolumeVariablesImplementation<TypeTag, GET_PROP_VALUE(TypeTag, EnableEnergyBalance)>;

/*!
 * \ingroup NavierStokesModel
 * \brief Volume variables for the isothermal single-phase Navier-Stokes model.
 */
template <class TypeTag>
class NavierStokesVolumeVariablesImplementation<TypeTag, false>
{
    using Scalar = typename GET_PROP_TYPE(TypeTag, Scalar);
    using Problem = typename GET_PROP_TYPE(TypeTag, Problem);
    using FVElementGeometry = typename GET_PROP_TYPE(TypeTag, FVElementGeometry);
    using SubControlVolume = typename FVElementGeometry::SubControlVolume;
    using ElementSolutionVector = typename GET_PROP_TYPE(TypeTag, ElementSolutionVector);
    using Indices = typename GET_PROP_TYPE(TypeTag, Indices);
    using FluidSystem = typename GET_PROP_TYPE(TypeTag, FluidSystem);
    using GridView = typename GET_PROP_TYPE(TypeTag, GridView);
    using Element = typename GridView::template Codim<0>::Entity;

    static const int defaultPhaseIdx = GET_PROP_VALUE(TypeTag, PhaseIdx);

public:

    using FluidState = typename GET_PROP_TYPE(TypeTag, FluidState);

    /*!
     * \brief Update all quantities for a given control volume
     *
     * \param elemSol A vector containing all primary variables connected to the element
     * \param problem The object specifying the problem which ought to
     *                be simulated
     * \param element An element which contains part of the control volume
     * \param scv The sub-control volume
     */
    void update(const ElementSolutionVector &elemSol,
                const Problem &problem,
                const Element &element,
                const SubControlVolume& scv)
    {
        extrusionFactor_ = problem.extrusionFactor(element, scv, elemSol);

        completeFluidState(elemSol, problem, element, scv, fluidState_);
    };

    /*!
     * \brief Returns the primary variables at the dof associated with a given scv.
     */
    static const auto& extractDofPriVars(const ElementSolutionVector& elemSol,
                                         const SubControlVolume& scv)
    { return elemSol[0]; }

    /*!
     * \brief Update the fluid state
     */
    static void completeFluidState(const ElementSolutionVector& elemSol,
                                   const Problem& problem,
                                   const Element& element,
                                   const SubControlVolume& scv,
                                   FluidState& fluidState)
    {
        const Scalar t = problem.temperatureAtPos(scv.dofPosition());
        fluidState.setTemperature(t);

        fluidState.setPressure(defaultPhaseIdx, elemSol[0][Indices::pressureIdx]);

        // saturation in a single phase is always 1 and thus redundant
        // to set. But since we use the fluid state shared by the
        // immiscible multi-phase models, so we have to set it here...
        fluidState.setSaturation(defaultPhaseIdx, 1.0);

        typename FluidSystem::ParameterCache paramCache;
        paramCache.updatePhase(fluidState, defaultPhaseIdx);

        Scalar value = FluidSystem::density(fluidState, paramCache, defaultPhaseIdx);
        fluidState.setDensity(defaultPhaseIdx, value);

        value = FluidSystem::viscosity(fluidState, paramCache, defaultPhaseIdx);
        fluidState.setViscosity(defaultPhaseIdx, value);
    }

    /*!
     * \brief Return how much the sub-control volume is extruded.
     *
     * This means the factor by which a lower-dimensional (1D or 2D)
     * entity needs to be expanded to get a full dimensional cell. The
     * default is 1.0 which means that 1D problems are actually
     * thought as pipes with a cross section of 1 m^2 and 2D problems
     * are assumed to extend 1 m to the back.
     */
    Scalar extrusionFactor() const
    { return extrusionFactor_; }

    /*!
     * \brief Return temperature \f$\mathrm{[K]}\f$ inside the sub-control volume.
     *
     * Note that we assume thermodynamic equilibrium, i.e. the
     * temperatures of the rock matrix and of all fluid phases are
     * identical.
     */
    Scalar temperature() const
    { return fluidState_.temperature(); }

    /*!
     * \brief Return the effective pressure \f$\mathrm{[Pa]}\f$ of a given phase within
     *        the control volume.
     */
    Scalar pressure(int phaseIdx = 0) const
    { return fluidState_.pressure(defaultPhaseIdx); }

    /*!
     * \brief Return the saturation
     */
    Scalar saturation(int phaseIdx = 0) const
    { return 1.0; }

    /*!
     * \brief Return the mass density \f$\mathrm{[kg/m^3]}\f$ of a given phase within the
     *        control volume.
     */
    Scalar density(int phaseIdx = 0) const
    { return fluidState_.density(defaultPhaseIdx); }

    /*!
    * \brief Returns the mass density of a given phase within the
    *        control volume.
    */
   Scalar molarDensity(int phaseIdx = 0) const
   {
       return fluidState_.molarDensity(defaultPhaseIdx);
   }

   /*!
   * \brief Returns the molar mass of a given phase within the
   *        control volume.
   */
  Scalar molarMass(int phaseIdx = 0) const
  {
      return fluidState_.averageMolarMass(defaultPhaseIdx);
  }

    /*!
     * \brief Return the dynamic viscosity \f$\mathrm{[Pa s]}\f$ of the fluid within the
     *        control volume.
     */
    Scalar viscosity(int phaseIdx = 0) const
    { return fluidState_.viscosity(defaultPhaseIdx); }

    /*!
     * \brief Return the fluid state of the control volume.
     */
    const FluidState &fluidState() const
    { return fluidState_; }

    //! The temperature is obtained from the problem as a constant for isothermal models
    static Scalar temperature(const ElementSolutionVector &elemSol,
                              const Problem& problem,
                              const Element &element,
                              const SubControlVolume &scv)
    {
        return problem.temperatureAtPos(scv.dofPosition());
    }

    //! The phase enthalpy is zero for isothermal models
    //! This is needed for completing the fluid state
    template<class FluidState, class ParameterCache>
    static Scalar enthalpy(const FluidState& fluidState,
                           const ParameterCache& paramCache,
                           const int defaultPhaseIdx)
    {
        return 0;
    }

protected:
    FluidState fluidState_;
    Scalar extrusionFactor_;
};

/*!
 * \ingroup NavierStokesModel
 * \brief Volume variables for the non-isothermal single-phase Navier-Stokes model.
 */
template <class TypeTag>
class NavierStokesVolumeVariablesImplementation<TypeTag, true>
: public NavierStokesVolumeVariablesImplementation<TypeTag, false>
{
    using ParentType = NavierStokesVolumeVariablesImplementation<TypeTag, false>;
    using Scalar = typename GET_PROP_TYPE(TypeTag, Scalar);
    using Problem = typename GET_PROP_TYPE(TypeTag, Problem);
    using FVElementGeometry = typename GET_PROP_TYPE(TypeTag, FVElementGeometry);
    using SubControlVolume = typename FVElementGeometry::SubControlVolume;
    using ElementSolutionVector = typename GET_PROP_TYPE(TypeTag, ElementSolutionVector);
    using GridView = typename GET_PROP_TYPE(TypeTag, GridView);
    using Element = typename GridView::template Codim<0>::Entity;
    using FluidSystem = typename GET_PROP_TYPE(TypeTag, FluidSystem);
    using Indices = typename GET_PROP_TYPE(TypeTag, Indices);

    static const int defaultPhaseIdx = GET_PROP_VALUE(TypeTag, PhaseIdx);
    static const int temperatureIdx = Indices::temperatureIdx;


public:

    using FluidState = typename GET_PROP_TYPE(TypeTag, FluidState);

    /*!
     * \brief Update all quantities for a given control volume
     *
     * \param elemSol A vector containing all primary variables connected to the element
     * \param problem The object specifying the problem which ought to
     *                be simulated
     * \param element An element which contains part of the control volume
     * \param scv The sub-control volume
     */
    void update(const ElementSolutionVector &elemSol,
                const Problem &problem,
                const Element &element,
                const SubControlVolume &scv)
    {
        this->extrusionFactor_ = problem.extrusionFactor(element, scv, elemSol);

        completeFluidState(elemSol, problem, element, scv, this->fluidState_);
    }

    /*!
     * \brief Update the fluid state
     */
    static void completeFluidState(const ElementSolutionVector& elemSol,
                                   const Problem& problem,
                                   const Element& element,
                                   const SubControlVolume& scv,
                                   FluidState& fluidState)
    {
        fluidState.setTemperature(elemSol[0][Indices::temperatureIdx]);
        fluidState.setPressure(defaultPhaseIdx, elemSol[0][Indices::pressureIdx]);

        // saturation in a single phase is always 1 and thus redundant
        // to set. But since we use the fluid state shared by the
        // immiscible multi-phase models, so we have to set it here...
        fluidState.setSaturation(defaultPhaseIdx, 1.0);

        typename FluidSystem::ParameterCache paramCache;
        paramCache.updatePhase(fluidState, defaultPhaseIdx);

        Scalar value = FluidSystem::density(fluidState, paramCache, defaultPhaseIdx);
        fluidState.setDensity(defaultPhaseIdx, value);

        value = FluidSystem::viscosity(fluidState, paramCache, defaultPhaseIdx);
        fluidState.setViscosity(defaultPhaseIdx, value);

        // compute and set the enthalpy
        value = FluidSystem::enthalpy(fluidState, paramCache, defaultPhaseIdx);
        fluidState.setEnthalpy(defaultPhaseIdx, value);
    }

    /*!
     * \brief Returns the total internal energy of the fluid phase in the
     *        sub-control volume.
     */
    Scalar internalEnergy() const
    { return this->fluidState_.internalEnergy(defaultPhaseIdx); }

    /*!
     * \brief Returns the total enthalpy of the fluid phase in the sub-control
     *        volume.
     */
    Scalar enthalpy() const
    { return this->fluidState_.enthalpy(defaultPhaseIdx); }

    /*!
     * \brief Return the specific isobaric heat capacity \f$\mathrm{[J/(kg*K)]}\f$
     *        in the sub-control volume.
     */
    Scalar heatCapacity() const
    { return FluidSystem::heatCapacity(this->fluidState_, defaultPhaseIdx); }

    /*!
     * \brief Returns the thermal conductivity \f$\mathrm{[W/(m*K)]}\f$
     *        of the fluid phase in the sub-control volume.
     */
    Scalar thermalConductivity() const
    { return FluidSystem::thermalConductivity(this->fluidState_, defaultPhaseIdx); }

    //! The temperature is a primary variable for non-isothermal models
    using ParentType::temperature;
    static Scalar temperature(const ElementSolutionVector &elemSol,
                              const Problem& problem,
                              const Element &element,
                              const SubControlVolume &scv)
    {
        return ParentType::extractDofPriVars(elemSol, scv)[temperatureIdx];
    }

    //! The phase enthalpy is zero for isothermal models
    //! This is needed for completing the fluid state
    template<class FluidState, class ParameterCache>
    static Scalar enthalpy(const FluidState& fluidState,
                           const ParameterCache& paramCache,
                           const int defaultPhaseIdx)
    {
        return FluidSystem::enthalpy(fluidState, paramCache, defaultPhaseIdx);
    }

};
}

#endif
