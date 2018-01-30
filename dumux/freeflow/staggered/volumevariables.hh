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
 * \brief Quantities required by the one-phase fully implicit model defined on a vertex.
 */
#ifndef DUMUX_NAVIERSTOKES_VOLUME_VARIABLES_HH
#define DUMUX_NAVIERSTOKES_VOLUME_VARIABLES_HH

#include "properties.hh"
#include <dumux/discretization/volumevariables.hh>

#include <dumux/material/fluidstates/immiscible.hh>

namespace Dumux
{

// forward declaration
template <class TypeTag, bool enableEnergyBalance>
class NavierStokesVolumeVariablesImplementation;

/*!
 * \ingroup ImplicitVolumeVariables
 * \brief Base class for the model specific class which provides
 *        access to all volume averaged quantities for the free flow.
 *        The volume variables base class
 *        is specialized for isothermal and non-isothermal models.
 */
template <class TypeTag>
using NavierStokesVolumeVariables = NavierStokesVolumeVariablesImplementation<TypeTag, GET_PROP_VALUE(TypeTag, EnableEnergyBalance)>;

/*!
 * \ingroup NavierStokesModel
 * \ingroup ImplicitVolumeVariables
 * \brief Isothermal base class, contains the quantities which are constant within a
 *        finite volume in the one-phase model.
 */
template <class TypeTag>
class NavierStokesVolumeVariablesImplementation<TypeTag, false>
{
    using Scalar = typename GET_PROP_TYPE(TypeTag, Scalar);
    using Problem = typename GET_PROP_TYPE(TypeTag, Problem);
    using SubControlVolume = typename GET_PROP_TYPE(TypeTag, SubControlVolume);
    using ElementSolutionVector = typename GET_PROP_TYPE(TypeTag, ElementSolutionVector);
    using Indices = typename GET_PROP_TYPE(TypeTag, Indices);
    using FluidSystem = typename GET_PROP_TYPE(TypeTag, FluidSystem);
    using GridView = typename GET_PROP_TYPE(TypeTag, GridView);
    using Element = typename GridView::template Codim<0>::Entity;
    using PrimaryVariables = typename GET_PROP_TYPE(TypeTag, PrimaryVariables);

    static constexpr bool isBox = GET_PROP_VALUE(TypeTag, ImplicitIsBox);
    static const int phaseIdx = Indices::phaseIdx;

public:

    using FluidState = typename GET_PROP_TYPE(TypeTag, FluidState);

    /*!
     * \copydoc ImplicitVolumeVariables::update
     */
    void update(const ElementSolutionVector &elemSol,
                const Problem &problem,
                const Element &element,
                const SubControlVolume& scv)
    {
        priVars_ = extractDofPriVars(elemSol, scv);
        extrusionFactor_ = problem.extrusionFactor(element, scv, elemSol);

        completeFluidState(elemSol, problem, element, scv, fluidState_);
    };

    /*!
     * \brief Return a component of primary variable vector
     *
     * \param pvIdx The index of the primary variable of interest
     */
    Scalar priVar(const int pvIdx) const
    { return priVars_[pvIdx]; } // TODO (method copied from dumux/discretization/volvars)

    /*!
     * \brief Returns the primary variables at the dof associated with a given scv.
     */
    static const PrimaryVariables& extractDofPriVars(const ElementSolutionVector& elemSol,
                                                     const SubControlVolume& scv)
    { return elemSol[0]; }

    /*!
     * \copydoc ImplicitModel::completeFluidState
     */
    static void completeFluidState(const ElementSolutionVector& elemSol,
                                   const Problem& problem,
                                   const Element& element,
                                   const SubControlVolume& scv,
                                   FluidState& fluidState)
    {
        Scalar t = problem.temperatureAtPos(scv.dofPosition());
        fluidState.setTemperature(t);

        fluidState.setPressure(phaseIdx, elemSol[0][Indices::pressureIdx]);

        // saturation in a single phase is always 1 and thus redundant
        // to set. But since we use the fluid state shared by the
        // immiscible multi-phase models, so we have to set it here...
        fluidState.setSaturation(phaseIdx, 1.0);

        typename FluidSystem::ParameterCache paramCache;
        paramCache.updatePhase(fluidState, phaseIdx);

        Scalar value = FluidSystem::density(fluidState, paramCache, phaseIdx);
        fluidState.setDensity(phaseIdx, value);

        value = FluidSystem::viscosity(fluidState, paramCache, phaseIdx);
        fluidState.setViscosity(phaseIdx, value);
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
    { return fluidState_.pressure(phaseIdx); }

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
    { return fluidState_.density(phaseIdx); }

    /*!
    * \brief Returns the mass density of a given phase within the
    *        control volume.
    *
    * \param phaseIdx The phase index
    */
   Scalar molarDensity() const
   {
       return fluidState_.molarDensity(phaseIdx);
   }

   /*!
   * \brief Returns the molar mass of a given phase within the
   *        control volume.
   *
   * \param phaseIdx The phase index
   */
  Scalar molarMass() const
  {
      return fluidState_.averageMolarMass(phaseIdx);
  }

    /*!
     * \brief Return the dynamic viscosity \f$\mathrm{[Pa s]}\f$ of the fluid within the
     *        control volume.
     */
    Scalar viscosity(int phaseIdx = 0) const
    { return fluidState_.viscosity(phaseIdx); }

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
                           const int phaseIdx)
    {
        return 0;
    }

protected:
    FluidState fluidState_;
    PrimaryVariables priVars_;
    Scalar extrusionFactor_;
};

/*!
 * \ingroup NavierStokesModel
 * \ingroup ImplicitVolumeVariables
 * \brief Non-isothermal base class
 */
template <class TypeTag>
class NavierStokesVolumeVariablesImplementation<TypeTag, true>
: public NavierStokesVolumeVariablesImplementation<TypeTag, false>
{
    using ParentType = NavierStokesVolumeVariablesImplementation<TypeTag, false>;
    using Scalar = typename GET_PROP_TYPE(TypeTag, Scalar);
    using Problem = typename GET_PROP_TYPE(TypeTag, Problem);
    using SubControlVolume = typename GET_PROP_TYPE(TypeTag, SubControlVolume);
    using ElementSolutionVector = typename GET_PROP_TYPE(TypeTag, ElementSolutionVector);
    using GridView = typename GET_PROP_TYPE(TypeTag, GridView);
    using Element = typename GridView::template Codim<0>::Entity;
    using FluidSystem = typename GET_PROP_TYPE(TypeTag, FluidSystem);
    using Indices = typename GET_PROP_TYPE(TypeTag, Indices);
    using PrimaryVariables = typename GET_PROP_TYPE(TypeTag, PrimaryVariables);

    static const int phaseIdx = Indices::phaseIdx;
    static const int temperatureIdx = Indices::temperatureIdx;

public:

    using FluidState = typename GET_PROP_TYPE(TypeTag, FluidState);

    /*!
     * \copydoc ImplicitVolumeVariables::update
     */
    void update(const ElementSolutionVector &elemSol,
                const Problem &problem,
                const Element &element,
                const SubControlVolume &scv)
    {
        this->priVars_ = ParentType::extractDofPriVars(elemSol, scv);
        this->extrusionFactor_ = problem.extrusionFactor(element, scv, elemSol);

        completeFluidState(elemSol, problem, element, scv, this->fluidState_);
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
        fluidState.setPressure(phaseIdx, elemSol[0][Indices::pressureIdx]);

        // saturation in a single phase is always 1 and thus redundant
        // to set. But since we use the fluid state shared by the
        // immiscible multi-phase models, so we have to set it here...
        fluidState.setSaturation(phaseIdx, 1.0);

        typename FluidSystem::ParameterCache paramCache;
        paramCache.updatePhase(fluidState, phaseIdx);

        Scalar value = FluidSystem::density(fluidState, paramCache, phaseIdx);
        fluidState.setDensity(phaseIdx, value);

        value = FluidSystem::viscosity(fluidState, paramCache, phaseIdx);
        fluidState.setViscosity(phaseIdx, value);

        // compute and set the enthalpy
        value = FluidSystem::enthalpy(fluidState, paramCache, phaseIdx);
        fluidState.setEnthalpy(phaseIdx, value);
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
                           const int phaseIdx)
    {
        return FluidSystem::enthalpy(fluidState, paramCache, phaseIdx);
    }

};
}

#endif
