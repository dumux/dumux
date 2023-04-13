// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup NavierStokesModel
 *
 * \copydoc Dumux::NavierStokesVolumeVariables
 */
#ifndef DUMUX_FREEFLOW_NAVIERSTOKES_MASS_1P_VOLUME_VARIABLES_HH
#define DUMUX_FREEFLOW_NAVIERSTOKES_MASS_1P_VOLUME_VARIABLES_HH

#include <dumux/freeflow/navierstokes/scalarvolumevariables.hh>
#include <dumux/freeflow/navierstokes/energy/volumevariables.hh>

namespace Dumux {

/*!
 * \ingroup NavierStokesModel
 * \brief Volume variables for the single-phase Navier-Stokes model.
 */
template <class Traits>
class NavierStokesMassOnePVolumeVariables
: public NavierStokesScalarConservationModelVolumeVariables<Traits>
, public NavierStokesEnergyVolumeVariables<Traits, NavierStokesMassOnePVolumeVariables<Traits>>
{
    using ParentType = NavierStokesScalarConservationModelVolumeVariables<Traits>;
    using EnergyVolumeVariables = NavierStokesEnergyVolumeVariables<Traits, NavierStokesMassOnePVolumeVariables<Traits>>;
    using Scalar = typename Traits::PrimaryVariables::value_type;

public:
    //! export the type used for the primary variables
    using PrimaryVariables = typename Traits::PrimaryVariables;
    //! export the indices type
    using Indices = typename Traits::ModelTraits::Indices;
    //! Export the underlying fluid system
    using FluidSystem = typename Traits::FluidSystem;
    //! Export the fluid state type
    using FluidState = typename Traits::FluidState;

    //! Return number of phases considered by the model
    static constexpr int numFluidPhases() { return Traits::ModelTraits::numFluidPhases(); }
    //! Return number of components considered by the model
    static constexpr int numFluidComponents() { return Traits::ModelTraits::numFluidComponents(); }

     /*!
     * \brief Update all quantities for a given control volume
     *
     * \param elemSol A vector containing all primary variables connected to the element
     * \param problem The object specifying the problem which ought to
     *                be simulated
     * \param element An element which contains part of the control volume
     * \param scv The sub-control volume
     */
    template<class ElementSolution, class Problem, class Element, class SubControlVolume>
    void update(const ElementSolution& elemSol,
                const Problem& problem,
                const Element& element,
                const SubControlVolume& scv)
    {
        ParentType::update(elemSol, problem, element, scv);
        completeFluidState(elemSol, problem, element,scv, fluidState_);
        EnergyVolumeVariables::updateEffectiveThermalConductivity();
    }

    /*!
     * \brief Sets complete fluid state
     *
     * \param elemSol A vector containing all primary variables connected to the element
     * \param problem The object specifying the problem which ought to
     *                be simulated
     * \param element An element which contains part of the control volume
     * \param scv The sub-control volume
     * \param fluidState A container with the current (physical) state of the fluid
     */
    template<class ElemSol, class Problem, class Element, class Scv>
    void completeFluidState(const ElemSol& elemSol,
                            const Problem& problem,
                            const Element& element,
                            const Scv& scv,
                            FluidState& fluidState)
    {
        fluidState.setTemperature(/*phaseIdx=*/0, EnergyVolumeVariables::getTemperature(elemSol, problem, element, scv));

        const auto& priVars = elemSol[scv.localDofIndex()];
        fluidState.setPressure(/*phaseIdx=*/0, priVars[Indices::pressureIdx]);

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
        value = EnergyVolumeVariables::enthalpy(fluidState, paramCache);
        fluidState.setEnthalpy(/*phaseIdx=*/0, value);
    }

     /*!
     * \brief Returns the effective pressure \f$\mathrm{[Pa]}\f$ of a given phase within
     *        the control volume.
     */
    Scalar pressure(int phaseIdx = 0) const
    { return fluidState_.pressure(phaseIdx); }

    /*!
     * \brief Returns the fluid state of the control volume.
     */
    const FluidState& fluidState() const
    { return fluidState_; }

    /*!
     * \brief Returns the dynamic viscosity \f$\mathrm{[Pa s]}\f$ of the fluid within the
     *        control volume.
     */
    Scalar viscosity(int phaseIdx = 0) const
    { return fluidState_.viscosity(phaseIdx); }

     /*!
     * \brief Returns the mass density \f$\mathrm{[kg/m^3]}\f$ of a given phase within the
     *        control volume.
     */
    Scalar density(int phaseIdx = 0) const
    { return fluidState_.density(phaseIdx); }

    /*!
     * \brief Returns the temperature \f$\mathrm{[K]}\f$ inside the sub-control volume.
     *
     * Note that we assume thermodynamic equilibrium, i.e. the
     * temperatures of the rock matrix and of all fluid phases are
     * identical.
     */
    Scalar temperature() const
    { return fluidState_.temperature(); }

protected:
    FluidState fluidState_;
};

} // end namespace Dumux

#endif
