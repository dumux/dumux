// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup NavierStokesModel
 * \brief Variables for the single-phase Navier-Stokes mass model.
 */
#ifndef DUMUX_NAVIERSTOKES_MASS_1P_VARIABLES_HH
#define DUMUX_NAVIERSTOKES_MASS_1P_VARIABLES_HH

#include <dumux/common/concepts/ipdata_.hh>

namespace Dumux {

/*!
 * \ingroup NavierStokesModel
 * \brief Variables for the single-phase Navier-Stokes mass model.
 *
 * The quantities are defined per local dof and updated at an interpolation point, so that
 * discretizations without sub-control volumes are served as well.
 */
template <class Traits>
class NavierStokesMassOnePVariables
{
    using Scalar = typename Traits::PrimaryVariables::value_type;

    static_assert(!Traits::ModelTraits::enableEnergyBalance(),
                  "The energy balance is currently not implemented for variables defined per local dof");

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
     * \brief Update all quantities for a local dof
     *
     * \param elemSol A vector containing all primary variables connected to the element
     * \param problem The object specifying the problem which ought to
     *                be simulated
     * \param elemDisc The element discretization
     * \param ipData The interpolation point data
     */
    template<class ElementSolution, class Problem, class ElementDiscretization, Concept::LocalDofIpData IpData>
    void update(const ElementSolution& elemSol,
                const Problem& problem,
                const ElementDiscretization& elemDisc,
                const IpData& ipData)
    {
        priVars_ = elemSol[ipData.localDofIndex()];
        extrusionFactor_ = problem.spatialParams().extrusionFactor(elemDisc, ipData, elemSol);

        completeFluidState(elemSol, problem, elemDisc, ipData, fluidState_);
    }

    /*!
     * \brief Sets the complete fluid state
     *
     * \param elemSol A vector containing all primary variables connected to the element
     * \param problem The object specifying the problem which ought to
     *                be simulated
     * \param elemDisc The element discretization
     * \param ipData The interpolation point data
     * \param fluidState A container with the current (physical) state of the fluid
     */
    template<class ElementSolution, class Problem, class ElementDiscretization, Concept::LocalDofIpData IpData>
    void completeFluidState(const ElementSolution& elemSol,
                            const Problem& problem,
                            const ElementDiscretization& elemDisc,
                            const IpData& ipData,
                            FluidState& fluidState) const
    {
        fluidState.setTemperature(/*phaseIdx=*/0, problem.spatialParams().temperature(elemDisc, ipData, elemSol));
        fluidState.setPressure(/*phaseIdx=*/0, elemSol[ipData.localDofIndex()][Indices::pressureIdx]);

        // saturation in a single phase is always 1 and thus redundant to set. But since we use
        // the fluid state shared by the immiscible multi-phase models, we have to set it here
        fluidState.setSaturation(/*phaseIdx=*/0, 1.0);

        typename FluidSystem::ParameterCache paramCache;
        paramCache.updatePhase(fluidState, /*phaseIdx=*/0);

        Scalar value = FluidSystem::density(fluidState, paramCache, /*phaseIdx=*/0);
        fluidState.setDensity(/*phaseIdx=*/0, value);

        value = FluidSystem::viscosity(fluidState, paramCache, /*phaseIdx=*/0);
        fluidState.setViscosity(/*phaseIdx=*/0, value);

        // the enthalpy is zero for isothermal models
        fluidState.setEnthalpy(/*phaseIdx=*/0, 0.0);
    }

    /*!
     * \brief Return how much the localDof is extruded.
     */
    Scalar extrusionFactor() const
    { return extrusionFactor_; }

    /*!
     * \brief Returns the effective pressure \f$\mathrm{[Pa]}\f$ of a given phase within
     *        the control volume.
     */
    Scalar pressure(int phaseIdx = 0) const
    { return fluidState_.pressure(phaseIdx); }

    /*!
     * \brief Returns the saturation.
     */
    Scalar saturation(int phaseIdx = 0) const
    { return 1.0; }

    /*!
     * \brief Returns the mass density \f$\mathrm{[kg/m^3]}\f$ of a given phase.
     */
    Scalar density(int phaseIdx = 0) const
    { return fluidState_.density(phaseIdx); }

    /*!
     * \brief Returns the molar density \f$\mathrm{[mol/m^3]}\f$ of a given phase.
     */
    Scalar molarDensity(int phaseIdx = 0) const
    { return fluidState_.molarDensity(phaseIdx); }

    /*!
     * \brief Returns the dynamic viscosity \f$\mathrm{[Pa s]}\f$ of the fluid.
     */
    Scalar viscosity(int phaseIdx = 0) const
    { return fluidState_.viscosity(phaseIdx); }

    /*!
     * \brief Returns the temperature \f$\mathrm{[K]}\f$ inside the control volume.
     */
    Scalar temperature() const
    { return fluidState_.temperature(); }

    /*!
     * \brief Returns the fluid state of the control volume.
     */
    const FluidState& fluidState() const
    { return fluidState_; }

    /*!
     * \brief Return a component of the primary variable vector
     * \param pvIdx The index of the primary variable of interest
     */
    Scalar priVar(const int pvIdx) const
    { return priVars_[pvIdx]; }

    /*!
     * \brief Return the primary variable vector
     */
    const PrimaryVariables& priVars() const
    { return priVars_; }

private:
    PrimaryVariables priVars_;
    FluidState fluidState_;
    Scalar extrusionFactor_;
};

} // end namespace Dumux

#endif
