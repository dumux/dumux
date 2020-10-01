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
 *
 * \copydoc Dumux::NavierStokesVolumeVariables
 */
#ifndef DUMUX_NAVIERSTOKES_VOLUME_VARIABLES_HH
#define DUMUX_NAVIERSTOKES_VOLUME_VARIABLES_HH

#include <dumux/freeflow/volumevariables.hh>

namespace Dumux {

/*!
 * \ingroup NavierStokesModel
 * \brief Volume variables for the single-phase Navier-Stokes model.
 */
template <class Traits>
class NavierStokesVolumeVariables : public FreeFlowVolumeVariables< Traits, NavierStokesVolumeVariables<Traits> >
{
    using ThisType = NavierStokesVolumeVariables<Traits>;
    using ParentType = FreeFlowVolumeVariables<Traits, ThisType>;

    using Scalar = typename Traits::PrimaryVariables::value_type;

public:
    //! export the underlying fluid system
    using FluidSystem = typename Traits::FluidSystem;
    //! export the fluid state type
    using FluidState = typename Traits::FluidState;
    //! export the indices type
    using Indices = typename Traits::ModelTraits::Indices;

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
        completeFluidState(elemSol, problem, element, scv, fluidState_);
    }

    /*!
     * \brief Update the fluid state
     */
    template<class ElementSolution, class Problem, class Element, class SubControlVolume>
    static void completeFluidState(const ElementSolution& elemSol,
                                   const Problem& problem,
                                   const Element& element,
                                   const SubControlVolume& scv,
                                   FluidState& fluidState)
    {
        const Scalar t = ParentType::temperature(elemSol, problem, element, scv);
        fluidState.setTemperature(t);

        fluidState.setPressure(0, elemSol[0][Indices::pressureIdx]);

        // saturation in a single phase is always 1 and thus redundant
        // to set. But since we use the fluid state shared by the
        // immiscible multi-phase models, so we have to set it here...
        fluidState.setSaturation(0, 1.0);

        typename FluidSystem::ParameterCache paramCache;
        paramCache.updateAll(fluidState);

        Scalar value = FluidSystem::density(fluidState, paramCache, 0);
        fluidState.setDensity(0, value);

        value = FluidSystem::viscosity(fluidState, paramCache, 0);
        fluidState.setViscosity(0, value);

        // compute and set the enthalpy
        value = ParentType::enthalpy(fluidState, paramCache);
        fluidState.setEnthalpy(0, value);
    }

    /*!
     * \brief Return the effective pressure \f$\mathrm{[Pa]}\f$ of a given phase within
     *        the control volume.
     */
    Scalar pressure(int phaseIdx = 0) const
    { return fluidState_.pressure(0); }

    /*!
     * \brief Return the mass density \f$\mathrm{[kg/m^3]}\f$ of a given phase within the
     *        control volume.
     */
    Scalar density(int phaseIdx = 0) const
    { return fluidState_.density(0); }

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
     * \brief Returns the molar mass of a given phase within the
     *        control volume.
     */
    Scalar molarMass(int phaseIdx = 0) const
    {
        return fluidState_.averageMolarMass(0);
    }

    /*!
     * \brief Return the dynamic viscosity \f$\mathrm{[Pa s]}\f$ of the fluid within the
     *        control volume.
     */
    Scalar viscosity(int phaseIdx = 0) const
    { return fluidState_.viscosity(0); }

    /*!
     * \brief Return the effective dynamic viscosity \f$\mathrm{[Pa s]}\f$ of the fluid within the
     *        control volume.
     */
    Scalar effectiveViscosity() const
    { return viscosity(); }

    /*!
     * \brief Return the fluid state of the control volume.
     */
    const FluidState& fluidState() const
    { return fluidState_; }

protected:
    FluidState fluidState_;
};

} // end namespace Dumux

#endif
