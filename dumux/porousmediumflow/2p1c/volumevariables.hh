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
 * \ingroup TwoPOneCModel
 * \copydoc Dumux::TwoPOneCVolumeVariables
 */
#ifndef DUMUX_2P1C_VOLUME_VARIABLES_HH
#define DUMUX_2P1C_VOLUME_VARIABLES_HH

#include <dumux/porousmediumflow/volumevariables.hh>
#include <dumux/material/fluidstates/compositional.hh>
#include "indices.hh"

namespace Dumux
{

/*!
 * \ingroup TwoPOneCModel
 * \brief The volume variables (i.e. secondary variables) for the two-phase one-component model.
 */
template <class TypeTag>
class TwoPOneCVolumeVariables : public PorousMediumFlowVolumeVariables<TypeTag>
{
    using ParentType = PorousMediumFlowVolumeVariables<TypeTag>;
    using Implementation = typename GET_PROP_TYPE(TypeTag, VolumeVariables);
    using Scalar = typename GET_PROP_TYPE(TypeTag, Scalar);
    using GridView = typename GET_PROP_TYPE(TypeTag, GridView);
    using Problem = typename GET_PROP_TYPE(TypeTag, Problem);
    using SpatialParams = typename GET_PROP_TYPE(TypeTag, SpatialParams);
    using PermeabilityType = typename SpatialParams::PermeabilityType;
    using SubControlVolume = typename GET_PROP_TYPE(TypeTag, SubControlVolume);
    using ElementSolutionVector = typename GET_PROP_TYPE(TypeTag, ElementSolutionVector);
    using FluidSystem = typename GET_PROP_TYPE(TypeTag, FluidSystem);
    using MaterialLaw = typename GET_PROP_TYPE(TypeTag, MaterialLaw);
    using Indices = typename GET_PROP_TYPE(TypeTag, Indices);

    enum {
        numPhases = GET_PROP_VALUE(TypeTag, NumPhases),

        wPhaseIdx = Indices::wPhaseIdx,
        nPhaseIdx = Indices::nPhaseIdx,

        switch1Idx = Indices::switch1Idx,
        pressureIdx = Indices::pressureIdx
    };

    // present phases
    enum {
        twoPhases = Indices::twoPhases,
        wPhaseOnly  = Indices::wPhaseOnly,
        nPhaseOnly  = Indices::nPhaseOnly,
    };

    using Element = typename GridView::template Codim<0>::Entity;

public:
    //! The type of the object returned by the fluidState() method
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
        ParentType::update(elemSol, problem, element, scv);

        completeFluidState(elemSol, problem, element, scv, fluidState_);

        /////////////
        // calculate the remaining quantities
        /////////////
        const auto& materialParams = problem.spatialParams().materialLawParams(element, scv, elemSol);

        // Second instance of a parameter cache.
        // Could be avoided if diffusion coefficients also
        // became part of the fluid state.
        typename FluidSystem::ParameterCache paramCache;
        paramCache.updateAll(fluidState_);
        for (int phaseIdx = 0; phaseIdx < numPhases; ++phaseIdx)
        {
            // relative permeabilities
            Scalar kr;
            if (phaseIdx == wPhaseIdx)
                kr = MaterialLaw::krw(materialParams, saturation(wPhaseIdx));
            else // ATTENTION: krn requires the wetting phase saturation
                // as parameter!
                kr = MaterialLaw::krn(materialParams, saturation(wPhaseIdx));
            relativePermeability_[phaseIdx] = kr;
            Valgrind::CheckDefined(relativePermeability_[phaseIdx]);
        }

        // porosity & permeability
        porosity_ = problem.spatialParams().porosity(element, scv, elemSol);
        permeability_ = problem.spatialParams().permeability(element, scv, elemSol);
    }

    //! Update the fluidstate
    static void completeFluidState(const ElementSolutionVector& elemSol,
                                   const Problem& problem,
                                   const Element& element,
                                   const SubControlVolume& scv,
                                   FluidState& fluidState)
    {

        // capillary pressure parameters
        const auto& materialParams = problem.spatialParams().materialLawParams(element, scv, elemSol);

        const auto& priVars = ParentType::extractDofPriVars(elemSol, scv);
        const auto phasePresence = priVars.state();

        // get saturations
        Scalar sw(0.0);
        Scalar sg(0.0);
        if (phasePresence == twoPhases)
        {
            sw = priVars[switch1Idx];
            sg = 1.0 - sw;
        }
        else if (phasePresence == wPhaseOnly)
        {
            sw = 1.0;
            sg = 0.0;
        }
        else if (phasePresence == nPhaseOnly)
        {
            sw = 0.0;
            sg = 1.0;
        }
        else DUNE_THROW(Dune::InvalidStateException, "phasePresence: " << phasePresence << " is invalid.");
        Valgrind::CheckDefined(sg);

        fluidState.setSaturation(wPhaseIdx, sw);
        fluidState.setSaturation(nPhaseIdx, sg);

        // get gas phase pressure
        const Scalar pg = priVars[pressureIdx];

        // calculate capillary pressure
        const Scalar pc = MaterialLaw::pc(materialParams, sw);

        // set wetting phase pressure
        const Scalar pw = pg - pc;

        //set pressures
        fluidState.setPressure(wPhaseIdx, pw);
        fluidState.setPressure(nPhaseIdx, pg);

        // get temperature
        Scalar temperature;
        if (phasePresence == wPhaseOnly || phasePresence == nPhaseOnly)
            temperature = priVars[switch1Idx];
        else if (phasePresence == twoPhases)
            temperature = FluidSystem::vaporTemperature(fluidState, wPhaseIdx);
        else
            DUNE_THROW(Dune::InvalidStateException, "phasePresence: " << phasePresence << " is invalid.");

        Valgrind::CheckDefined(temperature);

        fluidState.setTemperature(temperature);

        // set the densities
        const Scalar rhoW = FluidSystem::density(fluidState, wPhaseIdx);
        const Scalar rhoG = FluidSystem::density(fluidState, nPhaseIdx);

        fluidState.setDensity(wPhaseIdx, rhoW);
        fluidState.setDensity(nPhaseIdx, rhoG);

        //get the viscosity and mobility
        for (int phaseIdx = 0; phaseIdx < numPhases; ++phaseIdx)
        {
            // Mobilities
            const Scalar mu =
                FluidSystem::viscosity(fluidState,
                                       phaseIdx);
            fluidState.setViscosity(phaseIdx,mu);
        }

        // the enthalpies (internal energies are directly calculated in the fluidstate
        for (int phaseIdx = 0; phaseIdx < numPhases; ++phaseIdx)
        {
            const Scalar h = FluidSystem::enthalpy(fluidState, phaseIdx);
            fluidState.setEnthalpy(phaseIdx, h);
        }
    }

    /*!
     * \brief Returns the fluid state for the control-volume.
     */
    const FluidState &fluidState() const
    { return fluidState_; }

    /*!
     * \brief Returns the effective saturation of a given phase within
     *        the control volume.
     *
     * \param phaseIdx The phase index
     */
    Scalar saturation(const int phaseIdx) const
    { return fluidState_.saturation(phaseIdx); }

    /*!
     * \brief Returns the mass density of a given phase within the
     *        control volume.
     *
     * \param phaseIdx The phase index
     */
    Scalar density(const int phaseIdx) const
    { return fluidState_.density(phaseIdx); }

    /*!
     * \brief Returns the molar density of a given phase within the
     *        control volume.
     *
     * \param phaseIdx The phase index
     */
    Scalar molarDensity(const int phaseIdx) const
    { return fluidState_.density(phaseIdx) / fluidState_.averageMolarMass(phaseIdx); }

    /*!
     * \brief Returns the effective pressure of a given phase within
     *        the control volume.
     *
     * \param phaseIdx The phase index
     */
    Scalar pressure(const int phaseIdx) const
    { return fluidState_.pressure(phaseIdx); }

    /*!
     * \brief Returns temperature inside the sub-control volume.
     *
     * Note that we assume thermodynamic equilibrium, i.e. the
     * temperatures of the rock matrix and of all fluid phases are
     * identical.
     */
    Scalar temperature(const int phaseIdx = 0) const
    { return fluidState_.temperature(phaseIdx); }

    /*!
     * \brief Returns the effective mobility of a given phase within
     *        the control volume.
     *
     * \param phaseIdx The phase index
     */
    Scalar mobility(const int phaseIdx) const
    {
        return relativePermeability_[phaseIdx]/fluidState_.viscosity(phaseIdx);
    }

    /*!
     * \brief Returns the effective capillary pressure within the control volume
     *        in \f$[kg/(m*s^2)=N/m^2=Pa]\f$.
     */
    Scalar capillaryPressure() const
    { return fluidState_.pressure(nPhaseIdx) - fluidState_.pressure(wPhaseIdx); }

    /*!
     * \brief Returns the average porosity within the control volume.
     */
    Scalar porosity() const
    { return porosity_; }

    /*!
     * \brief Returns the average permeability within the control volume in \f$[m^2]\f$.
     */
    const PermeabilityType& permeability() const
    { return permeability_; }

    /*!
     * \brief Returns the vapor temperature \f$T_{vap}(p_n)\f$ of the fluid within the control volume.
     */
    Scalar vaporTemperature() const
    { return FluidSystem::vaporTemperature(fluidState_, wPhaseIdx);}

protected:
    FluidState fluidState_;

private:
    Scalar porosity_;               //!< Effective porosity within the control volume
    PermeabilityType permeability_; //!> Effective permeability within the control volume
    std::array<Scalar, numPhases> relativePermeability_;

    Implementation &asImp_()
    { return *static_cast<Implementation*>(this); }

    const Implementation &asImp_() const
    { return *static_cast<const Implementation*>(this); }
};

} // end namespace

#endif
