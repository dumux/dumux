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
 * \ingroup BlackOilModel
 * \brief Contains the quantities which are constant within a
 *        finite volume in the three-phase three-component model.
 */
#ifndef DUMUX_BLACKOIL_VOLUME_VARIABLES_HH
#define DUMUX_BLACKOIL_VOLUME_VARIABLES_HH

//#define WE_SEARCH_FOR_BUGS_BLACKOILVOLVARS

#include <dumux/material/constants.hh>
#include <dumux/material/fluidstates/compositional.hh>
#include <dumux/material/constraintsolvers/computefromreferencephase.hh>
#include <dumux/material/constraintsolvers/misciblemultiphasecomposition.hh>

#include <dumux/porousmediumflow/volumevariables.hh>
#include <dumux/porousmediumflow/nonisothermal/volumevariables.hh>
#include <dumux/material/solidstates/updatesolidvolumefractions.hh>
#include <dumux/common/optionalscalar.hh>


namespace Dumux {

/*!
 * \ingroup BlackOilModel
 * \brief Contains the quantities which are are constant within a
 *        finite volume in the three-phase three-component model.
 */
template <class Traits>
class BlackOilVolumeVariables
: public PorousMediumFlowVolumeVariables<Traits>
{
    using ParentType = PorousMediumFlowVolumeVariables<Traits>;

    using Scalar = typename Traits::PrimaryVariables::value_type;
    using PermeabilityType = typename Traits::PermeabilityType;

    using FS = typename Traits::FluidSystem;
    using MiscibleMultiPhaseComposition = Dumux::MiscibleMultiPhaseComposition<Scalar, FS>;
    using ComputeFromReferencePhase = Dumux::ComputeFromReferencePhase<Scalar, FS>;

    using ModelTraits = typename Traits::ModelTraits;
    using Idx = typename ModelTraits::Indices;
    static constexpr int numComponents = ParentType::numFluidComponents();
    static constexpr int numPhases = 3; //  ModelTraits::numFluidPhases();
    enum {
        wCompIdx = FS::wCompIdx,
        gCompIdx = FS::gCompIdx,
        oCompIdx = FS::oCompIdx,

        wPhaseIdx = FS::wPhaseIdx,
        gPhaseIdx = FS::gPhaseIdx,
        oPhaseIdx = FS::oPhaseIdx,

        // eine Wassersaettigung
        // eine Privar Gassaettigung
        switch1Idx = Idx::switch1Idx,
        switch2Idx = Idx::switch2Idx,
        pressureIdx = Idx::pressureIdx,

    };
    // Alternative names for the indices, TODO: Make code consistent!
    static constexpr int oSaturationIdx = switch1Idx;
    static constexpr int wSaturationIdx = switch2Idx;
public:
    //! export fluid state type
    using FluidState = typename Traits::FluidState;
    //! export fluid system type
    using FluidSystem = typename Traits::FluidSystem;
    //! export the indices
    using Indices = typename ModelTraits::Indices;
    //! export type of solid state
    using SolidState = typename Traits::SolidState;
    //! export type of solid system
    using SolidSystem = typename Traits::SolidSystem;


    template<class ElemSol, class Problem, class Element, class Scv>
    void update(const ElemSol &elemSol,
                const Problem &problem,
                const Element &element,
                const Scv& scv)
    {
        ParentType::update(elemSol, problem, element, scv);
        const auto& priVars = elemSol[scv.localDofIndex()];

        // TODO: Whatever to be done with this!! This is only for looking if something is happening!!!!!
        //fluidState_.setTemperature(problem_.temperature(element));
        fluidState_.setTemperature(311.0); // Note: This is just for now!!!!
        // Fluid matric interaction parameters
        const auto fluidMatrixInteraction = problem.spatialParams().fluidMatrixInteraction(element, scv, elemSol);

        // update the saturations
        Scalar sw_ = priVars[wSaturationIdx];
        Scalar so_ = priVars[oSaturationIdx];
        Scalar sg_ = 1-sw_-so_;
        fluidState_.setSaturation(wPhaseIdx, sw_);
        fluidState_.setSaturation(oPhaseIdx, so_);
        fluidState_.setSaturation(gPhaseIdx, sg_);

        /* now the pressures */
        pg_ = priVars[pressureIdx];
        Scalar pcgw = fluidMatrixInteraction.pcgw(sw_, so_);
        Scalar pcnw = fluidMatrixInteraction.pcnw(sw_, so_);
        Scalar pcgn = fluidMatrixInteraction.pcgn(sw_, so_);
        const Scalar pcAlpha = fluidMatrixInteraction.pcAlpha(sw_, so_);
        const Scalar pcNW1 = 0.0; // TODO: this should be possible to assign in the problem file
        po_ = pg_- pcAlpha * pcgn - (1.-pcAlpha)*(pcgw - pcNW1);
        pw_ = po_ - pcAlpha * pcnw - (1.-pcAlpha)*pcNW1;

        fluidState_.setPressure(wPhaseIdx, pw_);
        fluidState_.setPressure(gPhaseIdx, pg_);
        fluidState_.setPressure(oPhaseIdx, po_);

        typename FluidSystem::ParameterCache paramCache;

        // update phase compositions. first, set everything to 0, then
        // make the gas/water phases consist of only the gas/water
        // components and calculate the composition of the liquid oil
        // phase from the gas formation factor plus the gas/oil
        // formation volume factors and the reference densities
        for (int phaseIdx = 0; phaseIdx < numPhases; ++ phaseIdx)
            for (int compIdx = 0; compIdx < numComponents; ++compIdx)
                fluidState_.setMoleFraction(phaseIdx, compIdx, 0.0);
        // set composition of gas and water phases
        fluidState_.setMoleFraction(gPhaseIdx, gCompIdx, 1.0);
        fluidState_.setMoleFraction(wPhaseIdx, wCompIdx, 1.0);

        // retrieve the relevant black-oil parameters from the fluid
        // system.
        Scalar p = fluidState_.pressure(oPhaseIdx);
        Scalar Bg = FluidSystem::gasFormationVolumeFactor(p);
        Scalar Bo = FluidSystem::oilFormationVolumeFactor(p);
        Scalar Rs = FluidSystem::gasFormationFactor(p);
        Scalar rhoo = FluidSystem::surfaceDensity(oPhaseIdx)/Bo;
        Scalar rhorefg = FluidSystem::surfaceDensity(gPhaseIdx);
        Scalar rhog = rhorefg/Bg;
        Scalar MG = FluidSystem::molarMass(gPhaseIdx);
        Scalar MO = FluidSystem::molarMass(oPhaseIdx);

        // calculate composition of oil phase in terms of mass
        // fractions.
        Scalar XoG = Rs*rhorefg / rhoo;
        Scalar XoO = 1 - XoG;

        if (XoG < 0 || XoO < 0) {
            DUNE_THROW(NumericalProblem,
                       "Only positive values are allowed for the mass fractions "
                       "of the oil and the gas components in the oil phase");
        }

        // handle undersaturated oil. We interpret negative gas
        // saturations as the amount of gas that needs to get out pf
        // the oil. The saturation of the oil phase is then given by 1
        // minus the water saturation
        sg_ = fluidState_.saturation(gPhaseIdx);
        sw_ = fluidState_.saturation(wPhaseIdx);
        if (sg_ < 0) {
            so_ = 1 - sw_;

            if (so_ > 0) {
                // Calculate the total partial mass density of the gas and
                // oil components in the saturated oil phase
                Scalar rho_oGSat = so_*XoG*rhoo;

                // Calculate the amount of gas that cannot be subtracted
                // from the oil and thus needs to be accounted for in the
                // gas saturation
                Scalar rho_GInPhase = std::max(0.0, std::abs(sg_)*rhog - rho_oGSat);

                // calculate the composition of the undersaturated oil phase
                Scalar rho_oG = std::max(0.0, rho_oGSat - std::abs(sg_)*rhog);

                // convert to mass fractions
                XoG = rho_oG/(so_*rhoo);
                XoO = 1 - XoG;

                // calculate the bubble pressure of the oil with the new
                // composition
                Scalar pBubb = FluidSystem::oilSaturationPressure(XoG);

                // calculate the density of the oil at the saturation pressure
                rhoo = FluidSystem::surfaceDensity(oPhaseIdx)/FluidSystem::oilFormationVolumeFactor(p);
                // compress to the actual pressure of the system
                rhoo *= 1.0 + FluidSystem::oilCompressibility() * (fluidState_.pressure(oPhaseIdx) - pBubb);

                // convert the "residual gas" into a saturation
                sg_ = - rho_GInPhase/rhog;
            }

            // update the saturations. Gas phase is not present!
            fluidState_.setSaturation(gPhaseIdx, sg_);
            fluidState_.setSaturation(oPhaseIdx, so_);
        }

        // convert mass to mole fractions
        Scalar avgMolarMass = MO/(1 + XoG*(MO/MG - 1));
        Scalar xoG = XoG*avgMolarMass/MG;
        Scalar xoO = 1 - xoG;

        // set the oil-phase composition
        fluidState_.setMoleFraction(oPhaseIdx, gCompIdx, xoG);
        fluidState_.setMoleFraction(oPhaseIdx, oCompIdx, xoO);

        rhoo = FluidSystem::density(fluidState_, oPhaseIdx);

        // set the phase densities for the phases
        fluidState_.setDensity(wPhaseIdx, FluidSystem::density(fluidState_, wPhaseIdx));
        fluidState_.setDensity(gPhaseIdx, rhog);
        fluidState_.setDensity(oPhaseIdx, rhoo);

        for (int phaseIdx = 0; phaseIdx < numPhases; ++phaseIdx)
        {
            // compute and set the viscosity
            Scalar mu = FluidSystem::viscosity(fluidState_, phaseIdx);
            fluidState_.setViscosity(phaseIdx, mu);

             // mobilities
            const Scalar kr = fluidMatrixInteraction.kr(phaseIdx,
                                  fluidState_.saturation(wPhaseIdx),
                                  fluidState_.saturation(oPhaseIdx));

             mobility_[phaseIdx] = kr / mu;
        }

        // porosity & permeabilty
        updateSolidVolumeFractions(elemSol, problem, element, scv, solidState_, numComponents);
        permeability_ = problem.spatialParams().permeability(element, scv, elemSol);
    }
    /*!
     * \brief Returns the phase state for the control volume.
     */
    const FluidState &fluidState() const
    { return fluidState_; }

    /*!
     * \brief Returns the phase state for the control volume.
     */
    const SolidState &solidState() const
    { return solidState_; }

    /*!
     * \brief Returns the average molar mass \f$\mathrm{[kg/mol]}\f$ of the fluid phase.
     *
     * \param phaseIdx The phase index
     */
    Scalar averageMolarMass(int phaseIdx) const
    { return fluidState_.averageMolarMass(phaseIdx); }

    /*!
     * \brief Returns the effective saturation of a given phase within
     *        the control volume.
     *
     * \param phaseIdx The phase index
     */
    Scalar saturation(const int phaseIdx) const
    { return fluidState_.saturation(phaseIdx); }

    /*!
     * \brief Returns the mass fraction of a given component in a
     *        given phase within the control volume in \f$[-]\f$.
     *
     * \param phaseIdx The phase index
     * \param compIdx The component index
     */
    Scalar massFraction(const int phaseIdx, const int compIdx) const
    { return fluidState_.massFraction(phaseIdx, compIdx); }

    /*!
     * \brief Returns the mole fraction of a given component in a
     *        given phase within the control volume in \f$[-]\f$.
     *
     * \param phaseIdx The phase index
     * \param compIdx The component index
     */
    Scalar moleFraction(const int phaseIdx, const int compIdx) const
    { return fluidState_.moleFraction(phaseIdx, compIdx); }

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
    Scalar temperature() const
    { return fluidState_.temperature(/*phaseIdx=*/0); }

    /*!
     * \brief Returns the effective mobility of a given phase within
     *        the control volume.
     *
     * \param phaseIdx The phase index
     */
    Scalar mobility(const int phaseIdx) const
    { return mobility_[phaseIdx]; }

    /*!
     * \brief Returns the effective capillary pressure within the control volume.
     */
    Scalar capillaryPressure() const
    { return fluidState_.capillaryPressure(); }

    /*!
     * \brief Returns the average porosity within the control volume.
     */
    Scalar porosity() const
    { return solidState_.porosity(); }

    /*!
     * \brief Returns the adsorption information.
     */
    Scalar bulkDensTimesAdsorpCoeff() const
    {
        if (bulkDensTimesAdsorpCoeff_)
            return bulkDensTimesAdsorpCoeff_.value();
        else
            DUNE_THROW(Dune::NotImplemented, "Your spatialParams do not provide an adsorption model");
    }

    /*!
     * \brief Returns the average permeability within the control volume in \f$[m^2]\f$.
     */
    const PermeabilityType& permeability() const
    { return permeability_; }

    // Black-Oil model does not include diffusion coefficients!



protected:
    FluidState fluidState_;
    SolidState solidState_;


private:
    Scalar sw_, sg_, so_, pg_, pw_, po_;

    Scalar moleFrac_[ModelTraits::numFluidPhases()][ModelTraits::numFluidComponents()];
    Scalar massFrac_[ModelTraits::numFluidPhases()][ModelTraits::numFluidComponents()];

    PermeabilityType permeability_; //!< Effective permeability within the control volume
    Scalar mobility_[ModelTraits::numFluidPhases()];  //!< Effective mobility within the control volume
    OptionalScalar<Scalar> bulkDensTimesAdsorpCoeff_; //!< the basis for calculating adsorbed NAPL

};

} // end namespace Dumux

#endif
