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
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the           *
 *   GNU General Public License for more details.                            *
 *                                                                           *
 *   You should have received a copy of the GNU General Public License       *
 *   along with this program.  If not, see <http://www.gnu.org/licenses/>.   *
 *****************************************************************************/
/*!
 * \file
 *
 * \brief Contains the secondary variables (Quantities which are
 *        constant within a finite volume) of the M-phase, N-component
 *        model.
 */
#ifndef DUMUX_MPNC_VOLUME_VARIABLES_HH
#define DUMUX_MPNC_VOLUME_VARIABLES_HH

#include "diffusion/volumevariables.hh"
#include "energy/mpncvolumevariablesenergy.hh"
#include "mass/mpncvolumevariablesmass.hh"
#include "mpncvolumevariablesia.hh"

#include <dumux/boxmodels/common/boxvolumevariables.hh>
#include <dumux/material/constraintsolvers/ncpflash.hh>

namespace Dumux
{
/*!
 * \ingroup MPNCModel
 * \ingroup BoxVolumeVariables
 * \brief Contains the quantities which are are constant within a
 *        finite volume in the M-phase, N-component model.
 */
template <class TypeTag>
class MPNCVolumeVariables
    : public BoxVolumeVariables<TypeTag>
    , public MPNCVolumeVariablesIA<TypeTag, GET_PROP_VALUE(TypeTag, EnableKinetic), GET_PROP_VALUE(TypeTag, EnableKineticEnergy)>
    , public MPNCVolumeVariablesMass<TypeTag, GET_PROP_VALUE(TypeTag, EnableKinetic)>
    , public MPNCVolumeVariablesDiffusion<TypeTag, GET_PROP_VALUE(TypeTag, EnableDiffusion) || GET_PROP_VALUE(TypeTag, EnableKinetic)>
    , public MPNCVolumeVariablesEnergy<TypeTag, GET_PROP_VALUE(TypeTag, EnableEnergy), GET_PROP_VALUE(TypeTag, EnableKineticEnergy)>
{
    typedef BoxVolumeVariables<TypeTag> ParentType;
    typedef typename GET_PROP_TYPE(TypeTag, VolumeVariables) Implementation;
    typedef typename GET_PROP_TYPE(TypeTag, Problem) Problem;
    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
    typedef typename GET_PROP_TYPE(TypeTag, GridView) GridView;
    typedef typename GET_PROP_TYPE(TypeTag, FVElementGeometry) FVElementGeometry;
    typedef typename GET_PROP_TYPE(TypeTag, FluidSystem) FluidSystem;
    typedef typename GET_PROP_TYPE(TypeTag, FluidState) FluidState;
    typedef typename FluidSystem::ParameterCache ParameterCache;
    typedef typename GET_PROP_TYPE(TypeTag, MaterialLaw) MaterialLaw;
    typedef typename GET_PROP_TYPE(TypeTag, MaterialLawParams) MaterialLawParams;
    typedef typename GET_PROP_TYPE(TypeTag, PrimaryVariables) PrimaryVariables;
    typedef typename GET_PROP_TYPE(TypeTag, Indices) Indices;

    enum {numPhases = GET_PROP_VALUE(TypeTag, NumPhases)};
    enum {numComponents = GET_PROP_VALUE(TypeTag, NumComponents)};
    enum {enableEnergy = GET_PROP_VALUE(TypeTag, EnableEnergy)};
    enum {enableKinetic = GET_PROP_VALUE(TypeTag, EnableKinetic)};
    enum {enableKineticEnergy = GET_PROP_VALUE(TypeTag, EnableKineticEnergy)};
    enum {enableDiffusion = GET_PROP_VALUE(TypeTag, EnableDiffusion) || enableKinetic};
    enum {S0Idx = Indices::S0Idx};
    enum {p0Idx = Indices::p0Idx};

    typedef typename GridView::template Codim<0>::Entity Element;
    typedef MPNCVolumeVariablesMass<TypeTag, enableKinetic> MassVolumeVariables;
    typedef MPNCVolumeVariablesEnergy<TypeTag, enableEnergy, enableKineticEnergy> EnergyVolumeVariables;
    typedef MPNCVolumeVariablesIA<TypeTag, enableKinetic, enableKineticEnergy> IAVolumeVariables;
    typedef MPNCVolumeVariablesDiffusion<TypeTag, enableDiffusion> DiffusionVolumeVariables;

public:
    MPNCVolumeVariables()
    { hint_ = NULL; };

    /*!
     * \brief Set the volume variables which should be used as initial
     *        conditions for complex calculations.
     */
    void setHint(const Implementation *hint)
    {
        hint_ = hint;
    }

    /*!
     * \copydoc BoxVolumeVariables::update
     */
    void update(const PrimaryVariables &priVars,
                const Problem &problem,
                const Element &element,
                const FVElementGeometry &fvGeometry,
                const unsigned int scvIdx,
                const bool isOldSol)
    {
        Valgrind::CheckDefined(priVars);
        ParentType::update(priVars,
                           problem,
                           element,
                           fvGeometry,
                           scvIdx,
                           isOldSol);
        ParentType::checkDefined();

        ParameterCache paramCache;

        /////////////
        // set the phase saturations
        /////////////
        Scalar sumSat = 0;
        for (int phaseIdx = 0; phaseIdx < numPhases - 1; ++phaseIdx) {
            sumSat += priVars[S0Idx + phaseIdx];
            fluidState_.setSaturation(phaseIdx, priVars[S0Idx + phaseIdx]);
        }
        Valgrind::CheckDefined(sumSat);
        fluidState_.setSaturation(numPhases - 1, 1.0 - sumSat);

        /////////////
        // set the fluid phase temperatures
        /////////////
        EnergyVolumeVariables::updateTemperatures(fluidState_,
                                                  paramCache,
                                                  priVars,
                                                  element,
                                                  fvGeometry,
                                                  scvIdx,
                                                  problem);

        /////////////
        // set the phase pressures
        /////////////

        // capillary pressure parameters
        const MaterialLawParams &materialParams =
            problem.spatialParams().materialLawParams(element, fvGeometry, scvIdx);
        // capillary pressures
        Scalar capPress[numPhases];
        MaterialLaw::capillaryPressures(capPress, materialParams, fluidState_);
        // add to the pressure of the first fluid phase
        Scalar p0 = priVars[p0Idx];
        for (int phaseIdx = 0; phaseIdx < numPhases; ++ phaseIdx)
            fluidState_.setPressure(phaseIdx, p0 - capPress[0] + capPress[phaseIdx]);

        /////////////
        // set the fluid compositions
        /////////////
        MassVolumeVariables::update(fluidState_,
                                    paramCache,
                                    priVars,
                                    hint_,
                                    problem,
                                    element,
                                    fvGeometry,
                                    scvIdx);
        MassVolumeVariables::checkDefined();

        /////////////
        // Porosity
        /////////////

        // porosity
        porosity_ = problem.spatialParams().porosity(element,
                                                     fvGeometry,
                                                     scvIdx);
        Valgrind::CheckDefined(porosity_);

        /////////////
        // Phase mobilities
        /////////////

        // relative permeabilities
        MaterialLaw::relativePermeabilities(relativePermeability_,
                                            materialParams,
                                            fluidState_);

        // dynamic viscosities
        for (int phaseIdx = 0; phaseIdx < numPhases; ++phaseIdx) {
            // viscosities
            Scalar mu = FluidSystem::viscosity(fluidState_, paramCache, phaseIdx);
            fluidState_.setViscosity(phaseIdx, mu);
        }

        /////////////
        // diffusion
        /////////////

        // update the diffusion part of the volume data
        DiffusionVolumeVariables::update(fluidState_, paramCache, *this, problem);
        DiffusionVolumeVariables::checkDefined();

        /////////////
        // energy
        /////////////

        // update the remaining parts of the energy module
        EnergyVolumeVariables::update(fluidState_,
                                      paramCache,
                                      element,
                                      fvGeometry,
                                      scvIdx,
                                      problem);
        EnergyVolumeVariables::checkDefined();

        // make sure the quantities in the fluid state are well-defined
        fluidState_.checkDefined();

        // specific interfacial area,
        // well also all the dimensionless numbers :-)
        // well, also the mass transfer rate
        IAVolumeVariables::update(*this,
                                  fluidState_,
                                  paramCache,
                                  priVars,
                                  problem,
                                  element,
                                  fvGeometry,
                                  scvIdx);
        IAVolumeVariables::checkDefined();
        checkDefined();
    }

    /*!
     * \brief Return the fluid configuration at the given primary
     *        variables
     */
    const FluidState &fluidState() const
    { return fluidState_; }

    /*!
     * \brief Returns the effective mobility of a given phase within
     *        the control volume.
     */
    Scalar mobility(const unsigned int phaseIdx) const
    { return relativePermeability(phaseIdx)/fluidState_.viscosity(phaseIdx); }

    /*!
     * \brief Returns the viscosity of a given phase within
     *        the control volume.
     */
    Scalar viscosity(const unsigned int phaseIdx) const
    { return fluidState_.viscosity(phaseIdx); }

    /*!
     * \brief Returns the relative permeability of a given phase within
     *        the control volume.
     */
    Scalar relativePermeability(const unsigned int phaseIdx) const
    { return relativePermeability_[phaseIdx]; }

    /*!
     * \brief Returns the average porosity within the control volume.
     */
    Scalar porosity() const
    { return porosity_; }

    /*!
     * \brief Returns true iff the fluid state is in the active set
     *        for a phase,
     */
    bool isPhaseActive(const unsigned int phaseIdx) const
    {
        return
            phasePresentIneq(fluidState(), phaseIdx) -
            phaseNotPresentIneq(fluidState(), phaseIdx)
            >= 0;
    }

    /*!
     * \brief Returns the value of the NCP-function for a phase.
     */
    Scalar phaseNcp(const unsigned int phaseIdx) const
    {
        Scalar aEval = phaseNotPresentIneq(this->evalPoint().fluidState(), phaseIdx);
        Scalar bEval = phasePresentIneq(this->evalPoint().fluidState(), phaseIdx);
        if (aEval > bEval)
            return phasePresentIneq(fluidState(), phaseIdx);
        return phaseNotPresentIneq(fluidState(), phaseIdx);
    };

    /*!
     * \brief Returns the value of the inequality where a phase is
     *        present.
     */
    Scalar phasePresentIneq(const FluidState &fluidState,
                            const unsigned int phaseIdx) const
    { return fluidState.saturation(phaseIdx); }

    /*!
     * \brief Returns the value of the inequality where a phase is not
     *        present.
     */
    Scalar phaseNotPresentIneq(const FluidState &fluidState,
                               const unsigned int phaseIdx) const
    {
        // difference of sum of mole fractions in the phase from 100%
        Scalar a = 1;
        for (int compIdx = 0; compIdx < numComponents; ++compIdx)
            a -= fluidState.moleFraction(phaseIdx, compIdx);
        return a;
    }

    /*!
     * \brief If running in valgrind this makes sure that all
     *        quantities in the volume variables are defined.
     */
    void checkDefined() const
    {
#if !defined NDEBUG && HAVE_VALGRIND
        ParentType::checkDefined();

        Valgrind::CheckDefined(porosity_);
        Valgrind::CheckDefined(hint_);
        Valgrind::CheckDefined(relativePermeability_);

        fluidState_.checkDefined();
#endif
    }

protected:
    Scalar porosity_; //!< Effective porosity within the control volume
    Scalar relativePermeability_[numPhases]; //!< Effective relative permeability within the control volume

    const Implementation *hint_;

    //! Mass fractions of each component within each phase
    FluidState fluidState_;
};

} // end namepace

#endif
