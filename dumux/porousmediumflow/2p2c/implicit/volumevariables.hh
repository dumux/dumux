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
 * \brief Contains the quantities which are constant within a
 *        finite volume in the two-phase two-component model.
 */
#ifndef DUMUX_2P2C_VOLUME_VARIABLES_HH
#define DUMUX_2P2C_VOLUME_VARIABLES_HH

#include <dumux/implicit/model.hh>
#include <dumux/material/fluidstates/compositional.hh>
#include <dumux/material/constraintsolvers/computefromreferencephase.hh>
#include <dumux/material/constraintsolvers/misciblemultiphasecomposition.hh>
#include "properties.hh"
#include "indices.hh"

namespace Dumux
{

/*!
 * \ingroup TwoPTwoCModel
 * \ingroup ImplicitVolumeVariables
 * \brief Contains the quantities which are constant within a
 *        finite volume in the two-phase two-component model.
 */
template <class TypeTag>
class TwoPTwoCVolumeVariables : public ImplicitVolumeVariables<TypeTag>
{
    typedef ImplicitVolumeVariables<TypeTag> ParentType;

    typedef typename GET_PROP_TYPE(TypeTag, VolumeVariables) Implementation;
    typedef typename GET_PROP_TYPE(TypeTag, Problem) Problem;
    typedef typename GET_PROP_TYPE(TypeTag, FVElementGeometry) FVElementGeometry;
    typedef typename GET_PROP_TYPE(TypeTag, PrimaryVariables) PrimaryVariables;
    typedef typename GET_PROP_TYPE(TypeTag, MaterialLaw) MaterialLaw;
    enum {
        numPhases = GET_PROP_VALUE(TypeTag, NumPhases),
        numComponents = GET_PROP_VALUE(TypeTag, NumComponents)
    };

    typedef typename GET_PROP_TYPE(TypeTag, Indices) Indices;
    enum {
        wCompIdx = Indices::wCompIdx,
        nCompIdx = Indices::nCompIdx,
        wPhaseIdx = Indices::wPhaseIdx,
        nPhaseIdx = Indices::nPhaseIdx
    };

    // present phases
    enum {
        wPhaseOnly = Indices::wPhaseOnly,
        nPhaseOnly = Indices::nPhaseOnly,
        bothPhases = Indices::bothPhases
    };

    // formulations
    enum {
        formulation = GET_PROP_VALUE(TypeTag, Formulation),
        pwsn = TwoPTwoCFormulation::pwsn,
        pnsw = TwoPTwoCFormulation::pnsw
    };

    // primary variable indices
    enum {
        switchIdx = Indices::switchIdx,
        pressureIdx = Indices::pressureIdx
    };

    typedef typename GET_PROP_TYPE(TypeTag, GridView) GridView;
    typedef typename GridView::template Codim<0>::Entity Element;
    enum { dim = GridView::dimension};

    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
    typedef typename GET_PROP_TYPE(TypeTag, FluidSystem) FluidSystem;
    typedef Dumux::MiscibleMultiPhaseComposition<Scalar, FluidSystem> MiscibleMultiPhaseComposition;
    static const bool useMoles = GET_PROP_VALUE(TypeTag, UseMoles);
    static const bool useConstraintSolver = GET_PROP_VALUE(TypeTag, UseConstraintSolver);
    static_assert(useMoles || (!useMoles && useConstraintSolver),
                  "if UseMoles is set false, UseConstraintSolver has to be set to true");
    typedef Dumux::ComputeFromReferencePhase<Scalar, FluidSystem> ComputeFromReferencePhase;

    enum { isBox = GET_PROP_VALUE(TypeTag, ImplicitIsBox) };
    enum { dofCodim = isBox ? dim : 0 };

public:

    typedef typename GET_PROP_TYPE(TypeTag, FluidState) FluidState;

    /*!
     * \copydoc ImplicitVolumeVariables::update
     */
    void update(const PrimaryVariables &priVars,
                const Problem &problem,
                const Element &element,
                const FVElementGeometry &fvGeometry,
                const int scvIdx,
                const bool isOldSol)
    {
        ParentType::update(priVars,
                           problem,
                           element,
                           fvGeometry,
                           scvIdx,
                           isOldSol);

        int dofIdxGlobal = problem.model().dofMapper().subIndex(element, scvIdx, dofCodim);
        phasePresence_ = problem.model().phasePresence(dofIdxGlobal, isOldSol);

        completeFluidState(priVars, problem, element, fvGeometry, scvIdx, fluidState_, isOldSol);

        /////////////
        // calculate the remaining quantities
        /////////////
        const auto& materialParams =
            problem.spatialParams().materialLawParams(element, fvGeometry, scvIdx);

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

            // binary diffusion coefficients
            diffCoeff_[phaseIdx] =
            FluidSystem::binaryDiffusionCoefficient(fluidState_,
                                                    paramCache,
                                                    phaseIdx,
                                                    wCompIdx,
                                                    nCompIdx);
            Valgrind::CheckDefined(diffCoeff_[phaseIdx]);
        }

        // porosity
        porosity_ = problem.spatialParams().porosity(element,
                                                     fvGeometry,
                                                     scvIdx);
        Valgrind::CheckDefined(porosity_);

        // energy related quantities not contained in the fluid state
        asImp_().updateEnergy_(priVars, problem, element, fvGeometry, scvIdx, isOldSol);
    }

    /*!
     * \copydoc ImplicitModel::completeFluidState
     * \param isOldSol Specifies whether this is the previous solution or the current one
     */
    static void completeFluidState(const PrimaryVariables& priVars,
                                   const Problem& problem,
                                   const Element& element,
                                   const FVElementGeometry& fvGeometry,
                                   int scvIdx,
                                   FluidState& fluidState,
                                   bool isOldSol = false)
    {
        Scalar t = Implementation::temperature_(priVars, problem, element,
                                                fvGeometry, scvIdx);
        fluidState.setTemperature(t);
        int dofIdxGlobal = problem.model().dofMapper().subIndex(element, scvIdx, dofCodim);
        int phasePresence = problem.model().phasePresence(dofIdxGlobal, isOldSol);

        /////////////
        // set the saturations
        /////////////
        Scalar sn;
        if (phasePresence == nPhaseOnly)
            sn = 1.0;
        else if (phasePresence == wPhaseOnly) {
            sn = 0.0;
        }
        else if (phasePresence == bothPhases) {
            if (formulation == pwsn)
                sn = priVars[switchIdx];
            else if (formulation == pnsw)
                sn = 1.0 - priVars[switchIdx];
            else DUNE_THROW(Dune::InvalidStateException, "Formulation: " << formulation << " is invalid.");
        }
        else DUNE_THROW(Dune::InvalidStateException, "phasePresence: " << phasePresence << " is invalid.");
        fluidState.setSaturation(wPhaseIdx, 1 - sn);
        fluidState.setSaturation(nPhaseIdx, sn);

        /////////////
        // set the pressures of the fluid phases
        /////////////

        // calculate capillary pressure
        const auto& materialParams =
            problem.spatialParams().materialLawParams(element, fvGeometry, scvIdx);
        Scalar pc = MaterialLaw::pc(materialParams, 1 - sn);

        if (formulation == pwsn) {
            fluidState.setPressure(wPhaseIdx, priVars[pressureIdx]);
            fluidState.setPressure(nPhaseIdx, priVars[pressureIdx] + pc);
        }
        else if (formulation == pnsw) {
            fluidState.setPressure(nPhaseIdx, priVars[pressureIdx]);
            fluidState.setPressure(wPhaseIdx, priVars[pressureIdx] - pc);
        }
        else DUNE_THROW(Dune::InvalidStateException, "Formulation: " << formulation << " is invalid.");

        /////////////
        // calculate the phase compositions
        /////////////
        typename FluidSystem::ParameterCache paramCache;

        //get the phase pressures and set the fugacity coefficients here if constraintsolver is not used
        Scalar pn = 0;
        Scalar pw = 0;

        if(!useConstraintSolver) {
            if (formulation == pwsn) {
                pw = priVars[pressureIdx];
                pn = pw + pc;
            }
            else {
                pn = priVars[pressureIdx];
                pw = pn - pc;
            }

            for (int phaseIdx = 0; phaseIdx < numPhases; ++ phaseIdx) {
                assert(FluidSystem::isIdealMixture(phaseIdx));

                for (int compIdx = 0; compIdx < numComponents; ++ compIdx) {
                    Scalar phi = FluidSystem::fugacityCoefficient(fluidState, paramCache, phaseIdx, compIdx);
                    fluidState.setFugacityCoefficient(phaseIdx, compIdx, phi);
                }
            }
        }

        // now comes the tricky part: calculate phase compositions
        if (phasePresence == bothPhases) {
            // both phases are present, phase compositions are a
            // result of the the nonwetting <-> wetting equilibrium. This is
            // the job of the "MiscibleMultiPhaseComposition"
            // constraint solver
            if(useConstraintSolver) {
                MiscibleMultiPhaseComposition::solve(fluidState,
                                                     paramCache,
                                                     /*setViscosity=*/true,
                                                     /*setInternalEnergy=*/false);
            }
            // ... or calculated explicitly this way ...
            else {
                //get the partial pressure of the main component of the the wetting phase ("H20") within the nonwetting (gas) phase == vapor pressure due to equilibrium
                //note that in this case the fugacityCoefficient * pw is the vapor pressure (see implementation in respective fluidsystem)
                Scalar partPressH2O = FluidSystem::fugacityCoefficient(fluidState,
                                                                       wPhaseIdx,
                                                                       wCompIdx) * pw;

                // get the partial pressure of the main component of the the nonwetting (gas) phase ("Air")
                Scalar partPressAir = pn - partPressH2O;

                //calculate the mole fractions of the components within the nonwetting phase
                Scalar xnn = partPressAir/pn;
                Scalar xnw = partPressH2O/pn;

                // calculate the mole fractions of the components within the wetting phase
                //note that in this case the fugacityCoefficient * pw is the Henry Coefficient (see implementation in respective fluidsystem)
                Scalar xwn = partPressAir
                  / (FluidSystem::fugacityCoefficient(fluidState,
                                                      wPhaseIdx,nCompIdx)
                  * pw);

                Scalar xww = 1.0 -xwn;

                //set all mole fractions
                fluidState.setMoleFraction(wPhaseIdx, wCompIdx, xww);
                fluidState.setMoleFraction(wPhaseIdx, nCompIdx, xwn);
                fluidState.setMoleFraction(nPhaseIdx, wCompIdx, xnw);
                fluidState.setMoleFraction(nPhaseIdx, nCompIdx, xnn);

                paramCache.updateComposition(fluidState, wPhaseIdx);
                paramCache.updateComposition(fluidState, nPhaseIdx);

                //set the phase densities
                Scalar rhoW = FluidSystem::density(fluidState, paramCache, wPhaseIdx);
                Scalar rhoN = FluidSystem::density(fluidState, paramCache, nPhaseIdx);

                fluidState.setDensity(wPhaseIdx, rhoW);
                fluidState.setDensity(nPhaseIdx, rhoN);
            }
        }
        else if (phasePresence == nPhaseOnly) {
            // only the nonwetting phase is present, i.e. nonwetting phase
            // composition is stored explicitly.
            if(useMoles)
            {
                fluidState.setMoleFraction(nPhaseIdx, nCompIdx, 1 - priVars[switchIdx]);
                fluidState.setMoleFraction(nPhaseIdx, wCompIdx, priVars[switchIdx]);
            }
            else
            {
                // setMassFraction() has only to be called 1-numComponents times
                fluidState.setMassFraction(nPhaseIdx, wCompIdx, priVars[switchIdx]);
            }

            // calculate the composition of the remaining phases (as
            // well as the densities of all phases). This is the job
            // of the "ComputeFromReferencePhase" constraint solver
            if (useConstraintSolver) {
                ComputeFromReferencePhase::solve(fluidState,
                                                 paramCache,
                                                 nPhaseIdx,
                                                 /*setViscosity=*/true,
                                                 /*setInternalEnergy=*/false);
            }
            // ... or calculated explicitly this way ...
            else {
                // note that the water phase is actually not existing!
                // thus, this is used as phase switch criterion
                Scalar xnw = priVars[switchIdx];
                Scalar xnn = 1.0 -xnw;

                //first, xww:
                // xnw * pn = "actual" (hypothetical) vapor pressure
                // fugacityCoefficient * pw = vapor pressure given by thermodynamic conditions
                // Here, xww is not actually the mole fraction of water in the wetting phase
                // xww is only the ratio of "actual" vapor pressure / "thermodynamic" vapor pressure
                // If xww > 1 : gas is over-saturated with water vapor,
                // condensation takes place (see switch criterion in model)
                Scalar xww = xnw * pn
                  / (FluidSystem::fugacityCoefficient(fluidState,
                                                      wPhaseIdx,wCompIdx)
                     * pw);

                // now, xwn:
                //partialPressure / xwn = Henry
                //partialPressure = xnn * pn
                //xwn = xnn * pn / Henry
                // Henry = fugacityCoefficient * pw
                Scalar xwn = xnn * pn / (FluidSystem::fugacityCoefficient(fluidState,
                                                                          wPhaseIdx,nCompIdx)
                                         * pw);

                fluidState.setMoleFraction(wPhaseIdx, wCompIdx, xww);
                fluidState.setMoleFraction(wPhaseIdx, nCompIdx, xwn);

                paramCache.updateComposition(fluidState, wPhaseIdx);
                paramCache.updateComposition(fluidState, nPhaseIdx);

                Scalar rhoW = FluidSystem::density(fluidState, paramCache, wPhaseIdx);
                Scalar rhoN = FluidSystem::density(fluidState, paramCache, nPhaseIdx);

                fluidState.setDensity(wPhaseIdx, rhoW);
                fluidState.setDensity(nPhaseIdx, rhoN);
            }
        }
        else if (phasePresence == wPhaseOnly) {
            // only the wetting phase is present, i.e. wetting phase
            // composition is stored explicitly.
            if(useMoles) // mole-fraction formulation
            {
                fluidState.setMoleFraction(wPhaseIdx, wCompIdx, 1-priVars[switchIdx]);
                fluidState.setMoleFraction(wPhaseIdx, nCompIdx, priVars[switchIdx]);
            }
            else // mass-fraction formulation
            {
                // setMassFraction() has only to be called 1-numComponents times
                fluidState.setMassFraction(wPhaseIdx, nCompIdx, priVars[switchIdx]);
            }

            // calculate the composition of the remaining phases (as
            // well as the densities of all phases). This is the job
            // of the "ComputeFromReferencePhase" constraint solver
            if (useConstraintSolver) {
                ComputeFromReferencePhase::solve(fluidState,
                                                 paramCache,
                                                 wPhaseIdx,
                                                 /*setViscosity=*/true,
                                                 /*setInternalEnergy=*/false);
            }
            // ... or calculated explicitly this way ...
            else {
                // note that the gas phase is actually not existing!
                // thus, this is used as phase switch criterion
                Scalar xwn = priVars[switchIdx];

                //first, xnw:
                //psteam = xnw * pn = partial pressure of water in gas phase
                //psteam = fugacityCoefficient * pw
                Scalar xnw = (FluidSystem::fugacityCoefficient(fluidState,
                                                               wPhaseIdx,wCompIdx)
                              * pw) / pn ;

                //now, xnn:
                // xwn = partialPressure / Henry
                // partialPressure = pn * xnn
                // xwn = pn * xnn / Henry
                // xnn = xwn * Henry / pn
                // Henry = fugacityCoefficient * pw
                Scalar xnn = xwn * (FluidSystem::fugacityCoefficient(fluidState,
                                                                     wPhaseIdx,nCompIdx)
                                    * pw) / pn ;

                fluidState.setMoleFraction(nPhaseIdx, nCompIdx, xnn);
                fluidState.setMoleFraction(nPhaseIdx, wCompIdx, xnw);

                paramCache.updateComposition(fluidState, wPhaseIdx);
                paramCache.updateComposition(fluidState, nPhaseIdx);

                Scalar rhoW = FluidSystem::density(fluidState, paramCache, wPhaseIdx);
                Scalar rhoN = FluidSystem::density(fluidState, paramCache, nPhaseIdx);

                fluidState.setDensity(wPhaseIdx, rhoW);
                fluidState.setDensity(nPhaseIdx, rhoN);
            }
        }

        for (int phaseIdx = 0; phaseIdx < numPhases; ++phaseIdx) {
            //set the viscosity here if constraintsolver is not used
            if(!useConstraintSolver) {
                const Scalar mu =
                FluidSystem::viscosity(fluidState,
                                       paramCache,
                                       phaseIdx);
                fluidState.setViscosity(phaseIdx,mu);
            }
            // compute and set the enthalpy
            Scalar h = Implementation::enthalpy_(fluidState, paramCache, phaseIdx);
            fluidState.setEnthalpy(phaseIdx, h);
        }
   }


    /*!
     * \brief Returns the phase state within the control volume.
     */
    const FluidState &fluidState() const
    { return fluidState_; }

    /*!
     * \brief Returns the phase presence within the control volume.
     */
    const int phasePresence() const
    { return phasePresence_; }

    /*!
     * \brief Returns the saturation of a given phase within
     *        the control volume in \f$[-]\f$.
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
     *        control volume in \f$[kg/m^3]\f$.
     *
     * \param phaseIdx The phase index
     */
    Scalar density(const int phaseIdx) const
    { return fluidState_.density(phaseIdx); }

    /*!
     * \brief Returns the dynamic viscosity of the fluid within the
     *        control volume in \f$\mathrm{[Pa s]}\f$.
     *
     * \param phaseIdx The phase index
     */
    Scalar viscosity(const int phaseIdx) const
    { return fluidState_.viscosity(phaseIdx); }

    /*!
     * \brief Returns the mass density of a given phase within the
     *        control volume in \f$[mol/m^3]\f$.
     *
     * \param phaseIdx The phase index
     */
    Scalar molarDensity(const int phaseIdx) const
    { return fluidState_.density(phaseIdx) / fluidState_.averageMolarMass(phaseIdx); }

    /*!
     * \brief Returns the effective pressure of a given phase within
     *        the control volume in \f$[kg/(m*s^2)=N/m^2=Pa]\f$.
     *
     * \param phaseIdx The phase index
     */
    Scalar pressure(const int phaseIdx) const
    { return fluidState_.pressure(phaseIdx); }

    /*!
     * \brief Returns temperature within the control volume in \f$[K]\f$.
     *
     * Note that we assume thermodynamic equilibrium, i.e. the
     * temperature of the rock matrix and of all fluid phases are
     * identical.
     */
    Scalar temperature() const
    { return fluidState_.temperature(/*phaseIdx=*/0); }

    /*!
     * \brief Returns the relative permeability of a given phase within
     *        the control volume in \f$[-]\f$.
     *
     * \param phaseIdx The phase index
     */
    Scalar relativePermeability(const int phaseIdx) const
    {
        return relativePermeability_[phaseIdx];
    }

    /*!
     * \brief Returns the effective mobility of a given phase within
     *        the control volume in \f$[s*m/kg]\f$.
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
     * \brief Returns the average porosity within the control volume in \f$[-]\f$.
     */
    Scalar porosity() const
    { return porosity_; }

    /*!
     * \brief Returns the binary diffusion coefficients for a phase in \f$[m^2/s]\f$.
     */
    Scalar diffCoeff(const int phaseIdx) const
    { return diffCoeff_[phaseIdx]; }


protected:
    static Scalar temperature_(const PrimaryVariables &priVars,
                               const Problem& problem,
                               const Element &element,
                               const FVElementGeometry &fvGeometry,
                               int scvIdx)
    {
        return problem.temperatureAtPos(fvGeometry.subContVol[scvIdx].global);
    }

    template<class ParameterCache>
    static Scalar enthalpy_(const FluidState& fluidState,
                            const ParameterCache& paramCache,
                            const int phaseIdx)
    {
        return 0;
    }

    /*!
     * \brief Called by update() to compute the energy related quantities
     */
    void updateEnergy_(const PrimaryVariables &sol,
                       const Problem &problem,
                       const Element &element,
                       const FVElementGeometry &fvGeometry,
                       const int scvIdx,
                       bool isOldSol)
    { }

    Scalar porosity_; //!< Effective porosity within the control volume
    Scalar relativePermeability_[numPhases]; //!< Relative permeability within the control volume
    Scalar diffCoeff_[numPhases]; //!< Binary diffusion coefficients for the phases
    FluidState fluidState_;
    int phasePresence_;

private:
    Implementation &asImp_()
    { return *static_cast<Implementation*>(this); }

    const Implementation &asImp_() const
    { return *static_cast<const Implementation*>(this); }


};

} // end namespace

#endif
