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
 *        finite volume in the three-phase three-component model.
 */
#ifndef DUMUX_3P3C_VOLUME_VARIABLES_HH
#define DUMUX_3P3C_VOLUME_VARIABLES_HH

#include <dumux/implicit/model.hh>
#include <dumux/material/constants.hh>
#include <dumux/material/fluidstates/compositionalfluidstate.hh>
#include <dumux/material/constraintsolvers/computefromreferencephase.hh>
#include <dumux/material/constraintsolvers/misciblemultiphasecomposition.hh>
#include "properties.hh"

namespace Dumux
{

/*!
 * \ingroup ThreePThreeCModel
 * \ingroup ImplicitVolumeVariables
 * \brief Contains the quantities which are are constant within a
 *        finite volume in the three-phase three-component model.
 */
template <class TypeTag>
class ThreePThreeCVolumeVariables : public ImplicitVolumeVariables<TypeTag>
{
    typedef ImplicitVolumeVariables<TypeTag> ParentType;
    typedef typename GET_PROP_TYPE(TypeTag, VolumeVariables) Implementation;

    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
    typedef typename GET_PROP_TYPE(TypeTag, GridView) GridView;
    typedef typename GET_PROP_TYPE(TypeTag, Problem) Problem;
    typedef typename GET_PROP_TYPE(TypeTag, FVElementGeometry) FVElementGeometry;
    typedef typename GET_PROP_TYPE(TypeTag, PrimaryVariables) PrimaryVariables;
    typedef typename GET_PROP_TYPE(TypeTag, FluidSystem) FluidSystem;
    typedef typename GET_PROP_TYPE(TypeTag, MaterialLaw) MaterialLaw;
    typedef typename GET_PROP_TYPE(TypeTag, MaterialLawParams) MaterialLawParams;

    // constraint solvers
    typedef Dumux::MiscibleMultiPhaseComposition<Scalar, FluidSystem> MiscibleMultiPhaseComposition;
    typedef Dumux::ComputeFromReferencePhase<Scalar, FluidSystem> ComputeFromReferencePhase;

    typedef typename GET_PROP_TYPE(TypeTag, Indices) Indices;
    enum {
        dim = GridView::dimension,

        numPhases = GET_PROP_VALUE(TypeTag, NumPhases),
        numComponents = GET_PROP_VALUE(TypeTag, NumComponents),

        wCompIdx = Indices::wCompIdx,
        gCompIdx = Indices::gCompIdx,
        nCompIdx = Indices::nCompIdx,

        wPhaseIdx = Indices::wPhaseIdx,
        gPhaseIdx = Indices::gPhaseIdx,
        nPhaseIdx = Indices::nPhaseIdx,

        switch1Idx = Indices::switch1Idx,
        switch2Idx = Indices::switch2Idx,
        pressureIdx = Indices::pressureIdx
    };

    // present phases
    enum {
        threePhases = Indices::threePhases,
        wPhaseOnly  = Indices::wPhaseOnly,
        gnPhaseOnly = Indices::gnPhaseOnly,
        wnPhaseOnly = Indices::wnPhaseOnly,
        gPhaseOnly  = Indices::gPhaseOnly,
        wgPhaseOnly = Indices::wgPhaseOnly
    };

    typedef typename GridView::template Codim<0>::Entity Element;

    static const Scalar R; // universial gas constant

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
                bool isOldSol)
    {
        ParentType::update(priVars,
                           problem,
                           element,
                           fvGeometry,
                           scvIdx,
                           isOldSol);

        bool useConstraintSolver = GET_PROP_VALUE(TypeTag, UseConstraintSolver);

        // capillary pressure parameters
        const MaterialLawParams &materialParams =
            problem.spatialParams().materialLawParams(element, fvGeometry, scvIdx);

        int dofIdxGlobal = problem.model().dofMapper().subIndex(element, scvIdx, dofCodim);
        int phasePresence = problem.model().phasePresence(dofIdxGlobal, isOldSol);

        Scalar temp = Implementation::temperature_(priVars, problem, element, fvGeometry, scvIdx);
        fluidState_.setTemperature(temp);

        /* first the saturations */
        if (phasePresence == threePhases)
        {
            sw_ = priVars[switch1Idx];
            sn_ = priVars[switch2Idx];
            sg_ = 1. - sw_ - sn_;
        }
        else if (phasePresence == wPhaseOnly)
        {
            sw_ = 1.;
            sn_ = 0.;
            sg_ = 0.;
        }
        else if (phasePresence == gnPhaseOnly)
        {
            sw_ = 0.;
            sn_ = priVars[switch2Idx];
            sg_ = 1. - sn_;
        }
        else if (phasePresence == wnPhaseOnly)
        {
            sn_ = priVars[switch2Idx];
            sw_ = 1. - sn_;
            sg_ = 0.;
        }
        else if (phasePresence == gPhaseOnly)
        {
            sw_ = 0.;
            sn_ = 0.;
            sg_ = 1.;
        }
        else if (phasePresence == wgPhaseOnly)
        {
            sw_ = priVars[switch1Idx];
            sn_ = 0.;
            sg_ = 1. - sw_;
        }
        else DUNE_THROW(Dune::InvalidStateException, "phasePresence: " << phasePresence << " is invalid.");
        Valgrind::CheckDefined(sg_);

        fluidState_.setSaturation(wPhaseIdx, sw_);
        fluidState_.setSaturation(gPhaseIdx, sg_);
        fluidState_.setSaturation(nPhaseIdx, sn_);

        /* now the pressures */
        pg_ = priVars[pressureIdx];

        // calculate capillary pressures
        Scalar pcgw = MaterialLaw::pcgw(materialParams, sw_);
        Scalar pcnw = MaterialLaw::pcnw(materialParams, sw_);
        Scalar pcgn = MaterialLaw::pcgn(materialParams, sw_ + sn_);

        Scalar pcAlpha = MaterialLaw::pcAlpha(materialParams, sn_);
        Scalar pcNW1 = 0.0; // TODO: this should be possible to assign in the problem file

        pn_ = pg_- pcAlpha * pcgn - (1.-pcAlpha)*(pcgw - pcNW1);
        pw_ = pn_ - pcAlpha * pcnw - (1.-pcAlpha)*pcNW1;

        fluidState_.setPressure(wPhaseIdx, pw_);
        fluidState_.setPressure(gPhaseIdx, pg_);
        fluidState_.setPressure(nPhaseIdx, pn_);

        // calculate and set all fugacity coefficients. this is
        // possible because we require all phases to be an ideal
        // mixture, i.e. fugacity coefficients are not supposed to
        // depend on composition!
        typename FluidSystem::ParameterCache paramCache;
        // assert(FluidSystem::isIdealGas(gPhaseIdx));
        for (int phaseIdx = 0; phaseIdx < numPhases; ++ phaseIdx) {
            assert(FluidSystem::isIdealMixture(phaseIdx));

            for (int compIdx = 0; compIdx < numComponents; ++ compIdx) {
                Scalar phi = FluidSystem::fugacityCoefficient(fluidState_, paramCache, phaseIdx, compIdx);
                fluidState_.setFugacityCoefficient(phaseIdx, compIdx, phi);
            }
        }

        // now comes the tricky part: calculate phase composition
        if (phasePresence == threePhases) {
            // all phases are present, phase compositions are a
            // result of the the gas <-> liquid equilibrium. This is
            // the job of the "MiscibleMultiPhaseComposition"
            // constraint solver ...
            if (useConstraintSolver) {
                MiscibleMultiPhaseComposition::solve(fluidState_,
                                                     paramCache,
                                                     /*setViscosity=*/true,
                                                     /*setInternalEnergy=*/false);
            }
            // ... or calculated explicitly this way ...
            else {
                Scalar partPressH2O = FluidSystem::fugacityCoefficient(fluidState_,
                                                                      wPhaseIdx,
                                                                      wCompIdx) * pw_;
                Scalar partPressNAPL = FluidSystem::fugacityCoefficient(fluidState_,
                                                                       nPhaseIdx,
                                                                       nCompIdx) * pn_;
                Scalar partPressAir = pg_ - partPressH2O - partPressNAPL;

                Scalar xgn = partPressNAPL/pg_;
                Scalar xgw = partPressH2O/pg_;
                Scalar xgg = partPressAir/pg_;

                // actually, it's nothing else than Henry coefficient
                Scalar xwn = partPressNAPL
                             / (FluidSystem::fugacityCoefficient(fluidState_,
                                                                 wPhaseIdx,nCompIdx)
                                * pw_);
                Scalar xwg = partPressAir
                             / (FluidSystem::fugacityCoefficient(fluidState_,
                                                                 wPhaseIdx,gCompIdx)
                                * pw_);
                Scalar xww = 1.-xwg-xwn;

                Scalar xnn = 1.-2.e-10;
                Scalar xna = 1.e-10;
                Scalar xnw = 1.e-10;

                fluidState_.setMoleFraction(wPhaseIdx, wCompIdx, xww);
                fluidState_.setMoleFraction(wPhaseIdx, gCompIdx, xwg);
                fluidState_.setMoleFraction(wPhaseIdx, nCompIdx, xwn);
                fluidState_.setMoleFraction(gPhaseIdx, wCompIdx, xgw);
                fluidState_.setMoleFraction(gPhaseIdx, gCompIdx, xgg);
                fluidState_.setMoleFraction(gPhaseIdx, nCompIdx, xgn);
                fluidState_.setMoleFraction(nPhaseIdx, wCompIdx, xnw);
                fluidState_.setMoleFraction(nPhaseIdx, gCompIdx, xna);
                fluidState_.setMoleFraction(nPhaseIdx, nCompIdx, xnn);

                Scalar rhoW = FluidSystem::density(fluidState_, wPhaseIdx);
                Scalar rhoG = FluidSystem::density(fluidState_, gPhaseIdx);
                Scalar rhoN = FluidSystem::density(fluidState_, nPhaseIdx);

                fluidState_.setDensity(wPhaseIdx, rhoW);
                fluidState_.setDensity(gPhaseIdx, rhoG);
                fluidState_.setDensity(nPhaseIdx, rhoN);
            }
        }
        else if (phasePresence == wPhaseOnly) {
            // only the water phase is present, water phase composition is
            // stored explicitly.

            // extract mole fractions in the water phase
            Scalar xwg = priVars[switch1Idx];
            Scalar xwn = priVars[switch2Idx];
            Scalar xww = 1 - xwg - xwn;

            // write water mole fractions in the fluid state
            fluidState_.setMoleFraction(wPhaseIdx, wCompIdx, xww);
            fluidState_.setMoleFraction(wPhaseIdx, gCompIdx, xwg);
            fluidState_.setMoleFraction(wPhaseIdx, nCompIdx, xwn);

            // calculate the composition of the remaining phases (as
            // well as the densities of all phases). this is the job
            // of the "ComputeFromReferencePhase" constraint solver ...
            if (useConstraintSolver)
            {
                ComputeFromReferencePhase::solve(fluidState_,
                                                 paramCache,
                                                 wPhaseIdx,
                                                 /*setViscosity=*/true,
                                                 /*setInternalEnergy=*/false);
            }
            // ... or calculated explicitly this way ...
            else {
                // note that the gas phase is actually not existing!
                // thus, this is used as phase switch criterion
                Scalar xgg = xwg * FluidSystem::fugacityCoefficient(fluidState_,
                                                                    wPhaseIdx,gCompIdx)
                                   * pw_ / pg_;
                Scalar xgn = xwn * FluidSystem::fugacityCoefficient(fluidState_,
                                                                    wPhaseIdx,nCompIdx)
                                   * pw_ / pg_;
                Scalar xgw = FluidSystem::fugacityCoefficient(fluidState_,
                                                              wPhaseIdx,wCompIdx)
                                   * pw_ / pg_;


                // note that the gas phase is actually not existing!
                // thus, this is used as phase switch criterion
                Scalar xnn = xwn * FluidSystem::fugacityCoefficient(fluidState_,
                                                                    wPhaseIdx,nCompIdx)
                                   * pw_;
                Scalar xna = 1.e-10;
                Scalar xnw = 1.e-10;

                fluidState_.setMoleFraction(gPhaseIdx, wCompIdx, xgw);
                fluidState_.setMoleFraction(gPhaseIdx, gCompIdx, xgg);
                fluidState_.setMoleFraction(gPhaseIdx, nCompIdx, xgn);
                fluidState_.setMoleFraction(nPhaseIdx, wCompIdx, xnw);
                fluidState_.setMoleFraction(nPhaseIdx, gCompIdx, xna);
                fluidState_.setMoleFraction(nPhaseIdx, nCompIdx, xnn);

                Scalar rhoW = FluidSystem::density(fluidState_, wPhaseIdx);
                Scalar rhoG = FluidSystem::density(fluidState_, gPhaseIdx);
                Scalar rhoN = FluidSystem::density(fluidState_, nPhaseIdx);

                fluidState_.setDensity(wPhaseIdx, rhoW);
                fluidState_.setDensity(gPhaseIdx, rhoG);
                fluidState_.setDensity(nPhaseIdx, rhoN);
            }
        }
        else if (phasePresence == gnPhaseOnly) {
            // only gas and NAPL phases are present
            // we have all (partly hypothetical) phase pressures
            // and temperature and the mole fraction of water in
            // the gas phase

            // we have all (partly hypothetical) phase pressures
            // and temperature and the mole fraction of water in
            // the gas phase
            Scalar partPressNAPL = fluidState_.fugacityCoefficient(nPhaseIdx, nCompIdx)*pn_;

            Scalar xgw = priVars[switch1Idx];
            Scalar xgn = partPressNAPL/pg_;
            Scalar xgg = 1.-xgw-xgn;

            // write mole fractions in the fluid state
            fluidState_.setMoleFraction(gPhaseIdx, wCompIdx, xgw);
            fluidState_.setMoleFraction(gPhaseIdx, gCompIdx, xgg);
            fluidState_.setMoleFraction(gPhaseIdx, nCompIdx, xgn);

            // calculate the composition of the remaining phases (as
            // well as the densities of all phases). this is the job
            // of the "ComputeFromReferencePhase" constraint solver
            ComputeFromReferencePhase::solve(fluidState_,
                                             paramCache,
                                             gPhaseIdx,
                                             /*setViscosity=*/true,
                                             /*setInternalEnergy=*/false);
        }
        else if (phasePresence == wnPhaseOnly) {
            // only water and NAPL phases are present
            Scalar pPartialC = fluidState_.fugacityCoefficient(nPhaseIdx,nCompIdx)*pn_;
            Scalar henryC = fluidState_.fugacityCoefficient(wPhaseIdx,nCompIdx)*pw_;

            Scalar xwg = priVars[switch1Idx];
            Scalar xwn = pPartialC/henryC;
            Scalar xww = 1.-xwg-xwn;

            // write mole fractions in the fluid state
            fluidState_.setMoleFraction(wPhaseIdx, wCompIdx, xww);
            fluidState_.setMoleFraction(wPhaseIdx, gCompIdx, xwg);
            fluidState_.setMoleFraction(wPhaseIdx, nCompIdx, xwn);

            // calculate the composition of the remaining phases (as
            // well as the densities of all phases). this is the job
            // of the "ComputeFromReferencePhase" constraint solver
            ComputeFromReferencePhase::solve(fluidState_,
                                             paramCache,
                                             wPhaseIdx,
                                             /*setViscosity=*/true,
                                             /*setInternalEnergy=*/false);
        }
        else if (phasePresence == gPhaseOnly) {
            // only the gas phase is present, gas phase composition is
            // stored explicitly here below.

            const Scalar xgw = priVars[switch1Idx];
            const Scalar xgn = priVars[switch2Idx];
            Scalar xgg = 1 - xgw - xgn;

            // write mole fractions in the fluid state
            fluidState_.setMoleFraction(gPhaseIdx, wCompIdx, xgw);
            fluidState_.setMoleFraction(gPhaseIdx, gCompIdx, xgg);
            fluidState_.setMoleFraction(gPhaseIdx, nCompIdx, xgn);

            // calculate the composition of the remaining phases (as
            // well as the densities of all phases). this is the job
            // of the "ComputeFromReferencePhase" constraint solver ...
            if (useConstraintSolver)
            {
                ComputeFromReferencePhase::solve(fluidState_,
                                                 paramCache,
                                                 gPhaseIdx,
                                                 /*setViscosity=*/true,
                                                 /*setInternalEnergy=*/false);
            }
            // ... or calculated explicitly this way ...
            else {

                // note that the water phase is actually not existing!
                // thus, this is used as phase switch criterion
                Scalar xww = xgw * pg_
                             / (FluidSystem::fugacityCoefficient(fluidState_,
                                                                 wPhaseIdx,wCompIdx)
                                * pw_);
                Scalar xwn = 1.e-10;
                Scalar xwg = 1.e-10;

                // note that the NAPL phase is actually not existing!
                // thus, this is used as phase switch criterion
                Scalar xnn = xgn * pg_
                             / (FluidSystem::fugacityCoefficient(fluidState_,
                                                                 nPhaseIdx,nCompIdx)
                                * pn_);
                Scalar xna = 1.e-10;
                Scalar xnw = 1.e-10;

                fluidState_.setMoleFraction(wPhaseIdx, wCompIdx, xww);
                fluidState_.setMoleFraction(wPhaseIdx, gCompIdx, xwg);
                fluidState_.setMoleFraction(wPhaseIdx, nCompIdx, xwn);
                fluidState_.setMoleFraction(nPhaseIdx, wCompIdx, xnw);
                fluidState_.setMoleFraction(nPhaseIdx, gCompIdx, xna);
                fluidState_.setMoleFraction(nPhaseIdx, nCompIdx, xnn);

                Scalar rhoW = FluidSystem::density(fluidState_, wPhaseIdx);
                Scalar rhoG = FluidSystem::density(fluidState_, gPhaseIdx);
                Scalar rhoN = FluidSystem::density(fluidState_, nPhaseIdx);

                fluidState_.setDensity(wPhaseIdx, rhoW);
                fluidState_.setDensity(gPhaseIdx, rhoG);
                fluidState_.setDensity(nPhaseIdx, rhoN);
            }
        }
        else if (phasePresence == wgPhaseOnly) {
            // only water and gas phases are present
            Scalar xgn = priVars[switch2Idx];
            Scalar partPressH2O = fluidState_.fugacityCoefficient(wPhaseIdx, wCompIdx)*pw_;

            Scalar xgw = partPressH2O/pg_;
            Scalar xgg = 1.-xgn-xgw;

            // write mole fractions in the fluid state
            fluidState_.setMoleFraction(gPhaseIdx, wCompIdx, xgw);
            fluidState_.setMoleFraction(gPhaseIdx, gCompIdx, xgg);
            fluidState_.setMoleFraction(gPhaseIdx, nCompIdx, xgn);

            // calculate the composition of the remaining phases (as
            // well as the densities of all phases). this is the job
            // of the "ComputeFromReferencePhase" constraint solver ...
            if (useConstraintSolver)
            {
                ComputeFromReferencePhase::solve(fluidState_,
                                                 paramCache,
                                                 gPhaseIdx,
                                                 /*setViscosity=*/true,
                                                 /*setInternalEnergy=*/false);
            }
            // ... or calculated explicitly this way ...
            else {
                // actually, it's nothing else than Henry coefficient
                Scalar xwn = xgn * pg_
                             / (FluidSystem::fugacityCoefficient(fluidState_,
                                                                 wPhaseIdx,nCompIdx)
                                * pw_);
                Scalar xwg = xgg * pg_
                             / (FluidSystem::fugacityCoefficient(fluidState_,
                                                                 wPhaseIdx,gCompIdx)
                                * pw_);
                Scalar xww = 1.-xwg-xwn;

                // note that the NAPL phase is actually not existing!
                // thus, this is used as phase switch criterion
                Scalar xnn = xgn * pg_
                             / (FluidSystem::fugacityCoefficient(fluidState_,
                                                                 nPhaseIdx,nCompIdx)
                                * pn_);
                Scalar xna = 1.e-10;
                Scalar xnw = 1.e-10;

                fluidState_.setMoleFraction(wPhaseIdx, wCompIdx, xww);
                fluidState_.setMoleFraction(wPhaseIdx, gCompIdx, xwg);
                fluidState_.setMoleFraction(wPhaseIdx, nCompIdx, xwn);
                fluidState_.setMoleFraction(nPhaseIdx, wCompIdx, xnw);
                fluidState_.setMoleFraction(nPhaseIdx, gCompIdx, xna);
                fluidState_.setMoleFraction(nPhaseIdx, nCompIdx, xnn);

                Scalar rhoW = FluidSystem::density(fluidState_, wPhaseIdx);
                Scalar rhoG = FluidSystem::density(fluidState_, gPhaseIdx);
                Scalar rhoN = FluidSystem::density(fluidState_, nPhaseIdx);

                fluidState_.setDensity(wPhaseIdx, rhoW);
                fluidState_.setDensity(gPhaseIdx, rhoG);
                fluidState_.setDensity(nPhaseIdx, rhoN);
            }
        }
        else
            assert(false); // unhandled phase state

        for (int phaseIdx = 0; phaseIdx < numPhases; ++phaseIdx) {
            // Mobilities
            const Scalar mu =
                FluidSystem::viscosity(fluidState_,
                                       paramCache,
                                       phaseIdx);
            fluidState_.setViscosity(phaseIdx,mu);

            Scalar kr;
            kr = MaterialLaw::kr(materialParams, phaseIdx,
                                 fluidState_.saturation(wPhaseIdx),
                                 fluidState_.saturation(nPhaseIdx),
                                 fluidState_.saturation(gPhaseIdx));
            mobility_[phaseIdx] = kr / mu;
            Valgrind::CheckDefined(mobility_[phaseIdx]);
        }

        // material dependent parameters for NAPL adsorption
        bulkDensTimesAdsorpCoeff_ =
            MaterialLaw::bulkDensTimesAdsorpCoeff(materialParams);

        /* ATTENTION: The conversion to effective diffusion parameters
         *            for the porous media happens at another place!
         */

        // diffusivity coefficents
        diffusionCoefficient_[gPhaseIdx][wCompIdx] =
            FluidSystem::diffusionCoefficient(fluidState_,
                                              paramCache,
                                              gPhaseIdx,
                                              wCompIdx);
        diffusionCoefficient_[gPhaseIdx][nCompIdx] =
            FluidSystem::diffusionCoefficient(fluidState_,
                                              paramCache,
                                              gPhaseIdx,
                                              nCompIdx);
        diffusionCoefficient_[gPhaseIdx][gCompIdx] = 0.0; // dummy, should not be used !

        diffusionCoefficient_[wPhaseIdx][gCompIdx] =
            FluidSystem::diffusionCoefficient(fluidState_,
                                              paramCache,
                                              wPhaseIdx,
                                              gCompIdx);
        diffusionCoefficient_[wPhaseIdx][nCompIdx] =
            FluidSystem::diffusionCoefficient(fluidState_,
                                              paramCache,
                                              wPhaseIdx,
                                              nCompIdx);
        diffusionCoefficient_[wPhaseIdx][wCompIdx] = 0.0; // dummy, should not be used !

        /* no diffusion in NAPL phase considered  at the moment */
        diffusionCoefficient_[nPhaseIdx][nCompIdx] = 0.0;
        diffusionCoefficient_[nPhaseIdx][wCompIdx] = 0.0;
        diffusionCoefficient_[nPhaseIdx][gCompIdx] = 0.0;

        Valgrind::CheckDefined(diffusionCoefficient_);

        // porosity
        porosity_ = problem.spatialParams().porosity(element,
                                                         fvGeometry,
                                                         scvIdx);
        Valgrind::CheckDefined(porosity_);

        // energy related quantities not contained in the fluid state
        asImp_().updateEnergy_(priVars, problem, element, fvGeometry, scvIdx, isOldSol);
        // compute and set the enthalpy
        for (int phaseIdx = 0; phaseIdx < numPhases; ++phaseIdx)
        {
            Scalar h = Implementation::enthalpy_(fluidState_, paramCache, phaseIdx);
            fluidState_.setEnthalpy(phaseIdx, h);
        }
    }

    /*!
     * \brief Returns the phase state for the control volume.
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
    Scalar temperature() const
    { return fluidState_.temperature(/*phaseIdx=*/0); }

    /*!
     * \brief Returns the effective mobility of a given phase within
     *        the control volume.
     *
     * \param phaseIdx The phase index
     */
    Scalar mobility(const int phaseIdx) const
    {
        return mobility_[phaseIdx];
    }

    /*!
     * \brief Returns the effective capillary pressure within the control volume.
     */
    Scalar capillaryPressure() const
    { return fluidState_.capillaryPressure(); }

    /*!
     * \brief Returns the average porosity within the control volume.
     */
    Scalar porosity() const
    { return porosity_; }

    /*!
     * \brief Returns the diffusivity coefficient matrix.
     */
    Dune::FieldMatrix<Scalar, numPhases, numComponents> diffusionCoefficient() const
    { return diffusionCoefficient_; }

    /*!
     * \brief Returns the adsorption information.
     */
    Scalar bulkDensTimesAdsorpCoeff() const
    { return bulkDensTimesAdsorpCoeff_; }


protected:

    static Scalar temperature_(const PrimaryVariables &priVars,
                               const Problem &problem,
                               const Element &element,
                               const FVElementGeometry &fvGeometry,
                               const int scvIdx)
    {
        return problem.temperatureAtPos(fvGeometry.subContVol[scvIdx].global);
    }

    /*!
     * \brief Called by update() to compute the energy related quantities
     */
    void updateEnergy_(const PrimaryVariables &priVars,
                       const Problem &problem,
                       const Element &element,
                       const FVElementGeometry &fvGeometry,
                       const int scvIdx,
                       bool isOldSol)
    { }

    template<class ParameterCache>
    static Scalar enthalpy_(const FluidState& fluidState,
                            const ParameterCache& paramCache,
                            const int phaseIdx)
    {
        return 0;
    }

    Scalar sw_, sg_, sn_, pg_, pw_, pn_;

    Scalar moleFrac_[numPhases][numComponents];
    Scalar massFrac_[numPhases][numComponents];

    Scalar porosity_;        //!< Effective porosity within the control volume
    Scalar mobility_[numPhases];  //!< Effective mobility within the control volume
    Scalar bulkDensTimesAdsorpCoeff_; //!< the basis for calculating adsorbed NAPL
    /* We need a tensor here !! */
    //!< Binary diffusion coefficients of the 3 components in the phases
    Dune::FieldMatrix<Scalar, numPhases, numComponents> diffusionCoefficient_;
    FluidState fluidState_;

private:
    Implementation &asImp_()
    { return *static_cast<Implementation*>(this); }

    const Implementation &asImp_() const
    { return *static_cast<const Implementation*>(this); }
};

template <class TypeTag>
const typename ThreePThreeCVolumeVariables<TypeTag>::Scalar ThreePThreeCVolumeVariables<TypeTag>::R
                 = Constants<typename GET_PROP_TYPE(TypeTag, Scalar)>::R;

} // end namespace

#endif
