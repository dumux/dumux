// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*****************************************************************************
 *   Copyright (C) 2011-2012 by Holger Class                                 *
 *   Copyright (C) 2008-2009 by Klaus Mosthaf                                *
 *   Copyright (C) 2008-2009 by Andreas Lauser                               *
 *   Copyright (C) 2008-2009 by Bernd Flemisch                               *
 *   Institute of Hydraulic Engineering                                      *
 *   University of Stuttgart, Germany                                        *
 *   email: <givenname>.<name>@iws.uni-stuttgart.de                          *
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
 * \brief Contains the quantities which are constant within a
 *        finite volume in the two-phase, two-component model.
 */
#ifndef DUMUX_3P3C_VOLUME_VARIABLES_HH
#define DUMUX_3P3C_VOLUME_VARIABLES_HH

#include <dumux/boxmodels/common/boxmodel.hh>
#include <dumux/common/math.hh>

#include <dune/common/collectivecommunication.hh>
#include <vector>
#include <iostream>

#include "3p3cproperties.hh"

#include <dumux/material/constants.hh>
#include <dumux/material/fluidstates/compositionalfluidstate.hh>
#include <dumux/material/constraintsolvers/computefromreferencephase.hh>
#include <dumux/material/constraintsolvers/misciblemultiphasecomposition.hh>

namespace Dumux
{

/*!
 * \ingroup ThreePThreeCModel
 * \brief Contains the quantities which are are constant within a
 *        finite volume in the two-phase, two-component model.
 */
template <class TypeTag>
class ThreePThreeCVolumeVariables : public BoxVolumeVariables<TypeTag>
{
    typedef BoxVolumeVariables<TypeTag> ParentType;
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

    typedef typename GET_PROP_TYPE(TypeTag, ThreePThreeCIndices) Indices;
    enum {
        dim = GridView::dimension,

        numPhases = GET_PROP_VALUE(TypeTag, NumPhases),
        numComponents = GET_PROP_VALUE(TypeTag, NumComponents),

        wCompIdx = Indices::wCompIdx,
        aCompIdx = Indices::aCompIdx,
        cCompIdx = Indices::cCompIdx,

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

public:
    //! The type of the object returned by the fluidState() method
    typedef Dumux::CompositionalFluidState<Scalar, FluidSystem> FluidState;


    /*!
     * \brief Update all quantities for a given control volume.
     *
     * \param primaryVars The primary variables
     * \param problem The problem
     * \param element The element
     * \param elemGeom The finite-volume geometry in the box scheme
     * \param scvIdx The local index of the SCV (sub-control volume)
     * \param isOldSol Evaluate function with solution of current or previous time step
     */
    void update(const PrimaryVariables &primaryVars,
                const Problem &problem,
                const Element &element,
                const FVElementGeometry &elemGeom,
                int scvIdx,
                bool isOldSol)
    {
        ParentType::update(primaryVars,
                           problem,
                           element,
                           elemGeom,
                           scvIdx,
                           isOldSol);

        // capillary pressure parameters
        const MaterialLawParams &materialParams =
            problem.spatialParameters().materialLawParams(element, elemGeom, scvIdx);

        int globalVertIdx = problem.model().dofMapper().map(element, scvIdx, dim);
        int phasePresence = problem.model().phasePresence(globalVertIdx, isOldSol);

        Scalar temp = Implementation::temperature_(primaryVars, problem, element, elemGeom, scvIdx);
        fluidState_.setTemperature(temp);

        /* first the saturations */
        if (phasePresence == threePhases)
        {
            Sw_ = primaryVars[switch1Idx];
            Sn_ = primaryVars[switch2Idx];
            Sg_ = 1. - Sw_ - Sn_;
        }
        else if (phasePresence == wPhaseOnly)
        {
            Sw_ = 1.;
            Sn_ = 0.;
            Sg_ = 0.;
        }
        else if (phasePresence == gnPhaseOnly)
        {
            Sw_ = 0.;
            Sn_ = primaryVars[switch2Idx];
            Sg_ = 1. - Sn_;
        }
        else if (phasePresence == wnPhaseOnly)
        {
            Sn_ = primaryVars[switch2Idx];
            Sw_ = 1. - Sn_;
            Sg_ = 0.;
        }
        else if (phasePresence == gPhaseOnly)
        {
            Sw_ = 0.;
            Sn_ = 0.;
            Sg_ = 1.;
        }
        else if (phasePresence == wgPhaseOnly)
        {
            Sw_ = primaryVars[switch1Idx];
            Sn_ = 0.;
            Sg_ = 1. - Sw_;
        }
        else DUNE_THROW(Dune::InvalidStateException, "phasePresence: " << phasePresence << " is invalid.");
        Valgrind::CheckDefined(Sg_);

        fluidState_.setSaturation(wPhaseIdx, Sw_);
        fluidState_.setSaturation(gPhaseIdx, Sg_);
        fluidState_.setSaturation(nPhaseIdx, Sn_);

        /* now the pressures */
        pg_ = primaryVars[pressureIdx];

        // calculate capillary pressures
        Scalar pCGW = MaterialLaw::pCGW(materialParams, Sw_);
        Scalar pCNW = MaterialLaw::pCNW(materialParams, Sw_);
        Scalar pCGN = MaterialLaw::pCGN(materialParams, Sw_ + Sn_);

        Scalar pcAlpha = MaterialLaw::pCAlpha(materialParams, Sn_);
        Scalar pcNW1 = 0.0; // TODO: this should be possible to assign in the problem file

        pn_ = pg_- pcAlpha * pCGN - (1.-pcAlpha)*(pCGW - pcNW1);
        pw_ = pn_ - pcAlpha * pCNW - (1.-pcAlpha)*pcNW1;

        fluidState_.setPressure(wPhaseIdx, pw_);
        fluidState_.setPressure(gPhaseIdx, pg_);
        fluidState_.setPressure(nPhaseIdx, pn_);

        // calculate and set all fugacity coefficients. this is
        // possible because we require all phases to be an ideal
        // mixture, i.e. fugacity coefficients are not supposed to
        // depend on composition!
        typename FluidSystem::ParameterCache paramCache;
        assert(FluidSystem::isIdealGas(gPhaseIdx));
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
            // constraint solver
            MiscibleMultiPhaseComposition::solve(fluidState_,
                                                 paramCache,
                                                 /*setViscosity=*/true,
                                                 /*setInternalEnergy=*/false);
        }
        else if (phasePresence == wPhaseOnly) {
            // only the water phase is present, water phase composition is
            // stored explicitly.          

            // extract mole fractions in the water phase
            Scalar xwa = primaryVars[switch1Idx];
            Scalar xwc = primaryVars[switch2Idx];
            Scalar xww = 1 - xwa - xwc;

            // write water mole fractions in the fluid state
            fluidState_.setMoleFraction(wPhaseIdx, wCompIdx, xww);
            fluidState_.setMoleFraction(wPhaseIdx, aCompIdx, xwa);
            fluidState_.setMoleFraction(wPhaseIdx, cCompIdx, xwc);

            // calculate the composition of the remaining phases (as
            // well as the densities of all phases). this is the job
            // of the "ComputeFromReferencePhase" constraint solver
            ComputeFromReferencePhase::solve(fluidState_, 
                                             paramCache,
                                             wPhaseIdx,
                                             /*setViscosity=*/true,
                                             /*setInternalEnergy=*/false);
        }
        else if (phasePresence == gnPhaseOnly) {
            // only gas and NAPL phases are present
            // we have all (partly hypothetical) phase pressures
            // and temperature and the mole fraction of water in
            // the gas phase

            // we have all (partly hypothetical) phase pressures
            // and temperature and the mole fraction of water in
            // the gas phase
            Scalar partPressNAPL = fluidState_.fugacityCoefficient(nPhaseIdx, cCompIdx)*pn_;

            Scalar xgw = primaryVars[switch1Idx];
            Scalar xgc = partPressNAPL/pg_;
            Scalar xga = 1.-xgw-xgc;

            // write mole fractions in the fluid state
            fluidState_.setMoleFraction(gPhaseIdx, wCompIdx, xgw);
            fluidState_.setMoleFraction(gPhaseIdx, aCompIdx, xga);
            fluidState_.setMoleFraction(gPhaseIdx, cCompIdx, xgc);

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
            Scalar pPartialC = fluidState_.fugacityCoefficient(nPhaseIdx,cCompIdx)*pn_;
            Scalar henryC = fluidState_.fugacityCoefficient(wPhaseIdx,cCompIdx)*pw_;
            
            Scalar xwa = primaryVars[switch1Idx];
            Scalar xwc = pPartialC/henryC;
            Scalar xww = 1.-xwa-xwc;

            // write mole fractions in the fluid state
            fluidState_.setMoleFraction(wPhaseIdx, wCompIdx, xww);
            fluidState_.setMoleFraction(wPhaseIdx, aCompIdx, xwa);
            fluidState_.setMoleFraction(wPhaseIdx, cCompIdx, xwc);

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

            const Scalar xgw = primaryVars[switch1Idx];
            const Scalar xgc = primaryVars[switch2Idx];
            Scalar xga = 1 - xgw - xgc;

            // write mole fractions in the fluid state
            fluidState_.setMoleFraction(gPhaseIdx, wCompIdx, xgw);
            fluidState_.setMoleFraction(gPhaseIdx, aCompIdx, xga);
            fluidState_.setMoleFraction(gPhaseIdx, cCompIdx, xgc);

            // calculate the composition of the remaining phases (as
            // well as the densities of all phases). this is the job
            // of the "ComputeFromReferencePhase" constraint solver
            ComputeFromReferencePhase::solve(fluidState_,
                                             paramCache,
                                             gPhaseIdx,
                                             /*setViscosity=*/true,
                                             /*setInternalEnergy=*/false);
        }
        else if (phasePresence == wgPhaseOnly) {
            // only water and gas phases are present
            Scalar xgc = primaryVars[switch2Idx];
            Scalar partPressH2O = fluidState_.fugacityCoefficient(wPhaseIdx, wCompIdx)*pw_;

            Scalar xgw = partPressH2O/pg_;
            Scalar xga = 1.-xgc-xgw;

            // write mole fractions in the fluid state
            fluidState_.setMoleFraction(gPhaseIdx, wCompIdx, xgw);
            fluidState_.setMoleFraction(gPhaseIdx, aCompIdx, xga);
            fluidState_.setMoleFraction(gPhaseIdx, cCompIdx, xgc);

            // calculate the composition of the remaining phases (as
            // well as the densities of all phases). this is the job
            // of the "ComputeFromReferencePhase" constraint solver
            ComputeFromReferencePhase::solve(fluidState_,
                                             paramCache,
                                             gPhaseIdx,
                                             /*setViscosity=*/true,
                                             /*setInternalEnergy=*/false);
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
        diffusionCoefficient_[gPhaseIdx][cCompIdx] =
            FluidSystem::diffusionCoefficient(fluidState_,
                                              paramCache,
                                              gPhaseIdx,
                                              cCompIdx);
        diffusionCoefficient_[gPhaseIdx][aCompIdx] = 0.0; // dummy, should not be used !

        diffusionCoefficient_[wPhaseIdx][aCompIdx] =
            FluidSystem::diffusionCoefficient(fluidState_,
                                              paramCache,
                                              wPhaseIdx,
                                              aCompIdx);
        diffusionCoefficient_[wPhaseIdx][cCompIdx] =
            FluidSystem::diffusionCoefficient(fluidState_,
                                              paramCache,
                                              wPhaseIdx,
                                              cCompIdx);
        diffusionCoefficient_[wPhaseIdx][wCompIdx] = 0.0; // dummy, should not be used !

        /* no diffusion in NAPL phase considered  at the moment */
        diffusionCoefficient_[nPhaseIdx][cCompIdx] = 0.0;
        diffusionCoefficient_[nPhaseIdx][wCompIdx] = 0.0;
        diffusionCoefficient_[nPhaseIdx][aCompIdx] = 0.0;

        Valgrind::CheckDefined(diffusionCoefficient_);

        // porosity
        porosity_ = problem.spatialParameters().porosity(element,
                                                         elemGeom,
                                                         scvIdx);
        Valgrind::CheckDefined(porosity_);

        // permeability
        permeability_ = problem.spatialParameters().intrinsicPermeability(element,
                                                                          elemGeom,
                                                                          scvIdx);
        Valgrind::CheckDefined(permeability_);

        // energy related quantities not contained in the fluid state
        asImp_().updateEnergy_(primaryVars, problem, element, elemGeom, scvIdx, isOldSol);
    }

    /*!
     * \brief Returns the phase state for the control-volume.
     */
    const FluidState &fluidState() const
    { return fluidState_; }

    /*!
     * \brief Returns the effective saturation of a given phase within
     *        the control volume.
     *
     * \param phaseIdx The phase index
     */
    Scalar saturation(int phaseIdx) const
    { return fluidState_.saturation(phaseIdx); }

    /*!
     * \brief Returns the mass density of a given phase within the
     *        control volume.
     *
     * \param phaseIdx The phase index
     */
    Scalar density(int phaseIdx) const
    { return fluidState_.density(phaseIdx); }

    /*!
     * \brief Returns the molar density of a given phase within the
     *        control volume.
     *
     * \param phaseIdx The phase index
     */
    Scalar molarDensity(int phaseIdx) const
    { return fluidState_.density(phaseIdx) / fluidState_.averageMolarMass(phaseIdx); }

    /*!
     * \brief Returns the effective pressure of a given phase within
     *        the control volume.
     *
     * \param phaseIdx The phase index
     */
    Scalar pressure(int phaseIdx) const
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
    Scalar mobility(int phaseIdx) const
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
     * \brief Returns the permeability within the control volume.
     */
    Scalar permeability() const
    { return permeability_; }

    /*!
     * \brief Returns the diffusivity coefficient matrix
     */
    Dune::FieldMatrix<Scalar, numPhases, numComponents> diffusionCoefficient() const
    { return diffusionCoefficient_; }

    /*!
     * \brief Returns the adsorption information
     */
    Scalar bulkDensTimesAdsorpCoeff() const
    { return bulkDensTimesAdsorpCoeff_; }


protected:

    static Scalar temperature_(const PrimaryVariables &primaryVars,
                               const Problem &problem,
                               const Element &element,
                               const FVElementGeometry &elemGeom,
                               int scvIdx)
    {
        return problem.boxTemperature(element, elemGeom, scvIdx);
    }

    /*!
     * \brief Called by update() to compute the energy related quantities
     */
    void updateEnergy_(const PrimaryVariables &sol,
                       const Problem &problem,
                       const Element &element,
                       const FVElementGeometry &elemGeom,
                       int vertIdx,
                       bool isOldSol)
    { }

    Scalar Sw_, Sg_, Sn_, pg_, pw_, pn_;

    Scalar moleFrac_[numPhases][numComponents];
    Scalar massFrac_[numPhases][numComponents];

    Scalar porosity_;        //!< Effective porosity within the control volume
    Scalar permeability_;        //!< Effective porosity within the control volume
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
const typename ThreePThreeCVolumeVariables<TypeTag>::Scalar ThreePThreeCVolumeVariables<TypeTag>::R = Constants<typename GET_PROP_TYPE(TypeTag, Scalar)>::R;

} // end namepace

#endif
