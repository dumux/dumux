// $Id$
/*****************************************************************************
 *   Copyright (C) 2008,2009 by Klaus Mosthaf,                               *
 *                              Andreas Lauser,                              *
 *                              Bernd Flemisch                               *
 *   Institute of Hydraulic Engineering                                      *
 *   University of Stuttgart, Germany                                        *
 *   email: <givenname>.<name>@iws.uni-stuttgart.de                          *
 *                                                                           *
 *   This program is free software; you can redistribute it and/or modify    *
 *   it under the terms of the GNU General Public License as published by    *
 *   the Free Software Foundation; either version 2 of the License, or       *
 *   (at your option) any later version, as long as this copyright notice    *
 *   is included in its original form.                                       *
 *                                                                           *
 *   This program is distributed WITHOUT ANY WARRANTY.                       *
 *****************************************************************************/
/*!
 * \file
 *
 * \brief Contains the quantities which are constant within a
 *        finite volume in the two-phase, two-component model.
 */
#ifndef DUMUX_2P2C_VOLUME_VARIABLES_HH
#define DUMUX_2P2C_VOLUME_VARIABLES_HH

#include <dumux/boxmodels/common/boxmodel.hh>
#include <dumux/common/math.hh>

#include <dune/common/collectivecommunication.hh>

#include <dumux/material/fluidstates/equilibriumfluidstate.hh>
#include <dumux/material/constraintsolvers/compositionfromfugacities.hh>

#include <vector>
#include <iostream>

#include "2p2cproperties.hh"
#include "2p2cindices.hh"

namespace Dumux
{

/*!
 * \ingroup TwoPTwoCModel
 * \brief Contains the quantities which are are constant within a
 *        finite volume in the two-phase, two-component model.
 */
template <class TypeTag>
class TwoPTwoCVolumeVariables : public BoxVolumeVariables<TypeTag>
{
    typedef BoxVolumeVariables<TypeTag> ParentType;
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(VolumeVariables)) Implementation;
    
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(Scalar)) Scalar;
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(GridView)) GridView;
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(Problem)) Problem;
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(FVElementGeometry)) FVElementGeometry;
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(PrimaryVariables)) PrimaryVariables;
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(FluidSystem)) FluidSystem;
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(MaterialLaw)) MaterialLaw;
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(MaterialLawParams)) MaterialLawParams;

    typedef typename GET_PROP_TYPE(TypeTag, PTAG(TwoPTwoCIndices)) Indices;
    enum {
        dim = GridView::dimension,

        numPhases = GET_PROP_VALUE(TypeTag, PTAG(NumPhases)),
        formulation = GET_PROP_VALUE(TypeTag, PTAG(Formulation)),

        // component indices
        lCompIdx = Indices::lCompIdx,
        gCompIdx = Indices::gCompIdx,

        // phase indices
        lPhaseIdx = Indices::lPhaseIdx,
        gPhaseIdx = Indices::gPhaseIdx,

        // primary variable indices
        pressureIdx = Indices::pressureIdx,
        switchIdx = Indices::switchIdx,

        // phase states
        lPhaseOnly = Indices::lPhaseOnly,
        gPhaseOnly = Indices::gPhaseOnly,
        bothPhases = Indices::bothPhases,
        
        // formulations
        plSg = TwoPTwoCFormulation::plSg,
        pgSl = TwoPTwoCFormulation::pgSl,
    };

    typedef typename GridView::template Codim<0>::Entity Element;

public:
    //! The return type of the fluidState() method
    typedef Dumux::EquilibriumFluidState<Scalar, FluidSystem> FluidState;

    /*!
     * \brief Update all quantities for a given control volume.
     *
     * \param priVars The primary variables
     * \param problem The problem
     * \param element The element
     * \param elemGeom The finite-volume geometry in the box scheme
     * \param scvIdx The local index of the SCV (sub-control volume)
     * \param isOldSol Evaluate function with solution of current or previous time step
     */
    void update(const PrimaryVariables &priVars,
                const Problem &problem,
                const Element &element,
                const FVElementGeometry &elemGeom,
                int scvIdx,
                bool isOldSol)
    {
        ParentType::update(priVars,
                           problem,
                           element,
                           elemGeom,
                           scvIdx,
                           isOldSol);

        typename FluidSystem::MutableParameters mutParams;
        typename FluidSystem::MutableParameters::FluidState &fs 
            = mutParams.fluidState();      

        int globalVertIdx = problem.model().vertexMapper().map(element, scvIdx, dim);
        int phasePresence = problem.model().phasePresence(globalVertIdx, isOldSol);

        /////////////
        // set the phase saturations
        /////////////
        if (phasePresence == gPhaseOnly) {
            fs.setSaturation(lPhaseIdx, 0.0);
            fs.setSaturation(gPhaseIdx, 1.0);
        }
        else if (phasePresence == lPhaseOnly) {
            fs.setSaturation(lPhaseIdx, 1.0);
            fs.setSaturation(gPhaseIdx, 0.0);
        }
        else if (phasePresence == bothPhases) {
            Scalar Sg;
            if (formulation == plSg)
                Sg = priVars[switchIdx];
            else if (formulation == pgSl)
                Sg = 1.0 - priVars[switchIdx];
            
            fs.setSaturation(lPhaseIdx, 1 - Sg);
            fs.setSaturation(gPhaseIdx, Sg);
        }

        /////////////
        // set the phase temperatures
        /////////////
        // update the temperature part of the energy module
        Scalar T = asImp_().getTemperature(priVars,
                                             element,
                                             elemGeom,
                                             scvIdx,
                                             problem);
        Valgrind::CheckDefined(T);
        for (int i = 0; i < numPhases; ++i)
            fs.setTemperature(i, T);

        /////////////
        // set the phase pressures
        /////////////

        // capillary pressure parameters
        const MaterialLawParams &materialParams =
            problem.spatialParameters().materialLawParams(element, elemGeom, scvIdx);
        Scalar pC = MaterialLaw::pC(materialParams, fs.saturation(lPhaseIdx));
        if (formulation == plSg) {
            fs.setPressure(lPhaseIdx, priVars[pressureIdx]);
            fs.setPressure(gPhaseIdx, priVars[pressureIdx] + pC);
        }
        else if (formulation == pgSl){
            fs.setPressure(lPhaseIdx, priVars[pressureIdx] - pC);
            fs.setPressure(gPhaseIdx, priVars[pressureIdx]);
        }
        Valgrind::CheckDefined(fs.pressure(lPhaseIdx));
        Valgrind::CheckDefined(fs.pressure(gPhaseIdx));
        
        // update the mutable parameters for the pure components
        mutParams.updatePure(lPhaseIdx);
        mutParams.updatePure(gPhaseIdx);

        /////////////
        // set the phase compositions. 
        /////////////
        if (phasePresence == gPhaseOnly) {
            // mass fractions
            Scalar Xg1 = priVars[switchIdx];
            Scalar Xg2 = 1 - Xg1;

            // molar masses
            Scalar M1 = FluidSystem::molarMass(lCompIdx);
            Scalar M2 = FluidSystem::molarMass(gCompIdx);

            // convert mass to mole fractions
            Scalar meanM = M1*M2/(M2 + Xg2*(M1 - M2));
            fs.setMoleFrac(gPhaseIdx, lCompIdx, Xg1 * M1/meanM);
            fs.setMoleFrac(gPhaseIdx, gCompIdx, Xg2 * M2/meanM);
            mutParams.updateMix(gPhaseIdx);
            
            // calculate component fugacities in gas phase
            Scalar fug1 = FluidSystem::computeFugacity(mutParams, gPhaseIdx, lCompIdx);
            Scalar fug2 = FluidSystem::computeFugacity(mutParams, gPhaseIdx, gCompIdx);
            fs.setFugacity(gPhaseIdx, lCompIdx, fug1);
            fs.setFugacity(gPhaseIdx, gCompIdx, fug2);
            
            // use same fugacities in liquid phase
            fs.setFugacity(lPhaseIdx, lCompIdx, fug1);
            fs.setFugacity(lPhaseIdx, gCompIdx, fug2);

            // initial guess of liquid composition
            fs.setMoleFrac(lPhaseIdx, lCompIdx, 0.98);
            fs.setMoleFrac(lPhaseIdx, gCompIdx, 0.02);

            // calculate liquid composition from fugacities
            typedef Dumux::CompositionFromFugacities<Scalar, FluidSystem> CompFromFug;
            CompFromFug::run(mutParams, lPhaseIdx);

            Valgrind::CheckDefined(fs.moleFrac(gPhaseIdx, lCompIdx));
            Valgrind::CheckDefined(fs.moleFrac(gPhaseIdx, gCompIdx));
            Valgrind::CheckDefined(fs.moleFrac(lPhaseIdx, lCompIdx));
            Valgrind::CheckDefined(fs.moleFrac(lPhaseIdx, gCompIdx));

            // calculate molar volume of gas phase
            Scalar Vmg = FluidSystem::computeMolarVolume(mutParams, gPhaseIdx);
            fs.setMolarVolume(gPhaseIdx, Vmg);
        }
        else if (phasePresence == lPhaseOnly) {
            // mass fractions
            Scalar Xl2 = priVars[switchIdx];
            Scalar Xl1 = 1 - Xl2;

            // molar masses
            Scalar M1 = FluidSystem::molarMass(lCompIdx);
            Scalar M2 = FluidSystem::molarMass(gCompIdx);

            // convert mass to mole fractions
            Scalar meanM = M1*M2/(M2 + Xl2*(M1 - M2));
            fs.setMoleFrac(lPhaseIdx, lCompIdx, Xl1 * M1/meanM);
            fs.setMoleFrac(lPhaseIdx, gCompIdx, Xl2 * M2/meanM);
            mutParams.updateMix(lPhaseIdx);
            
            // calculate component fugacities in liquid phase
            Scalar fug1 = FluidSystem::computeFugacity(mutParams, lPhaseIdx, lCompIdx);
            Scalar fug2 = FluidSystem::computeFugacity(mutParams, lPhaseIdx, gCompIdx);
            fs.setFugacity(lPhaseIdx, lCompIdx, fug1);
            fs.setFugacity(lPhaseIdx, gCompIdx, fug2);          

            // use same fugacities in gas phase
            fs.setFugacity(gPhaseIdx, lCompIdx, fug1);
            fs.setFugacity(gPhaseIdx, gCompIdx, fug2);

            // initial guess of gas composition
            fs.setMoleFrac(gPhaseIdx, lCompIdx, 0.05);
            fs.setMoleFrac(gPhaseIdx, gCompIdx, 0.95);

            // calculate liquid composition from fugacities
            typedef Dumux::CompositionFromFugacities<Scalar, FluidSystem> CompFromFug;
            CompFromFug::run(mutParams, gPhaseIdx);
            Valgrind::CheckDefined(fs.moleFrac(gPhaseIdx, lCompIdx));
            Valgrind::CheckDefined(fs.moleFrac(gPhaseIdx, gCompIdx));
            Valgrind::CheckDefined(fs.moleFrac(lPhaseIdx, lCompIdx));
            Valgrind::CheckDefined(fs.moleFrac(lPhaseIdx, gCompIdx));

            // calculate molar volume of liquid phase
            Scalar Vml = FluidSystem::computeMolarVolume(mutParams, lPhaseIdx);
            fs.setMolarVolume(lPhaseIdx, Vml);
        }
        else if (phasePresence == bothPhases) {
            // HACK: assume both phases to be an ideal mixture,
            // i.e. the fugacity coefficents do not depend on the
            // composition
            Scalar phi_l1 = FluidSystem::computeFugacityCoeff(mutParams, lPhaseIdx, lCompIdx);
            Scalar phi_l2 = FluidSystem::computeFugacityCoeff(mutParams, lPhaseIdx, gCompIdx);
            Scalar phi_g1 = FluidSystem::computeFugacityCoeff(mutParams, gPhaseIdx, lCompIdx);
            Scalar phi_g2 = FluidSystem::computeFugacityCoeff(mutParams, gPhaseIdx, gCompIdx);
            Scalar pl = fs.pressure(lPhaseIdx);
            Scalar pg = fs.pressure(gPhaseIdx);
            Valgrind::CheckDefined(phi_l1);
            Valgrind::CheckDefined(phi_l2);
            Valgrind::CheckDefined(phi_g1);
            Valgrind::CheckDefined(phi_g2);

            Scalar xg2 = (phi_g2/phi_l1 - pl/pg) / (phi_g1/phi_l1 - phi_g2/phi_l2);
            Scalar xg1 = 1 - xg2;
            Scalar xl2 = (xg2*pg*phi_g2)/(pl*phi_l2);
            Scalar xl1 = 1 - xl2;

            fs.setMoleFrac(lPhaseIdx, lCompIdx, xl1);
            fs.setMoleFrac(lPhaseIdx, gCompIdx, xl2);
            fs.setMoleFrac(gPhaseIdx, lCompIdx, xg1);
            fs.setMoleFrac(gPhaseIdx, gCompIdx, xg2);

            mutParams.updateMix(lPhaseIdx);
            mutParams.updateMix(gPhaseIdx);

            Scalar Vml =  FluidSystem::computeMolarVolume(mutParams, lPhaseIdx);
            Scalar Vmg =  FluidSystem::computeMolarVolume(mutParams, gPhaseIdx);
            fs.setMolarVolume(lPhaseIdx, Vml);
            fs.setMolarVolume(gPhaseIdx, Vmg);
        }
       

        // Mobilities
        Scalar muL = FluidSystem::computeViscosity(mutParams, lPhaseIdx);
        Scalar krL = MaterialLaw::krw(materialParams, fs.saturation(lPhaseIdx));
        mobility_[lPhaseIdx] = krL / muL;
        Valgrind::CheckDefined(mobility_[lPhaseIdx]);

        // ATTENTION: krn requires the liquid saturation as parameter!
        Scalar muG = FluidSystem::computeViscosity(mutParams, gPhaseIdx);
        Scalar krG = MaterialLaw::krn(materialParams, fs.saturation(lPhaseIdx));
        mobility_[gPhaseIdx] = krG / muG;
        Valgrind::CheckDefined(mobility_[gPhaseIdx]);

#if 0
        // binary diffusion coefficents
        diffCoeff_[phaseIdx] =
            FluidSystem::diffCoeff(phaseIdx,
                                   lCompIdx,
                                   gCompIdx,
                                   fluidState_.temperature(),
                                   fluidState_.phasePressure(phaseIdx),
                                   fluidState_);
        Valgrind::CheckDefined(diffCoeff_[phaseIdx]);
#endif
        
        // porosity
        porosity_ = problem.spatialParameters().porosity(element,
                                                         elemGeom,
                                                         scvIdx);
        Valgrind::CheckDefined(porosity_);

        asImp_().updateEnergy(mutParams, priVars, element, elemGeom, scvIdx, problem);

        // assign the equilibrium fluid state from the generic one
        fluidState_.assign(fs);
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
     * \brief Returns the mass density of a given phase within the
     *        control volume.
     *
     * \param phaseIdx The phase index
     */
    Scalar molarDensity(int phaseIdx) const
    { return fluidState_.molarDensity(phaseIdx); }

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
     * temperature of the rock matrix and of all fluid phases are
     * identical.
     */
    Scalar temperature() const
    { return fluidState_.temperature(); }

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
    { return fluidState_.pressure(lPhaseIdx) - fluidState_.pressure(gPhaseIdx); }

    /*!
     * \brief Returns the average porosity within the control volume.
     */
    Scalar porosity() const
    { return porosity_; }

#if 0
    /*!
     * \brief Returns the binary diffusion coefficients for a phase
     */
    Scalar diffCoeff(int phaseIdx) const
    { return diffCoeff_[phaseIdx]; }
#endif

protected:

    Scalar getTemperature_(const PrimaryVariables &priVars,
                         const Element &element,
                         const FVElementGeometry &elemGeom,
                         int scvIdx,
                         const Problem &problem)
    {
        return problem.temperature(element, elemGeom, scvIdx);
    }

    Scalar porosity_;        //!< Effective porosity within the control volume
    Scalar mobility_[numPhases];  //!< Effective mobility within the control volume
//    Scalar diffCoeff_[numPhases]; //!< Binary diffusion coefficients for the phases
    FluidState fluidState_;

private:
    Implementation &asImp_()
    { return *static_cast<Implementation*>(this); }

    const Implementation &asImp_() const
    { return *static_cast<const Implementation*>(this); }
};

} // end namepace

#endif
