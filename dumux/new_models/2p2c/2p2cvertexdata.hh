/*****************************************************************************
 *   Copyright (C) 2008,2009 by Klaus Mosthaf,                               *
 *                              Andreas Lauser,                              *
 *                              Bernd Flemisch                               *
 *   Institute of Hydraulic Engineering                                      *
 *   University of Stuttgart, Germany                                        *
 *   email: klaus.mosthaf _at_ iws.uni-stuttgart.de                          *
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
 * \brief Contains the quantities which are are constant within a
 *        finite volume in the two-phase, two-component model.
 */
#ifndef DUMUX_2P2C_VERTEX_DATA_HH
#define DUMUX_2P2C_VERTEX_DATA_HH

#include <dumux/new_models/boxscheme/boxscheme.hh>
#include <dumux/new_models/boxscheme/p1boxtraits.hh>
#include <dumux/new_models/2p2c/2p2ctraits.hh>
#include <dumux/auxiliary/math.hh>

#include <dumux/material/multicomponentrelations.hh>
#include <dumux/auxiliary/apis.hh>
#include <dune/common/collectivecommunication.hh>
#include <vector>
#include <iostream>

namespace Dune
{

/*!
 * \brief Contains the quantities which are are constant within a
 *        finite volume in the two-phase, two-component model.
 */
template <class TwoPTwoCTraits, 
          class Problem>
class TwoPTwoCVertexData
{
    typedef TwoPTwoCTraits Tr;
    typedef typename Problem::DomainTraits::Scalar Scalar;
    typedef typename Problem::DomainTraits::Grid Grid;

    typedef typename Grid::template Codim<0>::Entity Element;
    
    typedef Dune::FieldVector<Scalar, Tr::numEq>      SolutionVector;
    typedef Dune::FieldVector<Scalar, Tr::numPhases>  PhasesVector;

    typedef Dune::FieldVector<Scalar, Grid::dimensionworld>  GlobalPosition;
    typedef Dune::FieldVector<Scalar, Grid::dimension>       LocalPosition;

public:
    /*!
     * \brief Update all quantities for a given control volume.
     */
    template <class JacobianImp>
    void update(const SolutionVector   &sol,
                const Element          &element,
                int                     vertIdx,
                bool                    isOldSol,
                JacobianImp            &jac) 
    {
        const GlobalPosition &global = element.geometry().corner(vertIdx);
        const LocalPosition   &local =
            Problem::DomainTraits::referenceElement(element.type()).position(vertIdx,
                                                                             Grid::dimension);
        int globalVertIdx = jac.problem().vertexIdx(element, vertIdx);
        int phaseState = jac.phaseState(globalVertIdx, isOldSol);
        Scalar temperature = jac.temperature(sol);

        // update saturations, pressures and mass fractions
        updateSaturations(sol, phaseState);
        Scalar pC = jac.problem().materialLaw().pC(saturation[Tr::wPhase],
                                                   global,
                                                   element,
                                                   local, 
                                                   temperature);
        updatePressures(sol, pC);
        updateMassFracs(sol, jac.problem().multicomp(), phaseState, temperature);
        
        // Densities
        density[Tr::wPhase] = jac.problem().wettingPhase().density(temperature,
                                                                   pressure[Tr::wPhase],
                                                                   massfrac[Tr::nComp][Tr::wPhase]);
        density[Tr::nPhase] = jac.problem().nonwettingPhase().density(temperature,
                                                                      pressure[Tr::nPhase],
                                                                      massfrac[Tr::wComp][Tr::nPhase]);
        
        // Mobilities
        mobility[Tr::wPhase] = jac.problem().materialLaw().mobW(saturation[Tr::wPhase],
                                                                global,
                                                                element,
                                                                local,
                                                                temperature,
                                                                pressure[Tr::wPhase]);
        mobility[Tr::nPhase] = jac.problem().materialLaw().mobN(saturation[Tr::nPhase],
                                                                global,
                                                                element,
                                                                local,
                                                                temperature,
                                                                pressure[Tr::nPhase]);
        
        // binary diffusion coefficents
        diffCoeff[Tr::wPhase] = 
            jac.problem().wettingPhase().diffCoeff(temperature, pressure[Tr::wPhase]);
        diffCoeff[Tr::nPhase] = 
            jac.problem().nonwettingPhase().diffCoeff(temperature, pressure[Tr::nPhase]);

        // porosity
        porosity = jac.problem().soil().porosity(global,
                                                 element,
                                                 local);
    }

    /*!
     * \brief Update saturations. 
     */
    void updateSaturations(const SolutionVector &sol, 
                           int phaseState)
    {
        // update saturations
        if (Tr::formulation == Tr::pWsN)
        {
            if (phaseState == Tr::bothPhases) 
                saturation[Tr::nPhase] = sol[Tr::switchIdx];
            else if (phaseState == Tr::wPhaseOnly)
                saturation[Tr::nPhase] = 0.0;
            else if (phaseState == Tr::nPhaseOnly) 
                saturation[Tr::nPhase] = 1.0;
            else
                DUNE_THROW(Dune::InvalidStateException, "Phase state " << phaseState << " is invalid.");

            saturation[Tr::wPhase] = 1.0 - saturation[Tr::nPhase];
        }
        else if (Tr::formulation == Tr::pNsW)
        {
            if (phaseState == Tr::bothPhases)
                saturation[Tr::wPhase] = sol[Tr::switchIdx];
            else if (phaseState == Tr::wPhaseOnly)
                saturation[Tr::wPhase] = 1.0;
            else if (phaseState == Tr::nPhaseOnly) 
                saturation[Tr::wPhase] = 0.0;
            else
                DUNE_THROW(Dune::InvalidStateException, "Phase state " << phaseState << " is invalid.");

            saturation[Tr::nPhase] = 1.0 - saturation[Tr::wPhase];
        }
        else DUNE_THROW(Dune::InvalidStateException, "Formulation: " << Tr::formulation << " is invalid.");
    }

    /*!
     * \brief Update phase pressures. 
     *
     * Requires up to date saturations.
     */
    void updatePressures(const SolutionVector &sol,
                         Scalar capillaryPressure)
    {
        // update pressures
        pC = capillaryPressure;

        if (Tr::formulation == Tr::pWsN) {
            pressure[Tr::wPhase] = sol[Tr::pressureIdx];
            pressure[Tr::nPhase] = pressure[Tr::wPhase] + pC;
        }
        else if (Tr::formulation == Tr::pNsW) {
            pressure[Tr::nPhase] = sol[Tr::pressureIdx];
            pressure[Tr::wPhase] = pressure[Tr::nPhase] - pC;
        }
        else DUNE_THROW(Dune::InvalidStateException, "Formulation: " << Tr::formulation << " is invalid.");
    }

    /*!
     * \brief Update mass fraction matrix.
     *
     * Requires up to date pressures and saturations.
     */
    void updateMassFracs(const SolutionVector &sol,
                         MultiComp &multicomp,
                         int phaseState,
                         Scalar temperature)
    {
        // Solubilities of components in phases
        if (phaseState == Tr::bothPhases) {
            massfrac[Tr::wComp][Tr::nPhase] = multicomp.xWN(pressure[Tr::nPhase], temperature);
            massfrac[Tr::nComp][Tr::wPhase] = multicomp.xAW(pressure[Tr::nPhase], temperature);
        }
        else if (phaseState == Tr::wPhaseOnly) {
            massfrac[Tr::wComp][Tr::nPhase] = 0.0;
            massfrac[Tr::nComp][Tr::wPhase] = sol[Tr::switchIdx];
        }
        else if (phaseState == Tr::nPhaseOnly) {
            massfrac[Tr::wComp][Tr::nPhase] = sol[Tr::switchIdx];
            massfrac[Tr::nComp][Tr::wPhase] = 0.0;
        }
        else DUNE_THROW(Dune::InvalidStateException, "Phase state " << phaseState << " is invalid.");

        massfrac[Tr::wComp][Tr::wPhase] = 1.0 - massfrac[Tr::nComp][Tr::wPhase];
        massfrac[Tr::nComp][Tr::nPhase] = 1.0 - massfrac[Tr::wComp][Tr::nPhase];
    };
    
public:
    PhasesVector saturation;//!< Effective phase saturations within the control volume
    PhasesVector pressure;  //!< Effective phase pressures within the control volume
    Scalar pC;              //!< Effective capillary pressure within the control volume
    Scalar porosity;        //!< Effective porosity within the control volume
    PhasesVector mobility;  //!< Effective mobility within the control volume
    PhasesVector density;   //!< Effective density within the control volume
    PhasesVector diffCoeff; //!< Binary diffusion coefficients for the phases

    //! Mass fractions of each component within each phase
    Dune::FieldMatrix<Scalar, Tr::numComponents, Tr::numPhases> massfrac; 
};

} // end namepace

#endif
