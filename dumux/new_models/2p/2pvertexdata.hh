/*****************************************************************************
 *   Copyright (C) 2008 by Klaus Mosthaf, Andreas Lauser, Bernd Flemisch                    *
 *   Institute of Hydraulic Engineering                                      *
 *   University of Stuttgart, Germany                                        *
 *   email: and _at_ poware.org                                              *
 *                                                                           *
 *   This program is free software; you can redistribute it and/or modify    *
 *   it under the terms of the GNU General Public License as published by    *
 *   the Free Software Foundation; either version 2 of the License, or       *
 *   (at your option) any later version, as long as this copyright notice    *
 *   is included in its original form.                                       *
 *                                                                           *
 *   This program is distributed WITHOUT ANY WARRANTY.                       *
 *****************************************************************************/
#ifndef DUMUX_VERTEX_DATA_HH
#define DUMUX_VERTEX_DATA_HH

#include <dune/common/fvector.hh>

namespace Dune
{
/*!
 * \brief Data which is attached to each vertex of an
 *        element. These quantities are coincidental with the
 *        averaged quantities inside a FV box.
 */
template <class TwoPTraits, 
          class Problem>
class TwoPVertexData
{
    typedef TwoPTraits Tr;
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
        // coordinates of the vertex
        const GlobalPosition &global = element.geometry().corner(vertIdx);
        const LocalPosition   &local =
            Problem::DomainTraits::referenceElement(element.type()).position(vertIdx,
                                                                             Grid::dimension);
            
        if (Tr::formulation == Tr::pWsN) {
            satN = sol[Tr::saturationIdx];
            satW = 1.0 - satN;
            pC =jac.problem().materialLaw().pC(satW,
                                               global,
                                               element,
                                               local);
            pressure[Tr::wPhase] = sol[Tr::pressureIdx];
            pressure[Tr::nPhase] = pressure[Tr::wPhase] + pC;
        }
        else if (Tr::formulation == Tr::pNsW) {
            satW = sol[Tr::saturationIdx];
            satN = 1.0 - satW;
            pC =jac.problem().materialLaw().pC(satW,
                                               global,
                                               element,
                                               local);
            pressure[Tr::nPhase] = sol[Tr::pressureIdx];
            pressure[Tr::wPhase] = pressure[Tr::nPhase] - pC;
        }
        
        density[Tr::wPhase] =jac.problem().wettingPhase().density(jac.temperature(sol),
                                                                  pressure[Tr::wPhase]);
        density[Tr::nPhase] =jac.problem().nonwettingPhase().density(jac.temperature(sol),
                                                                     pressure[Tr::nPhase]);
        
        mobility[Tr::wPhase] =jac.problem().materialLaw().mobW(satW,
                                                               global,
                                                               element,
                                                               local,
                                                               jac.temperature(sol),
                                                               pressure[Tr::wPhase]);
        mobility[Tr::nPhase] =jac.problem().materialLaw().mobN(satN,
                                                               global,
                                                               element,
                                                               local,
                                                               jac.temperature(sol),
                                                               pressure[Tr::nPhase]);

        // porosity
        porosity = jac.problem().soil().porosity(global,
                                                 element,
                                                 local);
    }
        
    Scalar satW;
    Scalar satN;
    Scalar pC;
    Scalar porosity;
        
    PhasesVector density;
    PhasesVector pressure;
    PhasesVector mobility;
};

}

#endif
