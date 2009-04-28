/*****************************************************************************
 *   Copyright (C) 2008 by Onur Dogan, Andreas Lauser, Bernd Flemisch                    *
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
template <class RichardsTraits,
          class Problem>
class RichardsVertexData
{
    typedef RichardsTraits Tr;
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

        /* pc = -pw || pc = 0 for computing Sw */

        pW = sol[Tr::pWIdx];
        if (pW >= 0)
            pC = 0.0;
        else
            pC = -pW;

        dSwdpC = jac.problem().materialLaw().dSdP(pC,
                                                  global,
                                                  element,
                                                  local);
        Sw = jac.problem().materialLaw().saturationW(pC,
                                                      global,
                                                      element,
                                                      local);
        mobilityW = jac.problem().materialLaw().mobW(Sw,
                                                     global,
                                                     element,
                                                     local,
                                                     jac.problem().temperature(),
                                                     pW + jac.problem().pNreference());
        densityW = jac.problem().wettingPhase().density(jac.problem().temperature(),
                                                                    pW + jac.problem().pNreference());
        porosity = jac.problem().soil().porosity(global,
                                                 element,
                                                 local);
    }

    Scalar pW;
    Scalar pC;
    Scalar Sw;
    Scalar dSwdpC;

    Scalar densityW;
    Scalar mobilityW;  //FieldVector with the number of phases
    Scalar porosity;
};

}

#endif
