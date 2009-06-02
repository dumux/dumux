//$Id:$
/*****************************************************************************
 *   Copyright (C) 2008 by Bernd Flemisch                                    *
 *   Copyright (C) 2008-2009 by Andreas Lauser                               *
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
 * \brief Quantities required by the twophase box model defined on a vertex.
 */
#ifndef DUMUX_2P_VERTEX_DATA_HH
#define DUMUX_2P_VERTEX_DATA_HH

#include <dumux/new_models/tags.hh>
#include "2pproperties.hh"

namespace Dune
{

/*!
 * \brief Contains the quantities which are are constant within a
 *        finite volume in the two-phase model.
 */
template <class TypeTag>
class TwoPVertexData
{
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(Scalar))   Scalar;
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(GridView)) GridView;

    typedef typename GridView::template Codim<0>::Entity Element;
    
    enum {
        numEq         = GET_PROP_VALUE(TypeTag, PTAG(NumEq)),
        numPhases     = GET_PROP_VALUE(TypeTag, PTAG(NumPhases)),

        formulation   = GET_PROP_VALUE(TypeTag, PTAG(Formulation)),

        dim           = GridView::dimension,
        dimWorld      = GridView::dimensionworld,
    };
    
    typedef typename GET_PROP(TypeTag, PTAG(SolutionTypes))     SolutionTypes;
    typedef typename GET_PROP(TypeTag, PTAG(ReferenceElements)) RefElemProp;
    typedef typename RefElemProp::Container                     ReferenceElements;

    typedef typename GET_PROP_TYPE(TypeTag, PTAG(TwoPIndices)) Indices;
    typedef typename SolutionTypes::PrimaryVarVector  PrimaryVarVector;
    typedef Dune::FieldVector<Scalar, numPhases>      PhasesVector;

    typedef Dune::FieldVector<Scalar, dimWorld>  GlobalPosition;
    typedef Dune::FieldVector<Scalar, dim>       LocalPosition;
    
public:
    /*!
     * \brief Update all quantities for a given control volume.
     */
    template <class JacobianImp>
    void update(const PrimaryVarVector &sol,
                const Element          &element,
                int                     vertIdx,
                bool                    isOldSol,
                JacobianImp            &jac) 
    {
        typedef Indices I;
        
        // coordinates of the vertex
        const GlobalPosition &global = element.geometry().corner(vertIdx);
        const LocalPosition   &local =
            ReferenceElements::general(element.type()).position(vertIdx,
                                                                dim);
        
        if (formulation == I::pWsN) {
            satN = sol[I::saturationIdx];
            satW = 1.0 - satN;
            pC =jac.problem().materialLaw().pC(satW,
                                               global,
                                               element,
                                               local);
            pressure[I::wPhase] = sol[I::pressureIdx];
            pressure[I::nPhase] = pressure[I::wPhase] + pC;
        }
        else if (formulation == I::pNsW) {
            satW = sol[I::saturationIdx];
            satN = 1.0 - satW;
            pC =jac.problem().materialLaw().pC(satW,
                                               global,
                                               element,
                                               local);
            pressure[I::nPhase] = sol[I::pressureIdx];
            pressure[I::wPhase] = pressure[I::nPhase] - pC;
        }
        
        density[I::wPhase] =jac.problem().wettingPhase().density(jac.temperature(sol),
                                                                  pressure[I::wPhase]);
        density[I::nPhase] =jac.problem().nonwettingPhase().density(jac.temperature(sol),
                                                                     pressure[I::nPhase]);
        
        mobility[I::wPhase] =jac.problem().materialLaw().mobW(satW,
                                                               global,
                                                               element,
                                                               local,
                                                               jac.temperature(sol),
                                                               pressure[I::wPhase]);
        mobility[I::nPhase] =jac.problem().materialLaw().mobN(satN,
                                                               global,
                                                               element,
                                                               local,
                                                               jac.temperature(sol),
                                                               pressure[I::nPhase]);

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
