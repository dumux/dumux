//$Id:$
/*****************************************************************************
 *   Copyright (C) 2009 by Onur Dogan                                        *
 *   Copyright (C) 2009 by Andreas Lauser                                    *
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
 * \brief Quantities required by the isothermal richards box model
 *        defined on a vertex.
 */
#ifndef DUMUX_RICHARDS_VERTEX_DATA_HH
#define DUMUX_RICHARDS_VERTEX_DATA_HH

#include <dumux/new_models/tags.hh>

#include "richardsproperties.hh"

namespace Dune
{

/*!
 * \brief Data which is attached to each vertex of the and can be
 *        shared between multiple calculations and should thus be
 *        cached in order to increase efficency.
 */
template <class TypeTag>
class RichardsVertexData
{
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(Scalar))   Scalar;
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(GridView)) GridView;

    typedef typename GridView::template Codim<0>::Entity Element;
    
    enum {
        dim           = GridView::dimension,
        dimWorld      = GridView::dimensionworld,
    };
    
    typedef typename GET_PROP(TypeTag, PTAG(ReferenceElements)) RefElemProp;
    typedef typename RefElemProp::Container                     ReferenceElements;

    typedef typename GET_PROP(TypeTag, PTAG(SolutionTypes))     SolutionTypes;
    typedef typename SolutionTypes::PrimaryVarVector            PrimaryVarVector;
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(RichardsIndices))  Indices;

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
                                                                GridView::dimension);

        /* pc = -pw || pc = 0 for computing Sw */

        pW = sol[I::pWIdx];
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
    Scalar mobilityW;
    Scalar porosity;
};

}

#endif
