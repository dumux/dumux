//$Id$
/*****************************************************************************
 *   Copyright (C) 2009 by Melanie Darcis                                    *
 *   Copyright (C) 2009 by Andreas Lauser                                    *
 *   Copyright (C) 2008 by Bernd Flemisch                                    *
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
 * \brief Quantities required by the linear-elasticity model which is
 *        defined on a vertex.
 */
#ifndef DUMUX_ELASTIC_VERTEX_DATA_HH
#define DUMUX_ELASTIC_VERTEX_DATA_HH

#include <dumux/new_models/tags.hh>

#include "elasticproperties.hh"

namespace Dune
{
/*!
 * \brief Contains the quantities which are are constant within a
 *        finite volume in the linear elasticity model.
 */
template <class TypeTag>
class ElasticVertexData
{
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(Scalar))   Scalar;
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(GridView)) GridView;

    typedef typename GridView::template Codim<0>::Entity Element;
    
    enum {
        numEq         = GET_PROP_VALUE(TypeTag, PTAG(NumEq)),

        dim           = GridView::dimension,
        dimWorld      = GridView::dimensionworld,
    };
    
    typedef typename GET_PROP(TypeTag, PTAG(SolutionTypes))     SolutionTypes;
    typedef typename GET_PROP(TypeTag, PTAG(ReferenceElements)) RefElemProp;
    typedef typename RefElemProp::Container                     ReferenceElements;

    typedef typename GET_PROP_TYPE(TypeTag, PTAG(ElasticIndices)) Indices;
    typedef typename SolutionTypes::PrimaryVarVector  PrimaryVarVector;

    typedef Dune::FieldVector<Scalar, dimWorld>  GlobalPosition;
    typedef Dune::FieldVector<Scalar, dim>       LocalPosition;
    typedef Dune::FieldVector<Scalar, dimWorld>  DisplacementVector;
    
public:
    /*!
     * \brief Update all quantities for a given control volume.
     */
    template <class JacobianImp>
    void update(const PrimaryVarVector &sol,
                const Element          &element,
                int                     vertIdx,
                bool                    isOldSol,
                JacobianImp            &jacobian) 
    {
        typedef Indices I;
        
        // coordinates of the vertex
        const GlobalPosition &global = element.geometry().corner(vertIdx);
        const LocalPosition   &local =
            ReferenceElements::general(element.type()).position(vertIdx,
                                                                dim);

        for (int i = 0; i < dim; ++i)
            displacement[i] = sol[I::disp2Primary(i)];
        
        const Dune::FieldVector<Scalar, 2> &lameParams =
            jacobian.problem().soil().lameParams(global, element, local);
        lambda = lameParams[0];
        mu = lameParams[1];
    }
    
    DisplacementVector displacement;
    Scalar lambda;
    Scalar mu;
};

}

#endif
