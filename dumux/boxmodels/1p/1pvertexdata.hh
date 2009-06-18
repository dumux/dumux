// $Id:$
/*****************************************************************************
 *   Copyright (C) 2008 by Onur Dogan                                        *
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
 * \brief Quantities required by the single-phase box model defined on a vertex.
 */
#ifndef DUMUX_1P_VERTEX_DATA_HH
#define DUMUX_1P_VERTEX_DATA_HH

#include <dumux/boxmodels/tags.hh>

#include "1pproperties.hh"

namespace Dune
{

/*!
 * \ingroup TwoPBoxModel
 * \brief Contains the quantities which are are constant within a
 *        finite volume in the one-phase model.
 */
template <class TypeTag>
class OnePVertexData
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
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(OnePIndices))  Indices;

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
        const GlobalPosition &global = element.geometry().corner(vertIdx);
        const LocalPosition   &local =
            ReferenceElements::general(element.type()).position(vertIdx,
                                                                GridView::dimension);

        pressure = sol[I::pressureIdx];
        density = jac.problem().fluid().density(jac.temperature(sol),
                                                pressure);
        viscosity = jac.problem().fluid().viscosity(jac.temperature(sol),
                                                    pressure);
        // porosity
        porosity = jac.problem().soil().porosity(global,
                                                 element,
                                                 local);
    };

    Scalar pressure;
    Scalar density;
    Scalar viscosity;
    Scalar porosity;
};

}

#endif
