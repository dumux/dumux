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
 * \brief Defines the properties required for the linear elasticity model.
 */

#ifndef DUMUX_ELASTIC_PROPERTIES_HH
#define DUMUX_ELASTIC_PROPERTIES_HH

#include <dumux/new_models/tags.hh>

namespace Dune
{
////////////////////////////////
// forward declarations
////////////////////////////////
template<class TypeTag>
class ElasticBoxModel;

template<class TypeTag>
class ElasticBoxJacobian;

template <class TypeTag>
class ElasticVertexData;

template <class TypeTag>
class ElasticElementData;

template <class TypeTag>
class ElasticFluxData;

/*!
 * \brief Defines the properties required for the linear elasticity model.
 *
 * \tparam PVOffset    The first index in a primary variable vector.
 */
template <int PVOffset = 0>
struct ElasticIndices
{ 
    /*!
     * \brief Convert an index in a displacement vector to an index in a
     *        vector of primary variables.
     */
    static int disp2Primary(int dispIdx)
    {
        return PVOffset + dispIdx;
    };
};

////////////////////////////////
// properties
////////////////////////////////
namespace Properties
{
SET_INT_PROP(BoxElastic, NumEq, GET_PROP_TYPE(TypeTag, PTAG(GridView))::dimension); //!< set the number of equations to the number of dimensions of the grid

//! Use the 2p local jacobian operator for the 2p model
SET_TYPE_PROP(BoxElastic, 
              LocalJacobian,
              ElasticBoxJacobian<TypeTag>);

//! the Model property
SET_TYPE_PROP(BoxElastic, Model, ElasticBoxModel<TypeTag>);

//! the VertexData property
SET_TYPE_PROP(BoxElastic, VertexData, ElasticVertexData<TypeTag>);

//! The indices required by the isothermal 2p model
SET_TYPE_PROP(BoxElastic, ElasticIndices, ElasticIndices<0>);
}

}

#endif

