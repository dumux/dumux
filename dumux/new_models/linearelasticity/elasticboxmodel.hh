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
#ifndef DUMUX_ELASTIC_BOX_MODEL_HH
#define DUMUX_ELASTIC_BOX_MODEL_HH

#include <dumux/new_models/linearelasticity/elasticboxjacobian.hh>
#include <dumux/new_models/boxscheme/boxscheme.hh>

#include "elasticproperties.hh"

namespace Dune
{
/*!
 * \brief Elastic deformation specific details needed to approximately
 *        calculate the local jacobian in the BOX scheme.
 */
template <class TypeTag >
class ElasticBoxModel : public BoxScheme<TypeTag,  ElasticBoxModel<TypeTag> >
{
    typedef ElasticBoxModel<TypeTag>      ThisType;
    typedef BoxScheme<TypeTag, ThisType>   ParentType;

    typedef typename GET_PROP_TYPE(TypeTag, PTAG(Problem))        Problem;
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(Scalar))         Scalar;
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(LocalJacobian))  LocalJacobian;

public:
    ElasticBoxModel(Problem &prob)
        : ParentType(prob, elasticLocalJacobian_),
          elasticLocalJacobian_(prob)
    {
    }

    /*!
     * \brief Add the mass fraction of air in water to VTK output of
     *        the current timestep.
     */
    template <class MultiWriter>
    void addVtkFields(MultiWriter &writer)
    {
        elasticLocalJacobian_.addVtkFields(writer, this->curSolFunction());
    }

private:
    // calculates the jacobian matrix at a given position
    LocalJacobian  elasticLocalJacobian_;
};
}

#endif
