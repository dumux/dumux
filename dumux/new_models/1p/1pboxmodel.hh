/*****************************************************************************
 *   Copyright (C) 2009 by Andreas Lauser                                    *
 *   Copyright (C) 2009 by Onur Dogan                                        *
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
#ifndef DUMUX_1P_BOX_MODEL_HH
#define DUMUX_1P_BOX_MODEL_HH

#include <dumux/new_models/boxscheme/boxscheme.hh>

#include "1pboxjacobian.hh"

namespace Dune
{
/*!
 * \brief Adaption of the BOX scheme to the single phase isothermal
 *        flow model.
 */
template<class TypeTag >
class OnePBoxModel : public BoxScheme<TypeTag,  OnePBoxModel<TypeTag> >
{
    typedef OnePBoxModel<TypeTag>          ThisType;
    typedef BoxScheme<TypeTag, ThisType>   ParentType;

    typedef typename GET_PROP_TYPE(TypeTag, PTAG(Problem))        Problem;
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(Scalar))         Scalar;
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(LocalJacobian))  LocalJacobian;

public:
    OnePBoxModel(Problem &prob)
        : ParentType(prob, onePLocalJacobian_),
          onePLocalJacobian_(prob)
    {
    }

    /*!
     * \brief All relevant primary and secondary of the current
     *        solution to an ouput writer.
     */
    template <class MultiWriter>
    void addVtkFields(MultiWriter &writer) 
    {
        onePLocalJacobian_.addVtkFields(writer, this->curSolFunction());
    }

private:
    // calculates the jacobian matrix at a given position
    LocalJacobian  onePLocalJacobian_;
};
}

#endif
