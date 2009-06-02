/*****************************************************************************
 *   Copyright (C) 2009 by Karin Erbertseder                                 *
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
#ifndef DUMUX_ONEP_TWOC_BOX_MODEL_HH
#define DUMUX_ONEP_TWOC_BOX_MODEL_HH

#include <dumux/new_models/1p2c/1p2cboxjacobian.hh>

namespace Dune
{
/*!
 * \brief Adaption of the BOX scheme to the single-phase, two-component model.
 */
template<class TypeTag >
class OnePTwoCBoxModel : public BoxScheme<TypeTag,  OnePTwoCBoxModel<TypeTag> >
{
    typedef OnePTwoCBoxModel<TypeTag>      ThisType;
    typedef BoxScheme<TypeTag, ThisType>   ParentType;

    typedef typename GET_PROP_TYPE(TypeTag, PTAG(Problem))        Problem;
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(Scalar))         Scalar;
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(LocalJacobian))  LocalJacobian;

public:
    OnePTwoCBoxModel(Problem &prob)
        : ParentType(prob, twoPLocalJacobian_),
          twoPLocalJacobian_(prob)
    {
    }

    /*!
     * \brief Add the mass fraction of air in water to VTK output of
     *        the current timestep.
     */
    template <class MultiWriter>
    void addVtkFields(MultiWriter &writer)
    {
        twoPLocalJacobian_.addVtkFields(writer, this->curSolFunction());
    }

    /*!
     * \brief Calculate the masses in the system for
     *        the current timestep.
     */
    void calculateMass(Dune::FieldVector<Scalar, 2> &mass)
    {
        twoPLocalJacobian_.calculateMass(this->curSolFunction(), mass);
    }



private:
    // calculates the jacobian matrix at a given position
    LocalJacobian  twoPLocalJacobian_;
};
}

#endif
