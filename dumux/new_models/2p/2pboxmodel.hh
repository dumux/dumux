//$Id:$
/*****************************************************************************
 *   Copyright (C) 2007 by Peter Bastian                                     *
 *   Institute of Parallel and Distributed System                            *
 *   Department Simulation of Large Systems                                  *
 *   University of Stuttgart, Germany                                        *
 *                                                                           *
 *   Copyright (C) 2008 by Andreas Lauser, Bernd Flemisch                    *
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
#ifndef DUMUX_TWOP_BOX_MODEL_HH
#define DUMUX_TWOP_BOX_MODEL_HH

#include <dumux/new_models/2p/2pboxjacobian.hh>

namespace Dune
{
/*!
 * \brief Adaption of the BOX scheme to the twophase flow model.
 *
 * You can pick the formulation by setting the "Formulation"
 * property. The default is pW-Sn.
 */
template<class TypeTag >
class TwoPBoxModel : public BoxScheme<TypeTag,  TwoPBoxModel<TypeTag> >
{
    typedef TwoPBoxModel<TypeTag>          ThisType;
    typedef BoxScheme<TypeTag, ThisType>   ParentType;

    typedef typename GET_PROP_TYPE(TypeTag, PTAG(Problem))        Problem;
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(Scalar))         Scalar;
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(LocalJacobian))  LocalJacobian;

public:
    TwoPBoxModel(Problem &prob)
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
