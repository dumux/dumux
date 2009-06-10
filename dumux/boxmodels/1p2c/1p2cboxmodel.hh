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

#include "1p2cboxjacobian.hh"
#include "1p2cboxproblem.hh"

namespace Dune
{

/*!
 * \ingroup BoxProblems
 * \defgroup OnePTwoCBoxProblems One-phase Two-component box problems
 */

/*!
 * \ingroup BoxModels
 * \defgroup OnePTwoCBoxModel One-phase Two-component box model
 */

/*!
 * \ingroup OnePTwoCBoxModel
 * \brief Adaption of the BOX scheme to the one-phase two-component flow model.
 *
 * This model implements an one-phase flow of an incompressible fluid, that consists of two components,
 * using a standard Darcy
 * approach (neglect of the gravitation) as the equation for the conservation of momentum:
 \f[
 v_{D} = - \frac{K}{\mu}
 \left(\text{grad} p  \right)
 \f]
 *
 * By inserting this into the continuity equation, one gets
 \f[
 - \text{div} \left\{
  \frac{K}{\mu}  \left(\text{grad} p \right)
 \right\} = q \;,
 \f]
 *
 * The transport of the components is described by the following equation:
 \f[
 \Phi \frac{ \partial x}{\partial t} - \text{div} \left( \frac{K x}{\mu} \text{grad} p
  + \tau \Phi D \text{grad} x \right) = q.
 \f]
 *
 * All equations are discretized by this model using the fully-coupled vertex
 * centered finite volume (box) scheme as spatial and
 * the implicit Euler method as time discretization.
 *
 * The primary variables are the pressure p and the mole fraction of dissolved component x.
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
