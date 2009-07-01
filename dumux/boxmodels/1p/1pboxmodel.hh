// $Id$
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

#include <dumux/boxmodels/boxscheme/boxscheme.hh>

#include "1pboxjacobian.hh"
#include "1pboxproblem.hh"

namespace Dune
{
/*!
 * \ingroup BoxProblems
 * \defgroup OnePBoxProblems One-phase box problems
 */

/*!
 * \ingroup BoxModels
 * \defgroup OnePBoxModel One-phase box model
 */

/*!
 * \ingroup OnePBoxModel
 * \brief Adaption of the BOX scheme to the single phase isothermal flow model.
 *
 * Single phase isothermal flow model is implemented for compressible flow.
 * \f{align*}
 * \phi \frac{\partial \varrho}{\partial t} + \vec{\nabla} \cdot (- \varrho \frac{\bar{\bar{K}}}{\mu} ( \nabla p -\varrho \vec{g})) = q
 * \f}
 * which is discretized by this model using vertex
 * centered finite volume (box) scheme as spatial and
 * the implicit Euler method as time discretization.
 * However, the model can also be used for incompressible single phase flow modeling, when in problem file a fluid with constant density is chosen.
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
