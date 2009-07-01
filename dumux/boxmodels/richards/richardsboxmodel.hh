// $Id$
/*****************************************************************************
 *   Copyright (C) 2007 by Peter Bastian                                     *
 *   Institute of Parallel and Distributed System                            *
 *   Department Simulation of Large Systems                                  *
 *   University of Stuttgart, Germany                                        *
 *                                                                           *
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
#ifndef DUMUX_RICHARDS_BOX_MODEL_HH
#define DUMUX_RICHARDS_BOX_MODEL_HH

#include <dumux/boxmodels/boxscheme/boxscheme.hh>

#include "richardsboxjacobian.hh"
#include "richardsboxproblem.hh"

namespace Dune
{
/*!
 * \ingroup BoxProblems
 * \defgroup RichardsBoxProblems Richards box problems
 */

/*!
 * \ingroup BoxModels
 * \defgroup RichardsBoxModel Richards box model
 */

/*!
 * \ingroup RichardsBoxModel
 * \brief Adaption of the BOX scheme to the isothermal Richards model.
 *
 *
 * In the unsaturated zone, Richards' equation can be used.
 * Gas has resistance against the water flow in porous media.
 * However, viscosity of air is about 1\% of the viscosity of water,
 * which makes it highly mobile compared to the water phase.
 * Therefore, in Richards` equation only water phase with capillary effects are considered,
 * where pressure of the gas phase is set to a reference pressure (\f${p_n}_{ref}\f$).
 *
 * \f{align*}
 * \varrho \hspace{1mm} \phi \hspace{1mm} \frac{\partial S_w}{\partial p_c} \frac{\partial p_c}{\partial t} - \nabla \cdot (\frac{kr_w}{\mu_w} \hspace{1mm} \varrho_w \hspace{1mm} K \hspace{1mm}
 * (\nabla p_w - \varrho_w \hspace{1mm} \vec{g})) \hspace{1mm} = \hspace{1mm} q,
 * \f}
 * where \f$p_w = {p_n}_{ref} - p_c\f$.
 * Here \f$ p_w \f$, \f$ p_c \f$, and \f$ {p_n}_{ref} \f$
 * denote water pressure, capillary pressure, and non-wetting phase reference pressure, repectively.
 *
 * To overcome convergence problems, \f$ \frac{\partial S_w}{\partial p_c} \f$ is taken from the old iteration step.
 *
 */
template<class TypeTag >
class RichardsBoxModel : public BoxScheme<TypeTag,  RichardsBoxModel<TypeTag> >
{
    typedef RichardsBoxModel<TypeTag>      ThisType;
    typedef BoxScheme<TypeTag, ThisType>   ParentType;

    typedef typename GET_PROP_TYPE(TypeTag, PTAG(Problem))        Problem;
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(Scalar))         Scalar;
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(LocalJacobian))  LocalJacobian;

public:
    RichardsBoxModel(Problem &prob)
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
