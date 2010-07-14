// $Id: 2pboxmodel.hh 3738 2010-06-15 14:01:09Z lauser $
/*****************************************************************************
 *   Copyright (C) 2007 by Peter Bastian                                     *
 *   Institute of Parallel and Distributed System                            *
 *   Department Simulation of Large Systems                                  *
 *   University of Stuttgart, Germany                                        *
 *                                                                           *
 *   Copyright (C) 2008-2009 by Andreas Lauser                               *
 *   Copyright (C) 2007-2007 by Bernd Flemisch                               *
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

#include "2pboxjacobian.hh"
#include "2pnewtoncontroller.hh"
#include "2pboxproblem.hh"

namespace Dumux
{

/*!
 * \ingroup BoxProblems
 * \defgroup TwoPBoxProblems Two-phase box problems
 */

/*!
 * \ingroup BoxModels
 * \defgroup TwoPBoxModel Two-phase box model
 */

/*!
 * \ingroup TwoPBoxModel
 * \brief Adaption of the BOX scheme to the twophase flow model.
 *
 * This model implements two-phase flow of two completely immiscible fluids
 * \f$\alpha \in \{ w, n \}\f$ using a standard multiphase Darcy
 * approach as the equation for the conservation of momentum:
 \f[
 v_\alpha = - \frac{k_{r\alpha}}{\mu_\alpha} K
 \left(\text{grad} p_\alpha - \varrho_{\alpha} \boldsymbol{g} \right)
 \f]
 *
 * By inserting this into the equation for the conservation of the
 * phase mass, one gets
 \f[
 \phi \frac{\partial \varrho_\alpha S_\alpha}{\partial t}
 -
 \text{div} \left\{
 \varrho_\alpha \frac{k_{r\alpha}}{\mu_\alpha} K \left(\text{grad} p_\alpha - \varrho_{\alpha} \boldsymbol{g} \right)
 \right\} = q_\alpha \;,
 \f]
 * discretized by a fully-coupled vertex
 * centered finite volume (box) scheme as spatial and
 * the implicit Euler method as time discretization.
 *
 * By using constitutive relations for the capillary pressure \f$p_c =
 * p_n - p_w\f$ and relative permability \f$k_{r\alpha}\f$ and taking
 * advantage of the fact that \f$S_w + S_n = 1\f$, the number of
 * unknowns can be reduced to two. Currently the model supports
 * choosing either \f$p_w\f$ and \f$S_n\f$ or \f$p_n\f$ and \f$S_w\f$
 * as primary variables. The formulation which ought to be used can be
 * specified by setting the <tt>Formulation</tt> property to either
 * <tt>TwoPCommonIndices::pWsN</tt> or <tt>TwoPCommonIndices::pNsW</tt>. By
 * default, the model uses \f$p_w\f$ and \f$S_n\f$.
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
        : ParentType(prob)
    {
    }
    
    /*!
     * \brief Append all quantities of interest which can be derived
     *        from the solution of the current time step to the VTK
     *        writer.
     */
    template <class MultiWriter>
    void addOutputVtkFields(MultiWriter &writer)
    {
        this->localJacobian().addOutputVtkFields(writer, this->curSol());
    }

    /*!
     * \brief Calculate the flux of the nonwetting phase across a given
     * layer for the current timestep
     */
    void calculateFluxAcrossLayer(Dune::FieldVector<Scalar, 2> &flux, int coord, Scalar coordVal)
    {
        this->localJacobian().calculateFluxAcrossLayer(this->curSolFunction(), flux, coord, coordVal);
    }
    
    /*!
     * \brief Calculate the phase masses in the system for the current
     *        timestep.
     */
    void calculateMass(Dune::FieldVector<Scalar, 2> &mass)
    {
        this->localJacobian().calculateMass(this->curSolFunction(), mass);
    }
};
}

#endif
