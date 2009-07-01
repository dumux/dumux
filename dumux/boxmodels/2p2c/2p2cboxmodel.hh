// $Id$
/*****************************************************************************
 *   Copyright (C) 2008 by Klaus Mosthaf, Andreas Lauser, Bernd Flemisch     *
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
#ifndef DUMUX_NEW_2P2C_BOX_MODEL_HH
#define DUMUX_NEW_2P2C_BOX_MODEL_HH

#include "2p2cboxjacobian.hh"
#include "2p2cboxproblem.hh"

namespace Dune
{
/*!
 * \ingroup BoxProblems
 * \defgroup TwoPTwoCBoxProblems Two-phase two-component box problems
 */

/*!
 * \ingroup BoxModels
 * \defgroup TwoPTwoCBoxModel Two-phase two-component box model
 */

/*!
 * \ingroup TwoPTwoCBoxModel
 * \brief Adaption of the BOX scheme to the two-phase two-component flow model.
 *
 * This model implements two-phase two-component flow of two compressible and
 * partially miscible fluids \f$\alpha \in \{ w, n \}\f$ composed of the two components
 * \f$\kappa \in \{ w, a \}\f$. The standard multiphase Darcy
 * approach is used as the equation for the conservation of momentum:
 * \f[
     v_\alpha = - \frac{k_{r\alpha}}{\mu_\alpha} K
     \left(\text{grad} p_\alpha - \varrho_{\alpha} \boldsymbol{g} \right)
 * \f]
 *
 * By inserting this into the equations for the conservation of the
 * components, one gets one transport equation for each component
 * \f{eqnarray*}
    &&  \phi \frac{\partial (\sum_\alpha \varrho_\alpha X_\alpha^\kappa S_\alpha )}
	{\partial t}
	- \sum_\alpha \nabla \cdot \left\{ \varrho_\alpha X_\alpha^\kappa
	\frac{k_{r\alpha}}{\mu_\alpha} \mbox{\bf K}
    ({\bf \nabla} p_\alpha - \varrho_{\alpha} \mbox{\bf g}) \right\}
	\nonumber \\ \nonumber \\
    &-& \sum_\alpha \nabla \cdot \left\{{\bf D_{pm}^\kappa} \varrho_{\alpha} {\bf \nabla} X^\kappa_{\alpha} \right\}
    - \sum_\alpha q_\alpha^\kappa = \quad 0 \qquad \kappa \in \{w, a\} \, ,
	\alpha \in \{w, g\}
    \f}
 *
 * This is discretized using a fully-coupled vertex
 * centered finite volume (box) scheme as spatial and
 * the implicit Euler method as temporal discretization.
 *
 * By using constitutive relations for the capillary pressure \f$p_c =
 * p_n - p_w\f$ and relative permeability \f$k_{r\alpha}\f$ and taking
 * advantage of the fact that \f$S_w + S_n = 1\f$ and \f$X^\kappa_w + X^\kappa_n = 1\f$, the number of
 * unknowns can be reduced to two.
 * The used primary variables are, like in the two-phase model, either \f$p_w\f$ and \f$S_n\f$
 * or \f$p_n\f$ and \f$S_w\f$. The formulation which ought to be used can be
 * specified by setting the <tt>Formulation</tt> property to either
 * TwoPTwoCIndices::pWsN or TwoPTwoCIndices::pNsW. By
 * default, the model uses \f$p_w\f$ and \f$S_n\f$.
 * Moreover, the second primary variable depends on the phase state, since a
 * primary variable switch is included. The phase state is stored for all nodes
 * of the system. Following cases can be distinguished:
 * <ul>
 *  <li> Both phases are present: The saturation is used (either \f$S_n\f$ or \f$S_w\f$, dependent on the chosen <tt>Formulation</tt>),
 *  	as long as \f$ 0 < S_\alpha < 1\f$</li>.
 *  <li> Only wetting phase is present: The mass fraction of, e.g., air in the wetting phase \f$X^a_w\f$ is used,
 *  	as long as the maximum mass fraction is not exceeded (\f$X^a_w<X^a_{w,max}\f$)</li>
 *  <li> Only non-wetting phase is present: The mass fraction of, e.g., water in the non-wetting phase, \f$X^w_n\f$, is used,
 *  	as long as the maximum mass fraction is not exceeded (\f$X^w_n<X^w_{n,max}\f$)</li>
 * </ul>
 */

template<class TypeTag, class Implementation >
class TwoPTwoCBoxModelBase
    : public BoxScheme<TypeTag,
                       // Implementation of the box scheme
                       Implementation >
{
    typedef TwoPTwoCBoxModelBase<TypeTag, Implementation>         ThisType;
    typedef BoxScheme<TypeTag, Implementation>                    ParentType;

    typedef typename GET_PROP_TYPE(TypeTag, PTAG(Scalar))        Scalar;
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(Problem))       Problem;
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(GridView))      GridView;
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(LocalJacobian)) LocalJacobian;


    enum {
        dim = GridView::dimension
    };
    typedef typename GridView::template Codim<0>::Entity     Element;
    typedef typename GridView::template Codim<dim>::Entity   Vertex;

public:
    TwoPTwoCBoxModelBase(Problem &prob)
        : ParentType(prob, twoPTwoCLocalJacobian_),
          twoPTwoCLocalJacobian_(prob)
    {
    }

    /*!
     * \brief Called by the update() method if applying the newton
     *         method was unsuccessful.
     */
    void updateFailedTry()
    {
        ParentType::updateFailedTry();

        twoPTwoCLocalJacobian_.setSwitched(false);
        twoPTwoCLocalJacobian_.resetPhaseState();
        twoPTwoCLocalJacobian_.updateStaticData(this->curSolFunction(),
                                                this->prevSolFunction());
    };

    /*!
     * \brief Called by the BoxScheme's update method.
     */
    void updateSuccessful()
    {
        ParentType::updateSuccessful();

        twoPTwoCLocalJacobian_.updateOldPhaseState();
        twoPTwoCLocalJacobian_.setSwitched(false);
    }

    /*!
     * \brief Append all quantities of interest which can be derived
     *        from the solution of the current time step to the VTK
     *        writer.
     */
    template <class MultiWriter>
    void addVtkFields(MultiWriter &writer)
    {
        twoPTwoCLocalJacobian_.addVtkFields(writer, this->curSolFunction());
    }

    /*!
     * \brief Returns true if there was a primary variable switch
     *        after the last time step.
     */
    bool switched() const
    { return twoPTwoCLocalJacobian_.switched(); }

    /*!
     * \brief Write the current solution to a restart file.
     */
    void serializeEntity(std::ostream &outStream,
                         const Vertex &vert)
    {
        // write primary variables
        ParentType::serializeEntity(outStream, vert);

        twoPTwoCLocalJacobian_.serializeEntity(outStream, vert);
    };

    /*!
     * \brief Reads the current solution for a vertex from a restart
     *        file.
     */
    void deserializeEntity(std::istream &inStream,
                           const Vertex &vert)
    {
        // read primary variables
        ParentType::deserializeEntity(inStream, vert);

        twoPTwoCLocalJacobian_.deserializeEntity(inStream, vert);
    };


private:
    // calculates the jacobian matrix at a given position
    LocalJacobian twoPTwoCLocalJacobian_;
};

/**
 * \brief Isothermal two-phase two-component model.
 *
 * This implements an isothermal two phase two component
 * model. This class is just a simple wrapper for \ref TwoPTwoCBoxModelBase .
 */
template<class TypeTag >
class TwoPTwoCBoxModel
    : public TwoPTwoCBoxModelBase<TypeTag, TwoPTwoCBoxModel<TypeTag> >
{
public:
    typedef TwoPTwoCBoxModel<TypeTag>                           ThisType;
    typedef TwoPTwoCBoxModelBase<TypeTag, ThisType>             ParentType;

    typedef typename GET_PROP_TYPE(TypeTag, PTAG(Problem))      Problem;

public:
    TwoPTwoCBoxModel(Problem &prob)
        : ParentType(prob)
    {
    }
};

}

#endif
