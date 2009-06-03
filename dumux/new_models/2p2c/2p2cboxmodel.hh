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

#include <dumux/new_models/2p2c/2p2cboxjacobian.hh>

namespace Dune
{
/**
 * \brief Isothermal two-phase two-component model.
 *
 * This implements an isothermal two phase two component
 * model. 
 *
 * Depending on the value of the "Formulation" property, the primary
 * variables are either $p_w$ and $S_n;X$ or $p_n$ or $S_w;X$. By
 * default they are $p_w$ and $S_n$
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
     * \brief Add the mass fraction of air in water to VTK output of
     *        the current timestep.
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
