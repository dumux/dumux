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
 *   email: and _at_ poware.org                                              *
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

#include <dumux/new_models/2p/2pboxjacobianbase.hh>

namespace Dune
{
/*!
 * \brief Calculate the local Jacobian for the two phase model in the
 *        BOX scheme.
 */
template<class ProblemT, class BoxTraitsT, class TwoPTraitsT>
class TwoPBoxJacobian : public TwoPBoxJacobianBase<ProblemT,
                                           BoxTraitsT,
                                           TwoPTraitsT,
                                           TwoPVertexData<TwoPTraitsT,
                                                              ProblemT>,
                                           TwoPFluxData<TwoPTraitsT,
                                                              ProblemT,
                                                              TwoPVertexData<TwoPTraitsT,
                                                              ProblemT> >,
                                           TwoPBoxJacobian<ProblemT,
                                                              BoxTraitsT,
                                                              TwoPTraitsT> >
{
    typedef TwoPBoxJacobian<ProblemT,
                            BoxTraitsT,
                            TwoPTraitsT>            ThisType;
    typedef TwoPVertexData<TwoPTraitsT, ProblemT>   VertexData;
    typedef TwoPFluxData<TwoPTraitsT,
                             ProblemT,
                             VertexData >          FluxData;
    typedef TwoPBoxJacobianBase<ProblemT,
                                BoxTraitsT,
                                TwoPTraitsT,
                                VertexData,
                                FluxData,
                                ThisType>     ParentType;

public:
    TwoPBoxJacobian(ProblemT &problem)
        : ParentType(problem)
    {
    };
};


///////////////////////////////////////////////////////////////////////////
// TwoPBoxModel (The actual numerical model.)
///////////////////////////////////////////////////////////////////////////
/*!
 * \brief Adaption of the BOX scheme to the pW-Sn twophase flow model.
 */
template<class ProblemT,
         class TwoPTraitsT = PwSnTwoPTraits<typename ProblemT::DomainTraits::Scalar> >
class TwoPBoxModel : public BoxScheme<TwoPBoxModel<ProblemT>, // Implementation

                                      // The Traits for the BOX method
                                      P1BoxTraits<typename ProblemT::DomainTraits::Scalar,
                                                  typename ProblemT::DomainTraits::Grid,
                                                  TwoPTraitsT::numEq>,

                                      // The actual problem we would like to solve
                                      ProblemT,

                                      // The local jacobian operator
                                      TwoPBoxJacobian<ProblemT,
                                                      P1BoxTraits<typename ProblemT::DomainTraits::Scalar,
                                                                  typename ProblemT::DomainTraits::Grid,
                                                                  TwoPTraitsT::numEq>,
                                                      TwoPTraitsT> >
{
    typedef typename ProblemT::DomainTraits::Grid   Grid;
    typedef typename ProblemT::DomainTraits::Scalar Scalar;
    typedef typename ProblemT::DomainTraits::Vertex Vertex;
    typedef TwoPBoxModel<ProblemT>                  ThisType;

public:
    typedef TwoPTraitsT                                  TwoPTraits;
    typedef P1BoxTraits<Scalar, Grid, TwoPTraits::numEq> BoxTraits;

private:
    typedef TwoPBoxJacobian<ProblemT, BoxTraits, TwoPTraits>  TwoPLocalJacobian;
    typedef BoxScheme<ThisType,
                      BoxTraits,
                      ProblemT,
                      TwoPLocalJacobian>  ParentType;

public:
    typedef NewNewtonMethod<ThisType> NewtonMethod;

    TwoPBoxModel(ProblemT &prob)
        : ParentType(prob, twoPLocalJacobian_),
          twoPLocalJacobian_(prob)
    {
        Api::require<Api::BasicDomainTraits, typename ProblemT::DomainTraits>();
    }

    /*!
     * \brief Add the mass fraction of air in water to VTK output of
     *        the current timestep.
     */
    template <class MultiWriter>
    void addVtkFields(MultiWriter &writer)
    {
        twoPLocalJacobian_.addVtkFields(writer, this->currentSolution());
    }

    /*!
     * \brief Calculate the masses in the system for
     *        the current timestep.
     */
    void calculateMass(Dune::FieldVector<Scalar, 2> &mass)
    {
        twoPLocalJacobian_.calculateMass(this->currentSolution(), mass);
    }

        /*!
     * \brief Write the current solution to a restart file.
     */
    void serializeEntity(std::ostream &outStream,
                         const Vertex &vert)
    {
        // write primary variables
        ParentType::serializeEntity(outStream, vert);
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
    };

private:
    // calculates the jacobian matrix at a given position
    TwoPLocalJacobian  twoPLocalJacobian_;
};
}

#endif
