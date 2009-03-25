/*****************************************************************************
 *   Copyright (C) 2008 by Klaus Mosthaf, Andreas Lauser, Bernd Flemisch     *
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
#ifndef DUMUX_NEW_2P2C_BOX_MODEL_HH
#define DUMUX_NEW_2P2C_BOX_MODEL_HH

#include <dumux/new_models/2p2c/2p2cboxjacobianbase.hh>

namespace Dune
{

/** \todo Please doc me! */

template<class ProblemT,
         class BoxTraitsT,
         class TwoPTwoCTraitsT>
class TwoPTwoCBoxJacobian : public TwoPTwoCBoxJacobianBase<ProblemT,
                                                           BoxTraitsT,
                                                           TwoPTwoCTraitsT,
                                                           TwoPTwoCBoxJacobian<ProblemT,
                                                                               BoxTraitsT,
                                                                               TwoPTwoCTraitsT> >
{
    typedef TwoPTwoCBoxJacobian<ProblemT,
                                BoxTraitsT,
                                TwoPTwoCTraitsT> ThisType;
    typedef TwoPTwoCBoxJacobianBase<ProblemT,
                                    BoxTraitsT,
                                    TwoPTwoCTraitsT,
                                    ThisType>      ParentType;

    typedef ProblemT                               Problem;

    typedef typename Problem::DomainTraits         DomTraits;
    typedef TwoPTwoCTraitsT                        TwoPTwoCTraits;
    typedef BoxTraitsT                             BoxTraits;

    typedef typename DomTraits::Scalar             Scalar;
    typedef typename DomTraits::Vertex             Vertex;

    typedef typename DomTraits::GlobalPosition     GlobalPosition;
    typedef typename DomTraits::LocalPosition      LocalPosition;

    typedef typename BoxTraits::SolutionVector     SolutionVector;
    typedef typename BoxTraits::FVElementGeometry  FVElementGeometry;

    typedef typename ParentType::VariableVertexData  VariableVertexData;
    typedef typename ParentType::LocalFunction       LocalFunction;
    typedef typename ParentType::ElementData         ElementData;

    typedef typename TwoPTwoCTraits::PhasesVector  PhasesVector;


public:
    TwoPTwoCBoxJacobian(ProblemT &problem)
        : ParentType(problem)
    {
    };


    /*!
     * \brief The storage term of heat
     */
    void heatStorage(SolutionVector &result,
                     int scvId,
                     const LocalFunction &sol,
                     const ElementData &elementCache) const
    {
        // only relevant for the non-isothermal model!
    }

    /*!
     * \brief Update the temperature gradient at a face of a FV
     *        element.
     */
    void updateTempGrad(GlobalPosition &tempGrad,
                        const GlobalPosition &feGrad,
                        const LocalFunction &sol,
                        int vertexIdx) const
    {
        // only relevant for the non-isothermal model!
    }

    /*!
     * \brief Sets the temperature term of the flux vector to the
     *        heat flux due to advection of the fluids.
     */
    void advectiveHeatFlux(SolutionVector &flux,
                           const PhasesVector &darcyOut,
                           Scalar alpha, // upwind parameter
                           const VariableVertexData *upW, // up/downstream verts
                           const VariableVertexData *dnW,
                           const VariableVertexData *upN,
                           const VariableVertexData *dnN) const
    {
        // only relevant for the non-isothermal model!
    }

    /*!
     * \brief Adds the diffusive heat flux to the flux vector over
     *        the face of a sub-control volume.
     */
    void diffusiveHeatFlux(SolutionVector &flux,
                           int faceIdx,
                           const GlobalPosition &tempGrad) const
    {
        // only relevant for the non-isothermal model!
    }

    // internal method!
    static Scalar temperature_(const SolutionVector &sol)
    { return 283.15; /* -> 10 °C */ }
};

///////////////////////////////////////////////////////////////////////////
// TwoPTwoCBoxModel (The actual numerical model.)
///////////////////////////////////////////////////////////////////////////
/**
 * \brief Isothermal two phase two component model with Pw and
 *        Sn/X as primary unknowns.
 *
 * This implements an isothermal two phase two component model
 * with Pw and Sn/X as primary unknowns.
 */
template<class ProblemT, class TwoPTwoCTraitsT = TwoPTwoCPwSnTraits<typename ProblemT::DomainTraits::Scalar> >
class TwoPTwoCBoxModel
    : public BoxScheme<TwoPTwoCBoxModel<ProblemT, TwoPTwoCTraitsT>, // Implementation of the box scheme

                       // The Traits for the BOX method
                       P1BoxTraits<typename ProblemT::DomainTraits::Scalar,
                                   typename ProblemT::DomainTraits::Grid,
                                   TwoPTwoCTraitsT::numEq>,

                       // The actual problem we would like to solve
                       ProblemT,
                       // The local jacobian operator
                       TwoPTwoCBoxJacobian<ProblemT,
                                           P1BoxTraits<typename ProblemT::DomainTraits::Scalar,
                                                       typename ProblemT::DomainTraits::Grid,
                                                       TwoPTwoCTraitsT::numEq>, TwoPTwoCTraitsT > >
{
    typedef typename ProblemT::DomainTraits::Grid       Grid;
    typedef typename ProblemT::DomainTraits::Scalar     Scalar;
    typedef TwoPTwoCBoxModel<ProblemT,TwoPTwoCTraitsT>  ThisType;

public:
    typedef TwoPTwoCTraitsT                                  TwoPTwoCTraits;
    typedef P1BoxTraits<Scalar, Grid, TwoPTwoCTraits::numEq> BoxTraits;

private:
    typedef TwoPTwoCBoxJacobian<ProblemT, BoxTraits, TwoPTwoCTraits>  TwoPTwoCLocalJacobian;
    typedef BoxScheme<ThisType,
                      BoxTraits,
                      ProblemT,
                      TwoPTwoCLocalJacobian>        ParentType;

    typedef typename ProblemT::DomainTraits              DomTraits;
    typedef typename DomTraits::Element                  Element;
    typedef typename DomTraits::Vertex                   Vertex;
    typedef typename DomTraits::ElementIterator          ElementIterator;
    typedef typename DomTraits::LocalPosition            LocalPosition;
    typedef typename DomTraits::GlobalPosition           GlobalPosition;

    enum {
        dim              = DomTraits::dim,
        dimWorld         = DomTraits::dimWorld
    };

public:
    typedef NewNewtonMethod<ThisType> NewtonMethod;

    TwoPTwoCBoxModel(ProblemT &prob)
        : ParentType(prob, twoPTwoCLocalJacobian_),
          twoPTwoCLocalJacobian_(prob)
    {
        Api::require<Api::BasicDomainTraits, typename ProblemT::DomainTraits>();
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
        twoPTwoCLocalJacobian_.updateStaticData(this->currentSolution(),
                                                this->previousSolution());
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
        twoPTwoCLocalJacobian_.addVtkFields(writer, this->currentSolution());
    }

    /*!
     * \brief Calculate the masses in the system for
     *        the current timestep.
     */
    void calculateMass(Dune::FieldVector<Scalar, 4> &mass)
    {
        twoPTwoCLocalJacobian_.calculateMass(this->currentSolution(), mass);
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

        int vertIdx = this->problem().vertexIdx(vert);

        // write phase state
        if (!outStream.good()) {
            DUNE_THROW(IOError,
                       "Could not serialize vertex "
                       << vertIdx);
        }
        outStream << twoPTwoCLocalJacobian_.staticVertexDat_[vertIdx].phaseState
                  << " ";
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

        int vertIdx = this->problem().vertexIdx(vert);

        // read phase state
        if (!inStream.good()) {
            DUNE_THROW(IOError,
                       "Could not deserialize vertex "
                       << vertIdx);
        }

        inStream >> twoPTwoCLocalJacobian_.staticVertexDat_[vertIdx].phaseState;
        twoPTwoCLocalJacobian_.staticVertexDat_[vertIdx].oldPhaseState
            = twoPTwoCLocalJacobian_.staticVertexDat_[vertIdx].phaseState;
    };


private:
    // calculates the jacobian matrix at a given position
    TwoPTwoCLocalJacobian  twoPTwoCLocalJacobian_;
};
}

#endif
