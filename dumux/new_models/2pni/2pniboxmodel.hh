//$Id$
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
#ifndef DUMUX_NEW_2PNI_BOX_MODEL_HH
#define DUMUX_NEW_2PNI_BOX_MODEL_HH

#include <dumux/new_models/2p/2pboxjacobianbase.hh>

#include <dumux/new_models/2pni/2pnivertexdata.hh>
#include <dumux/new_models/2pni/2pnifluxdata.hh>
#include <dumux/new_models/2pni/2pnitraits.hh>

namespace Dune
{

///////////////////////////////////////////////////////////////////////////
// TwoPNIBoxJacobian (evaluate the local jacobian for the newton method.)
///////////////////////////////////////////////////////////////////////////

/** \todo Please doc me! */

template<class ProblemT,
         class BoxTraitsT,
         class TwoPNITraitsT>
class TwoPNIBoxJacobian : public TwoPBoxJacobianBase<ProblemT,
                                                     BoxTraitsT,
                                                     TwoPNITraitsT,
                                                     TwoPNIVertexData<TwoPNITraitsT,
                                                     ProblemT>,
                                                     TwoPNIFluxData<TwoPNITraitsT,
                                                     ProblemT,TwoPNIVertexData<TwoPNITraitsT,
                                                     ProblemT> >,
                                                     TwoPNIBoxJacobian<ProblemT,
                                                     BoxTraitsT,
                                                     TwoPNITraitsT> >
{
    typedef TwoPNITraitsT TwoPNITraits;
    typedef TwoPNIBoxJacobian<ProblemT,
                                  BoxTraitsT,
                                  TwoPNITraits>  ThisType;
    typedef TwoPNIVertexData<TwoPNITraitsT,
                                 ProblemT>           VertexData;
    typedef TwoPNIFluxData<TwoPNITraitsT,
                               ProblemT,
                               VertexData >          FluxData;
    typedef TwoPBoxJacobianBase<ProblemT,
                                    BoxTraitsT,
                                    TwoPNITraitsT,
                                    VertexData,
                                    FluxData,
                                    ThisType>     ParentType;


    typedef ProblemT                       Problem;
    typedef typename Problem::DomainTraits DomTraits;
    typedef typename DomTraits::Scalar     Scalar;
    typedef typename DomTraits::Element    Element;

    enum {
        dim            = DomTraits::dim,
        dimWorld       = DomTraits::dimWorld,

        numPhases      = TwoPNITraits::numPhases,

        pressureIdx    = TwoPNITraits::pressureIdx,
        saturationIdx  = TwoPNITraits::saturationIdx,
        temperatureIdx = TwoPNITraits::temperatureIdx,

        wMassIdx       = TwoPNITraits::wMassIdx,
        nMassIdx       = TwoPNITraits::nMassIdx,

        wPhase         = TwoPNITraits::wPhase,
        nPhase         = TwoPNITraits::nPhase,
     };

    typedef BoxTraitsT                              BoxTraits;
    typedef typename BoxTraits::SolutionVector      SolutionVector;

    typedef typename ParentType::LocalFunction       LocalFunction;

    typedef typename ParentType::VertexDataArray     VertexDataArray;

    typedef typename BoxTraits::FVElementGeometry    FVElementGeometry;

    typedef typename DomTraits::GlobalPosition        GlobalPosition;
    typedef typename DomTraits::LocalPosition         LocalPosition;

    typedef FieldMatrix<Scalar, dim, dim>  Tensor;

public:
    TwoPNIBoxJacobian(ProblemT &problem)
        : ParentType(problem)
    {
    };

    /*!
     * \brief Evaluate the amount all conservation quantites
     *        (e.g. phase mass) within a sub-control volume.
     *
     * The result should be averaged over the volume (e.g. phase mass
     * inside a sub control volume divided by the volume)
     */
    void computeStorage(SolutionVector &result, int scvIdx, bool usePrevSol) const
    {
        // compute the storage term for phase mass
        ParentType::computeStorage(result, scvIdx, usePrevSol);

        // if flag usePrevSol is set, the solution from the previous
        // time step is used, otherwise the current solution is
        // used. The secondary variables are used accordingly.  This
        // is required to compute the derivative of the storage term
        // using the implicit euler method.
        const VertexDataArray &elemDat = usePrevSol ? this->prevElemDat_  : this->curElemDat_;
        const VertexData  &vertDat = elemDat[scvIdx];

        result[temperatureIdx] =
            vertDat.porosity*(vertDat.density[wPhase] *
                              vertDat.intEnergy[wPhase] *
                              vertDat.satW
                              +
                              vertDat.density[nPhase] *
                              vertDat.intEnergy[nPhase] *
                              vertDat.satN)
            +
            vertDat.temperature *
            this->problem_.soil().heatCap(this->curElementGeom_.elementGlobal,
                                          this->curElement_(),
                                          this->curElementGeom_.elementLocal);
    }

    /*!
     * \brief Sets the temperature term of the flux vector to the
     *        heat flux due to advection of the fluids.
     */
    void computeFluxdvectiveFlux(SolutionVector &flux,
                              const FluxData &fluxData) const
    {
        // advective mass flux
        ParentType::computeFluxdvectiveFlux(flux, fluxData);

        // advective heat flux in all phases
        const Scalar alpha = TwoPNITraits::upwindAlpha;
        flux[temperatureIdx] = 0;
        for (int phase = 0; phase < numPhases; ++phase) {
            // vertex data of the upstream and the downstream vertices
            const VertexData &up = this->curElemDat_[fluxData.upstreamIdx[phase]];
            const VertexData &dn = this->curElemDat_[fluxData.downstreamIdx[phase]];

            flux[temperatureIdx] +=
                fluxData.vDarcyNormal[phase] * (
                    alpha * // upstream vertex
                    (  up.density[phase] *
                       up.mobility[phase] *
                       up.enthalpy[phase])
                    +
                    (1-alpha) * // downstream vertex
                    (  dn.density[phase] *
                       dn.mobility[phase] *
                       dn.enthalpy[phase]) );
        }
    }

    /*!
     * \brief Adds the diffusive heat flux to the flux vector over
     *        the face of a sub-control volume.
     */
    void computeDiffusiveFlux(SolutionVector &flux,
                              const FluxData &fluxData) const
    {
        // diffusive mass flux
        ParentType::computeDiffusiveFlux(flux, fluxData);

        // diffusive heat flux
        flux[temperatureIdx] += (fluxData.temperatureGrad*fluxData.face->normal)*fluxData.heatCondAtIp;
    }

    // internal method!
    template <class SolutionVector>
    Scalar temperature(const SolutionVector &sol)
    { return sol[temperatureIdx]; }
};

///////////////////////////////////////////////////////////////////////////
// TwoPNIBoxModel (The actual numerical model.)
///////////////////////////////////////////////////////////////////////////
/**
 * \brief Non-isothermal two phase two component model with Pw and
 *        Sn/X as primary unknowns.
 *
 * This implements a non-isothermal two-phase two-component model
 * with Pw and Sn/X as primary unknowns.
 */
template<class ProblemT,
         class TwoPNITraitsT = TwoPNITraits<typename ProblemT::DomainTraits::Scalar,
                                                    PwSnTwoPTraits<typename ProblemT::DomainTraits::Scalar> > >
class TwoPNIBoxModel
    : public BoxScheme<TwoPNIBoxModel<ProblemT,TwoPNITraitsT>, // Implementation of the box scheme

                       // The Traits for the BOX method
                       P1BoxTraits<typename ProblemT::DomainTraits::Scalar,
                                   typename ProblemT::DomainTraits::Grid,
                                   TwoPNITraitsT::numEq>,

                       // The actual problem we would like to solve
                       ProblemT,

                       // The local jacobian operator
                       TwoPNIBoxJacobian<ProblemT,
                                             P1BoxTraits<typename ProblemT::DomainTraits::Scalar,
                                                         typename ProblemT::DomainTraits::Grid,
                                                         TwoPNITraitsT::numEq>,
                                             TwoPNITraitsT
                                             >
                       >
{
    typedef typename ProblemT::DomainTraits::Grid           Grid;
    typedef typename ProblemT::DomainTraits::Scalar         Scalar;
    typedef TwoPNIBoxModel<ProblemT, TwoPNITraitsT> ThisType;

public:
    typedef TwoPNITraitsT                                   TwoPNITraits;
    typedef P1BoxTraits<Scalar, Grid, TwoPNITraits::numEq>  BoxTraits;

private:
    typedef TwoPNIBoxJacobian<ProblemT, BoxTraits, TwoPNITraits>  TwoPNILocalJacobian;
    typedef BoxScheme<ThisType,
                      BoxTraits,
                      ProblemT,
                      TwoPNILocalJacobian>        ParentType;

    typedef typename ProblemT::DomainTraits           DomTraits;
    typedef typename DomTraits::Element               Element;
    typedef typename DomTraits::ElementIterator       ElementIterator;
    typedef typename DomTraits::LocalPosition         LocalPosition;
    typedef typename DomTraits::GlobalPosition        GlobalPosition;

    enum {
        dim          = DomTraits::dim,
        dimWorld     = DomTraits::dimWorld
    };

public:
    typedef NewNewtonMethod<ThisType> NewtonMethod;

    TwoPNIBoxModel(ProblemT &prob)
        : ParentType(prob, twoPLocalJacobian_),
          twoPLocalJacobian_(prob)
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
    };

    /*!
     * \brief Called by the BoxScheme's update method.
     */
    void updateSuccessful()
    {
        ParentType::updateSuccessful();
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
    void calculateMass(Dune::FieldVector<Scalar, 4> &mass)
    {
        twoPLocalJacobian_.calculateMass(this->currentSolution(), mass);
    }

private:
    // calculates the jacobian matrix at a given position
    TwoPNILocalJacobian  twoPLocalJacobian_;
};
}

#endif
