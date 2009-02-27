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
#ifndef DUMUX_NEW_2P2CNI_BOX_MODEL_HH
#define DUMUX_NEW_2P2CNI_BOX_MODEL_HH

#include <dumux/new_models/2p2c/2p2cboxmodel.hh>

#include <dumux/new_models/2p2c/2p2ctraits.hh>

namespace Dune
{

    ///////////////////////////////////////////////////////////////////////////
    // TwoPTwoCNIBoxJacobian (evaluate the local jacobian for the newton method.)
    ///////////////////////////////////////////////////////////////////////////

/** \todo Please doc me! */

    template<class ProblemT,
             class BoxTraitsT,
             class TwoPTwoCNITraitsT>
    class TwoPTwoCNIBoxJacobian : public TwoPTwoCBoxJacobianBase<ProblemT,
                                                                 BoxTraitsT,
                                                                 TwoPTwoCNITraitsT,
                                                                 TwoPTwoCNIBoxJacobian<ProblemT,
                                                                                       BoxTraitsT,
                                                                                       TwoPTwoCNITraitsT>
                                                                 >
    {
        typedef TwoPTwoCNITraitsT TwoPTwoCNITraits;
        typedef TwoPTwoCNIBoxJacobian<ProblemT,
                                      BoxTraitsT,
                                      TwoPTwoCNITraits> ThisType;
        typedef TwoPTwoCBoxJacobianBase<ProblemT,
                                        BoxTraitsT,
                                        TwoPTwoCNITraits,
                                        ThisType> ParentType;


        typedef ProblemT                       Problem;
        typedef typename Problem::DomainTraits DomTraits;
        typedef typename DomTraits::Scalar     Scalar;
        typedef typename DomTraits::Element       Element;

        enum {
            dim  = DomTraits::dim,
            dimWorld = DomTraits::dimWorld,

            pressureIdx    = TwoPTwoCNITraits::pressureIdx,
            switchIdx      = TwoPTwoCNITraits::switchIdx,
            temperatureIdx = TwoPTwoCNITraits::temperatureIdx,

            numPhases   = TwoPTwoCNITraits::numPhases,
            wPhase      = TwoPTwoCNITraits::wPhase,
            nPhase      = TwoPTwoCNITraits::nPhase,

            wComp      = TwoPTwoCNITraits::wComp,
            nComp      = TwoPTwoCNITraits::nComp,

            wPhaseOnly      = TwoPTwoCNITraits::wPhaseOnly,
            nPhaseOnly      = TwoPTwoCNITraits::nPhaseOnly,
            bothPhases      = TwoPTwoCNITraits::bothPhases,

        };

        typedef BoxTraitsT                              BoxTraits;
        typedef typename BoxTraits::SolutionVector      SolutionVector;
        typedef typename TwoPTwoCNITraits::PhasesVector PhasesVector;

        typedef typename ParentType::LocalFunction       LocalFunction;
        typedef typename ParentType::ElementData         ElementData;

        typedef typename ParentType::VariableVertexData  VariableVertexData;
        typedef typename BoxTraits::FVElementGeometry    FVElementGeometry;

        typedef typename DomTraits::GlobalPosition        GlobalPosition;
        typedef typename DomTraits::LocalPosition         LocalPosition;

        typedef FieldMatrix<Scalar, dim, dim>  Tensor;

    public:
        TwoPTwoCNIBoxJacobian(ProblemT &problem)
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
                const VariableVertexData &vertDat = elementCache.vertex[scvId];

                Scalar satN = vertDat.satN;
                Scalar satW = vertDat.satW;

                // assume porosity defined at verts
                Scalar porosity =
                    this->problem_.porosity(this->curElement_(), scvId);

                result[temperatureIdx] =
                    porosity*(vertDat.density[wPhase] *
                              vertDat.intenergy[wPhase] *
                              satW
                              +
                              vertDat.density[nPhase] *
                              vertDat.intenergy[nPhase] *
                              satN)
                    +
                    this->problem_.soil().heatCap(this->curElementGeom_.elementGlobal,
                                                  this->curElement_(),
                                                  this->curElementGeom_.elementLocal)*
                    this->temperature_(sol[scvId]);
            }

        /*!
         * \brief Update the temperature gradient at a face of an FV
         *        element.
         */
        void updateTempGrad(GlobalPosition &tempGrad,
                            const GlobalPosition &feGrad,
                            const LocalFunction &sol,
                            int vertIdx) const
            {
                GlobalPosition tmp = feGrad;
                tmp *= sol[vertIdx][temperatureIdx];
                tempGrad += tmp;
            }

        /*!
         * \brief Sets the temperature term of the flux vector to the
         *        heat flux due to advection of the fluids.
         */
        void advectiveHeatFlux(SolutionVector &flux,
                               const PhasesVector &darcyOut,
                               Scalar alpha, // upwind parameter
                               const VariableVertexData *upW, // up/downstream vertices
                               const VariableVertexData *dnW,
                               const VariableVertexData *upN,
                               const VariableVertexData *dnN) const
            {
                // heat flux in non-wetting phase
                flux[temperatureIdx]  =  darcyOut[nPhase] * (
                    alpha * // upstream vertex
                    (  upN->density[nPhase] *
                       upN->mobility[nPhase] *
                       upN->enthalpy[nPhase])
                    +
                    (1-alpha) * // downstream vertex
                    (  dnN->density[nPhase] *
                       dnN->mobility[nPhase] *
                       dnN->enthalpy[nPhase]) );

                // advective heat flux in wetting phase
                flux[temperatureIdx] +=  darcyOut[wPhase] * (
                    alpha * // upstream vertex
                    (  upW->density[wPhase] *
                       upW->mobility[wPhase] *
                       upW->enthalpy[wPhase] )
                    +
                    (1-alpha) * // downstream vertex
                    (  dnW->density[wPhase] *
                       dnW->mobility[wPhase] *
                       dnW->enthalpy[wPhase]) );
            }

        /*!
         * \brief Adds the diffusive heat flux to the flux vector over
         *        the face of a sub-control volume.
         */
        void diffusiveHeatFlux(SolutionVector &flux,
                               int faceIdx,
                               const GlobalPosition &tempGrad) const
            {
                const ElementData &eDat = this->curElemDat_;
                const FVElementGeometry &fvGeom = this->curElementGeom_;
                int i = fvGeom.subContVolFace[faceIdx].i;
                int j = fvGeom.subContVolFace[faceIdx].j;

                // unit normal vector of the face
                const GlobalPosition &normal = fvGeom.subContVolFace[faceIdx].normal;

                // Harmonic mean of the heat conducitivities of the
                // sub control volumes adjacent to the face
                Scalar lambda =  ParentType::harmonicMean_(eDat.vertex[i].lambda,
                                                           eDat.vertex[j].lambda);

                // heat conduction
                flux[temperatureIdx] += (tempGrad*normal) * lambda;;
            }

    public:
        // internal method!
        void updateVarVertexData_(VariableVertexData &vertDat,
                                  const SolutionVector &vertSol,
                                  int phaseState,
                                  const Element &element,
                                  int localIdx,
                                  Problem &problem,
                                  Scalar temperature) const
            {
                // update data for the isothermal stuff
                ParentType::updateVarVertexData_(vertDat,
                                               vertSol,
                                               phaseState,
                                               element,
                                               localIdx,
                                               problem,
                                               temperature);

                const LocalPosition &local =
                    DomTraits::referenceElement(element.type()).position(localIdx,
                                                                         dim);
                const GlobalPosition &global =
                    element.geometry().corner(localIdx);

                // update data for the energy equation
                vertDat.lambda = problem.soil().heatCond(global, element, local, vertDat.satW);
                vertDat.enthalpy[pressureIdx] = problem.wettingPhase().enthalpy(temperature,
                                                                                vertDat.pW,
                                                                                vertDat.massfrac[nComp][wPhase]);
                vertDat.enthalpy[switchIdx] = problem.nonwettingPhase().enthalpy(temperature,
                                                                              vertDat.pN,
                                                                              vertDat.massfrac[wComp][nPhase]);
                vertDat.intenergy[pressureIdx] = problem.wettingPhase().intEnergy(temperature,
                                                                          vertDat.pW,
                                                                          vertDat.massfrac[nComp][wPhase]);
                vertDat.intenergy[switchIdx] = problem.nonwettingPhase().intEnergy(temperature,
                                                                          vertDat.pN,
                                                                          vertDat.massfrac[wComp][nPhase]);
            }


        // internal method!
        static Scalar temperature_(const SolutionVector &sol)
            { return sol[temperatureIdx]; }
    };

    ///////////////////////////////////////////////////////////////////////////
    // TwoPTwoCNIBoxModel (The actual numerical model.)
    ///////////////////////////////////////////////////////////////////////////
    /**
     * \brief Non-isothermal two phase two component model with Pw and
     *        Sn/X as primary unknowns.
     *
     * This implements a non-isothermal two-phase two-component model
     * with Pw and Sn/X as primary unknowns.
     */
    template<class ProblemT,
             class TwoPTwoCNITraitsT = TwoPTwoCNITraits<typename ProblemT::DomainTraits::Scalar,
                                                        TwoPTwoCPwSnTraits<typename ProblemT::DomainTraits::Scalar> > >
    class TwoPTwoCNIBoxModel
        : public BoxScheme<TwoPTwoCNIBoxModel<ProblemT,TwoPTwoCNITraitsT>, // Implementation of the box scheme

                           // The Traits for the BOX method
                           P1BoxTraits<typename ProblemT::DomainTraits::Scalar,
                                       typename ProblemT::DomainTraits::Grid,
                                       TwoPTwoCNITraitsT::numEq>,

                           // The actual problem we would like to solve
                           ProblemT,

                           // The local jacobian operator
                           TwoPTwoCNIBoxJacobian<ProblemT,
                                                 P1BoxTraits<typename ProblemT::DomainTraits::Scalar,
                                                             typename ProblemT::DomainTraits::Grid,
                                                             TwoPTwoCNITraitsT::numEq>,
                                                 TwoPTwoCNITraitsT
                                                 >
                           >
    {
        typedef typename ProblemT::DomainTraits::Grid           Grid;
        typedef typename ProblemT::DomainTraits::Scalar         Scalar;
        typedef TwoPTwoCNIBoxModel<ProblemT, TwoPTwoCNITraitsT> ThisType;

    public:
        typedef TwoPTwoCNITraitsT                                   TwoPTwoCNITraits;
        typedef P1BoxTraits<Scalar, Grid, TwoPTwoCNITraits::numEq>  BoxTraits;

    private:
        typedef TwoPTwoCNIBoxJacobian<ProblemT, BoxTraits, TwoPTwoCNITraits>  TwoPTwoCNILocalJacobian;
        typedef BoxScheme<ThisType,
                          BoxTraits,
                          ProblemT,
                          TwoPTwoCNILocalJacobian>        ParentType;

        typedef typename ProblemT::DomainTraits           DomTraits;
        typedef typename DomTraits::Element                  Element;
        typedef typename DomTraits::ElementIterator          ElementIterator;
        typedef typename DomTraits::LocalPosition            LocalPosition;
        typedef typename DomTraits::GlobalPosition            GlobalPosition;

        enum {
            dim          = DomTraits::dim,
            dimWorld     = DomTraits::dimWorld
        };

    public:
        typedef NewNewtonMethod<ThisType> NewtonMethod;

        TwoPTwoCNIBoxModel(ProblemT &prob)
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
         * \brief Returns true if there was a primary variable switch
         *        after the last time step.
         */
        bool switched() const
            { return twoPTwoCLocalJacobian_.switched(); }


    private:
        // calculates the jacobian matrix at a given position
        TwoPTwoCNILocalJacobian  twoPTwoCLocalJacobian_;
    };
}

#endif
