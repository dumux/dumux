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

namespace Dune
{
    ///////////////////////////////////////////////////////////////////////////
    // two-phase two-component traits (central place for names and
    // indices required by the TwoPTwoCNIBoxJacobian and TwoPTwoCNIBoxModel)
    ///////////////////////////////////////////////////////////////////////////
    /*!
     * \brief The 2P-2C specific traits.
     */
    template <class DT>
    class TwoPTwoCNITraits
    {
    public:
        enum {
            numEq         = 3,  //!< Number of primary variables
            numPhases     = 2,  //!< Number of fluid phases
            numComponents = 2   //!< Number of fluid components within a phase
        };
        enum { // Primary variable indices
            pWIndex = 0,           //!< Index for the wetting phase pressure in a field vector
            switchIndex = 1,       //!< Index for the non-wetting phase quantity
            temperatureIndex = 2   //!< Index for the temperature
        };
        enum { // Phase Indices
            wPhase = 0,  //!< Index of the wetting phase
            nPhase = 1   //!< Index of the non-wetting phase
        };
        enum { // Component indices
            wComp = 0,  //!< Index of the wetting component
            nComp = 1   //!< Index of the non-wetting component
        };
        enum { // present phases
            nPhaseOnly = 0, //!< Only the non-wetting phase is present
            wPhaseOnly = 1, //!< Only the wetting phase is present
            bothPhases = 2 //!< Both phases are present
        };

        typedef FieldVector<DT, numPhases>         PhasesVector;

        struct VariableNodeData : public TwoPTwoCTraits<DT>::VariableNodeData
        {
            PhasesVector intenergy;
            PhasesVector enthalpy;
            DT       lambda;
        };
    };


    ///////////////////////////////////////////////////////////////////////////
    // TwoPTwoCNIBoxJacobian (evaluate the local jacobian for the newton method.)
    ///////////////////////////////////////////////////////////////////////////
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
        typedef typename DomTraits::Scalar     DT;
        typedef typename DomTraits::Cell       Cell;

        enum {
            dim  = DomTraits::GridDim,
            WorldDim = DomTraits::WorldDim,

            pWIndex          = TwoPTwoCNITraits::pWIndex,
            switchIndex      = TwoPTwoCNITraits::switchIndex,
            temperatureIndex = TwoPTwoCNITraits::temperatureIndex,

            NumPhases        = TwoPTwoCNITraits::numPhases,
            wPhase      = TwoPTwoCNITraits::wPhase,
            nPhase      = TwoPTwoCNITraits::nPhase,

            wComp      = TwoPTwoCNITraits::wComp,
            nComp      = TwoPTwoCNITraits::nComp,

            wPhaseOnly      = TwoPTwoCNITraits::wPhaseOnly,
            nPhaseOnly      = TwoPTwoCNITraits::nPhaseOnly,
            bothPhases      = TwoPTwoCNITraits::bothPhases,

        };

        typedef BoxTraitsT                              BoxTraits;
        typedef typename BoxTraits::UnknownsVector      UnknownsVector;
        typedef typename TwoPTwoCNITraits::PhasesVector PhasesVector;

        typedef typename ParentType::LocalFunction     LocalFunction;
        typedef typename ParentType::CellCache         CellCache;

        typedef typename ParentType::VariableNodeData  VariableNodeData;
        typedef typename BoxTraits::FVElementGeometry  FVElementGeometry;

        typedef typename DomTraits::WorldCoord         WorldCoord;
        typedef typename DomTraits::LocalCoord         LocalCoord;

        typedef FieldMatrix<DT, dim, dim>  Tensor;

    public:
        TwoPTwoCNIBoxJacobian(ProblemT &problem)
            : ParentType(problem)
            {
            };

        /*!
         * \brief The storage term of heat
         */
        void heatStorage(UnknownsVector &result,
                         int scvId,
                         const LocalFunction &sol,
                         const CellCache &cellCache) const
            {
                const VariableNodeData &nodeDat = cellCache.atSCV[scvId];

                DT satN = nodeDat.satN;
                DT satW = nodeDat.satW;

                // assume porosity defined at nodes
                DT porosity =
                    this->problem_.porosity(this->curCell_(), scvId);

                result[temperatureIndex] =
                    porosity*(nodeDat.density[wPhase] *
                              nodeDat.intenergy[wPhase] *
                              satW
                              +
                              nodeDat.density[nPhase] *
                              nodeDat.intenergy[nPhase] *
                              satN)
                    +
                    this->problem_.soil().heatCap(this->curCellGeom_.cellGlobal,
                                                  this->curCell_(),
                                                  this->curCellGeom_.cellLocal)*
                    this->temperature_(sol[scvId]);
            }

        /*!
         * \brief Update the temperature gradient at a face of an FV
         *        cell.
         */
        void updateTempGrad(WorldCoord &tempGrad,
                            const WorldCoord &feGrad,
                            const LocalFunction &sol,
                            int nodeIdx) const
            {
                WorldCoord tmp = feGrad;
                tmp *= sol[nodeIdx][temperatureIndex];
                tempGrad += tmp;
            }

        /*!
         * \brief Sets the temperature term of the flux vector to the
         *        heat flux due to advection of the fluids.
         */
        void advectiveHeatFlux(UnknownsVector &flux,
                               const PhasesVector &darcyOut,
                               DT alpha, // upwind parameter
                               const VariableNodeData *upW, // up/downstream nodes
                               const VariableNodeData *dnW,
                               const VariableNodeData *upN,
                               const VariableNodeData *dnN) const
            {
                // heat flux in non-wetting phase
                flux[temperatureIndex]  =  darcyOut[nPhase] * (
                    alpha * // upstream node
                    (  upN->density[nPhase] *
                       upN->mobility[nPhase] *
                       upN->enthalpy[nPhase])
                    +
                    (1-alpha) * // downstream node
                    (  dnN->density[nPhase] *
                       dnN->mobility[nPhase] *
                       dnN->enthalpy[nPhase]) );

                // advective heat flux in wetting phase
                flux[temperatureIndex] +=  darcyOut[wPhase] * (
                    alpha * // upstream node
                    (  upW->density[wPhase] *
                       upW->mobility[wPhase] *
                       upW->enthalpy[wPhase] )
                    +
                    (1-alpha) * // downstream node
                    (  dnW->density[wPhase] *
                       dnW->mobility[wPhase] *
                       dnW->enthalpy[wPhase]) );
            }

        /*!
         * \brief Adds the diffusive heat flux to the flux vector over
         *        the face of a sub-control volume.
         */
        void diffusiveHeatFlux(UnknownsVector &flux,
                               int faceIdx,
                               const WorldCoord &tempGrad) const
            {
                const CellCache &cc = this->curSolCache_;
                const FVElementGeometry &fvGeom = this->curCellGeom_;
                int i = fvGeom.subContVolFace[faceIdx].i;
                int j = fvGeom.subContVolFace[faceIdx].j;

                // unit normal vector of the face
                const WorldCoord &normal = fvGeom.subContVolFace[faceIdx].normal;

                // Harmonic mean of the heat conducitivities of the
                // sub control volumes adjacent to the face
                DT lambda =  2./((1. / cc.atSCV[i].lambda) + (1. / cc.atSCV[j].lambda));

                // heat conduction
                flux[temperatureIndex] += (tempGrad*normal) * lambda;;
            }

    public:
        // internal method!
        void updateVarNodeData_(VariableNodeData &d,
                                const UnknownsVector &nodeSol,
                                int phaseState,
                                const Cell &cell,
                                int localIdx,
                                Problem &problem,
                                DT temperature) const
            {
                // update data for the isothermal stuff
                ParentType::updateVarNodeData_(d,
                                               nodeSol,
                                               phaseState,
                                               cell,
                                               localIdx,
                                               problem,
                                               temperature);

                const LocalCoord &local =
                    DomTraits::referenceElement(cell.type()).position(localIdx,
                                                                      dim);
                const WorldCoord &global =
                    cell.geometry()[localIdx];

                // update data for the energy equation
                d.lambda = problem.soil().heatCond(global, cell, local, d.satW);
                d.enthalpy[pWIndex] = problem.wettingPhase().enthalpy(temperature, d.pW);
                d.enthalpy[switchIndex] = problem.nonwettingPhase().enthalpy(temperature, d.pN);

                d.intenergy[pWIndex] = problem.wettingPhase().intEnergy(temperature,
                                                                        d.pW);
                d.intenergy[switchIndex] = problem.nonwettingPhase().intEnergy(temperature,
                                                                               d.pN);
            }


        // internal method!
        static DT temperature_(const UnknownsVector &sol)
            { return sol[temperatureIndex]; }
    };

    ///////////////////////////////////////////////////////////////////////////
    // TwoPTwoCNIBoxModel (The actual numerical model.)
    ///////////////////////////////////////////////////////////////////////////
    /**
     * \brief Isothermal two phase two component model with Pw and
     *        Sn/X as primary unknowns.
     *
     * This implements an isothermal two phase two component model
     * with Pw and Sn/X as primary unknowns.
     */
    template<class ProblemT>
    class TwoPTwoCNIBoxModel
        : public BoxScheme<TwoPTwoCNIBoxModel<ProblemT>, // Implementation of the box scheme

                           // The Traits for the BOX method
                           P1BoxTraits<typename ProblemT::DomainTraits::Scalar,
                                       typename ProblemT::DomainTraits::Grid,
                                       TwoPTwoCNITraits<typename ProblemT::DomainTraits::Scalar>::numEq>,

                           // The actual problem we would like to solve
                           ProblemT,

                           // The local jacobian operator
                           TwoPTwoCNIBoxJacobian<ProblemT,
                                                 P1BoxTraits<typename ProblemT::DomainTraits::Scalar,
                                                             typename ProblemT::DomainTraits::Grid,
                                                             TwoPTwoCNITraits<typename ProblemT::DomainTraits::Scalar>::numEq>,
                                                 TwoPTwoCNITraits<typename ProblemT::DomainTraits::Scalar>
                                                 >
                           >
    {
        typedef typename ProblemT::DomainTraits::Grid   Grid;
        typedef typename ProblemT::DomainTraits::Scalar DT;
        typedef TwoPTwoCNIBoxModel<ProblemT>              ThisType;

    public:
        typedef Dune::TwoPTwoCNITraits<DT>                      TwoPTwoCNITraits;
        typedef P1BoxTraits<DT, Grid, TwoPTwoCNITraits::numEq>  BoxTraits;

    private:
        typedef TwoPTwoCNIBoxJacobian<ProblemT, BoxTraits, TwoPTwoCNITraits>  TwoPTwoCNILocalJacobian;
        typedef BoxScheme<ThisType,
                          BoxTraits,
                          ProblemT,
                          TwoPTwoCNILocalJacobian>        ParentType;

        typedef typename ProblemT::DomainTraits           DomTraits;
        typedef typename DomTraits::Cell                  Cell;
        typedef typename DomTraits::CellIterator          CellIterator;
        typedef typename DomTraits::LocalCoord            LocalCoord;
        typedef typename DomTraits::WorldCoord            WorldCoord;

        enum {
            dim          = DomTraits::GridDim,
            WorldDim         = DomTraits::WorldDim
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
        void addVtkFields(MultiWriter &writer) const
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
