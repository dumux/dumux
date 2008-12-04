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
#ifndef DUMUX_NEW_2P2C_BOX_MODEL_HH
#define DUMUX_NEW_2P2C_BOX_MODEL_HH

#include <dumux/new_models/boxscheme/boxscheme.hh>
#include <dumux/new_models/boxscheme/p1boxtraits.hh>

#include <dumux/auxiliary/apis.hh>

#include <vector>

namespace Dune
{
    ///////////////////////////////////////////////////////////////////////////
    // two-phase two-component traits (central place for names and
    // indices required by the TwoPTwoCBoxJacobian and TwoPTwoCBoxModel)
    ///////////////////////////////////////////////////////////////////////////
    /*!
     * \brief The 2P-2C specific traits.
     */
    template <class DT>
    class TwoPTwoCTraits
    {
    public:
        enum {
            numEq		  = 2,  //!< Number of primary variables / equations
            numPhases     = 2,  //!< Number of fluid phases
            numComponents = 2   //!< Number of fluid components within a phase
        };
        enum { // Primary variable indices
            pWIndex = 0,           //!< Index for the wetting phase pressure in a field vector
            switchIndex = 1,       //!< Index for the non-wetting phase quantity
        };
        enum { // Phase Indices
            wPhase = 0,  //!< Index of the wetting phase
            nPhase = 1   //!< Index of the non-wetting phase
        };
        enum { // Component indices
            wComp = 0,  //!< Index of the major component of the wetting phase
            nComp = 1   //!< Index of the major component of the non-wetting phase
        };
        enum { // present phases
            nPhaseOnly = 0, //!< Only the non-wetting phase is present
            wPhaseOnly = 1, //!< Only the wetting phase is present
            bothPhases = 2  //!< Both phases are present
        };

        typedef FieldVector<DT, numPhases>         PhasesVector;

        /*!
         * \brief Data which is attached to each node of the and can
         *        be shared between multiple calculations and should
         *        thus be cached in order to increase efficency.
         */
        struct VariableNodeData
        {
            DT satN;
            DT satW;
            DT pW;
            DT pC;
            DT pN;
            PhasesVector mobility;  //Vector with the number of phases
            PhasesVector density;
            FieldMatrix<DT, numComponents, numPhases> massfrac;
            int phaseState;
        };

    };


    ///////////////////////////////////////////////////////////////////////////
    // TwoPTwoCBoxJacobian (evaluate the local jacobian for the newton method.)
    ///////////////////////////////////////////////////////////////////////////
    /*!
     * \brief 2P-2C specific details needed to approximately calculate
     *        the local jacobian in the BOX scheme.
     *
     * This class is used to fill the gaps in BoxJacobian for the 2P-2C twophase flow.
     */
    template<class ProblemT,
             class BoxTraitsT,
             class TwoPTwoCTraitsT,
             class Implementation>
    class TwoPTwoCBoxJacobianBase : public BoxJacobian<ProblemT,
                                                       BoxTraitsT,
                                                       Implementation >
    {
    protected:
        typedef TwoPTwoCBoxJacobianBase<ProblemT,
                                        BoxTraitsT,
                                        TwoPTwoCTraitsT,
                                        Implementation>              ThisType;
        typedef BoxJacobian<ProblemT, BoxTraitsT, Implementation>    ParentType;

        typedef ProblemT                                Problem;
        typedef typename Problem::DomainTraits          DomTraits;
        typedef BoxTraitsT                              BoxTraits;
        typedef TwoPTwoCTraitsT                         TwoPTwoCTraits;

        enum {
            dim              = DomTraits::GridDim,
            dimWorld         = DomTraits::WorldDim,

            PrimaryVariables = BoxTraits::PrimaryVariables,
            numPhases        = TwoPTwoCTraits::numPhases,
            numComponents    = TwoPTwoCTraits::numComponents,

            pWIndex          = TwoPTwoCTraits::pWIndex,
            switchIndex      = TwoPTwoCTraits::switchIndex,

            wPhase           = TwoPTwoCTraits::wPhase,
            nPhase           = TwoPTwoCTraits::nPhase,

            wComp            = TwoPTwoCTraits::wComp,
            nComp            = TwoPTwoCTraits::nComp,

            wPhaseOnly       = TwoPTwoCTraits::wPhaseOnly,
            nPhaseOnly       = TwoPTwoCTraits::nPhaseOnly,
            bothPhases       = TwoPTwoCTraits::bothPhases
        };


        typedef typename DomTraits::Scalar                DT;
        typedef typename DomTraits::CoordScalar           CoordScalar;
        typedef typename DomTraits::Grid                  Grid;
        typedef typename DomTraits::Cell                  Cell;
        typedef typename DomTraits::CellIterator          CellIterator;
        typedef typename Cell::EntityPointer              CellPointer;
        typedef typename DomTraits::LocalCoord            LocalCoord;
        typedef typename DomTraits::WorldCoord            WorldCoord;
        typedef typename DomTraits::NodeIterator          NodeIterator;

        typedef typename BoxTraits::UnknownsVector      UnknownsVector;
        typedef typename BoxTraits::FVElementGeometry   FVElementGeometry;
        typedef typename BoxTraits::SpatialFunction     SpatialFunction;
        typedef typename BoxTraits::LocalFunction       LocalFunction;

        typedef typename TwoPTwoCTraits::PhasesVector      PhasesVector;
        typedef typename TwoPTwoCTraits::VariableNodeData  VariableNodeData;
        typedef FieldMatrix<DT, dim, dim>  Tensor;

        /*!
         * \brief Cached data for the each node of the cell.
         */
        struct CellCache
        {
            VariableNodeData atSCV[BoxTraits::ShapeFunctionSetContainer::maxsize];
        };


        /*!
         * \brief Data which is attached to each node and is not only
         *        locally.
         */
        struct StaticNodeData {
            int phaseState;
            int oldPhaseState;
        };

        void updateVarNodeData_(VariableNodeData &d,
                                const UnknownsVector &nodeSol,
                                int phaseState,
                                const Cell &cell,
                                int localIdx,
                                Problem &problem,
                                DT temperature) const
                {
                    const WorldCoord &global = cell.geometry()[localIdx];
                    const LocalCoord &local =
                        DomTraits::referenceElement(cell.type()).position(localIdx,
                                                                          dim);
                    d.pW = nodeSol[pWIndex];
                    if (phaseState == bothPhases) d.satN = nodeSol[switchIndex];
                    else if (phaseState == wPhaseOnly) d.satN = 0.0;
                    else if (phaseState == nPhaseOnly) d.satN = 1.0;
                    else DUNE_THROW(Dune::InvalidStateException, "Phase state " << phaseState << " is invalid.");

                    d.satW = 1.0 - d.satN;
                    d.pC = problem.materialLaw().pC(d.satW,
                                                    global,
                                                    cell,
                                                    local);
                    d.pN = d.pW + d.pC;

                    // Solubilities of components in phases
                    if (phaseState == bothPhases) {
                        d.massfrac[nComp][wPhase] = problem.multicomp().xAW(d.pN, temperature);
                        d.massfrac[wComp][nPhase] = problem.multicomp().xWN(d.pN, temperature);
                    }
                    else if (phaseState == wPhaseOnly) {
                        d.massfrac[wComp][nPhase] = 0.0;
                        d.massfrac[nComp][wPhase] =  nodeSol[switchIndex];
                    }
                    else if (phaseState == nPhaseOnly){
                        d.massfrac[wComp][nPhase] = nodeSol[switchIndex];
                        d.massfrac[nComp][wPhase] = 0.0;
                    }
                    else DUNE_THROW(Dune::InvalidStateException, "Phase state " << phaseState << " is invalid.");

                    d.massfrac[wComp][wPhase] = 1.0 - d.massfrac[nComp][wPhase];
                    d.massfrac[nComp][nPhase] = 1.0 - d.massfrac[wComp][nPhase];
                    d.phaseState = phaseState;

                    // Density of Water is set constant here!
                    d.density[wPhase] = problem.wettingPhase().density(temperature,
                                                                            d.pW,
                                                                            d.massfrac[nComp][wPhase]);
                    d.density[nPhase] = problem.nonwettingPhase().density(temperature,
                                                                               d.pN,
                                                                               d.massfrac[wComp][nPhase]);

                    // Mobilities & densities
                    d.mobility[wPhase] = problem.materialLaw().mobW(d.satW, global, cell, local, temperature, d.pW);
#if 1
                    d.mobility[nPhase] = problem.materialLaw().mobN(d.satN, global, cell, local, temperature, d.pN);
#else
#warning "CO2 specific"
                    DT viscosityCO2 = problem.nonwettingPhase().viscosityCO2(temperature,
                                                                                 d.pN,
                                                                                 d.density[nPhase]);
                    DT krCO2 = problem.materialLaw().krn(d.satN,
                                                             this->curCellGeom_.subContVol[localIdx].global,
                                                             cell,
                                                             this->curCellGeom_.subContVol[localIdx].local);

                    d.mobility[nPhase] = krCO2 / viscosityCO2;
#endif
                }

    public:
        TwoPTwoCBoxJacobianBase(ProblemT &problem)
            : ParentType(problem),
              staticNodeDat_(problem.numNodes())
            {
                switchFlag_ = false;
            };

        /*!
         * \brief Set the current grid cell.
         */
        void setCurrentCell(const Cell &cell)
            {
                ParentType::setCurrentCell_(cell);
            };

        /*!
         * \brief Set the parameters for the calls to the remaining
         *        members.
         */
        void setParams(const Cell &cell, LocalFunction &curSol, LocalFunction &prevSol)
            {
                setCurrentCell(cell);

                // TODO: scheme which allows not to copy curSol and
                // prevSol all the time
                curSol_ = curSol;
                updateCellCache_(curSolCache_, curSol_, false);
                curSolDeflected_ = false;

                prevSol_ = prevSol;
                updateCellCache_(prevSolCache_, prevSol_, true);
            };

        /*!
         * \brief Vary a single component of a single node of the
         *        local solution for the current cell.
         *
         * This method is a optimization, since if varying a single
         * component at a degree of freedom not the whole cell cache
         * needs to be recalculated. (Updating the cell cache is very
         * expensive since material laws need to be evaluated.)
         */
        void deflectCurSolution(int node, int component, DT value)
            {
                // make sure that the orignal state can be restored
                if (!curSolDeflected_) {
                    curSolDeflected_ = true;

                    curSolOrigValue_ = curSol_[node][component];
                    curSolOrigVarData_ = curSolCache_.atSCV[node];
                }

                int globalIdx = ParentType::problem_.nodeIndex(ParentType::curCell_(),
                                                               node);

                curSol_[node][component] = value;
                asImp_()->updateVarNodeData_(curSolCache_.atSCV[node],
                                             curSol_[node],
                                             staticNodeDat_[globalIdx].phaseState,
                                             this->curCell_(),
                                             node,
                                             this->problem_,
                                             Implementation::temperature_(curSol_[node]));
            }

        /*!
         * \brief Restore the local jacobian to the state before
         *        deflectCurSolution() was called.
         *
         * This only works if deflectSolution was only called with
         * (node, component) as arguments.
         */
        void restoreCurSolution(int node, int component)
            {
                curSolDeflected_ = false;
                curSol_[node][component] = curSolOrigValue_;
                curSolCache_.atSCV[node] = curSolOrigVarData_;
            };

        /*!
         * \brief Evaluate the rate of change of all conservation
         *        quantites (e.g. phase mass) within a sub control
         *        volume of a finite volume cell in the 2P-2C
         *        model.
         *
         * This function should not include the source and sink terms.
         */
        void localRate(UnknownsVector &result, int scvId, bool usePrevSol) const
            {
                result = DT(0);

                const LocalFunction &sol   = usePrevSol ? this->prevSol_ : this->curSol_;
                const CellCache &cellCache = usePrevSol ? prevSolCache_  : curSolCache_;

                DT satN = cellCache.atSCV[scvId].satN;
                DT satW = cellCache.atSCV[scvId].satW;

                // assume porosity defined at nodes
                DT porosity =
                    this->problem_.porosity(this->curCell_(), scvId);

                // storage of component water
                result[pWIndex] =
                    porosity*(cellCache.atSCV[scvId].density[wPhase]*
                              satW*
                              cellCache.atSCV[scvId].massfrac[wComp][wPhase]
                              + cellCache.atSCV[scvId].density[nPhase]*
                              satN*
                              cellCache.atSCV[scvId].massfrac[wComp][nPhase]);

                // storage of component air
                result[switchIndex] =
                    porosity*(cellCache.atSCV[scvId].density[nPhase]*
                              satN*
                              cellCache.atSCV[scvId].massfrac[nComp][nPhase]
                              + cellCache.atSCV[scvId].density[wPhase]*
                              satW*
                              cellCache.atSCV[scvId].massfrac[nComp][wPhase]);

                // storage of energy
                asImp_()->heatStorage(result, scvId, sol, cellCache);
            }

        /*!
         * \brief Evaluates the mass flux over a face of a subcontrol
         *        volume.
         */
        void fluxRate(UnknownsVector &flux, int faceId) const
            {
                // set flux vector to zero
                int i = this->curCellGeom_.subContVolFace[faceId].i;
                int j = this->curCellGeom_.subContVolFace[faceId].j;

                // normal vector, value of the area of the scvf
                const WorldCoord &normal(this->curCellGeom_.subContVolFace[faceId].normal);

                // get global coordinates of nodes i,j
                const WorldCoord &global_i = this->curCellGeom_.subContVol[i].global;
                const WorldCoord &global_j = this->curCellGeom_.subContVol[j].global;

                // get local coordinates of nodes i,j
                const LocalCoord &local_i = this->curCellGeom_.subContVol[i].local;
                const LocalCoord &local_j = this->curCellGeom_.subContVol[j].local;

                WorldCoord pGrad[numPhases];
                WorldCoord xGrad[numPhases];
                WorldCoord tempGrad(0.0);
                for (int k = 0; k < numPhases; ++k) {
                    pGrad[k] = DT(0);
                    xGrad[k] = DT(0);
                }

                WorldCoord tmp(0.0);
                PhasesVector pressure(0.0), massfrac(0.0);

                // calculate FE gradient (grad p for each phase)
                for (int k = 0; k < this->curCellGeom_.nNodes; k++) // loop over adjacent nodes
                {
                    // FEGradient at node k
                    const LocalCoord &feGrad = this->curCellGeom_.subContVolFace[faceId].grad[k];

                    pressure[wPhase] = this->curSolCache_.atSCV[k].pW;
                    pressure[nPhase] = this->curSolCache_.atSCV[k].pN;

                    // compute sum of pressure gradients for each phase
                    for (int phase = 0; phase < numPhases; phase++)
                    {
                        tmp = feGrad;
                        tmp *= pressure[phase];
                        pGrad[phase] += tmp;
                    }

                    // for diffusion of air in wetting phase
                    tmp = feGrad;
                    tmp *= this->curSolCache_.atSCV[k].massfrac[nComp][wPhase];
                    xGrad[wPhase] += tmp;

                    // for diffusion of water in nonwetting phase
                    tmp = feGrad;
                    tmp *= this->curSolCache_.atSCV[k].massfrac[wComp][nPhase];
                    xGrad[nPhase] += tmp;

                    // temperature gradient
                    asImp_()->updateTempGrad(tempGrad, feGrad, this->curSol_, k);
                }

                // correct the pressure gradients by the hydrostatic
                // pressure due to gravity
                for (int phase=0; phase < numPhases; phase++)
                {
                    tmp = this->problem_.gravity();
                    tmp *= this->curSolCache_.atSCV[i].density[phase];
                    pGrad[phase] -= tmp;
                }


                // calculate the permeability tensor
                Tensor K         = this->problem_.soil().K(global_i, ParentType::curCell_(), local_i);
                const Tensor &Kj = this->problem_.soil().K(global_j, ParentType::curCell_(), local_j);
                harmonicMeanK_(K, Kj);

                // magnitute of darcy velocity of each phase projected
                // on the normal of the sub-control volume's face
                PhasesVector vDarcyOut;
                // temporary vector for the Darcy velocity
                WorldCoord vDarcy;
                for (int phase=0; phase < numPhases; phase++)
                {
                    K.mv(pGrad[phase], vDarcy);  // vDarcy = K * grad p
                    vDarcyOut[phase] = vDarcy*normal;
                }

                // find upsteam and downstream nodes
                const VariableNodeData *upW = &(this->curSolCache_.atSCV[i]);
                const VariableNodeData *dnW = &(this->curSolCache_.atSCV[j]);
                const VariableNodeData *upN = &(this->curSolCache_.atSCV[i]);
                const VariableNodeData *dnN = &(this->curSolCache_.atSCV[j]);
                if (vDarcyOut[wPhase] > 0) {
                    std::swap(upW, dnW);
                };
                if (vDarcyOut[nPhase] > 0)  {
                    std::swap(upN, dnN);
                };

                // Upwind parameter
                DT alpha = 1.0; // -> use only the upstream node

                ////////
                // advective flux of the wetting component
                ////////

                // flux in the wetting phase
                flux[pWIndex] =  vDarcyOut[wPhase] * (
                    alpha* // upstream nodes
                    (  upW->density[wPhase] *
                       upW->mobility[wPhase] *
                       upW->massfrac[wComp][wPhase])
                    +
                    (1-alpha)* // downstream nodes
                    (  dnW->density[wPhase] *
                       dnW->mobility[wPhase] *
                       dnW->massfrac[wComp][wPhase]));

                // flux in the non-wetting phase
                flux[pWIndex] += vDarcyOut[nPhase] * (
                    alpha* // upstream node
                    (  upN->density[nPhase] *
                       upN->mobility[nPhase] *
                       upN->massfrac[wComp][nPhase])
                    +
                    (1-alpha)* // downstream node
                    (  dnN->density[nPhase] *
                       dnN->mobility[nPhase] *
                       dnN->massfrac[wComp][nPhase]) );

                ////////
                // advective flux of the non-wetting component
                ////////

                // flux in the wetting phase
                flux[switchIndex]   = vDarcyOut[nPhase] * (
                    alpha * // upstream nodes
                    (  upN->density[nPhase] *
                       upN->mobility[nPhase] *
                       upN->massfrac[nComp][nPhase])
                    +
                    (1-alpha) * // downstream node
                    (  dnN->density[nPhase] *
                       dnN->mobility[nPhase] *
                       dnN->massfrac[nComp][nPhase]) );

                // flux in the non-wetting phase
                flux[switchIndex]  += vDarcyOut[wPhase] * (
                    alpha * // upstream node
                    (  upW->density[wPhase] *
                       upW->mobility[wPhase] *
                       upW->massfrac[nComp][wPhase])
                    +
                    (1-alpha) * // downstream node
                    (  dnW->density[wPhase] *
                       dnW->mobility[wPhase] *
                       dnW->massfrac[nComp][wPhase]) );

                ////////
                // advective flux of energy
                ////////
                asImp_()->advectiveHeatFlux(flux, vDarcyOut, alpha, upW, dnW, upN, dnN);

                asImp_()->diffusiveHeatFlux(flux, faceId, tempGrad);

                return;

                // DIFFUSION
                UnknownsVector normDiffGrad;

                // get local to global id map
                int state_i = this->curSolCache_.atSCV[i].phaseState;
                int state_j = this->curSolCache_.atSCV[j].phaseState;

                DT diffusionWW(0.0), diffusionWN(0.0); // diffusion of water
                DT diffusionAW(0.0), diffusionAN(0.0); // diffusion of air
                UnknownsVector avgDensity, avgDpm;

                // Diffusion coefficent
                // TODO: needs to be continuously dependend on the phase saturations
                avgDpm[wPhase]=2e-9;
                avgDpm[nPhase]=2.25e-5;
                if (state_i == nPhaseOnly || state_j == nPhaseOnly)
                {
                    // only the nonwetting phase is present in at
                    // least one cell -> no diffusion within the
                    // wetting phase
                    avgDpm[wPhase] = 0;
                }
                if (state_i == wPhaseOnly || state_j == wPhaseOnly)
                {
                    // only the wetting phase is present in at least
                    // one cell -> no diffusion within the non wetting
                    // phase
                    avgDpm[nPhase] = 0;
                }

                // length of the diffusion gradient
                normDiffGrad[wPhase] = xGrad[wPhase]*normal;
                normDiffGrad[nPhase] = xGrad[nPhase]*normal;

                // calculate the arithmetic mean of densities
                avgDensity[wPhase] = 0.5*(this->curSolCache_.atSCV[i].density[wPhase] + this->curSolCache_.atSCV[j].density[wPhase]);
                avgDensity[nPhase] = 0.5*(this->curSolCache_.atSCV[i].density[nPhase] + this->curSolCache_.atSCV[j].density[nPhase]);

                diffusionAW = avgDpm[wPhase] * avgDensity[wPhase] * normDiffGrad[wPhase];
                diffusionWW = - diffusionAW;
                diffusionWN = avgDpm[nPhase] * avgDensity[nPhase] * normDiffGrad[nPhase];
                diffusionAN = - diffusionWN;

                // add diffusion of water to water flux
                flux[wComp] += (diffusionWW + diffusionWN);

                // add diffusion of air to air flux
                flux[nComp] += (diffusionAN + diffusionAW);

/*
                // set flux vector to zero
                int i = this->curCellGeom_.subContVolFace[faceId].i;
                int j = this->curCellGeom_.subContVolFace[faceId].j;

                // normal vector, value of the area of the scvf
                const WorldCoord &normal(this->curCellGeom_.subContVolFace[faceId].normal);

                // get global coordinates of nodes i,j
                const WorldCoord &global_i = this->curCellGeom_.subContVol[i].global;
                const WorldCoord &global_j = this->curCellGeom_.subContVol[j].global;

                const LocalCoord &local_i = this->curCellGeom_.subContVol[i].local;
                const LocalCoord &local_j = this->curCellGeom_.subContVol[j].local;

                WorldCoord pGrad[numPhases];
                WorldCoord xGrad[numPhases];
                for (int k = 0; k < numPhases; ++k) {
                    pGrad[k] = Scalar(0);
                    xGrad[k] = Scalar(0);
                }

                WorldCoord tmp(0.0);
                PhasesVector pressure(0.0), massfrac(0.0);

                // calculate harmonic mean of permeabilities of nodes i and j
                Tensor K         = this->problem_.soil().K(global_i, ParentType::curCell_(), local_i);
                const Tensor &Kj = this->problem_.soil().K(global_j, ParentType::curCell_(), local_j);
                harmonicMeanK_(K, Kj);

                // calculate FE gradient (grad p for each phase)
                for (int k = 0; k < this->curCellGeom_.nNodes; k++) // loop over adjacent nodes
                {
                    // FEGradient at node k
                    const LocalCoord &feGrad = this->curCellGeom_.subContVolFace[faceId].grad[k];

                    pressure[wPhase] = curSolCache_.atSCV[k].pW;
                    pressure[nPhase] = curSolCache_.atSCV[k].pN;

                    // compute sum of pressure gradients for each phase
                    for (int phase = 0; phase < numPhases; phase++)
                    {
                        tmp = feGrad;

                        tmp *= pressure[phase];

                        pGrad[phase] += tmp;
                    }

                    // for diffusion of air in wetting phase
                    tmp = feGrad;
                    tmp *= curSolCache_.atSCV[k].massfrac[nComp][wPhase];
                    xGrad[wPhase] += tmp;

                    // for diffusion of water in nonwetting phase
                    tmp = feGrad;
                    tmp *= curSolCache_.atSCV[k].massfrac[wComp][nPhase];
                    xGrad[nPhase] += tmp;
                }

                // deduce gravity*density of each phase
                WorldCoord contribComp[numPhases];
                for (int phase=0; phase < numPhases; phase++)
                {
                    contribComp[phase] = this->problem_.gravity();
                    contribComp[phase] *= curSolCache_.atSCV[i].density[phase];
                    pGrad[phase] -= contribComp[phase]; // grad p - rho*g
                }

                // calculate the advective flux using upwind: K*n(grad p -rho*g)
                PhasesVector outward;  // Darcy velocity of each phase
                WorldCoord v_tilde(0);
                for (int phase=0; phase < numPhases; phase++)
                {
                    K.mv(pGrad[phase], v_tilde);  // v_tilde=K*gradP
                    outward[phase] = v_tilde*normal;
                }

                // evaluate upwind nodes
                int up_w, dn_w, up_n, dn_n;
                if (outward[wPhase] <= 0) {
                    up_w = i; dn_w = j;
                }
                else {
                    up_w = j; dn_w = i;
                };

                if (outward[nPhase] <= 0) {
                    up_n = i; dn_n = j;
                }
                else {
                    up_n = j; dn_n = i;
                };

                DT alpha = 1.0;  // Upwind parameter

                // water conservation
                flux[wComp] =   (alpha* curSolCache_.atSCV[up_w].density[wPhase]*curSolCache_.atSCV[up_w].mobility[wPhase]
                                      * curSolCache_.atSCV[up_w].massfrac[wComp][wPhase]
                                      + (1-alpha)* curSolCache_.atSCV[dn_w].density[wPhase]*curSolCache_.atSCV[dn_w].mobility[wPhase]
                                      * curSolCache_.atSCV[dn_w].massfrac[wComp][wPhase])
                    * outward[wPhase];
                flux[wComp] +=  (alpha* curSolCache_.atSCV[up_n].density[nPhase]*curSolCache_.atSCV[up_n].mobility[nPhase]
                                      * curSolCache_.atSCV[up_n].massfrac[wComp][nPhase]
                                      + (1-alpha)* curSolCache_.atSCV[dn_n].density[nPhase]*curSolCache_.atSCV[dn_n].mobility[nPhase]
                                      * curSolCache_.atSCV[dn_n].massfrac[wComp][nPhase])
                    * outward[nPhase];


                // air conservation
                flux[nComp]   = (alpha* curSolCache_.atSCV[up_n].density[nPhase]*curSolCache_.atSCV[up_n].mobility[nPhase]
                                      * curSolCache_.atSCV[up_n].massfrac[nComp][nPhase]
                                      + (1-alpha)* curSolCache_.atSCV[dn_n].density[nPhase]*curSolCache_.atSCV[dn_n].mobility[nPhase]
                                      * curSolCache_.atSCV[dn_n].massfrac[nComp][nPhase])
                    * outward[nPhase];

                flux[nComp]  +=   (alpha* curSolCache_.atSCV[up_w].density[wPhase]*curSolCache_.atSCV[up_w].mobility[wPhase]
                                        * curSolCache_.atSCV[up_w].massfrac[nComp][wPhase]
                                        + (1-alpha)* curSolCache_.atSCV[dn_w].density[wPhase]*curSolCache_.atSCV[dn_w].mobility[wPhase]
                                        * curSolCache_.atSCV[dn_w].massfrac[nComp][wPhase])
                    * outward[wPhase];


                return;

                // DIFFUSION
                UnknownsVector normDiffGrad;

                // get local to global id map
                int state_i = curSolCache_.atSCV[i].phaseState;
                int state_j = curSolCache_.atSCV[j].phaseState;

                DT diffusionWW(0.0), diffusionWN(0.0); // diffusion of water
                DT diffusionAW(0.0), diffusionAN(0.0); // diffusion of air
                UnknownsVector avgDensity, avgDpm;

                // Diffusion coefficent
                // TODO: needs to be continuously dependend on the phase saturations
                avgDpm[wPhase]=2e-9;
                avgDpm[nPhase]=2.25e-5;
                if (state_i == nPhaseOnly || state_j == nPhaseOnly)
                {
                    // only the nonwetting phase is present in at
                    // least one cell -> no diffusion within the
                    // wetting phase
                    avgDpm[wPhase] = 0;
                }
                if (state_i == wPhaseOnly || state_j == wPhaseOnly)
                {
                    // only the wetting phase is present in at least
                    // one cell -> no diffusion within the non wetting
                    // phase
                    avgDpm[nPhase] = 0;
                }

                // length of the diffusion gradient
                normDiffGrad[wPhase] = xGrad[wPhase]*normal;
                normDiffGrad[nPhase] = xGrad[nPhase]*normal;

                // calculate the arithmetic mean of densities
                avgDensity[wPhase] = 0.5*(curSolCache_.atSCV[i].density[wPhase] + curSolCache_.atSCV[j].density[wPhase]);
                avgDensity[nPhase] = 0.5*(curSolCache_.atSCV[i].density[nPhase] + curSolCache_.atSCV[j].density[nPhase]);

                diffusionAW = avgDpm[wPhase] * avgDensity[wPhase] * normDiffGrad[wPhase];
                diffusionWW = - diffusionAW;
                diffusionWN = avgDpm[nPhase] * avgDensity[nPhase] * normDiffGrad[nPhase];
                diffusionAN = - diffusionWN;

                // add diffusion of water to water flux
                flux[wComp] += (diffusionWW + diffusionWN);

                // add diffusion of air to air flux
                flux[nComp] += (diffusionAN + diffusionAW);
*/
            }


        /*!
         * \brief Initialize the static data with the initial solution.
         *
         * Called by TwoPTwoCBoxModel::initial()
         */
        void initStaticData()
            {
                setSwitched(false);

                NodeIterator it = this->problem_.nodeBegin();
                NodeIterator endit = this->problem_.nodeEnd();
                for (; it != endit; ++it)
                {
                    int globalIdx = this->problem_.nodeIndex(*it);
                    const WorldCoord &globalPos = it->geometry()[0];

                    // initialize phase state
                    staticNodeDat_[globalIdx].phaseState =
                        this->problem_.initialPhaseState(*it, globalIdx, globalPos);
                    staticNodeDat_[globalIdx].oldPhaseState =
                        staticNodeDat_[globalIdx].phaseState;
                }
            }

        /*!
         * \brief Update the static data of a single node and do a
         *        variable switch if necessary.
         */
        void updateStaticData(SpatialFunction &curSol, SpatialFunction &oldSol)
            {
                bool wasSwitched = false;

                NodeIterator it = this->problem_.nodeBegin();
                for (; it != this->problem_.nodeEnd(); ++it)
                {
                    int globalIdx = this->problem_.nodeIndex(*it);
                    const WorldCoord &global = it->geometry()[0];

                    wasSwitched = primaryVarSwitch_(curSol,
                                                    globalIdx,
                                                    global)
                        || wasSwitched;
                }

                setSwitched(wasSwitched);
            }

        /*!
         * \brief Set the old phase of all nodes state to the current one.
         */
        void updateOldPhaseState()
            {
                int nNodes = this->problem_.numNodes();
                for (int i = 0; i < nNodes; ++i)
                    staticNodeDat_[i].oldPhaseState = staticNodeDat_[i].phaseState;
            }

        /*!
         * \brief Reset the current phase state of all nodes to the old one after an update failed
         */
        void resetPhaseState()
            {
                int nNodes = this->problem_.numNodes();
                for (int i = 0; i < nNodes; ++i)
                    staticNodeDat_[i].phaseState = staticNodeDat_[i].oldPhaseState;
            }

        /*!
         * \brief Return true if the primary variables were switched
         *        after the last timestep.
         */
        bool switched() const
            {
                return switchFlag_;
            }

        /*!
         * \brief Set whether there was a primary variable switch after in the last
         *        timestep.
         */
        void setSwitched(bool yesno)
            {
                switchFlag_ = yesno;
            }

        /*!
         * \brief Add the mass fraction of air in water to VTK output of
         *        the current timestep.
         */
        template <class MultiWriter>
        void addVtkFields(MultiWriter &writer, const SpatialFunction &globalSol) const
            {
                typedef Dune::BlockVector<Dune::FieldVector<DT, 1> > ScalarField;

                // create the required scalar fields
                unsigned nNodes = this->problem_.numNodes();
                ScalarField *pW =           writer.template createField<DT, 1>(nNodes);
                ScalarField *pN =           writer.template createField<DT, 1>(nNodes);
                ScalarField *pC =           writer.template createField<DT, 1>(nNodes);
                ScalarField *Sw =           writer.template createField<DT, 1>(nNodes);
                ScalarField *Sn =           writer.template createField<DT, 1>(nNodes);
                ScalarField *mobW =         writer.template createField<DT, 1>(nNodes);
                ScalarField *mobN =         writer.template createField<DT, 1>(nNodes);
                ScalarField *massfracAinW = writer.template createField<DT, 1>(nNodes);
                ScalarField *massfracAinA = writer.template createField<DT, 1>(nNodes);
                ScalarField *massfracWinW = writer.template createField<DT, 1>(nNodes);
                ScalarField *massfracWinA = writer.template createField<DT, 1>(nNodes);
                ScalarField *temperature  = writer.template createField<DT, 1>(nNodes);
                ScalarField *phaseState   = writer.template createField<DT, 1>(nNodes);

                VariableNodeData tmp;
                CellIterator it = this->problem_.cellBegin();
                CellIterator endit = this->problem_.cellEnd();
                for (; it != endit; ++it) {
                    for (int i = 0; i < it->template count<dim>(); ++i) {
                        int globalI = this->problem_.nodeIndex(*it, i);
                        asImp_()->updateVarNodeData_(tmp,
                                                     (*globalSol)[globalI],
                                                     staticNodeDat_[globalI].phaseState,
                                                     *it,
                                                     i,
                                                     this->problem_,
                                                     Implementation::temperature_((*globalSol)[globalI]));

                        (*pW)[globalI] = tmp.pW;
                        (*pN)[globalI] = tmp.pN;
                        (*pC)[globalI] = tmp.pC;
                        (*Sw)[globalI] = tmp.satW;
                        (*Sn)[globalI] = tmp.satN;
                        (*mobW)[globalI] = tmp.mobility[wPhase];
                        (*mobN)[globalI] = tmp.mobility[nPhase];
                        (*massfracAinW)[globalI] = tmp.massfrac[nComp][wPhase];
                        (*massfracAinA)[globalI] = tmp.massfrac[nComp][nPhase];
                        (*massfracWinW)[globalI] = tmp.massfrac[wComp][wPhase];
                        (*massfracWinA)[globalI] = tmp.massfrac[wComp][nPhase];
                        (*temperature)[globalI] = Implementation::temperature_((*globalSol)[globalI]);
                        (*phaseState)[globalI] = staticNodeDat_[globalI].phaseState;
                    };
                }


                writer.addVertexData(pW, "pW");
                writer.addVertexData(pN, "pN");
                writer.addVertexData(pC, "pC");
                writer.addVertexData(Sw, "Sw");
                writer.addVertexData(Sn, "Sn");
                writer.addVertexData(mobW, "mobW");
                writer.addVertexData(mobN, "mobN");
                writer.addVertexData(massfracAinW, "Xaw");
                writer.addVertexData(massfracAinA, "Xaa");
                writer.addVertexData(massfracWinW, "Xww");
                writer.addVertexData(massfracWinA, "Xwa");
                writer.addVertexData(temperature, "T");
                writer.addVertexData(phaseState, "phase state");
            }


    protected:
        Implementation *asImp_()
        { return static_cast<Implementation *>(this); }
        const Implementation *asImp_() const
        { return static_cast<const Implementation *>(this); }


        void updateCellCache_(CellCache &dest, const LocalFunction &sol, bool isOldSol)
            {
                int phaseState;
                int nNodes = this->curCell_().template count<dim>();
                for (int i = 0; i < nNodes; i++) {
                    int iGlobal = ParentType::problem_.nodeIndex(ParentType::curCell_(), i);
                    phaseState = isOldSol?staticNodeDat_[iGlobal].oldPhaseState:staticNodeDat_[iGlobal].phaseState;
                    asImp_()->updateVarNodeData_(dest.atSCV[i],
                                                 sol[i],
                                                 phaseState,
                                                 this->curCell_(),
                                                 i,
                                                 this->problem_,
                                                 Implementation::temperature_(sol[i]));
                }
            }


        //  perform variable switch at a node. Retrurns true iff a
        //  variable switch was performed.
        bool primaryVarSwitch_(SpatialFunction &sol,
                               int globalIdx,
                               const WorldCoord &globalPos)
            {
                // evaluate primary variable switch
                int phaseState    = staticNodeDat_[globalIdx].phaseState;
                int newPhaseState = phaseState;

                // Evaluate saturation and pressures
                DT pW = (*sol)[globalIdx][pWIndex];
                DT satW = 0.0;
                DT temperature = Implementation::temperature_((*sol)[globalIdx]);
                if      (phaseState == bothPhases) satW = 1.0 - (*sol)[globalIdx][switchIndex];
                else if (phaseState == wPhaseOnly) satW = 1.0;
                else if (phaseState == nPhaseOnly) satW = 0.0;

                DT pC = this->problem_.pC(satW, globalIdx, globalPos);
                DT pN = pW + pC;

                if (phaseState == nPhaseOnly)
                {
                    DT xWN = (*sol)[globalIdx][switchIndex];
                    DT xWNmax = this->problem_.multicomp().xWN(pN, temperature);

                    if (xWN > xWNmax*(1 + 2e-5))
                    {
                        // wetting phase appears
                        std::cout << "wetting phase appears at node " << globalIdx
                                  << ", coordinates: " << globalPos << std::endl;
                        newPhaseState = bothPhases;
                        (*sol)[globalIdx][switchIndex] = 1.0 - 2e-5;
                    };
                }
                else if (phaseState == wPhaseOnly)
                {
                    DT xAW = (*sol)[globalIdx][switchIndex];
                    DT xAWmax = this->problem_.multicomp().xAW(pN, temperature);

                    if (xAW > xAWmax*(1 + 2e-5))
                    {
                        // non-wetting phase appears
                        std::cout << "Non-wetting phase appears at node " << globalIdx
                                  << ", coordinates: " << globalPos << std::endl;
                        (*sol)[globalIdx][switchIndex] = 2e-5;
                        newPhaseState = bothPhases;
                    }
                }
                else if (phaseState == bothPhases) {
                    DT satN = 1 - satW;

                    if (satN < -1e-5) {
                        // non-wetting phase disappears
                        std::cout << "Non-wetting phase disappears at node " << globalIdx
                                  << ", coordinates: " << globalPos << std::endl;
                        (*sol)[globalIdx][switchIndex]
                            = this->problem_.multicomp().xAW(pN, temperature);
                        newPhaseState = wPhaseOnly;
                    }
                    else if (satW < -1e-5) {
                        // wetting phase disappears
                        std::cout << "Wetting phase disappears at node " << globalIdx
                                  << ", coordinates: " << globalPos << std::endl;
                        (*sol)[globalIdx][switchIndex]
                            = this->problem_.multicomp().xWN(pN, temperature);
                        newPhaseState = nPhaseOnly;
                    }
                }

                staticNodeDat_[globalIdx].phaseState = newPhaseState;

                return phaseState != newPhaseState;
            }

        // harmonic mean of the permeability computed directly.  the
        // first parameter is used to store the result.
        static void harmonicMeanK_(Tensor &Ki, const Tensor &Kj)
            {
                double eps = 1e-20;

                for (int kx=0; kx < Tensor::rows; kx++){
                    for (int ky=0; ky< Tensor::cols; ky++){
                        if (Ki[kx][ky] != Kj[kx][ky]) {
                            Ki[kx][ky] = 2 / (1/(Ki[kx][ky]+eps) + (1/(Kj[kx][ky]+eps)));
                        }
                    }
                }
            }

        // parameters given in constructor
        std::vector<StaticNodeData> staticNodeDat_;
        bool                        switchFlag_;

        // current solution
        LocalFunction    curSol_;
        CellCache        curSolCache_;

        // needed for restoreCurSolution()
        bool             curSolDeflected_;
        DT           curSolOrigValue_;
        VariableNodeData curSolOrigVarData_;

        // previous solution
        LocalFunction   prevSol_;
        CellCache       prevSolCache_;
    };


    template<class ProblemT,
             class BoxTraitsT,
             class TwoPTwoCTraitsT>
    class TwoPTwoCBoxJacobian : public TwoPTwoCBoxJacobianBase<ProblemT,
                                                               BoxTraitsT,
                                                               TwoPTwoCTraitsT,
                                                               TwoPTwoCBoxJacobian<ProblemT,
                                                                                   BoxTraitsT,
                                                                                   TwoPTwoCTraitsT>
                                                               >
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

        typedef typename DomTraits::Scalar             DT;

        typedef typename DomTraits::WorldCoord         WorldCoord;
        typedef typename DomTraits::LocalCoord         LocalCoord;

        typedef typename BoxTraits::UnknownsVector     UnknownsVector;
        typedef typename BoxTraits::FVElementGeometry  FVElementGeometry;

        typedef typename ParentType::VariableNodeData  VariableNodeData;
        typedef typename ParentType::LocalFunction     LocalFunction;
        typedef typename ParentType::CellCache         CellCache;

        typedef typename TwoPTwoCTraits::PhasesVector PhasesVector;


    public:
        TwoPTwoCBoxJacobian(ProblemT &problem)
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
                // only relevant for the non-isothermal model!
            }

        /*!
         * \brief Update the temperature gradient at a face of a FV
         *        cell.
         */
        void updateTempGrad(WorldCoord &tempGrad,
                            const WorldCoord &feGrad,
                            const LocalFunction &sol,
                            int nodeIdx) const
            {
                // only relevant for the non-isothermal model!
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
                // only relevant for the non-isothermal model!
            }

        /*!
         * \brief Adds the diffusive heat flux to the flux vector over
         *        the face of a sub-control volume.
         */
        void diffusiveHeatFlux(UnknownsVector &flux,
                               int faceIdx,
                               const WorldCoord &tempGrad) const
            {
                // only relevant for the non-isothermal model!
            }

        // internal method!
        static DT temperature_(const UnknownsVector &sol)
            { return 283.15; }
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
    template<class ProblemT>
    class TwoPTwoCBoxModel
        : public BoxScheme<TwoPTwoCBoxModel<ProblemT>, // Implementation of the box scheme

                           // The Traits for the BOX method
                           P1BoxTraits<typename ProblemT::DomainTraits::Scalar,
                                       typename ProblemT::DomainTraits::Grid,
                                       TwoPTwoCTraits<typename ProblemT::DomainTraits::Scalar>::numEq>,

                           // The actual problem we would like to solve
                           ProblemT,

                           // The local jacobian operator
                           TwoPTwoCBoxJacobian<ProblemT,
                                               P1BoxTraits<typename ProblemT::DomainTraits::Scalar,
                                                           typename ProblemT::DomainTraits::Grid,
                                                           TwoPTwoCTraits<typename ProblemT::DomainTraits::Scalar>::numEq>,
                                               TwoPTwoCTraits<typename ProblemT::DomainTraits::Scalar> > >
    {
        typedef typename ProblemT::DomainTraits::Grid   Grid;
        typedef typename ProblemT::DomainTraits::Scalar DT;
        typedef TwoPTwoCBoxModel<ProblemT>              ThisType;

    public:
        typedef Dune::TwoPTwoCTraits<DT>                           TwoPTwoCTraits;
        typedef P1BoxTraits<DT, Grid, TwoPTwoCTraits::numEq> BoxTraits;

    private:
        typedef TwoPTwoCBoxJacobian<ProblemT, BoxTraits, TwoPTwoCTraits>  TwoPTwoCLocalJacobian;
        typedef BoxScheme<ThisType,
                          BoxTraits,
                          ProblemT,
                          TwoPTwoCLocalJacobian>        ParentType;

        typedef typename ProblemT::DomainTraits           DomTraits;
        typedef typename DomTraits::Cell                  Cell;
        typedef typename DomTraits::CellIterator          CellIterator;
        typedef typename DomTraits::LocalCoord            LocalCoord;
        typedef typename DomTraits::WorldCoord            WorldCoord;

        enum {
            dim              = DomTraits::GridDim,
            dimWorld         = DomTraits::WorldDim
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
        TwoPTwoCLocalJacobian  twoPTwoCLocalJacobian_;
    };
}

#endif
