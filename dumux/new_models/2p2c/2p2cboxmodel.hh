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
    class TwoPTwoCTraits
    {
    public:
        enum {
            PrimaryVars   = 2,  //!< Number of primary variables
            NumPhases     = 2,  //!< Number of fluid phases
            NumComponents = 2   //!< Number of fluid components within a phase
        };
        enum { // Solution vector indices
            PwIndex = 0,      //!< Index for the wetting phase pressure in a field vector
            SwitchIndex = 1   //!< Index for the non-wetting phase quantity
        };
        enum { // Phase Indices
            WPhaseIndex = 0,  //!< Index of the wetting phase
            NPhaseIndex = 1   //!< Index of the non-wetting phase
        };
        enum { // Component indices
            WCompIndex = 0,  //!< Index of the wetting component
            NCompIndex = 1   //!< Index of the non-wetting component
        };
        enum { // present phases 
            NPhaseOnly = 0, //!< Only the non-wetting phase is present
            WPhaseOnly = 1, //!< Only the wetting phase is present
            BothPhases = 2 //!< Both phases are present
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
    template<class ProblemT, class BoxTraitsT, class TwoPTwoCTraitsT>
    class TwoPTwoCBoxJacobian : public BoxJacobian<ProblemT,
                                                   BoxTraitsT, 
                                                   TwoPTwoCBoxJacobian<ProblemT, 
                                                                       BoxTraitsT, 
                                                                       TwoPTwoCTraitsT> >
    {
    private:
        typedef TwoPTwoCBoxJacobian<ProblemT, BoxTraitsT, TwoPTwoCTraitsT>  ThisType;
        typedef BoxJacobian<ProblemT, BoxTraitsT, ThisType>                 ParentType;

        typedef ProblemT                                Problem;
        typedef typename Problem::DomainTraits          DomTraits;
        typedef BoxTraitsT                              BoxTraits;
        typedef TwoPTwoCTraitsT                         TwoPTwoCTraits;

        enum {
            GridDim          = DomTraits::GridDim,
            WorldDim         = DomTraits::WorldDim,

            PrimaryVariables = BoxTraits::PrimaryVariables,
            NumPhases        = TwoPTwoCTraits::NumPhases,
            NumComponents    = TwoPTwoCTraits::NumComponents,

            PwIndex          = TwoPTwoCTraits::PwIndex,
            SwitchIndex      = TwoPTwoCTraits::SwitchIndex,

            WPhaseIndex      = TwoPTwoCTraits::WPhaseIndex,
            NPhaseIndex      = TwoPTwoCTraits::NPhaseIndex,

            WCompIndex       = TwoPTwoCTraits::WCompIndex,
            NCompIndex       = TwoPTwoCTraits::NCompIndex,

            WPhaseOnly       = TwoPTwoCTraits::WPhaseOnly,
            NPhaseOnly       = TwoPTwoCTraits::NPhaseOnly,
            BothPhases       = TwoPTwoCTraits::BothPhases
        };

        
        typedef typename DomTraits::Scalar                Scalar;
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

        typedef FieldMatrix<Scalar, GridDim, GridDim>  Tensor;

        /*!
         * \brief Data which is attached to each node of the and can
         *        be shared between multiple calculations and should
         *        thus be cached in order to increase efficency.
         */
        struct VariableNodeData
        {
            Scalar satN;
            Scalar satW;
            Scalar pW;
            Scalar pC;
            Scalar pN;
            UnknownsVector mobility;  //Vector with the number of phases
            UnknownsVector density;
            FieldMatrix<Scalar, NumComponents, NumPhases> massfrac;
            int phasestate;

            void update(const UnknownsVector &nodeSol, 
                        int phaseState,
                        const Cell &cell, 
                        int localIdx,
                        Problem &problem,
                        Scalar temperature)
                {
                    const WorldCoord &global = cell.geometry()[localIdx];
                    const LocalCoord &local =
                        DomTraits::referenceElement(cell.type()).position(localIdx,
                                                                          GridDim);
                    
                    pW = nodeSol[PwIndex];
                    if (phaseState == BothPhases) satN = nodeSol[SwitchIndex];
                    else if (phaseState == WPhaseOnly) satN = 0.0;
                    else if (phaseState == NPhaseOnly) satN = 1.0;
                    else DUNE_THROW(Dune::InvalidStateException, "Phase state " << phaseState << " is invalid.");
                    
                    satW = 1.0 - satN;
                    pC = problem.materialLaw().pC(satW, 
                                                  global,
                                                  cell, 
                                                  local);
                    pN = pW + pC;
                    
                    // Solubilities of components in phases
                    if (phaseState == BothPhases) {
                        massfrac[NCompIndex][WPhaseIndex] = problem.multicomp().xAW(pN, temperature);
                        massfrac[WCompIndex][NPhaseIndex] = problem.multicomp().xWN(pN, temperature);
                    }
                    else if (phaseState == WPhaseOnly) {
                        massfrac[WCompIndex][NPhaseIndex] = 0.0;
                        massfrac[NCompIndex][WPhaseIndex] =  nodeSol[SwitchIndex];
                    }
                    else if (phaseState == NPhaseOnly){
                        massfrac[WCompIndex][NPhaseIndex] = nodeSol[SwitchIndex];
                        massfrac[NCompIndex][WPhaseIndex] = 0.0;
                    }
                    else DUNE_THROW(Dune::InvalidStateException, "Phase state " << phaseState << " is invalid.");

                    massfrac[WCompIndex][WPhaseIndex] = 1.0 - massfrac[NCompIndex][WPhaseIndex];
                    massfrac[NCompIndex][NPhaseIndex] = 1.0 - massfrac[WCompIndex][NPhaseIndex];
                    phasestate = phaseState;

                    // Mobilities & densities
                    mobility[WPhaseIndex] = problem.materialLaw().mobW(satW, global, cell, local, temperature, pW);
                    mobility[NPhaseIndex] = problem.materialLaw().mobN(satN, global, cell, local, temperature, pN);
                    
                    // Density of Water is set constant here!
                    density[WPhaseIndex] = 1000;//this->problem_.wettingPhase().density(temperature, pN);
                    density[NPhaseIndex] = problem.nonwettingPhase().density(temperature, 
                                                                             pN,
                                                                             massfrac[WCompIndex][NPhaseIndex]);
                }
        };
        
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
            bool visited;
            int numSwitches;
            int phaseState;
            int oldPhaseState;
            Scalar porosity;
        };

    public:
        TwoPTwoCBoxJacobian(ProblemT &problem) 
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
        void deflectCurSolution(int node, int component, Scalar value)
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
                curSolCache_.atSCV[node].update(curSol_[node],
                                                staticNodeDat_[globalIdx].phaseState,
                                                this->curCell_(),
                                                node,
                                                this->problem_,
                                                temperature);
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
                result = Scalar(0);

//                const LocalFunction &sol   = usePrevSol ? *prevSol_     : *curSol_;
                const CellCache &cellCache = usePrevSol ? prevSolCache_ : curSolCache_;
                
                int globalIdx = this->problem_.nodeIndex(this->curCell_(), scvId);
                
                Scalar satN = cellCache.atSCV[scvId].satN;
                Scalar satW = cellCache.atSCV[scvId].satW;
                
                Scalar porosity = staticNodeDat_[globalIdx].porosity;

                // storage of component water
                result[PwIndex] =
                    porosity*(cellCache.atSCV[scvId].density[WPhaseIndex]*
                              satW*
                              cellCache.atSCV[scvId].massfrac[WCompIndex][WPhaseIndex]
                              + cellCache.atSCV[scvId].density[NPhaseIndex]*
                                satN*
                                cellCache.atSCV[scvId].massfrac[WCompIndex][NPhaseIndex]);

                // storage of component air
                result[SwitchIndex] =
                    porosity*(cellCache.atSCV[scvId].density[NPhaseIndex]*
                              satN*
                              cellCache.atSCV[scvId].massfrac[NCompIndex][NPhaseIndex]
                              + cellCache.atSCV[scvId].density[WPhaseIndex]*
                                satW*
                                cellCache.atSCV[scvId].massfrac[NCompIndex][WPhaseIndex]);
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

                const LocalCoord &local_i = this->curCellGeom_.subContVol[i].local;
                const LocalCoord &local_j = this->curCellGeom_.subContVol[j].local;
               
                WorldCoord pGrad[PrimaryVariables];
                WorldCoord xGrad[PrimaryVariables];
                for (int k = 0; k < PrimaryVariables; ++k) {
                    pGrad[k] = Scalar(0);
                    xGrad[k] = Scalar(0);
                }

                UnknownsVector tmp(0.0);
                UnknownsVector pressure(0.0), massfrac(0.0);

                // calculate harmonic mean of permeabilities of nodes i and j
                Tensor K         = this->problem_.soil().K(global_i, ParentType::curCell_(), local_i);
                const Tensor &Kj = this->problem_.soil().K(global_j, ParentType::curCell_(), local_j);
                harmonicMeanK_(K, Kj);

                // calculate FE gradient (grad p for each phase)
                for (int k = 0; k < this->curCellGeom_.nNodes; k++) // loop over adjacent nodes
                {
                    // FEGradient at node k
                    const LocalCoord &feGrad = this->curCellGeom_.subContVolFace[faceId].grad[k];

                    pressure[WPhaseIndex] = curSolCache_.atSCV[k].pW;
                    pressure[NPhaseIndex] = curSolCache_.atSCV[k].pN;

                    // compute sum of pressure gradients for each phase
                    for (int phase = 0; phase < PrimaryVariables; phase++)
                    {
                        tmp = feGrad;
                        
                        tmp *= pressure[phase];

                        pGrad[phase] += tmp;
                    }

                    // for diffusion of air in wetting phase
                    tmp = feGrad;
                    tmp *= curSolCache_.atSCV[k].massfrac[NCompIndex][WPhaseIndex];
                    xGrad[WPhaseIndex] += tmp;

                    // for diffusion of water in nonwetting phase
                    tmp = feGrad;
                    tmp *= curSolCache_.atSCV[k].massfrac[WCompIndex][NPhaseIndex];
                    xGrad[NPhaseIndex] += tmp;
                }

                // deduce gravity*density of each phase
                UnknownsVector contribComp[NumPhases];
                for (int phase=0; phase < NumPhases; phase++)
                {
                    contribComp[phase] = this->problem_.gravity();
                    contribComp[phase] *= curSolCache_.atSCV[i].density[phase];
                    pGrad[phase] -= contribComp[phase]; // grad p - rho*g
                }

                UnknownsVector outward;  // Darcy velocity of each phase

                // calculate the advective flux using upwind: K*n(grad p -rho*g)
                for (int phase=0; phase < NumPhases; phase++)
                {
                    LocalCoord v_tilde(0);
                    K.mv(pGrad[phase], v_tilde);  // v_tilde=K*gradP
                    outward[phase] = v_tilde*normal;
                }

                // evaluate upwind nodes
                int up_w, dn_w, up_n, dn_n;
                if (outward[WPhaseIndex] <= 0) {
                    up_w = i; dn_w = j;
                }
                else {
                    up_w = j; dn_w = i;
                };
                
                if (outward[NPhaseIndex] <= 0) {
                    up_n = i; dn_n = j;
                }
                else {
                    up_n = j; dn_n = i;
                };

                Scalar alpha = 1.0;  // Upwind parameter

                // water conservation
                flux[WCompIndex] =   (alpha* curSolCache_.atSCV[up_w].density[WPhaseIndex]*curSolCache_.atSCV[up_w].mobility[WPhaseIndex]
                                      * curSolCache_.atSCV[up_w].massfrac[WCompIndex][WPhaseIndex]
                                      + (1-alpha)* curSolCache_.atSCV[dn_w].density[WPhaseIndex]*curSolCache_.atSCV[dn_w].mobility[WPhaseIndex]
                                      * curSolCache_.atSCV[dn_w].massfrac[WCompIndex][WPhaseIndex])
                    * outward[WPhaseIndex];
                flux[WCompIndex] +=  (alpha* curSolCache_.atSCV[up_n].density[NPhaseIndex]*curSolCache_.atSCV[up_n].mobility[NPhaseIndex]
                                      * curSolCache_.atSCV[up_n].massfrac[WCompIndex][NPhaseIndex]
                                      + (1-alpha)* curSolCache_.atSCV[dn_n].density[NPhaseIndex]*curSolCache_.atSCV[dn_n].mobility[NPhaseIndex]
                                      * curSolCache_.atSCV[dn_n].massfrac[WCompIndex][NPhaseIndex])
                    * outward[NPhaseIndex];
                

                // air conservation
                flux[NCompIndex]   = (alpha* curSolCache_.atSCV[up_n].density[NPhaseIndex]*curSolCache_.atSCV[up_n].mobility[NPhaseIndex]
                                      * curSolCache_.atSCV[up_n].massfrac[NCompIndex][NPhaseIndex]
                                      + (1-alpha)* curSolCache_.atSCV[dn_n].density[NPhaseIndex]*curSolCache_.atSCV[dn_n].mobility[NPhaseIndex]
                                      * curSolCache_.atSCV[dn_n].massfrac[NCompIndex][NPhaseIndex])
                    * outward[NPhaseIndex];

                flux[NCompIndex]  +=   (alpha* curSolCache_.atSCV[up_w].density[WPhaseIndex]*curSolCache_.atSCV[up_w].mobility[WPhaseIndex]
                                        * curSolCache_.atSCV[up_w].massfrac[NCompIndex][WPhaseIndex]
                                        + (1-alpha)* curSolCache_.atSCV[dn_w].density[WPhaseIndex]*curSolCache_.atSCV[dn_w].mobility[WPhaseIndex]
                                        * curSolCache_.atSCV[dn_w].massfrac[NCompIndex][WPhaseIndex]) 
                    * outward[WPhaseIndex];
                

                return;
                
                // DIFFUSION
                UnknownsVector normDiffGrad;

                // get local to global id map
                int state_i = curSolCache_.atSCV[i].phasestate;
                int state_j = curSolCache_.atSCV[j].phasestate;

                Scalar diffusionWW(0.0), diffusionWN(0.0); // diffusion of water
                Scalar diffusionAW(0.0), diffusionAN(0.0); // diffusion of air
                UnknownsVector avgDensity, avgDpm;

                // Diffusion coefficent
                // TODO: needs to be continuously dependend on the phase saturations
                avgDpm[WPhaseIndex]=2e-9; 
                avgDpm[NPhaseIndex]=2.25e-5;
                if (state_i == NPhaseOnly || state_j == NPhaseOnly)
                {
                    // only the nonwetting phase is present in at
                    // least one cell -> no diffusion within the
                    // wetting phase
                    avgDpm[WPhaseIndex] = 0;
                }
                if (state_i == WPhaseOnly || state_j == WPhaseOnly)
                {
                    // only the wetting phase is present in at least
                    // one cell -> no diffusion within the non wetting
                    // phase
                    avgDpm[NPhaseIndex] = 0;
                }
                
                // length of the diffusion gradient
                normDiffGrad[WPhaseIndex] = xGrad[WPhaseIndex]*normal;
                normDiffGrad[NPhaseIndex] = xGrad[NPhaseIndex]*normal;

                // calculate the arithmetic mean of densities
                avgDensity[WPhaseIndex] = 0.5*(curSolCache_.atSCV[i].density[WPhaseIndex] + curSolCache_.atSCV[j].density[WPhaseIndex]);
                avgDensity[NPhaseIndex] = 0.5*(curSolCache_.atSCV[i].density[NPhaseIndex] + curSolCache_.atSCV[j].density[NPhaseIndex]);

                diffusionAW = avgDpm[WPhaseIndex] * avgDensity[WPhaseIndex] * normDiffGrad[WPhaseIndex];
                diffusionWW = - diffusionAW;
                diffusionWN = avgDpm[NPhaseIndex] * avgDensity[NPhaseIndex] * normDiffGrad[NPhaseIndex];
                diffusionAN = - diffusionWN;
                
                // add diffusion of water to water flux
                flux[WCompIndex] += (diffusionWW + diffusionWN);

                // add diffusion of air to air flux
                flux[NCompIndex] += (diffusionAN + diffusionAW);
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

                    // ASSUME porosity defined at nodes
                    staticNodeDat_[globalIdx].porosity = 
                        this->problem_.porosity(*it, globalIdx, globalPos);

                    // number of primary variable switches at the node
                    staticNodeDat_[globalIdx].numSwitches = 0;

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

        enum VtkFieldBits {
            WettingPressureBit      = 0x001,
            NonwettingPressureBit   = 0x002,
            CapillaryPressureBit    = 0x004,
            
            WettingSaturationBit    = 0x008,
            NonwettingSaturationBit = 0x010,
            
            WettingMobilityBit      = 0x020,
            NonwettingMobilityBit   = 0x040,
            
            MassfracAinWBit         = 0x080,
            MassfracAinABit         = 0x100,
            MassfracWinWBit         = 0x200,
            MassfracWinABit         = 0x400,
            
            PhaseStateBit           = 0x800,
            
            VtkAllFields = 0xfffff
        };
                    
        /*!
         * \brief Add the mass fraction of air in water to VTK output of
         *        the current timestep.
         */
        template <class MultiWriter>
        void addVtkFields(MultiWriter &writer, const SpatialFunction &globalSol, int bitmask) const
            {
                typedef Dune::BlockVector<Dune::FieldVector<Scalar, 1> > ScalarField;
                
                // create the required scalar fields
                unsigned nNodes = this->problem_.numNodes();
                ScalarField *pW =           (bitmask&WettingPressureBit)      ? writer.template createField<Scalar, 1>(nNodes):NULL;
                ScalarField *pN =           (bitmask&NonwettingPressureBit)   ? writer.template createField<Scalar, 1>(nNodes):NULL;
                ScalarField *pC =           (bitmask&CapillaryPressureBit)    ? writer.template createField<Scalar, 1>(nNodes):NULL;
                ScalarField *Sw =           (bitmask&WettingSaturationBit)    ? writer.template createField<Scalar, 1>(nNodes):NULL;
                ScalarField *Sn =           (bitmask&NonwettingSaturationBit) ? writer.template createField<Scalar, 1>(nNodes):NULL;
                ScalarField *mobW =         (bitmask&WettingMobilityBit)      ? writer.template createField<Scalar, 1>(nNodes):NULL;
                ScalarField *mobN =         (bitmask&NonwettingMobilityBit)   ? writer.template createField<Scalar, 1>(nNodes):NULL;
                ScalarField *massfracAinW = (bitmask&MassfracAinWBit)         ? writer.template createField<Scalar, 1>(nNodes):NULL;
                ScalarField *massfracAinA = (bitmask&MassfracAinABit)         ? writer.template createField<Scalar, 1>(nNodes):NULL;
                ScalarField *massfracWinW = (bitmask&MassfracWinWBit)         ? writer.template createField<Scalar, 1>(nNodes):NULL;
                ScalarField *massfracWinA = (bitmask&MassfracWinABit)         ? writer.template createField<Scalar, 1>(nNodes):NULL;
                ScalarField *phaseState   = (bitmask&PhaseStateBit)           ? writer.template createField<Scalar, 1>(nNodes):NULL;
                
                VariableNodeData tmp;
                CellIterator it = this->problem_.cellBegin();
                CellIterator endit = this->problem_.cellEnd();
                for (; it != endit; ++it) {
                    for (int i = 0; i < it->template count<GridDim>(); ++i) {
                        int globalI = this->problem_.nodeIndex(*it, i);
                        tmp.update((*globalSol)[globalI],
                                   staticNodeDat_[globalI].phaseState,
                                   *it,
                                   i,
                                   this->problem_,
                                   temperature);
                        
                        if (bitmask&WettingPressureBit)
                            (*pW)[globalI] = tmp.pW;
                        if (bitmask&NonwettingPressureBit)
                            (*pN)[globalI] = tmp.pN;
                        if (bitmask&CapillaryPressureBit)
                            (*pC)[globalI] = tmp.pC;
                        if (bitmask&WettingSaturationBit)
                            (*Sw)[globalI] = tmp.satW;
                        if (bitmask&NonwettingSaturationBit)
                            (*Sn)[globalI] = tmp.satN;
                        if (bitmask&WettingMobilityBit)
                            (*mobW)[globalI] = tmp.mobility[PwIndex];
                        if (bitmask&NonwettingMobilityBit)
                            (*mobN)[globalI] = tmp.mobility[SwitchIndex];
                        if (bitmask&MassfracAinWBit)
                            (*massfracAinW)[globalI] = tmp.massfrac[NCompIndex][WPhaseIndex];
                        if (bitmask&MassfracAinABit)
                            (*massfracAinA)[globalI] = tmp.massfrac[NCompIndex][NPhaseIndex];
                        if (bitmask&MassfracWinWBit)
                            (*massfracWinW)[globalI] = tmp.massfrac[WCompIndex][WPhaseIndex];
                        if (bitmask&MassfracWinABit)
                            (*massfracWinA)[globalI] = tmp.massfrac[WCompIndex][NPhaseIndex];
                        if (bitmask&PhaseStateBit)
                            (*phaseState)[globalI] = staticNodeDat_[globalI].phaseState;
                    };
                }

                
                if (bitmask&WettingPressureBit)
                    writer.addVertexData(pW, "pW");
                if (bitmask&NonwettingPressureBit)
                    writer.addVertexData(pN, "pN");
                if (bitmask&CapillaryPressureBit)
                    writer.addVertexData(pC, "pC");
                if (bitmask&WettingSaturationBit)
                    writer.addVertexData(Sw, "Sw");
                if (bitmask&NonwettingSaturationBit)
                    writer.addVertexData(Sn, "Sn");
                if (bitmask&WettingMobilityBit)
                    writer.addVertexData(mobW, "mobW");
                if (bitmask&NonwettingMobilityBit)
                    writer.addVertexData(mobN, "mobN");
                if (bitmask&MassfracAinWBit)
                    writer.addVertexData(massfracAinW, "massfrac Air in Water");
                if (bitmask&MassfracAinABit)
                    writer.addVertexData(massfracAinA, "massfrac Air in Air");
                if (bitmask&MassfracWinWBit)
                    writer.addVertexData(massfracWinW, "massfrac Water in Water");
                if (bitmask&MassfracWinABit)
                    writer.addVertexData(massfracWinA, "massfrac Water in Air");
                if (bitmask&PhaseStateBit)
                    writer.addVertexData(phaseState, "phaseState");
            }

        
    private:
        void updateCellCache_(CellCache &dest, const LocalFunction &sol, bool isOldSol)
            {
                int phaseState;
                int nNodes = this->curCell_().template count<GridDim>();
                for (int i = 0; i < nNodes; i++) {
                    int iGlobal = ParentType::problem_.nodeIndex(ParentType::curCell_(), i);
                    phaseState = isOldSol?staticNodeDat_[iGlobal].oldPhaseState:staticNodeDat_[iGlobal].phaseState;
                    
                    dest.atSCV[i].update(sol[i], 
                                         phaseState,
                                         this->curCell_(), 
                                         i,
                                         this->problem_,
                                         temperature);
                }
            }

        
        //  perform variable switch at a node. Retrurns true iff a
        //  variable switch was performed.
   	bool primaryVarSwitch_(SpatialFunction &sol,
                               int globalIdx,
                               const WorldCoord &globalPos)
            {
                // evaluate primary variable switch
                int  phaseState = staticNodeDat_[globalIdx].phaseState;

                // Evaluate saturation and pressures
                Scalar pW = (*sol)[globalIdx][PwIndex];
                Scalar satW = 0.0;
                if (phaseState == BothPhases) satW = 1.0 - (*sol)[globalIdx][SwitchIndex];
  		if (phaseState == WPhaseOnly) satW = 1.0;
  		if (phaseState == NPhaseOnly) satW = 0.0;

                Scalar pC = this->problem_.pC(satW, globalIdx, globalPos);
  		Scalar pN = pW + pC;

                if (phaseState == NPhaseOnly) {
                    // auxiliary variables
                    Scalar xWNmass, xWNmolar, pwn, pWSat;
                    
                    xWNmass = (*sol)[globalIdx][SwitchIndex];
                    xWNmolar = this->problem_.multicomp().xWNmolar(pN, temperature);
                    pwn = xWNmolar * pN;
                    pWSat = this->problem_.multicomp().vaporPressure(temperature);
                    
                    if (pwn > (1 + 1e-5)*pWSat)
                    {
                        // appearance of water phase
                        std::cout << "Water appears at node " << globalIdx << "  Coordinates: " << globalPos << std::endl;
                        staticNodeDat_[globalIdx].phaseState = BothPhases;
                        (*sol)[globalIdx][SwitchIndex] = 1.0 - 2e-5; // initialize solution vector
                        staticNodeDat_[globalIdx].numSwitches += 1;
                    }
                }
                else if (phaseState == WPhaseOnly) 
                {
                    Scalar pbub, henryInv, pWSat;
                    
                    Scalar xAWmass = (*sol)[globalIdx][SwitchIndex];
                    Scalar xAWmolar = this->problem_.multicomp().convertMassToMoleFraction(xAWmass, WPhaseIndex);

//                    Scalar xAWmolarMax = this->problem_.multicomp().xAWmolar(pN, temperature);

                    henryInv = this->problem_.multicomp().henry(temperature);
                    pWSat = this->problem_.multicomp().vaporPressure(temperature);
                    pbub = pWSat + xAWmolar/henryInv; // pWSat + pAW

                    if (pN < (1 - 1e-5)*pbub)
//                        if (xAWmolar > (1 + 1e-2)*xAWmolarMax)
                    {
                        // appearance of gas phase
                        std::cout << "Gas appears at node " << globalIdx << ",  Coordinates: " << globalPos << std::endl;
//                        std::cerr << "xAWmolarMax: " << xAWmolarMax << " henry: " << henryInv << " temperature: " << temperature << " pbub: " << pbub << " pN: " << pN << " pWSat: " << pWSat << "\n";
                        staticNodeDat_[globalIdx].phaseState = BothPhases;
                        (*sol)[globalIdx][SwitchIndex] = 2e-5; // initialize solution vector
                        staticNodeDat_[globalIdx].numSwitches += 1;
                    }
                }
                else if (phaseState == BothPhases) {
                    Scalar satN = (*sol)[globalIdx][SwitchIndex];
                    
                    if (satN < -1e-5)
                    {
                        // disappearance of gas phase
                        std::cout << "Gas disappears at node " << globalIdx << "  Coordinates: " << globalPos << std::endl;
                        staticNodeDat_[globalIdx].phaseState = WPhaseOnly;
                        (*sol)[globalIdx][SwitchIndex] = this->problem_.multicomp().xAW(pN); // initialize solution vector
                        staticNodeDat_[globalIdx].numSwitches += 1;
                    }
                    else if (satW < -1e-5)
                    {
                        // disappearance of water phase
                        std::cout << "Water disappears at node " << globalIdx << "  Coordinates: " << globalPos << std::endl;
                        staticNodeDat_[globalIdx].phaseState = NPhaseOnly;
                        (*sol)[globalIdx][SwitchIndex] = this->problem_.multicomp().xWN(pN); // initialize solution vector
                        staticNodeDat_[globalIdx].numSwitches += 1;
                    }
                }
                else 
                    DUNE_THROW(Dune::InvalidStateException, "Phase state " << phaseState << " is invalid.");
               
                return phaseState != staticNodeDat_[globalIdx].phaseState;
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

        static const Scalar temperature = 283.15;

        // current solution
        LocalFunction    curSol_;
        CellCache        curSolCache_;

        // needed for restoreCurSolution()
        bool             curSolDeflected_;
        Scalar           curSolOrigValue_;
        VariableNodeData curSolOrigVarData_;

        // previous solution
        LocalFunction   prevSol_;
        CellCache       prevSolCache_;
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
                                       TwoPTwoCTraits::PrimaryVars>,
        
                           // The actual problem we would like to solve
                           ProblemT, 
                           
                           // The local jacobian operator
                           TwoPTwoCBoxJacobian<ProblemT, 
                                               P1BoxTraits<typename ProblemT::DomainTraits::Scalar,
                                                           typename ProblemT::DomainTraits::Grid,
                                                           TwoPTwoCTraits::PrimaryVars>,
                                               TwoPTwoCTraits> >
    {
        typedef typename ProblemT::DomainTraits::Grid   Grid;
        typedef typename ProblemT::DomainTraits::Scalar Scalar;
        typedef TwoPTwoCBoxModel<ProblemT>              ThisType;
        
    public:
        typedef P1BoxTraits<Scalar, Grid, TwoPTwoCTraits::PrimaryVars> BoxTraits;
        typedef Dune::TwoPTwoCTraits                                   TwoPTwoCTraits;
        
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
            GridDim          = DomTraits::GridDim,
            WorldDim         = DomTraits::WorldDim
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
        void addVtkFields(MultiWriter &writer, 
                       int bitmask = 0xffffff) const
            {
                twoPTwoCLocalJacobian_.addVtkFields(writer, this->currentSolution(), bitmask);
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
