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

#include <dumux/new_models/box/boxscheme.hh>
#include <dumux/new_models/box/p1boxtraits.hh>

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
            WPhaseOnly = 0, //!< Only the wetting phase is present
            NPhaseOnly = 1, //!< Only the non-wetting phase is present
            BothPhases  = 2 //!< Both phases are present
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
    class TwoPTwoCBoxJacobian : public BoxJacobian<ProblemT, BoxTraitsT, TwoPTwoCBoxJacobian<ProblemT, BoxTraitsT, TwoPTwoCTraitsT> >
    {
    private:
        typedef TwoPTwoCBoxJacobian<ProblemT, BoxTraitsT, TwoPTwoCTraitsT>  ThisType;
        typedef BoxJacobian<ProblemT, BoxTraitsT, ThisType >                ParentType;

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
        typedef typename Cell::EntityPointer              CellPointer;
        typedef typename DomTraits::LocalCoord            LocalCoord;
        typedef typename DomTraits::WorldCoord            WorldCoord;

        typedef typename BoxTraits::UnknownsVector      UnknownsVector;
        typedef typename BoxTraits::FVElementGeometry   FVElementGeometry;
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
            Scalar temperature;
            UnknownsVector mobility;  //Vector with the number of phases
            UnknownsVector density;
            FieldMatrix<Scalar, NumComponents, NumPhases> massfrac;
            int phasestate;
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
            int switched;
            int phaseState;
            int oldPhaseState;
            Scalar cellVolume;
            Scalar porosity;
            Tensor K;
        };

    public:
        TwoPTwoCBoxJacobian(ProblemT &problem,
                            bool levelBoundaryAsDirichlet = false,
                            bool procBoundaryAsDirichlet = true) 
            : ParentType(problem,
                         levelBoundaryAsDirichlet,
                         procBoundaryAsDirichlet),
              _staticNodeDat(problem.numNodes())
            {};

        /*!
         * \brief Set the current grid cell.
         */
        void setCurrentCell(const Cell &cell) 
            {
                if (ParentType::_setCurrentCell(cell)) {
//                    _curCellPorosity = ParentType::_problem.cellPorosity(ParentType::_curCell());
                }
            };

        /*!
         * \brief Set the parameters for the calls to the remaining
         *        members.
         */
        void setParams(const Cell &cell, LocalFunction &curSol, LocalFunction &prevSol)
            {
                setCurrentCell(cell);
                
                _curSol = &curSol;
                _updateCellCache(_curSolCache, *_curSol, false);
                _curSolDeflected = false;
                
                _prevSol = &prevSol;
                _updateCellCache(_prevSolCache, *_prevSol, true);
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
                if (!_curSolDeflected) {
                    _curSolDeflected = true;

                    _curSolOrigValue = (*_curSol)[node][component];
                    _curSolOrigVarData = _curSolCache.atSCV[node];
                }
                
                int globalIdx = ParentType::_problem.nodeIndex(ParentType::_curCell(), 
                                                               node);
                int state = _staticNodeDat[globalIdx].phaseState;

                (*_curSol)[node][component] = value;
                _partialCellCacheUpdate(_curSolCache,
                                        ParentType::_problem.cellIndex(ParentType::_curCell()),
                                        *_curSol,
                                        node,
                                        globalIdx,
                                        state); 

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
                _curSolDeflected = false;
                (*_curSol)[node][component] = _curSolOrigValue;
                _curSolCache.atSCV[node] = _curSolOrigVarData;
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
//                const LocalFunction &sol   = usePrevSol ? *_prevSol     : *_curSol;
                const CellCache &cellCache = usePrevSol ? _prevSolCache : _curSolCache;
                
                int globalIdx = this->_problem.nodeIndex(this->_curCell(), scvId);
                
                Scalar satN = cellCache.atSCV[scvId].satN;
                Scalar satW = cellCache.atSCV[scvId].satW;
                
                Scalar porosity = _staticNodeDat[globalIdx].porosity;
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
                int i = this->_curCellGeom.subContVolFace[faceId].i;
                int j = this->_curCellGeom.subContVolFace[faceId].j;

                // normal vector, value of the area of the scvf
                const UnknownsVector &normal(this->_curCellGeom.subContVolFace[faceId].normal);

                // get global coordinates of nodes i,j
                const LocalCoord &global_i = this->_curCellGeom.subContVol[i].global;
                const LocalCoord &global_j = this->_curCellGeom.subContVol[j].global;
                const LocalCoord &local_i = this->_curCellGeom.subContVol[i].local;
                const LocalCoord &local_j = this->_curCellGeom.subContVol[j].local;

                LocalCoord      pGrad[PrimaryVariables];
                LocalCoord      xGrad[PrimaryVariables];
                for (int k = 0; k < PrimaryVariables; ++k) {
                    pGrad[k] = 0;
                    xGrad[k] = 0;
                }

                UnknownsVector tmp(0.);
                UnknownsVector pressure(0.0), massfrac(0.0);

                // calculate harmonic mean of permeabilities of nodes i and j
                Tensor K         = this->_problem.soil().K(global_i, ParentType::_curCell(), local_i);
                const Tensor &Kj = this->_problem.soil().K(global_j, ParentType::_curCell(), local_j);
                _harmonicMeanK(K, Kj);

                // calculate FE gradient (grad p for each phase)
                for (int k = 0; k < this->_curCellGeom.nNodes; k++) // loop over adjacent nodes
                {
                    // FEGradient at node k
                    const LocalCoord &feGrad = this->_curCellGeom.subContVolFace[faceId].grad[k];

                    pressure[WPhaseIndex] = _curSolCache.atSCV[k].pW;
                    pressure[NPhaseIndex] = _curSolCache.atSCV[k].pN;

                    // compute sum of pressure gradients for each phase
                    for (int phase = 0; phase < PrimaryVariables; phase++)
                    {
                        tmp = feGrad;
                        tmp *= pressure[phase];
                        pGrad[phase] += tmp;
                    }
                    // for diffusion of air in wetting phase
                    tmp = feGrad;
                    tmp *= _curSolCache.atSCV[k].massfrac[NCompIndex][WPhaseIndex];
                    xGrad[WPhaseIndex] += tmp;

                    // for diffusion of water in nonwetting phase
                    tmp = feGrad;
                    tmp *= _curSolCache.atSCV[k].massfrac[WCompIndex][NPhaseIndex];
                    xGrad[NPhaseIndex] += tmp;
                }

                // deduce gravity*density of each phase
                UnknownsVector contribComp[NumPhases];
                for (int phase=0; phase < NumPhases; phase++)
                {
                    contribComp[phase] = this->_problem.gravity();
                    contribComp[phase] *= _curSolCache.atSCV[i].density[phase];
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
                flux[WCompIndex] =   (alpha* _curSolCache.atSCV[up_w].density[WPhaseIndex]*_curSolCache.atSCV[up_w].mobility[WPhaseIndex]
                                 * _curSolCache.atSCV[up_w].massfrac[WCompIndex][WPhaseIndex]
                                 + (1-alpha)* _curSolCache.atSCV[dn_w].density[WPhaseIndex]*_curSolCache.atSCV[dn_w].mobility[WPhaseIndex]
                                 * _curSolCache.atSCV[dn_w].massfrac[WCompIndex][WPhaseIndex])
                    * outward[WPhaseIndex];
                flux[WCompIndex] +=  (alpha* _curSolCache.atSCV[up_n].density[NPhaseIndex]*_curSolCache.atSCV[up_n].mobility[NPhaseIndex]
                                 * _curSolCache.atSCV[up_n].massfrac[WCompIndex][NPhaseIndex]
                                 + (1-alpha)* _curSolCache.atSCV[dn_n].density[NPhaseIndex]*_curSolCache.atSCV[dn_n].mobility[NPhaseIndex]
                                 * _curSolCache.atSCV[dn_n].massfrac[WCompIndex][NPhaseIndex])
                    * outward[NPhaseIndex];

                // air conservation
                flux[NCompIndex]   =   (alpha* _curSolCache.atSCV[up_n].density[NPhaseIndex]*_curSolCache.atSCV[up_n].mobility[NPhaseIndex]
                                 * _curSolCache.atSCV[up_n].massfrac[NCompIndex][NPhaseIndex]
                                 + (1-alpha)* _curSolCache.atSCV[dn_n].density[NPhaseIndex]*_curSolCache.atSCV[dn_n].mobility[NPhaseIndex]
                                 * _curSolCache.atSCV[dn_n].massfrac[NCompIndex][NPhaseIndex])
                    * outward[NPhaseIndex];
                flux[NCompIndex]  +=   (alpha* _curSolCache.atSCV[up_w].density[WPhaseIndex]*_curSolCache.atSCV[up_w].mobility[WPhaseIndex]
                                 * _curSolCache.atSCV[up_w].massfrac[NCompIndex][WPhaseIndex]
                                 + (1-alpha)* _curSolCache.atSCV[dn_w].density[WPhaseIndex]*_curSolCache.atSCV[dn_w].mobility[WPhaseIndex]
                                 * _curSolCache.atSCV[dn_w].massfrac[NCompIndex][WPhaseIndex])
                    * outward[WPhaseIndex];

                // DIFFUSION
                UnknownsVector normDiffGrad;

                // get local to global id map
                int state_i = _curSolCache.atSCV[i].phasestate;
                int state_j = _curSolCache.atSCV[j].phasestate;

                Scalar diffusionWW(0.0), diffusionWN(0.0); // diffusion of water
                Scalar diffusionAW(0.0), diffusionAN(0.0); // diffusion of air
                UnknownsVector avgDensity, avgDpm;

                avgDpm[WPhaseIndex]=2e-9; // TODO: needs to be changed !!!
                avgDpm[NPhaseIndex]=2.25e-5; // water in the gasphase

                normDiffGrad[WPhaseIndex] = xGrad[WPhaseIndex]*normal;
                normDiffGrad[NPhaseIndex] = xGrad[NPhaseIndex]*normal;

                // calculate the arithmetic mean of densities
                avgDensity[WPhaseIndex] = 0.5*(_curSolCache.atSCV[i].density[WPhaseIndex] + _curSolCache.atSCV[j].density[WPhaseIndex]);
                avgDensity[NPhaseIndex] = 0.5*(_curSolCache.atSCV[i].density[NPhaseIndex] + _curSolCache.atSCV[j].density[NPhaseIndex]);

                if (state_i==2 && state_j==2)
                {
                    diffusionAW = avgDpm[WPhaseIndex] * avgDensity[WPhaseIndex] * normDiffGrad[WPhaseIndex];
                    diffusionWW = - diffusionAW;
                    diffusionWN = avgDpm[NPhaseIndex] * avgDensity[NPhaseIndex] * normDiffGrad[NPhaseIndex];
                    diffusionAN = - diffusionWN;
                }
                else if ((state_i == 1 || state_j == 1) || (state_i == 1 && state_j == 1))
                {
                    diffusionAW = avgDpm[WPhaseIndex] * avgDensity[WPhaseIndex] * normDiffGrad[WPhaseIndex];
                    diffusionWW = - diffusionAW;
                }
                else if ((state_i == 0 || state_j == 0) || (state_i == 0 && state_j == 0))
                {
                    diffusionWN = avgDpm[NPhaseIndex] * avgDensity[NPhaseIndex] * normDiffGrad[NPhaseIndex];
                    diffusionAN = - diffusionWN;
                }

                // add diffusion of water to flux
                flux[WCompIndex] += (diffusionWW + diffusionWN);
                //	std::cout << "Water Flux: " << flux[WCompIndex] << std::endl;

                // add diffusion of air to flux
                flux[NCompIndex] += (diffusionAN + diffusionAW);
                // std::cout << "Air Flux: " << flux[NCompIndex] << std::endl;
            }
        
        void clearVisited ()
            {
                for (int i = 0; i < this->_problem.numNodes(); i++){
                    _staticNodeDat[i].visited = false;
//                    sNDat[i].switched = false;
                }
            }

        /*!
         * \brief Adapt the static data of the current cell to the current solution
         */
        void updateStaticData(const Cell& cell)
            {
                int n = cell.template count<GridDim>();
                for (int i = 0; i < n; ++i) {
                    int globalNodeIdx = this->_problem.nodeIndex(cell, i);
                    
                    // evaluate primary variable switch
                    _primaryVarSwitch(cell, i, globalNodeIdx);
                }
            }
        
        /*!
         * \brief Initialize the static data of a single node
         */
        void initStaticNodeData(const Cell& cell,
                                int globalNodeIdx,
                                const WorldCoord &global,
                                const LocalCoord &local)
            {
                // if nodes are not already visited
                if (_staticNodeDat[globalNodeIdx].visited)
                    return;

                // initialize phase state
                _staticNodeDat[globalNodeIdx].phaseState =
                    this->_problem.initialPhaseState(cell, global, local);
                
                _staticNodeDat[globalNodeIdx].oldPhaseState =
                    this->_problem.initialPhaseState(cell, global, local);
                
                // ASSUME porosity defined at nodes
                _staticNodeDat[globalNodeIdx].porosity = this->_problem.soil().porosity(this->_curCellGeom.cellGlobal, 
                                                                                        cell,
                                                                                        this->_curCellGeom.cellLocal);
                
                // set counter for variable switch to zero
                _staticNodeDat[globalNodeIdx].switched = 0;
                
                // mark elements that were already visited
                _staticNodeDat[globalNodeIdx].visited = true;
            }

        /*!
         * \brief Return true if the primary variables were switched 
         *        after the last timestep.
         */
        bool checkSwitched()
            { bool tmp = _switchFlag; _switchFlag = false; return tmp; }

    private:
        /*!
         * \brief Pre-compute the cell cache data.
         *
         * This method is called by BoxJacobian (which in turn is
         * called by the operator assembler) every time the current
         * cell changes.
         */
        void _updateCellCache(CellCache &dest, const LocalFunction &sol, bool oldSol)
            {
                int phaseState;
                int cellIndex   = ParentType::_problem.cellIndex(ParentType::_curCell());
                for (int i = 0; i < ParentType::_curCellGeom.nNodes; i++) {
                    int iGlobal = ParentType::_problem.nodeIndex(ParentType::_curCell(), i);
                    phaseState = oldSol?_staticNodeDat[iGlobal].oldPhaseState:_staticNodeDat[iGlobal].phaseState;
                    _partialCellCacheUpdate(dest,
                                            cellIndex,
                                            sol,
                                            i,        // index of sub volume to update,
                                            iGlobal,  // global node index of the sub volume's grid node
                                            phaseState);
                }
            }

        
        void _partialCellCacheUpdate(CellCache           &dest,
                                     int                  cellIndex,
                                     const LocalFunction &sol,
                                     int                  i, // index of the subvolume/grid node
                                     int                  iGlobal,
                                     int                  state) // global index of the sub-volume's node
            {
//                const int globalIdx = this->_problem.nodeIndex(cell, i);
                const WorldCoord &global = this->_curCellGeom.subContVol[i].global;
                const WorldCoord &local  = this->_curCellGeom.subContVol[i].local;
                                
                dest.atSCV[i].pW = sol[i][PwIndex];
                if (state == BothPhases) dest.atSCV[i].satN = sol[i][SwitchIndex];
                if (state == WPhaseOnly) dest.atSCV[i].satN = 0.0;
                if (state == NPhaseOnly) dest.atSCV[i].satN = 1.0;

                dest.atSCV[i].satW = 1.0 - dest.atSCV[i].satN;

                dest.atSCV[i].pC = this->_problem.materialLaw().pC(dest.atSCV[i].satW, 
                                                                   global,
                                                                   this->_curCell(), 
                                                                   local);
                dest.atSCV[i].pN = dest.atSCV[i].pW + dest.atSCV[i].pC;
                dest.atSCV[i].temperature = 283.15; // in [K], constant

                // Solubilities of components in phases
                if (state == BothPhases) {
                    dest.atSCV[i].massfrac[NCompIndex][WPhaseIndex] = this->_problem.multicomp().xAW(dest.atSCV[i].pN, dest.atSCV[i].temperature);
                    dest.atSCV[i].massfrac[WCompIndex][NPhaseIndex] = this->_problem.multicomp().xWN(dest.atSCV[i].pN, dest.atSCV[i].temperature);
                }
                if (state == WPhaseOnly) {
                    dest.atSCV[i].massfrac[WCompIndex][NPhaseIndex] = 0.0;
                    dest.atSCV[i].massfrac[NCompIndex][WPhaseIndex] =  sol[i][SwitchIndex];
                }
                if (state == NPhaseOnly){
                    dest.atSCV[i].massfrac[WCompIndex][NPhaseIndex] = sol[i][SwitchIndex];
                    dest.atSCV[i].massfrac[NCompIndex][WPhaseIndex] = 0.0;
                }
                dest.atSCV[i].massfrac[WCompIndex][WPhaseIndex] = 1.0 - dest.atSCV[i].massfrac[NCompIndex][WPhaseIndex];
                dest.atSCV[i].massfrac[NCompIndex][NPhaseIndex] = 1.0 - dest.atSCV[i].massfrac[WCompIndex][NPhaseIndex];
                dest.atSCV[i].phasestate = state;

                // Mobilities & densities
                dest.atSCV[i].mobility[WPhaseIndex] = this->_problem.materialLaw().mobW(dest.atSCV[i].satW, global, this->_curCell(), local, dest.atSCV[i].temperature, dest.atSCV[i].pW);
                dest.atSCV[i].mobility[NPhaseIndex] = this->_problem.materialLaw().mobN(dest.atSCV[i].satN, global, this->_curCell(), local, dest.atSCV[i].temperature, dest.atSCV[i].pN);
                // Density of Water is set constant here!
                dest.atSCV[i].density[WPhaseIndex] = 1000;//this->_problem.wettingPhase().density(dest.atSCV[i].temperature, dest.atSCV[i].pN);
                dest.atSCV[i].density[NPhaseIndex] = this->_problem.nonwettingPhase().density(dest.atSCV[i].temperature, 
                                                                            dest.atSCV[i].pN,
                                                                            dest.atSCV[i].massfrac[WCompIndex][NPhaseIndex]);
                
                // CONSTANT solubility (for comparison with twophase)
                // dest.atSCV[i].massfrac[NCompIndex][WPhaseIndex] = 0.0; dest.atSCV[i].massfrac[WCompIndex][WPhaseIndex] = 1.0;
                // dest.atSCV[i].massfrac[WCompIndex][NPhaseIndex] = 0.0; dest.atSCV[i].massfrac[NCompIndex][NPhaseIndex] = 1.0;

                //std::cout << "water in gasphase: " << dest.atSCV[i].massfrac[WCompIndex][NPhaseIndex] << std::endl;
                //std::cout << "air in waterphase: " << dest.atSCV[i].massfrac[NCompIndex][WPhaseIndex] << std::endl;

/*
                // for output
                (*outPressureN)[globalIdx] = dest.atSCV[i].pN;
                (*outCapillaryP)[globalIdx] = dest.atSCV[i].pC;
                (*outSaturationW)[globalIdx] = dest.atSCV[i].satW;
                (*outSaturationN)[globalIdx] = dest.atSCV[i].satN;
                (*outMassFracAir)[globalIdx] = dest.atSCV[i].massfrac[NCompIndex][WPhaseIndex];
                (*outMassFracWater)[globalIdx] = dest.atSCV[i].massfrac[WCompIndex][NPhaseIndex];
                (*outDensityW)[globalIdx] = dest.atSCV[i].density[WPhaseIndex];
                (*outDensityN)[globalIdx] = dest.atSCV[i].density[NPhaseIndex];
                (*outMobilityW)[globalIdx] = dest.atSCV[i].mobility[WPhaseIndex];
                (*outMobilityN)[globalIdx] = dest.atSCV[i].mobility[NPhaseIndex];
                (*outPhaseState)[globalIdx] = dest.atSCV[i].phasestate;
*/
            }

  	/** @brief perform variable switch
  	 *  @param global global node id
	 *  @param sol solution vector
	 *  @param local local node id
	 */
   	void _primaryVarSwitch (const Cell& cell, int localIdx, int globalIdx)
            {
                bool switched = false;
                int switch_counter = _staticNodeDat[globalIdx].switched;
                int state = _staticNodeDat[globalIdx].phaseState;

                const LocalFunction &sol = *this->_curSol;
                Scalar pW = sol[localIdx][PwIndex];
                // if (pW < 0.) pW = (-1)*pW;

                Scalar satW = 0.0;
                if (state == BothPhases) satW = 1.0-sol[localIdx][SwitchIndex];
  		if (state == WPhaseOnly) satW = 1.0;
  		if (state == NPhaseOnly) satW = 0.0;

  		const WorldCoord &global = this->_curCellGeom.subContVol[localIdx].global;
  		const LocalCoord  &local = this->_curCellGeom.subContVol[localIdx].local;

                Scalar pC = this->_problem.materialLaw().pC(satW, global, cell, local);
  		Scalar pN = pW + pC;

                switch(state)
                {
                    case NPhaseOnly:
                        Scalar xWNmass, xWNmolar, pwn, pWSat;
                        xWNmass = sol[localIdx][SwitchIndex];
                        xWNmolar = this->_problem.multicomp().convertMassToMoleFraction(xWNmass, NPhaseOnly);
                        pwn = xWNmolar * pN;
                        pWSat = this->_problem.multicomp().vaporPressure(_curSolCache.atSCV[localIdx].temperature);

                        if (pwn > 1.01*pWSat && switched == false)// && switch_counter < 3)
                        {
                            // appearance of water phase
                            std::cout << "Water appears at node " << globalIdx << "  Coordinates: " << global << std::endl;
                            _staticNodeDat[globalIdx].phaseState = BothPhases;
                            sol[localIdx][SwitchIndex] = 1.0 - 1e-6; // initialize solution vector
                            _staticNodeDat[globalIdx].switched += 1;
                            switched = true;
                        }
                        break;

                    case WPhaseOnly:
                        Scalar pbub, xAWmass, xAWmolar, henryInv;
                        xAWmass = sol[localIdx][SwitchIndex];
                        xAWmolar = this->_problem.multicomp().convertMassToMoleFraction(xAWmass, WPhaseIndex);
                        henryInv = this->_problem.multicomp().henry(_curSolCache.atSCV[localIdx].temperature);
                        pWSat = this->_problem.multicomp().vaporPressure(_curSolCache.atSCV[localIdx].temperature);
                        pbub = pWSat + xAWmolar/henryInv; // pWSat + pAW

                        if (pN < pbub && switched == false)// && switch_counter < 3)
                        {
                            // appearance of gas phase
                            std::cout << "Gas appears at node " << globalIdx << ",  Coordinates: " << global << std::endl;
                            _staticNodeDat[globalIdx].phaseState = BothPhases;
                            sol[localIdx][SwitchIndex] = 1e-5; // initialize solution vector
                            _staticNodeDat[globalIdx].switched += 1;
                            switched = true;
                        }
                        break;

                    case BothPhases:
                        Scalar satN = sol[localIdx][SwitchIndex];

                        if (satN < 0.0  && switched == false)// && switch_counter < 3)
                        {
                            // disappearance of gas phase
                            std::cout << "Gas disappears at node " << globalIdx << "  Coordinates: " << global << std::endl;
                            _staticNodeDat[globalIdx].phaseState = WPhaseOnly;
                            sol[localIdx][SwitchIndex] = this->_problem.multicomp().xAW(pN); // initialize solution vector
                            _staticNodeDat[globalIdx].switched += 1;
                            switched = true;
                        }
                        else if (satW < 0.0  && switched == false)// && switch_counter < 3)
                        {
                            // disappearance of water phase
                            std::cout << "Water disappears at node " << globalIdx << "  Coordinates: " << global << std::endl;
                            _staticNodeDat[globalIdx].phaseState = NPhaseOnly;
                            sol[localIdx][SwitchIndex] = this->_problem.multicomp().xWN(pN); // initialize solution vector
                            _staticNodeDat[globalIdx].switched += 1;
                            switched = true;
                        }
                        break;

                }
                
                if (switched){
                    updateVariableData(cell, sol, localIdx, _curSolCache, _staticNodeDat[globalIdx].phaseState);
                    _switchFlag = true;
                }
            }

        // harmonic mean of the permeability computed directly.  the
        // first parameter is used to store the result.
        static void _harmonicMeanK(Tensor &Ki, const Tensor &Kj)
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
        CWaterAir                   _multicomp;
        std::vector<StaticNodeData> _staticNodeDat;
	bool                        _switchFlag;

        // the porosity of the current cell
        Scalar           _curCellPorosity;

        // current solution
        LocalFunction   *_curSol;
        CellCache        _curSolCache;

        // needed for restoreCurSolution()
        bool             _curSolDeflected;
        Scalar           _curSolOrigValue;
        VariableNodeData _curSolOrigVarData;

        // previous solution
        LocalFunction  *_prevSol;
        CellCache       _prevSolCache;
    };
    

    ///////////////////////////////////////////////////////////////////////////
    // TwoPTwoCBoxModel (The actual numerical model.)
    ///////////////////////////////////////////////////////////////////////////
    /*!
     * \brief Adaption of the BOX scheme to the 2P-2C twophase flow model.
     */
    template<class ProblemT>
    class TwoPTwoCBoxModel : public BoxScheme< // The implementation of the model
                                           TwoPTwoCBoxModel<ProblemT>,
        
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
        typedef typename DomTraits::CellReferenceElements CellReferenceElements;
        typedef typename DomTraits::LocalCoord            LocalCoord;
        typedef typename DomTraits::WorldCoord            WorldCoord;

        enum {
            GridDim          = DomTraits::GridDim,
            WorldDim         = DomTraits::WorldDim
        };
        
    public:
        typedef NewNewtonMethod<ThisType> NewtonMethod;

        TwoPTwoCBoxModel(ProblemT &prob)
            : ParentType(prob, _twoPTwoCLocalJacobian),
              _twoPTwoCLocalJacobian(prob, false)
            {
                Api::require<Api::BasicDomainTraits, typename ProblemT::DomainTraits>();
            }

        void initial() {
            // TODO: make sure the calls to clearVisited and updateStaticData are correct!
            this->localJacobian().clearVisited();

            // iterate through leaf grid an evaluate c0 at cell center
            const CellIterator &endit = this->_problem.cellEnd();

            for (CellIterator it = this->_problem.cellBegin(); 
                 it!= endit;
                 ++it)
            {
                this->_localJacobian.setCurrentCell(*it);

                int nNodes = it->template count<GridDim>();
                for (int i = 0; i < nNodes; i++)
                {
                    // get cell center in reference element
                    
                    const LocalCoord &local = CellReferenceElements::general(it->type()).position(i, GridDim);
                    
                    // get global coordinate of cell center
                    WorldCoord global = it->geometry().global(local);

                    int globalId = this->_problem.nodeIndex(*it, i);
                    
                    // initialize cell concentration
                    this->_problem.initial((*(this->currentSolution()))[globalId],
                                           *it,
                                           global,
                                           local);
                    
                    this->localJacobian().initStaticNodeData(*it, globalId, global, local);
                }
            }

/*
            // set Dirichlet boundary conditions
            for (Iterator it = gridview.template begin<0>();
                 it != eendit; ++it)
            {
                // get geometry type
                Dune::GeometryType gt = it->geometry().type();

                // get entity
                const Entity& entity = *it;

                const typename Dune::LagrangeShapeFunctionSetContainer<DT,RT,dim>::value_type
                    &sfs=Dune::LagrangeShapeFunctions<DT, RT, dim>::general(gt,1);
                int size = sfs.size();

                // set type of boundary conditions
                this->localJacobian.template assembleBC<LeafTag>(entity);

                IntersectionIterator
                    endit = IntersectionIteratorGetter<G, LeafTag>::end(entity);

                for (IntersectionIterator is = IntersectionIteratorGetter<G,
                         LeafTag>::begin(entity); is!=endit; ++is)
                    if (is->boundary())
                    {
                        for (int i = 0; i < size; i++)
                            // handle subentities of this face
                            for (int j = 0; j < ReferenceElements<DT,dim>::general(gt).size(is->numberInSelf(), 1, sfs[i].codim()); j++)
                                if (sfs[i].entity()
                                    == ReferenceElements<DT,dim>::general(gt).subEntity(is->numberInSelf(), 1,
											j, sfs[i].codim()))
                                {
                                    for (int equationNumber = 0; equationNumber<m; equationNumber++)
                                    {
                                        if (this->localJacobian.bc(i)[equationNumber]
                                            == BoundaryConditions::dirichlet)
                                        {
                                            // get cell center in reference element
                                            Dune::FieldVector<DT,dim>
                                                local = sfs[i].position();

                                            // get global coordinate of cell center
                                            Dune::FieldVector<DT,dimworld>
                                                global = it->geometry().global(local);

                                            int globalId = this->vertexmapper.template map<dim>(entity, sfs[i].entity());
                                            FieldVector<int,m> dirichletIndex;
                                            FieldVector<BoundaryConditions::Flags, m>
                                                bctype = this->problem.bctype(global, entity, is, local);
                                            this->problem.dirichletIndex(global, entity, is,
                                                                         local, dirichletIndex);

                                            if (bctype[equationNumber] == BoundaryConditions::dirichlet)
                                            {
                                                FieldVector<RT,m>
                                                    ghelp = this->problem.g(
                                                        global, entity, is,
                                                        local);
                                                (*(this->u))[globalId][dirichletIndex[equationNumber]]
                                                    = ghelp[dirichletIndex[equationNumber]];
                                            }
                                        }
                                    }
                                }
                    }
		this->localJacobian.setLocalSolution(entity);
		for (int i = 0; i < size; i++)
                    this->localJacobian.updateVariableData(entity, this->localJacobian.u, i, false);

            }
*/
            *(this->previousSolution()) = *(this->currentSolution());
	}

    private:
        // calculates the jacobian matrix at a given position
        TwoPTwoCLocalJacobian  _twoPTwoCLocalJacobian;
    };
}

#endif
