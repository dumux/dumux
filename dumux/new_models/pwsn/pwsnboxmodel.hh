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
#ifndef DUMUX_PWSN_BOX_MODEL_HH
#define DUMUX_PWSN_BOX_MODEL_HH

#include <dumux/new_models/boxscheme/boxscheme.hh>
#include <dumux/new_models/boxscheme/p1boxtraits.hh>

#include <dumux/auxiliary/apis.hh>

namespace Dune
{
    ///////////////////////////////////////////////////////////////////////////
    // PwSn traits (central place for names and indices required by the 
    // PwSnBoxJacobian and PwSnBoxModel)
    ///////////////////////////////////////////////////////////////////////////
    /*!
     * \brief The pw-Sn specific traits.
     */
    class PwSnTraits
    {
    public:
        enum {
            numEq = 2 //!< Number of primary variables
        };
        enum {
            pWIndex = 0,  //!< Index for the wetting phase pressure in a field vector
            snIndex = 1   //!< Index for the non-wetting phase saturation in a field vector        
        };
    };

    ///////////////////////////////////////////////////////////////////////////
    // PwSnBoxJacobian (evaluate the local jacobian for the newton method.)
    ///////////////////////////////////////////////////////////////////////////
    /*!
     * \brief pw-Sn formulation specific details needed to
     *        approximately calculate the local jacobian in the BOX
     *        scheme.
     *
     * This class is used to fill the gaps in BoxJacobian for the Pw-Sn twophase flow.
     */
    template<class ProblemT, class BoxTraitsT, class PwSnTraitsT>
    class PwSnBoxJacobian : public BoxJacobian<ProblemT, BoxTraitsT, PwSnBoxJacobian<ProblemT, BoxTraitsT, PwSnTraitsT> >
    {
    private:
        typedef PwSnBoxJacobian<ProblemT, BoxTraitsT, PwSnTraitsT>  ThisType;
        typedef BoxJacobian<ProblemT, BoxTraitsT, ThisType >        ParentType;

        typedef ProblemT                                Problem;
        typedef typename Problem::DomainTraits          DomTraits;
        typedef BoxTraitsT                              BoxTraits;
        typedef PwSnTraitsT                             PwSnTraits;

        enum {
            GridDim          = DomTraits::GridDim,
            WorldDim         = DomTraits::WorldDim,
            numEq = BoxTraits::numEq,
            pWIndex          = PwSnTraits::pWIndex,
            snIndex          = PwSnTraits::snIndex
        };
        
        typedef typename DomTraits::Scalar              Scalar;
        typedef typename DomTraits::CoordScalar         CoordScalar;
        typedef typename DomTraits::Grid                Grid;
        typedef typename DomTraits::Cell                Cell;
        typedef typename Cell::EntityPointer            CellPointer;
        typedef typename DomTraits::LocalCoord          LocalCoord;

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
            Scalar Sw;
            Scalar pC;
            Scalar pN;
            
            UnknownsVector mobility;  //Vector with the number of phases
        };
        
        /*!
         * \brief Cached data for the each node of the cell.
         */
        struct CellCache
        {
            VariableNodeData  atSCV[BoxTraits::ShapeFunctionSetContainer::maxsize];
        };
        
    public:
        PwSnBoxJacobian(ProblemT &problem) 
            : ParentType(problem)
            {};

        /*!
         * \brief Set the current grid cell.
         */
        void setCurrentCell(const Cell &cell) 
            {
                if (ParentType::setCurrentCell_(cell)) {
                    curCellPorosity_ = ParentType::problem_.porosity(ParentType::curCell_());
                }
            };

        /*!
         * \brief Set the parameters for the calls to the remaining
         *        members.
         */
        void setParams(const Cell &cell, LocalFunction &curSol, LocalFunction &prevSol)
            {
                setCurrentCell(cell);
                
                curSol_ = &curSol;
                updateCellCache_(curSolCache_, *curSol_);
                curSolDeflected_ = false;
                
                prevSol_ = &prevSol;
                updateCellCache_(prevSolCache_, *prevSol_);
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

                    curSolOrigValue_ = (*curSol_)[node][component];
                    curSolOrigVarData_ = curSolCache_.atSCV[node];
                }
                
                (*curSol_)[node][component] = value;
                partialCellCacheUpdate_(curSolCache_,
                                        ParentType::problem_.cellIndex(ParentType::curCell_()),
                                        *curSol_,
                                        node,
                                        ParentType::problem_.nodeIndex(ParentType::curCell_(), 
                                                                         node)); 

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
                (*curSol_)[node][component] = curSolOrigValue_;
                curSolCache_.atSCV[node] = curSolOrigVarData_;
            };
        
        /*!
         * \brief Evaluate the rate of change of all conservation
         *        quantites (e.g. phase mass) within a sub control
         *        volume of a finite volume cell in the pw-Sn
         *        formulation.
         * 
         * This function should not include the source and sink terms.
         */
        void localRate(UnknownsVector &result, int scvId, bool usePrevSol) const
            {
                LocalFunction *sol = usePrevSol?prevSol_:curSol_;

                // partial time derivative of the wetting phase mass
                result[pWIndex] = -ParentType::problem_.densityW()
                                   * curCellPorosity_
                                   * (*sol)[scvId][snIndex];
                // partial time derivative of the non-wetting phase mass
                result[snIndex] = ParentType::problem_.densityN()
                                  * curCellPorosity_
                                  * (*sol)[scvId][snIndex];

            }

        
        /*!
         * \brief Evaluates the mass flux over a face of a subcontrol
         *        volume.
         */
        void fluxRate(UnknownsVector &flux, int faceId) const
            {
                Api::require<Api::BasicDomainTraits, typename ProblemT::DomainTraits>();
                assert(numEq == 2);

                const typename FVElementGeometry::SubControlVolumeFace
                    &face = ParentType::curCellGeom_.subContVolFace[faceId];
                const int i = face.i;
                const int j = face.j;

                LocalCoord Kij(0);
                
                // Kij = K*normal
                ParentType::problem_.applyPermeabilityTensor(Kij,
                                                             ParentType::curCell_(), 
                                                             face.normal);

                for (int phase = 0; phase < numEq; phase++) {
                    // calculate FE gradient
                    LocalCoord pGrad(0);
                    for (int k = 0; k < ParentType::curCellGeom_.nNodes; k++) {
                        LocalCoord grad(face.grad[k]);
                        if (phase == snIndex)
                            grad *= curSolCache_.atSCV[k].pN;
                        else
                            grad *= (*curSol_)[k][pWIndex];

                        pGrad += grad;
                    }

                    // adjust pressure gradient by gravity force
                    Scalar phaseDensity = ParentType::problem_.density(phase);
                    LocalCoord gravity = ParentType::problem_.gravity();
                    gravity *= phaseDensity;
                    pGrad   -= gravity;
                    
                    // calculate the flux using upwind
                    Scalar outward = pGrad*Kij;
                    if (outward < 0)
                        flux[phase] = phaseDensity*curSolCache_.atSCV[i].mobility[phase]*outward;
                    else
                        flux[phase] = phaseDensity*curSolCache_.atSCV[j].mobility[phase]*outward;
                }
            }

    private:
        /*!
         * \brief Pre-compute the cell cache data.
         *
         * This method is called by BoxJacobian (which in turn is
         * called by the operator assembler) every time the current
         * cell changes.
         */
        void updateCellCache_(CellCache &dest, const LocalFunction &sol)
            {
                assert(numEq == 2);

                int cellIndex   = ParentType::problem_.cellIndex(ParentType::curCell_());
                for (int i = 0; i < ParentType::curCellGeom_.nNodes; i++) {
                    int iGlobal = ParentType::problem_.nodeIndex(ParentType::curCell_(), i);
                    partialCellCacheUpdate_(dest,
                                            cellIndex,
                                            sol,
                                            i,        // index of sub volume to update,
                                            iGlobal); // global node index of the sub volume's grid node
                }
            }

        
        void partialCellCacheUpdate_(CellCache           &dest,
                                     int                  cellIndex,
                                     const LocalFunction &sol,
                                     int                  i, // index of the subvolume/grid node
                                     int                  iGlobal) // global index of the sub-volume's node
            {
                // Current cache at sub-controlvolume
                VariableNodeData &scvCache = dest.atSCV[i];
                // Current solution for sub-controlvolume
                const UnknownsVector &scvSol = sol[i];

                scvCache.Sw = 1.0 - scvSol[snIndex];
                scvCache.pC = ParentType::problem_.pC(ParentType::curCell_(),
                                                      cellIndex,
                                                      i,
                                                      iGlobal,
                                                      scvCache.Sw);
                scvCache.pN = scvSol[pWIndex] + scvCache.pC;
                scvCache.mobility[pWIndex] = ParentType::problem_.mobilityW(ParentType::curCell_(),
                                                                            cellIndex,
                                                                            i,
                                                                            iGlobal,
                                                                            scvCache.Sw);
                scvCache.mobility[snIndex] = ParentType::problem_.mobilityN(ParentType::curCell_(),
                                                                            cellIndex,
                                                                            i,
                                                                            iGlobal,
                                                                            scvSol[snIndex]);
            }

//        CellPointer     curCell_;
        Scalar          curCellPorosity_;
        Tensor         *curCellPermeability_;

        LocalFunction   *curSol_;
        CellCache        curSolCache_;

        bool             curSolDeflected_;
        Scalar           curSolOrigValue_;
        VariableNodeData curSolOrigVarData_;

        LocalFunction  *prevSol_;
        CellCache       prevSolCache_;
    };
    

    ///////////////////////////////////////////////////////////////////////////
    // PwSnBoxModel (The actual numerical model.)
    ///////////////////////////////////////////////////////////////////////////
    /*!
     * \brief Adaption of the BOX scheme to the pW-Sn twophase flow model.
     */
    template<class ProblemT>
    class PwSnBoxModel : public BoxScheme< // The implementation of the model
                                           PwSnBoxModel<ProblemT>,
        
                                           // The Traits for the BOX method
                                           P1BoxTraits<typename ProblemT::DomainTraits::Scalar,
                                                       typename ProblemT::DomainTraits::Grid,
                                                       PwSnTraits::numEq>,

                                           // The actual problem we would like to solve
                                           ProblemT, 
        
                                           // The local jacobian operator
                                           PwSnBoxJacobian<ProblemT, 
                                                           P1BoxTraits<typename ProblemT::DomainTraits::Scalar,
                                                                       typename ProblemT::DomainTraits::Grid,
                                                                       PwSnTraits::numEq>,
                                                           PwSnTraits> >
    {
        typedef typename ProblemT::DomainTraits::Grid   Grid;
        typedef typename ProblemT::DomainTraits::Scalar Scalar;
        typedef PwSnBoxModel<ProblemT>                  ThisType;
        
    public:
        typedef P1BoxTraits<Scalar, Grid, PwSnTraits::numEq> BoxTraits;
        typedef Dune::PwSnTraits                                   PwSnTraits;
        
    private:
        typedef PwSnBoxJacobian<ProblemT, BoxTraits, PwSnTraits>  PwSnLocalJacobian;
        typedef BoxScheme<ThisType,
                          BoxTraits,
                          ProblemT, 
                          PwSnLocalJacobian>  ParentType;
        
    public:
        typedef NewNewtonMethod<ThisType> NewtonMethod;

        PwSnBoxModel(ProblemT &prob)
            : ParentType(prob, pwSnLocalJacobian_),
              pwSnLocalJacobian_(prob)
            {
                Api::require<Api::BasicDomainTraits, typename ProblemT::DomainTraits>();
            }

    private:
        // calculates the jacobian matrix at a given position
        PwSnLocalJacobian  pwSnLocalJacobian_;
    };
}

#endif
