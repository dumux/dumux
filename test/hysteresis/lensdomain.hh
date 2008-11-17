/*****************************************************************************
 *   Copyright (C) 2008 by Bernd Flemisch, Andreas Lauser                    *
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
/*!
 * \file
 * \brief Manages the domain for the lens problem. The domain consists of a
 *        grid and the states for the given entities.
 */
#ifndef DUMUX_LENS_DOMAIN_HH
#define DUMUX_LENS_DOMAIN_HH

#include "lensstate.hh"

#include <dumux/auxiliary/basicdomain.hh>

#include <dumux/new_models/pwsn/pwsnboxmodel.hh>

#include <dumux/new_material/regularizedvangenuchten.hh>
#include <dumux/new_material/parkerlenhard.hh>

#include <dumux/material/properties.hh>

#include <dune/grid/sgrid.hh>

namespace Dune
{
namespace Lens
{
    typedef Dune::SGrid<2, 2> LensGrid;
//      typedef Dune::YaspGrid<2, 2> LensGrid;
//      typedef Dune::UGGrid<2> LensGrid;

    /*!
     * \brief The domain for the Pw-Sn lens problem
     */
    template<class ScalarT>
    class PwSnLensDomain : public BasicDomain<LensGrid, ScalarT>
    {
        typedef PwSnLensDomain<ScalarT>                               ThisType;
        typedef BasicDomain<LensGrid, ScalarT>                        ParentType;
        typedef ScalarT                                               Scalar;
        
    public:
        // traits for the domain
        typedef typename ParentType::DomainTraits DomainTraits;

        // traits of the material relations
        struct MaterialTraits
        {
            // The basic capillary pressure model used for the
            // hystersis model.
            typedef RegularizedVanGenuchtenState<ScalarT>      VanGenuchtenState;
            typedef RegularizedVanGenuchten<VanGenuchtenState> VanGenuchten;

            // States for globally constant parameters, parameters
            // specfic for a medium, and cell specific parameters
            typedef LensGlobalState<ScalarT>                   GlobalState;
            typedef LensMediumState<VanGenuchtenState>         MediumState;
            typedef LensNodeState<MediumState>                 NodeState;
            typedef LensCellState<MediumState>                 CellState;
            
            // The parker-lenhard hysteresis model
#if !USE_NODE_PARAMETERS
            typedef Dune::ParkerLenhard<CellState, VanGenuchten, USE_SPLINES>     ParkerLenhard;
#else // USE_NODE_PARAMETERS
            typedef Dune::ParkerLenhard<NodeState, VanGenuchten, USE_SPLINES>   ParkerLenhard;
#endif
        };


    private:

        typedef typename DomainTraits::Grid                  Grid;
        typedef typename DomainTraits::Cell                  Cell;
        typedef typename DomainTraits::CellReferenceElement  CellReferenceElement;
        typedef typename DomainTraits::CellReferenceElements CellReferenceElements;

        typedef typename DomainTraits::Node                  Node;

        typedef typename DomainTraits::CellIterator         CellIterator;
        typedef typename DomainTraits::NodeIterator       NodeIterator;

        typedef typename DomainTraits::IntersectionIterator       IntersectionIterator;
        typedef typename DomainTraits::IntersectionIteratorGetter IntersectionIteratorGetter;

        typedef typename DomainTraits::LocalCoord    LocalCoord;
        typedef typename DomainTraits::WorldCoord    WorldCoord;
        typedef typename DomainTraits::Vector        Vector;
        typedef typename DomainTraits::Matrix        Matrix;

        typedef typename MaterialTraits::GlobalState   GlobalState;
        typedef typename MaterialTraits::MediumState   MediumState;
        typedef typename MaterialTraits::CellState     CellState;
        typedef typename MaterialTraits::NodeState   NodeState;

        typedef typename MaterialTraits::VanGenuchtenState VanGenuchtenState;
        typedef typename MaterialTraits::VanGenuchten      VanGenuchten;
        typedef typename MaterialTraits::ParkerLenhard     ParkerLenhard;

        enum {
            GridDim = DomainTraits::GridDim
        };

    public:
        PwSnLensDomain()
            {
                Api::require<Api::BasicDomainTraits, DomainTraits>();

                initGrid_();

                initGlobalState_();
                initMediaStates_();
                initCellStates_();
                initNodeStates_();
            };

        ~PwSnLensDomain()
            {
                delete outerMedium_;
                delete lensMedium_;
                delete globalState_;
            }

        //! returns the body force vector within a cell
        const Vector &gravity(const Cell &cell,
                              int cellIdx) const
            {
                return globalState_->gravity();
            }

        const Vector &gravity() const
            {
                return globalState_->gravity();
            }

        //! return the permeability tensor for a cell
        void applyPermeabilityTensor(Vector &dest, const Cell &cell, const Vector &pressureGradient) const
            {
                cellState(cellIndex(cell)).permeability().umv(pressureGradient, dest);
            }

        //! Return the capillary pressure for a given node of a cell
        Scalar pC(const Cell &cell,
                  int cellIdx,
                  int localVertIdx,
                  int globalVertIdx,
                  Scalar Sw) const
            {
#if USE_NODE_PARAMETERS
                return ParkerLenhard::pC(nodeState(globalVertIdx), Sw);
#else // !USE_NODE_PARAMETERS
#if USE_INTERFACE_CONDITION
                const NodeState &vertState = nodeState(globalVertIdx);
                if (!vertState.isOnInterface())
                    return ParkerLenhard::pC(cellState(cellIdx),
                                             Sw);

                // find the minimum capilarry pressure around the
                // node if the node is on an interface
                Scalar minPc = 1e100;
                int nNeighbors = vertState.numNeighbourCells();
                for (int i = 0; i < nNeighbors; ++i) {
                    const CellState &ncState = cellState(vertState.neighbourCellIdx(i));
                    minPc = std::min(minPc,
                                     ParkerLenhard::pC(ncState, Sw));
                }

                return minPc;
#else // ! USE_INTERFACE_CONDITION
                return ParkerLenhard::pC(cellState(cellIdx),
                                         Sw);
#endif
#endif
            }

        // return the derivative of the capillary pressure regarding
        // the wetting phase saturation
        Scalar dpC_dSw(const Cell &cell,
                       int cellIdx,
                       int localVertIdx,
                       int globalVertIdx,
                       Scalar Sw) const
            {
                /* TODO
                const CellState &state = cellState(cellIdx);
                Scalar a = ParkerLenhard::pC(state, Sw + 5e-4);
                Scalar b = ParkerLenhard::pC(state, Sw - 5e-4);
                return (a - b) / 1e-3;
                */
                return 0;
            }

        // return the capillary pressure for a given cell
        Scalar porosity(const Cell &cell) const
            {
                return cellState(cell).porosity();
            }

        // return the density of the wetting phase
        Scalar densityW() const
            {
                return globalState_->densityW();
            }

        // return the density of the non-wetting phase
        Scalar densityN() const
            {
                return globalState_->densityN();
            }

        // return the density of a phase given by an index. (lower
        // indices means that the phase is more wetting)
        Scalar density(int phase) const
            {
                return (phase == 0)? densityW() : densityN();
            }

        // return the viscosity of the wetting phase
        Scalar viscosityW() const
            {
                return globalState_->viscosityW();
            }

        // return the viscosity of the non-wetting phase
        Scalar viscosityN() const
            {
                return globalState_->viscosityN();
            }

        // return the viscosity of a phase given by an index. (lower
        // indices means that the phase is more wetting)
        Scalar viscosity(int phase) const
            {
                return (phase == 0)? viscosityW() : viscosityN();
            }

        // return the mobility of the wetting phase at a node
        Scalar mobilityW(const Cell &cell,
                         int cellIdx,
                         int localVertIdx,
                         int globalVertIdx,
                         Scalar Sw) const
            {
#if USE_NODE_PARAMETERS
                return ParkerLenhard::krw(nodeState(globalVertIdx), Sw) / viscosityW();
#else // !USE_NODE_PARAMETERS
#if USE_INTERFACE_CONDITION
                const NodeState &vertState = nodeState(globalVertIdx);
                if (!vertState.isOnInterface())
                    return ParkerLenhard::krw(cellState(cellIdx), Sw) / viscosityW();

                // find the minimum wetting phase mobility around the
                // node if the node is on an interface
                Scalar minMob = 1e100;
                int nNeighbors = vertState.numNeighbourCells();
                for (int i = 0; i < nNeighbors; ++i) {
                    const CellState &ncState = cellState(vertState.neighbourCellIdx(i));
                    minMob = std::min(minMob,
                                      ParkerLenhard::krw(ncState, Sw) / viscosityW());
                };

                return minMob;
#else // ! USE_INTERFACE_CONDITION
                return ParkerLenhard::krw(cellState(cellIdx), Sw) / viscosityW();
#endif
#endif
            }

        // return the mobility of the non-wetting phase at a node
        Scalar mobilityN(const Cell &cell,
                         int cellIdx,
                         int localVertIdx,
                         int globalVertIdx,
                         Scalar Sn) const
            {
                // the capillary pressure models expect Sw as the
                // independend variable even for the non-wetting phase
                // mobility!
                Scalar Sw = 1 - Sn;

#if USE_NODE_PARAMETERS
                return ParkerLenhard::krn(nodeState(globalVertIdx), Sw) / viscosityN();
#else // !USE_NODE_PARAMETERS
#if USE_INTERFACE_CONDITION
                const NodeState &vertState = nodeState(globalVertIdx);
                if (!vertState.isOnInterface())
                    return ParkerLenhard::krn(cellState(cellIdx), Sw) / viscosityN();

                // find the minimum wetting phase mobility around the
                // node if the node is on an interface
                Scalar minMob = 1e100;
                int nNeighbors = vertState.numNeighbourCells();
                for (int i = 0; i < nNeighbors; ++i) {
                    const CellState &ncState = cellState(vertState.neighbourCellIdx(i));
                    minMob = std::min(minMob,
                                      ParkerLenhard::krn(ncState, Sw) / viscosityN());
                };

                return minMob;
#else // ! USE_INTERFACE_CONDITION
                return ParkerLenhard::krn(cellState(cellIdx), Sw) / viscosityN();
#endif
#endif
            }

        // given a cell index, return the corresponding cell state
        CellState &cellState(int index)
            { return cellStates_[index]; }
        const CellState &cellState(int index) const
            { return cellStates_[index]; }

        CellState &cellState(const Cell &cell)
            { return cellStates_[cellIndex(cell)]; }
        const CellState &cellState(const Cell &cell) const
            { return cellStates_[cellIndex(cell)]; }

        // given a global node index, return the corresponding state
        NodeState &nodeState(int index)
            { return nodeStates_[index]; }
        const NodeState &nodeState(int index) const
            { return nodeStates_[index]; }

        // given a cell and a local node index, return the corresponding state
        NodeState &nodeState(const Cell &cell, int i)
            { return nodeStates_[ParentType::nodeIndex(cell, i)]; }
        const NodeState &nodeState(const Cell &cell, int i) const
            { return nodeStates_[ParentType::nodeIndex(cell, i)]; }

        // given a node, return it's state object
        NodeState &nodeState(const Node &vert)
            { return nodeStates_[ParentType::nodeIndex(vert)]; }
        const NodeState &nodeState(const Node &vert) const
            { return nodeStates_[ParentType::nodeIndex(vert)]; }

        const WorldCoord &lowerLeft() const
            { return gridLowerLeft_; }
        const WorldCoord &upperRight() const
            { return gridUpperRight_; }

        Scalar width() const
            { return gridUpperRight_[0] - gridLowerLeft_[0]; }
        Scalar height() const
            { return gridUpperRight_[1] - gridLowerLeft_[1]; }

        bool onUpperBoundary(const WorldCoord &pos) const
            { return pos[1] > gridUpperRight_[1] - eps_; }
        bool onLowerBoundary(const WorldCoord &pos) const
            { return pos[1] < gridLowerLeft_[1] + eps_; }
        bool onLeftBoundary(const WorldCoord &pos) const
            { return pos[0] < gridLowerLeft_[0] + eps_; }
        bool onRightBoundary(const WorldCoord &pos) const
            { return pos[0] > gridUpperRight_[0] - eps_; }

        // returns true iff a world coordinate is within the fine
        // sand lens.
        bool isInLens(const WorldCoord &coord)
            {
                for (int i = 0; i < DomainTraits::WorldDim; ++i) {
                    if (lensLowerLeft_[i] > coord[i] ||
                        lensUpperRight_[i] < coord[i])
                    {
                        return false;
                    }
                }
                return true;
            }

/*
        // returns true iff the cell corrosponding to the state is
        // within the fine sand lens.
        bool isInLens(const NodeState &state)
            {
                return state.mediumState() == lensMedium_;
            }
*/

    private:
        void initGrid_()
            {
                // specify the position and size of the grid in world
                // coordinates
                gridLowerLeft_[0] = 0;
                gridLowerLeft_[1] = 0;
                gridUpperRight_[0] = 6;
                gridUpperRight_[1] = 4;

                // the epsilon constant
                eps_ = 1e-8 * width();

                // create the grid
                Dune::FieldVector<int, GridDim> cellRes;
#if USE_ORIG_PROB
                cellRes[0] = 48;
                cellRes[1] = 32;
#else
                cellRes[0] = CELLRES_X;
                cellRes[1] = CELLRES_Y;
#endif

                Grid *grid = new Grid(cellRes,
                                      gridLowerLeft_,
                                      gridUpperRight_);
                ParentType::setGrid(grid);

                // specify the location and the size of the fine sand
                // lens in world coordinates
                lensLowerLeft_[0] = width()*(1./6);
                lensLowerLeft_[1] = height()*(2./4);
                lensUpperRight_[0] = width()*(4./6);
                lensUpperRight_[1] = height()*(3./4);
                // make sure the global coordniates of the lens are
                // aligned to the grid's cell boundaries in order to
                // make sure the interface occures at the same
                // position for sized lenses for node and cell based
                // material parameter storage. (this is not strictly
                // necessary, but makes it easier to compare the
                // different approaches.)
                lensLowerLeft_[0] = width()/cellRes[0]
                                  * (int) (lensLowerLeft_[0]*cellRes[0]/width())
                                  + (0.5 - eps_)/cellRes[0];
                lensLowerLeft_[1] = height()/cellRes[1]
                                  * (int) (lensLowerLeft_[1]*cellRes[1]/height())
                                  + (0.5 - eps_)/cellRes[1];
                lensUpperRight_[0] = width()/cellRes[0]
                                   * (int) (lensUpperRight_[0]*cellRes[0]/width())
                                   + (-0.5 + eps_)/cellRes[0];
                lensUpperRight_[1] = height()/cellRes[1]
                                   * (int) (lensUpperRight_[1]*cellRes[1]/height())
                                   + (-0.5 + eps_)/cellRes[1];
            }

        void initGlobalState_()
            {
                // initialize the globally uniform parameters
                globalState_ = new GlobalState();

                Water water;
                DNAPL dnapl;
                globalState_->setDensityW(water.density());
                globalState_->setDensityN(dnapl.density());

                globalState_->setViscosityW(water.viscosity());
                globalState_->setViscosityN(dnapl.viscosity());
                Vector gravity(0); gravity[1] = -9.81;
                globalState_->setGravity(gravity);
            }

        void initMediaStates_()
            {
                // initializing the state of the lens' medium
                lensMedium_ =  new MediumState();
                lensMedium_->setGlobalState(globalState_);
                lensMedium_->setSwr(0.18);
#if USE_ORIG_PROB
                lensMedium_->setSnr(0.0);
#else
                lensMedium_->setSnr(0.25);
#endif
                lensMedium_->setPorosity(0.4);
                lensMedium_->setMdcParams(VanGenuchtenState(0.00045, 7.3));
#if USE_DIFFERENT_MAIN_CURVES
                lensMedium_->setMicParams(VanGenuchtenState(0.00045*2, 7.5));
#else
                lensMedium_->setMicParams(VanGenuchtenState(0.00045, 7.3));
#endif
                lensMedium_->micParams().setVgMaxPC(lensMedium_->mdcParams().vgMaxPC());
                lensMedium_->setPermeability(9.05e-13);

                // initializing the state of the outer medium
                outerMedium_ =  new MediumState();
                outerMedium_->setGlobalState(globalState_);
                outerMedium_->setSwr(0.05);
#if USE_ORIG_PROB
                outerMedium_->setSnr(0.0);
#else
                outerMedium_->setSnr(0.20);
#endif
                outerMedium_->setPorosity(0.4);
                outerMedium_->setMdcParams(VanGenuchtenState(0.0037, 4.7));
#if USE_DIFFERENT_MAIN_CURVES
                outerMedium_->setMicParams(VanGenuchtenState(0.0037*2, 4.8));
#else
                outerMedium_->setMicParams(VanGenuchtenState(0.0037, 4.7));
#endif

                outerMedium_->micParams().setVgMaxPC(outerMedium_->mdcParams().vgMaxPC());
                outerMedium_->setPermeability(4.6e-10);

#if 0
                lensMedium_->setSwr(outerMedium_->Swr());
                lensMedium_->setSnr(outerMedium_->Snr());
                lensMedium_->setPorosity(outerMedium_->porosity());
                lensMedium_->setMicParams(outerMedium_->micParams());
                lensMedium_->setMdcParams(outerMedium_->mdcParams());
                lensMedium_->setPermeability(outerMedium_->permeability());
#endif

            }

        void initCellStates_()
            {
                cellStates_.resize(ParentType::numCells());

                // initialize the cell state objects depending on
                // wether the cell is inside or outside of the lens of
                // fine sand.
                WorldCoord cellCenter;
                CellIterator it = ParentType::cellBegin();
                CellIterator endit = ParentType::cellEnd();
                for (; it != endit; ++it) {
                    ParentType::cellCenter(*it, cellCenter);
                    int arrayPos = ParentType::cellIndex(*it);
                    if (isInLens(cellCenter))
                        cellStates_[arrayPos].setMediumState(lensMedium_);
                    else
                        cellStates_[arrayPos].setMediumState(outerMedium_);
                }
            }

        void initNodeStates_()
            {
                nodeStates_.resize(ParentType::numNodes());
#if USE_NODE_PARAMETERS
                NodeIterator vertIt = ParentType::nodeBegin();
                const NodeIterator &endVertIt = ParentType::nodeEnd();
                for (; vertIt != endVertIt; ++vertIt) {
                    WorldCoord pos;
                    ParentType::nodePosition(pos, *vertIt);
                    NodeState &vertState = nodeState(*vertIt);
                    if (isInLens(pos))
                        vertState.setMediumState(lensMedium_);
                    else
                        vertState.setMediumState(outerMedium_);
                }

#else // !USE_NODE_PARAMETERS
                // initialize the node state objects

                // loop over all cells
                CellIterator cellIt = ParentType::cellBegin();
                const CellIterator &endCellIt = ParentType::cellEnd();
                for (; cellIt != endCellIt; ++cellIt) {
                    // loop over all faces of the current cell
                    IntersectionIterator faceIt = IntersectionIteratorGetter::begin(*cellIt);
                    const IntersectionIterator &endFaceIt = IntersectionIteratorGetter::end(*cellIt);
                    for (; faceIt != endFaceIt; ++ faceIt) {
                        // make sure the face not on the boundary of
                        // the grid
                        if (!faceIt.neighbor())
                            continue;

                        // check whether both cells share the same
                        // medium. if yes we don't need to do anything
                        // about the current face.
                        const CellState &inState = cellState(*faceIt.inside());
                        const CellState &outState = cellState(*faceIt.outside());
                        if (inState.mediumState() != outState.mediumState())
                            continue;

                        // alright, the face is on a inhomogenity, so
                        // we have to mark all its vertices
                        markInterfaceVertices_(faceIt);
                    }
                }
#endif
            }

#if !USE_NODE_PARAMETERS
        //  mark all vertices on a face as as belonging to an
        //  inhomogenity
        void markInterfaceVertices_(IntersectionIterator &interfaceIt)
            {
                const Cell &cell = *interfaceIt.inside();
                const CellReferenceElement &refElem =
                    CellReferenceElements::general(cell.geometry().type());

                int inIdx = cellIndex(*interfaceIt.inside());
                int outIdx = cellIndex(*interfaceIt.outside());

                int faceIdx = interfaceIt.numberInSelf();
                int nVerticesOfFace = refElem.size(faceIdx, 1, GridDim);
                for (int nodeInFace = 0;
                     nodeInFace < nVerticesOfFace;
                     nodeInFace++)
                {
                    int nodeIdxInElement = refElem.subEntity(faceIdx, 1, nodeInFace, GridDim);

                    nodeState(*interfaceIt.inside(), nodeIdxInElement).addNeighbourCellIdx(inIdx);
                    nodeState(*interfaceIt.inside(), nodeIdxInElement).addNeighbourCellIdx(outIdx);
                }
            }
#endif // !USE_NODE_PARAMETERS


        // the actual type which stores the cell states.
        typedef std::vector<NodeState>  NodeStateArray_;

        // the actual type which stores the cell states.
        typedef std::vector<CellState>    CellStateArray_;

        // the lower left and upper right coordinates of the complete
        // grid
        WorldCoord gridLowerLeft_;
        WorldCoord gridUpperRight_;

        // the lower left and upper right coordinates of the fine sand
        // lens
        WorldCoord lensLowerLeft_;
        WorldCoord lensUpperRight_;

        // global state and states of the media
        GlobalState *globalState_;
        MediumState *outerMedium_;
        MediumState *lensMedium_;

        // stores a state for each node
        NodeStateArray_ nodeStates_;

        // stores a state for each node
        CellStateArray_  cellStates_;

        // A small epsilon value and the current simulated
        // time. FIXME/TODO: should probably not be here..
        Scalar eps_;
    };
}
}


#endif
