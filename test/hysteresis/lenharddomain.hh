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
 * \brief Manages the domain for the lenhard problem. The domain consists of a
 *        grid and the states for the given entities.
 *
 * This has been modeled to match as closely as possible to the one
 * described at:
 *
 * Sheta, Hussam: "Simulation von Mehrphasenvorgaengen in poroesen
 *     Medien unter Einbeziehung von Hystereseeffekten", PhD theses,
 *     Braunschweig 1999, pp. 112
 */
#ifndef DUMUX_LENHARD_DOMAIN_HH
#define DUMUX_LENHARD_DOMAIN_HH

#include "lenhardstate.hh"

#include <dumux/auxiliary/basicdomain.hh>

#include <dumux/new_material/regularizedvangenuchten.hh>
#include <dumux/new_material/parkerlenhard.hh>

#include <dumux/material/properties.hh>

#include <dune/grid/sgrid.hh>
#include <dune/grid/onedgrid.hh>
//#include <dune/grid/yaspgrid.hh>
//#include <dumux/onedinndgrid/onedinndgrid.hh>
//#include <dumux/onedinndgrid/onedinndgrid.cc>

namespace Dune
{
namespace Lenhard
{
    typedef Dune::SGrid<1, 1> LenhardGrid;
//    typedef Dune::YaspGrid<1, 1> LenhardGrid;
//    typedef Dune::OneDGrid       LenhardGrid;

    /*!
     * \brief The domain for the Pw-Sn lenhard problem
     */
    template<class ScalarT>
    class PwSnLenhardDomain : public BasicDomain<LenhardGrid, ScalarT>
    {
        typedef PwSnLenhardDomain<ScalarT>                   ThisType;
        typedef BasicDomain<LenhardGrid, ScalarT>            ParentType;
        typedef ScalarT                                      Scalar;
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
            typedef LenhardGlobalState<ScalarT>                   GlobalState;
            typedef LenhardMediumState<VanGenuchtenState>         MediumState;
            typedef LenhardNodeState<MediumState>                 NodeState;
            typedef LenhardCellState<MediumState>                 CellState;
            
            // The parker-lenhard hysteresis model
#if !defined USE_NODE_PARAMETERS
            typedef Dune::ParkerLenhard<CellState, VanGenuchten>     ParkerLenhard;
#else // defined USE_NODE_PARAMETERS
            typedef Dune::ParkerLenhard<NodeState, VanGenuchten>   ParkerLenhard;
#endif
        };

    private:
        typedef typename DomainTraits::Grid                  Grid;
        typedef typename DomainTraits::Cell                  Cell;
        typedef typename DomainTraits::ReferenceElement      ReferenceElement;

        typedef typename DomainTraits::Node                Node;

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
        PwSnLenhardDomain()
            {
                Api::require<Api::BasicDomainTraits, DomainTraits>();

                initGrid_();

                initGlobalState_();
                initMediaStates_();
                initCellStates_();
                initNodeStates_();
            };

        ~PwSnLenhardDomain()
            {
                delete coarseSand_;
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
#if defined USE_NODE_PARAMETERS
                return ParkerLenhard::pC(nodeState(globalVertIdx), Sw);
#else // !defined USE_NODE_PARAMETERS
#if defined USE_INTERFACE_CONDITION
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
/*                printf("krw(%f): %f, viscosity: %f\n",
                       Sw,
                       ParkerLenhard::krw(cellState(cellIdx), Sw),
                       viscosityW());
*/

#if defined USE_NODE_PARAMETERS
                return ParkerLenhard::krw(nodeState(globalVertIdx), Sw) / viscosityW();
#else // !defined USE_NODE_PARAMETERS
#ifdef USE_INTERFACE_CONDITION
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
/*                printf("krn(%f): %f, viscosity: %f\n",
                  1 - Sn,
                  ParkerLenhard::krn(cellState(cellIdx), 1 - Sn),
                  viscosityW());
*/

                // the capillary pressure models expect Sw as the
                // independend variable even for the non-wetting phase
                // mobility!
                Scalar Sw = 1 - Sn;

#if defined USE_NODE_PARAMETERS
                return ParkerLenhard::krn(nodeState(globalVertIdx), Sw) / viscosityN();
#else // !defined USE_NODE_PARAMETERS
#ifdef USE_INTERFACE_CONDITION
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

        Scalar height() const
            { return gridUpperRight_[0] - gridLowerLeft_[0]; }

        bool onUpperBoundary(const WorldCoord &pos) const
            { return pos[0] > gridUpperRight_[0] - eps_; }
        bool onLowerBoundary(const WorldCoord &pos) const
            { return pos[0] < gridLowerLeft_[0] + eps_; }

    private:
        void initGrid_()
            {
                // specify the position and size of the grid in world
                // coordinates
                gridLowerLeft_[0] = 0;

                // column was 72 cm high, but we need some additional
                // space so that the boundary condition at the upper
                // boundary doesn't destroy the experiment
                gridUpperRight_[0] = 0.80;


                // the epsilon constant
                eps_ = 1e-8 * height();

                // create the grid
                Dune::FieldVector<int, GridDim> cellRes;
                cellRes[0] = 160;

                Grid *grid = new Grid(cellRes,
                                      gridLowerLeft_,
                                      gridUpperRight_);
                ParentType::setGrid(grid);
            }

        void initGlobalState_()
            {
                // initialize the globally uniform parameters
                globalState_ = new GlobalState();

                Water wettingFluid;
                Air   nonwettingFluid;
                globalState_->setDensityW(wettingFluid.density());
                globalState_->setDensityN(nonwettingFluid.density());

                globalState_->setViscosityW(wettingFluid.viscosity());
                globalState_->setViscosityN(nonwettingFluid.viscosity());

                Vector gravity(0);
                gravity[0] = -9.81;
                globalState_->setGravity(gravity);
            }

        void initMediaStates_()
            {
                // initializing the state of the porous medium
                coarseSand_ =  new MediumState();
                coarseSand_->setGlobalState(globalState_);
#if LENHARD_EXPERIMENT == 1
                coarseSand_->setSwr(0.0);
                coarseSand_->setSnr(0.25);
                coarseSand_->setPorosity(0.36);
                coarseSand_->setMicParams(VanGenuchtenState(0.00052, 4.63));
                coarseSand_->setMdcParams(VanGenuchtenState(0.00052, 4.63));
                coarseSand_->micParams().setVgMaxPC(coarseSand_->mdcParams().vgMaxPC());
                coarseSand_->setPermeability(3.37e-11);
#elif LENHARD_EXPERIMENT == 2
                coarseSand_->setSwr(0.17);
                coarseSand_->setSnr(0.25);
                coarseSand_->setPorosity(0.36);
                coarseSand_->setMdcParams(VanGenuchtenState(0.00042, 5.25));
//                coarseSand_->setMicParams(VanGenuchtenState(0.00042, 5.25));
                coarseSand_->setMicParams(VanGenuchtenState(0.00042*2, 5.25));
                coarseSand_->micParams().setVgMaxPC(coarseSand_->mdcParams().vgMaxPC());
                coarseSand_->setPermeability(3.37e-11);
#endif
            }

        void initCellStates_()
            {
                cellStates_.resize(ParentType::numCells());

                // initialize the cell state objects depending on
                // wether the cell is inside or outside of the lenhard of
                // fine sand.
                WorldCoord cellCenter;
                CellIterator it = ParentType::cellBegin();
                CellIterator endit = ParentType::cellEnd();
                for (; it != endit; ++it) {
                    ParentType::cellCenter(*it, cellCenter);
                    int arrayPos = ParentType::cellIndex(*it);
                    cellStates_[arrayPos].setMediumState(coarseSand_);
                }

            }

        void initNodeStates_()
            {
                nodeStates_.resize(ParentType::numNodes());

#if defined USE_NODE_PARAMETERS
                NodeIterator vertIt = ParentType::nodeBegin();
                const NodeIterator &endVertIt = ParentType::nodeEnd();
                for (; vertIt != endVertIt; ++vertIt) {
                    WorldCoord pos;
                    ParentType::nodePosition(pos, *vertIt);
                    NodeState &vertState = nodeState(*vertIt);

                    vertState.setMediumState(coarseSand_);
                }

#else // ! defined USE_NODE_PARAMETERS
                // initialize the node state objects
                /*
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
                */
#endif
            }

#if !defined USE_NODE_PARAMETERS
        //  mark all vertices on a face as as belonging to an
        //  inhomogenity
        /*
        void markInterfaceVertices_(IntersectionIterator &interfaceIt)
            {
                const Cell &cell = *interfaceIt.inside();
                const ReferenceElement &refElem =
                    DomainTraits::referenceElement(cell.geometry().type());

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
        */
#endif // !defined USE_NODE_PARAMETERS


        // the actual type which stores the cell states.
        typedef std::vector<NodeState>  NodeStateArray_;

        // the actual type which stores the cell states.
        typedef std::vector<CellState>    CellStateArray_;

        // the lower left and upper right coordinates of the complete
        // grid
        WorldCoord gridLowerLeft_;
        WorldCoord gridUpperRight_;

        // the lower left and upper right coordinates of the fine sand
        // lenhard
        WorldCoord lenhardLowerLeft_;
        WorldCoord lenhardUpperRight_;

        // global state and states of the media
        GlobalState *globalState_;
        MediumState *coarseSand_;

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
