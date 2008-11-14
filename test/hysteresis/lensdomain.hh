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

#include <dumux/new_models/box/pwsn/pwsnboxmodel.hh>

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

                _initGrid();

                _initGlobalState();
                _initMediaStates();
                _initCellStates();
                _initNodeStates();
            };

        ~PwSnLensDomain()
            {
                delete _outerMedium;
                delete _lensMedium;
                delete _globalState;
            }

        //! returns the body force vector within a cell
        const Vector &gravity(const Cell &cell,
                              int cellIdx) const
            {
                return _globalState->gravity();
            }

        const Vector &gravity() const
            {
                return _globalState->gravity();
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
                return _globalState->densityW();
            }

        // return the density of the non-wetting phase
        Scalar densityN() const
            {
                return _globalState->densityN();
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
                return _globalState->viscosityW();
            }

        // return the viscosity of the non-wetting phase
        Scalar viscosityN() const
            {
                return _globalState->viscosityN();
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
            { return _cellStates[index]; }
        const CellState &cellState(int index) const
            { return _cellStates[index]; }

        CellState &cellState(const Cell &cell)
            { return _cellStates[cellIndex(cell)]; }
        const CellState &cellState(const Cell &cell) const
            { return _cellStates[cellIndex(cell)]; }

        // given a global node index, return the corresponding state
        NodeState &nodeState(int index)
            { return _nodeStates[index]; }
        const NodeState &nodeState(int index) const
            { return _nodeStates[index]; }

        // given a cell and a local node index, return the corresponding state
        NodeState &nodeState(const Cell &cell, int i)
            { return _nodeStates[ParentType::nodeIndex(cell, i)]; }
        const NodeState &nodeState(const Cell &cell, int i) const
            { return _nodeStates[ParentType::nodeIndex(cell, i)]; }

        // given a node, return it's state object
        NodeState &nodeState(const Node &vert)
            { return _nodeStates[ParentType::nodeIndex(vert)]; }
        const NodeState &nodeState(const Node &vert) const
            { return _nodeStates[ParentType::nodeIndex(vert)]; }

        const WorldCoord &lowerLeft() const
            { return _gridLowerLeft; }
        const WorldCoord &upperRight() const
            { return _gridUpperRight; }

        Scalar width() const
            { return _gridUpperRight[0] - _gridLowerLeft[0]; }
        Scalar height() const
            { return _gridUpperRight[1] - _gridLowerLeft[1]; }

        bool onUpperBoundary(const WorldCoord &pos) const
            { return pos[1] > _gridUpperRight[1] - _eps; }
        bool onLowerBoundary(const WorldCoord &pos) const
            { return pos[1] < _gridLowerLeft[1] + _eps; }
        bool onLeftBoundary(const WorldCoord &pos) const
            { return pos[0] < _gridLowerLeft[0] + _eps; }
        bool onRightBoundary(const WorldCoord &pos) const
            { return pos[0] > _gridUpperRight[0] - _eps; }

        // returns true iff a world coordinate is within the fine
        // sand lens.
        bool isInLens(const WorldCoord &coord)
            {
                for (int i = 0; i < DomainTraits::WorldDim; ++i) {
                    if (_lensLowerLeft[i] > coord[i] ||
                        _lensUpperRight[i] < coord[i])
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
                return state.mediumState() == _lensMedium;
            }
*/

    private:
        void _initGrid()
            {
                // specify the position and size of the grid in world
                // coordinates
                _gridLowerLeft[0] = 0;
                _gridLowerLeft[1] = 0;
                _gridUpperRight[0] = 6;
                _gridUpperRight[1] = 4;

                // the epsilon constant
                _eps = 1e-8 * width();

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
                                      _gridLowerLeft,
                                      _gridUpperRight);
                ParentType::setGrid(grid);

                // specify the location and the size of the fine sand
                // lens in world coordinates
                _lensLowerLeft[0] = width()*(1./6);
                _lensLowerLeft[1] = height()*(2./4);
                _lensUpperRight[0] = width()*(4./6);
                _lensUpperRight[1] = height()*(3./4);
                // make sure the global coordniates of the lens are
                // aligned to the grid's cell boundaries in order to
                // make sure the interface occures at the same
                // position for sized lenses for node and cell based
                // material parameter storage. (this is not strictly
                // necessary, but makes it easier to compare the
                // different approaches.)
                _lensLowerLeft[0] = width()/cellRes[0]
                                  * (int) (_lensLowerLeft[0]*cellRes[0]/width())
                                  + (0.5 - _eps)/cellRes[0];
                _lensLowerLeft[1] = height()/cellRes[1]
                                  * (int) (_lensLowerLeft[1]*cellRes[1]/height())
                                  + (0.5 - _eps)/cellRes[1];
                _lensUpperRight[0] = width()/cellRes[0]
                                   * (int) (_lensUpperRight[0]*cellRes[0]/width())
                                   + (-0.5 + _eps)/cellRes[0];
                _lensUpperRight[1] = height()/cellRes[1]
                                   * (int) (_lensUpperRight[1]*cellRes[1]/height())
                                   + (-0.5 + _eps)/cellRes[1];
            }

        void _initGlobalState()
            {
                // initialize the globally uniform parameters
                _globalState = new GlobalState();

                Water water;
                DNAPL dnapl;
                _globalState->setDensityW(water.density());
                _globalState->setDensityN(dnapl.density());

                _globalState->setViscosityW(water.viscosity());
                _globalState->setViscosityN(dnapl.viscosity());
                Vector gravity(0); gravity[1] = -9.81;
                _globalState->setGravity(gravity);
            }

        void _initMediaStates()
            {
                // initializing the state of the lens' medium
                _lensMedium =  new MediumState();
                _lensMedium->setGlobalState(_globalState);
                _lensMedium->setSwr(0.18);
#if USE_ORIG_PROB
                _lensMedium->setSnr(0.0);
#else
                _lensMedium->setSnr(0.25);
#endif
                _lensMedium->setPorosity(0.4);
                _lensMedium->setMdcParams(VanGenuchtenState(0.00045, 7.3));
#if USE_DIFFERENT_MAIN_CURVES
                _lensMedium->setMicParams(VanGenuchtenState(0.00045*2, 7.5));
#else
                _lensMedium->setMicParams(VanGenuchtenState(0.00045, 7.3));
#endif
                _lensMedium->micParams().setVgMaxPC(_lensMedium->mdcParams().vgMaxPC());
                _lensMedium->setPermeability(9.05e-13);

                // initializing the state of the outer medium
                _outerMedium =  new MediumState();
                _outerMedium->setGlobalState(_globalState);
                _outerMedium->setSwr(0.05);
#if USE_ORIG_PROB
                _outerMedium->setSnr(0.0);
#else
                _outerMedium->setSnr(0.20);
#endif
                _outerMedium->setPorosity(0.4);
                _outerMedium->setMdcParams(VanGenuchtenState(0.0037, 4.7));
#if USE_DIFFERENT_MAIN_CURVES
                _outerMedium->setMicParams(VanGenuchtenState(0.0037*2, 4.8));
#else
                _outerMedium->setMicParams(VanGenuchtenState(0.0037, 4.7));
#endif

                _outerMedium->micParams().setVgMaxPC(_outerMedium->mdcParams().vgMaxPC());
                _outerMedium->setPermeability(4.6e-10);

#if 0
                _lensMedium->setSwr(_outerMedium->Swr());
                _lensMedium->setSnr(_outerMedium->Snr());
                _lensMedium->setPorosity(_outerMedium->porosity());
                _lensMedium->setMicParams(_outerMedium->micParams());
                _lensMedium->setMdcParams(_outerMedium->mdcParams());
                _lensMedium->setPermeability(_outerMedium->permeability());
#endif

            }

        void _initCellStates()
            {
                _cellStates.resize(ParentType::numCells());

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
                        _cellStates[arrayPos].setMediumState(_lensMedium);
                    else
                        _cellStates[arrayPos].setMediumState(_outerMedium);
                }
            }

        void _initNodeStates()
            {
                _nodeStates.resize(ParentType::numVertices());
#if USE_NODE_PARAMETERS
                NodeIterator vertIt = ParentType::nodeBegin();
                const NodeIterator &endVertIt = ParentType::nodeEnd();
                for (; vertIt != endVertIt; ++vertIt) {
                    WorldCoord pos;
                    ParentType::nodePosition(pos, *vertIt);
                    NodeState &vertState = nodeState(*vertIt);
                    if (isInLens(pos))
                        vertState.setMediumState(_lensMedium);
                    else
                        vertState.setMediumState(_outerMedium);
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
                        _markInterfaceVertices(faceIt);
                    }
                }
#endif
            }

#if !USE_NODE_PARAMETERS
        //  mark all vertices on a face as as belonging to an
        //  inhomogenity
        void _markInterfaceVertices(IntersectionIterator &interfaceIt)
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
        typedef std::vector<NodeState>  _NodeStateArray;

        // the actual type which stores the cell states.
        typedef std::vector<CellState>    _CellStateArray;

        // the lower left and upper right coordinates of the complete
        // grid
        WorldCoord _gridLowerLeft;
        WorldCoord _gridUpperRight;

        // the lower left and upper right coordinates of the fine sand
        // lens
        WorldCoord _lensLowerLeft;
        WorldCoord _lensUpperRight;

        // global state and states of the media
        GlobalState *_globalState;
        MediumState *_outerMedium;
        MediumState *_lensMedium;

        // stores a state for each node
        _NodeStateArray _nodeStates;

        // stores a state for each node
        _CellStateArray  _cellStates;

        // A small epsilon value and the current simulated
        // time. FIXME/TODO: should probably not be here..
        Scalar _eps;
    };
}
}


#endif
