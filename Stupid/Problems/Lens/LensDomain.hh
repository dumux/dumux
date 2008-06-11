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
#ifndef STUPID_LENS_DOMAIN_HH
#define STUPID_LENS_DOMAIN_HH

#include "LensState.hh"

#include <Stupid/Auxilary/BasicDomain.hh>

#include <Stupid/Models/Box/PwSn/PwSnBoxTraits.hh>

#include <Stupid/Material/RegularizedVanGenuchten.hh>
#include <Stupid/Material/ParkerLenhard.hh>

#include <dumux/material/properties.hh>

#include <dune/grid/sgrid.hh>

namespace Stupid
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
            typedef LensVertexState<MediumState>               VertexState;
            typedef LensCellState<MediumState>                 CellState;
            
            // The parker-lenhard hysteresis model
#if !defined USE_VERTEX_PARAMETERS
            typedef Stupid::ParkerLenhard<CellState, VanGenuchten>     ParkerLenhard;
#else // defined USE_VERTEX_PARAMETERS
            typedef Stupid::ParkerLenhard<VertexState, VanGenuchten>   ParkerLenhard;
#endif
        };


    private:

        typedef typename DomainTraits::Grid                  Grid;
        typedef typename DomainTraits::Cell                  Cell;
        typedef typename DomainTraits::CellReferenceElement  CellReferenceElement;
        typedef typename DomainTraits::CellReferenceElements CellReferenceElements;

        typedef typename DomainTraits::Vertex                Vertex;

        typedef typename DomainTraits::CellIterator         CellIterator;
        typedef typename DomainTraits::VertexIterator       VertexIterator;

        typedef typename DomainTraits::IntersectionIterator       IntersectionIterator;
        typedef typename DomainTraits::IntersectionIteratorGetter IntersectionIteratorGetter;

        typedef typename DomainTraits::LocalCoord    LocalCoord;
        typedef typename DomainTraits::WorldCoord    WorldCoord;
        typedef typename DomainTraits::Vector        Vector;
        typedef typename DomainTraits::Matrix        Matrix;

        typedef typename MaterialTraits::GlobalState   GlobalState;
        typedef typename MaterialTraits::MediumState   MediumState;
        typedef typename MaterialTraits::CellState     CellState;
        typedef typename MaterialTraits::VertexState   VertexState;

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
                _initVertexStates();
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
        const Matrix &permeability(const Cell &cell) const
            {
                return cellState(cellIndex(cell)).permeability();
            }

        //! Return the capillary pressure for a given vertex of a cell
        Scalar pC(const Cell &cell,
                  int cellIdx,
                  int localVertIdx,
                  int globalVertIdx,
                  Scalar Sw) const
            {
#if defined USE_VERTEX_PARAMETERS
                return ParkerLenhard::pC(vertexState(globalVertIdx), Sw);
#else // !defined USE_VERTEX_PARAMETERS
#if defined USE_INTERFACE_CONDITION
                const VertexState &vertState = vertexState(globalVertIdx);
                if (!vertState.isOnInterface())
                    return ParkerLenhard::pC(cellState(cellIdx),
                                             Sw);

                // find the minimum capilarry pressure around the
                // vertex if the vertex is on an interface
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

        // return the mobility of the wetting phase at a vertex
        Scalar mobilityW(const Cell &cell,
                         int cellIdx,
                         int localVertIdx,
                         int globalVertIdx,
                         Scalar Sw) const
            {
#if defined USE_VERTEX_PARAMETERS
                return ParkerLenhard::krw(vertexState(globalVertIdx), Sw) / viscosityW();
#else // !defined USE_VERTEX_PARAMETERS
#ifdef USE_INTERFACE_CONDITION
                const VertexState &vertState = vertexState(globalVertIdx);
                if (!vertState.isOnInterface())
                    return ParkerLenhard::krw(cellState(cellIdx), Sw) / viscosityW();

                // find the minimum wetting phase mobility around the
                // vertex if the vertex is on an interface
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

        // return the mobility of the non-wetting phase at a vertex
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

#if defined USE_VERTEX_PARAMETERS
                return ParkerLenhard::krn(vertexState(globalVertIdx), Sw) / viscosityN();
#else // !defined USE_VERTEX_PARAMETERS
#ifdef USE_INTERFACE_CONDITION
                const VertexState &vertState = vertexState(globalVertIdx);
                if (!vertState.isOnInterface())
                    return ParkerLenhard::krn(cellState(cellIdx), Sw) / viscosityN();

                // find the minimum wetting phase mobility around the
                // vertex if the vertex is on an interface
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

        // given a global vertex index, return the corresponding state
        VertexState &vertexState(int index)
            { return _vertexStates[index]; }
        const VertexState &vertexState(int index) const
            { return _vertexStates[index]; }

        // given a cell and a local vertex index, return the corresponding state
        VertexState &vertexState(const Cell &cell, int i)
            { return _vertexStates[ParentType::vertexIndex(cell, i)]; }
        const VertexState &vertexState(const Cell &cell, int i) const
            { return _vertexStates[ParentType::vertexIndex(cell, i)]; }

        // given a vertex, return it's state object
        VertexState &vertexState(const Vertex &vert)
            { return _vertexStates[ParentType::vertexIndex(vert)]; }
        const VertexState &vertexState(const Vertex &vert) const
            { return _vertexStates[ParentType::vertexIndex(vert)]; }

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
        bool isInLens(const VertexState &state)
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
#ifndef USE_ORIG_PROB
                cellRes[0] = 8;
                cellRes[1] = 24;
//                cellRes[0] = 48;
//                cellRes[1] = 32;
#else
                cellRes[0] = 48;
                cellRes[1] = 32;
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
                // position for sized lenses for vertex and cell based
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
#ifdef USE_ORIG_PROB
                _lensMedium->setSnr(0.0);
#else
                _lensMedium->setSnr(0.25);
#endif
                _lensMedium->setPorosity(0.4);
                _lensMedium->setMicParams(VanGenuchtenState(0.00045, 7.3));
                _lensMedium->setMdcParams(VanGenuchtenState(0.00045, 7.3));
                _lensMedium->micParams().setVgMaxPC(_lensMedium->mdcParams().vgMaxPC());
                _lensMedium->setPermeability(9.05e-13);

                // initializing the state of the outer medium
                _outerMedium =  new MediumState();
                _outerMedium->setGlobalState(_globalState);
                _outerMedium->setSwr(0.05);
#ifdef USE_ORIG_PROB
                _outerMedium->setSnr(0.0);
#else
                _outerMedium->setSnr(0.20);
#endif
                _outerMedium->setPorosity(0.4);
                _outerMedium->setMicParams(VanGenuchtenState(0.0037, 4.7));
                _outerMedium->setMdcParams(VanGenuchtenState(0.0037, 4.7));
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

        void _initVertexStates()
            {
                _vertexStates.resize(ParentType::numVertices());
#if defined USE_VERTEX_PARAMETERS
                VertexIterator vertIt = ParentType::vertexBegin();
                const VertexIterator &endVertIt = ParentType::vertexEnd();
                for (; vertIt != endVertIt; ++vertIt) {
                    WorldCoord pos;
                    ParentType::vertexPosition(pos, *vertIt);
                    VertexState &vertState = vertexState(*vertIt);
                    if (isInLens(pos))
                        vertState.setMediumState(_lensMedium);
                    else
                        vertState.setMediumState(_outerMedium);
                }

#else // ! defined USE_VERTEX_PARAMETERS
                // initialize the vertex state objects

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

#if !defined USE_VERTEX_PARAMETERS
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
                for (int vertexInFace = 0;
                     vertexInFace < nVerticesOfFace;
                     vertexInFace++)
                {
                    int vertexIdxInElement = refElem.subEntity(faceIdx, 1, vertexInFace, GridDim);

                    vertexState(*interfaceIt.inside(), vertexIdxInElement).addNeighbourCellIdx(inIdx);
                    vertexState(*interfaceIt.inside(), vertexIdxInElement).addNeighbourCellIdx(outIdx);
                }
            }
#endif // !defined USE_VERTEX_PARAMETERS


        // the actual type which stores the cell states.
        typedef std::vector<VertexState>  _VertexStateArray;

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

        // stores a state for each vertex
        _VertexStateArray _vertexStates;

        // stores a state for each vertex
        _CellStateArray  _cellStates;

        // A small epsilon value and the current simulated
        // time. FIXME/TODO: should probably not be here..
        Scalar _eps;
    };
}
}


#endif
