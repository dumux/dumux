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
#include <dumux/new_models/box/pwsn/pwsnboxtraits.hh>

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
            typedef LenhardVertexState<MediumState>               VertexState;
            typedef LenhardCellState<MediumState>                 CellState;
            
            // The parker-lenhard hysteresis model
#if !defined USE_VERTEX_PARAMETERS
            typedef Dune::ParkerLenhard<CellState, VanGenuchten>     ParkerLenhard;
#else // defined USE_VERTEX_PARAMETERS
            typedef Dune::ParkerLenhard<VertexState, VanGenuchten>   ParkerLenhard;
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
        PwSnLenhardDomain()
            {
                Api::require<Api::BasicDomainTraits, DomainTraits>();

                _initGrid();

                _initGlobalState();
                _initMediaStates();
                _initCellStates();
                _initVertexStates();
            };

        ~PwSnLenhardDomain()
            {
                delete _coarseSand;
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
/*                printf("krw(%f): %f, viscosity: %f\n",
                       Sw,
                       ParkerLenhard::krw(cellState(cellIdx), Sw),
                       viscosityW());
*/

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
/*                printf("krn(%f): %f, viscosity: %f\n",
                  1 - Sn,
                  ParkerLenhard::krn(cellState(cellIdx), 1 - Sn),
                  viscosityW());
*/

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

        Scalar height() const
            { return _gridUpperRight[0] - _gridLowerLeft[0]; }

        bool onUpperBoundary(const WorldCoord &pos) const
            { return pos[0] > _gridUpperRight[0] - _eps; }
        bool onLowerBoundary(const WorldCoord &pos) const
            { return pos[0] < _gridLowerLeft[0] + _eps; }

    private:
        void _initGrid()
            {
                // specify the position and size of the grid in world
                // coordinates
                _gridLowerLeft[0] = 0;

                // column was 72 cm high, but we need some additional
                // space so that the boundary condition at the upper
                // boundary doesn't destroy the experiment
                _gridUpperRight[0] = 0.80;


                // the epsilon constant
                _eps = 1e-8 * height();

                // create the grid
                Dune::FieldVector<int, GridDim> cellRes;
                cellRes[0] = 160;

                Grid *grid = new Grid(cellRes,
                                      _gridLowerLeft,
                                      _gridUpperRight);
                ParentType::setGrid(grid);
            }

        void _initGlobalState()
            {
                // initialize the globally uniform parameters
                _globalState = new GlobalState();

                Water wettingFluid;
                Air   nonwettingFluid;
                _globalState->setDensityW(wettingFluid.density());
                _globalState->setDensityN(nonwettingFluid.density());

                _globalState->setViscosityW(wettingFluid.viscosity());
                _globalState->setViscosityN(nonwettingFluid.viscosity());

                Vector gravity(0);
                gravity[0] = -9.81;
                _globalState->setGravity(gravity);
            }

        void _initMediaStates()
            {
                // initializing the state of the porous medium
                _coarseSand =  new MediumState();
                _coarseSand->setGlobalState(_globalState);
#if LENHARD_EXPERIMENT == 1
                _coarseSand->setSwr(0.0);
                _coarseSand->setSnr(0.25);
                _coarseSand->setPorosity(0.36);
                _coarseSand->setMicParams(VanGenuchtenState(0.00052, 4.63));
                _coarseSand->setMdcParams(VanGenuchtenState(0.00052, 4.63));
                _coarseSand->micParams().setVgMaxPC(_coarseSand->mdcParams().vgMaxPC());
                _coarseSand->setPermeability(3.37e-11);
#elif LENHARD_EXPERIMENT == 2
                _coarseSand->setSwr(0.17);
                _coarseSand->setSnr(0.25);
                _coarseSand->setPorosity(0.36);
                _coarseSand->setMdcParams(VanGenuchtenState(0.00042, 5.25));
//                _coarseSand->setMicParams(VanGenuchtenState(0.00042, 5.25));
                _coarseSand->setMicParams(VanGenuchtenState(0.00042*2, 5.25));
                _coarseSand->micParams().setVgMaxPC(_coarseSand->mdcParams().vgMaxPC());
                _coarseSand->setPermeability(3.37e-11);
#endif
            }

        void _initCellStates()
            {
                _cellStates.resize(ParentType::numCells());

                // initialize the cell state objects depending on
                // wether the cell is inside or outside of the lenhard of
                // fine sand.
                WorldCoord cellCenter;
                CellIterator it = ParentType::cellBegin();
                CellIterator endit = ParentType::cellEnd();
                for (; it != endit; ++it) {
                    ParentType::cellCenter(*it, cellCenter);
                    int arrayPos = ParentType::cellIndex(*it);
                    _cellStates[arrayPos].setMediumState(_coarseSand);
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

                    vertState.setMediumState(_coarseSand);
                }

#else // ! defined USE_VERTEX_PARAMETERS
                // initialize the vertex state objects
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
                        _markInterfaceVertices(faceIt);
                    }
                }
                */
#endif
            }

#if !defined USE_VERTEX_PARAMETERS
        //  mark all vertices on a face as as belonging to an
        //  inhomogenity
        /*
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
        */
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
        // lenhard
        WorldCoord _lenhardLowerLeft;
        WorldCoord _lenhardUpperRight;

        // global state and states of the media
        GlobalState *_globalState;
        MediumState *_coarseSand;

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
