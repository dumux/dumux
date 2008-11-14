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
 * \brief Provides some convenient interfaces to cells and vertices of grids.
 */
#ifndef DUMUX_BASIC_DOMAIN_HH
#define DUMUX_BASIC_DOMAIN_HH

#include <dumux/auxiliary/apis.hh>

#include <dune/common/fvector.hh>
#include <dune/common/fmatrix.hh>

#include <dune/grid/common/referenceelements.hh>
#include <dune/grid/common/intersectioniterator.hh>
#include <dune/grid/io/file/dgfparser.hh>
#include <dune/grid/utility/intersectiongetter.hh>
#include <dune/grid/common/mcmgmapper.hh>

#include <dune/disc/operators/boundaryconditions.hh>

namespace Dune
{
    namespace Api {
        BEGIN_API_DEF(BasicDomainTraits)
        {
            typedef Implementation I;

            // make sure the required constants exist
            int i;
            i = I::GridDim;
            i = I::WorldDim;

            // make sure the required types exist
            {typename I::Grid                       *x;x=NULL;}
            {typename I::Scalar                     *x;x=NULL;}
            {typename I::CoordScalar                *x;x=NULL;}
            {typename I::LocalCoord                 *x;x=NULL;}
            {typename I::WorldCoord                 *x;x=NULL;}
            {typename I::Cell                       *x;x=NULL;}
            {typename I::CellIterator               *x;x=NULL;}
            {typename I::Node                       *x;x=NULL;}
            {typename I::NodeIterator               *x;x=NULL;}
            {typename I::CellReferenceElement       *x;x=NULL;}
            {typename I::CellReferenceElements      *x;x=NULL;}
            {typename I::IntersectionIteratorGetter *x;x=NULL;}
            {typename I::IntersectionIterator       *x;x=NULL;}
            {typename I::Vector                     *x;x=NULL;}
            {typename I::Matrix                     *x;x=NULL;}
        }
        END_API_DEF
    };

    /*!
     * \brief Provides some convenient interfaces to cells and vertices of grids.
     *
     * This class doesn't handle any user data itself, it just
     * provides a way to map form cells and vertices to their
     * respective indices.
     */
    template<class GridT, class ScalarT>
    class BasicDomain
    {
    public:
        /*!
         * \brief Provides a few default types to interface cell or node
         *        based grids.
         */
        class DomainTraits
        {
        public:
            typedef GridT       Grid;
            typedef ScalarT     Scalar;
            
            enum {
                GridDim = GridT::dimension,        //!< Dimension of the grid.
                WorldDim = GridT::dimensionworld   //!< Dimension of the world the grid is embedded into.
            };
            
            // coordinate stuff
            typedef typename Grid::ctype                     CoordScalar;
            typedef Dune::FieldVector<CoordScalar, GridDim>  LocalCoord;
            typedef Dune::FieldVector<CoordScalar, WorldDim> WorldCoord;
            
            // Dune grid related stuff: iterators, cells, vertices, traits...
            typedef typename Grid::Traits::template Codim<0>                  CellTraits;
            typedef typename CellTraits::Entity                               Cell;
            typedef typename CellTraits::LeafIterator                         CellIterator;
            
            typedef typename Grid::Traits::template Codim<GridDim>            NodeTraits;
            typedef typename NodeTraits::Entity                               Node;
            typedef typename NodeTraits::LeafIterator                         NodeIterator;
            
            typedef Dune::ReferenceElement<CoordScalar, GridDim>              CellReferenceElement;
            typedef Dune::ReferenceElements<CoordScalar, GridDim>             CellReferenceElements;
            typedef Dune::IntersectionIteratorGetter<Grid,Dune::LeafTag>      IntersectionIteratorGetter;
            typedef typename IntersectionIteratorGetter::IntersectionIterator IntersectionIterator;
            
            // grid-space vector and matrix types
            typedef Dune::FieldVector<Scalar, GridDim>         Vector;
            typedef Dune::FieldMatrix<Scalar, GridDim,GridDim> Matrix;
        };

    private:
        // some types from the traits for convenience
        typedef typename DomainTraits::Grid                  Grid;
        typedef typename DomainTraits::Cell                  Cell;
        typedef typename DomainTraits::CellIterator          CellIterator;
        typedef typename DomainTraits::Node                  Node;
        typedef typename DomainTraits::NodeIterator          NodeIterator;
        typedef typename DomainTraits::LocalCoord            LocalCoord;
        typedef typename DomainTraits::WorldCoord            WorldCoord;
        typedef typename DomainTraits::CellReferenceElement  CellReferenceElement;
        typedef typename DomainTraits::CellReferenceElements CellReferenceElements;

        typedef Dune::GridPtr<Grid>                          GridPointer;

        enum {
            GridDim = DomainTraits::GridDim,
            WorldDim = DomainTraits::WorldDim
        };

        // helper class for DUNE. This class tells the element mapper
        // whether data should be attached to cells of a given
        // geometry type. Here we only attach data data to any
        // geometry type as long as the entity represents a call
        // (i.e. leaf in DUNE-lingo).
        template<int dim>
        struct _CellLayout
        { bool contains(Dune::GeometryType geoType) { return geoType.dim() == dim; } };

        // mapper from grid cells to the actual array indices
        // of the cell state vector.
        typedef Dune::LeafMultipleCodimMultipleGeomTypeMapper<Grid, _CellLayout> CellMap;

        // the set of indices for the grid's cells
        typedef typename Grid::Traits::LeafIndexSet _CellIndexSet;

        // mapper: one data element per node
        template<int dim>
        struct _NodeLayout
        { bool contains (Dune::GeometryType gt) { return gt.dim() == 0; } };
        typedef Dune::MultipleCodimMultipleGeomTypeMapper<Grid,_CellIndexSet,_NodeLayout> NodeMap;
        
    public:
        BasicDomain()
            {
                Api::require<Api::BasicDomainTraits, DomainTraits>();

                _cellMap = NULL;
                _nodeMap = NULL;
            };

        BasicDomain(GridPointer grid)
            {
                Api::require<Api::BasicDomainTraits, DomainTraits>();

                _cellMap = NULL;
                _nodeMap = NULL;
                setGrid(grid);
            };

        ~BasicDomain()
            {
                delete _cellMap;
                delete _nodeMap;
            }

        /*!
         * \brief Returns the current grid.
         */
        Grid &grid()
            { return *_grid; }

        /*!
         * \brief Returns the current grid.
         */
        const Grid &grid() const
            { return *_grid; }

        /*!
         * \brief Returns the current grid's number of cells.
         */
        int numCells() const
            { return _cellMap->size(); }

        /*!
         * \brief Given a cell, return it's respective index.
         */
        int cellIndex(const Cell &e) const
            { return _cellMap->map(e); }

        /*!
         * \brief Returns the iterator pointing to the first cell of
         *        the grid.
         */
        CellIterator cellBegin() {
            return _grid->template leafbegin<0>();
        }

        /*!
         * \brief Returns the iterator pointing to the next to last
         *        cell of the grid.
         */
        CellIterator cellEnd() {
            return _grid->template leafend<0>();
        }

        /*!
         * \brief Calculate the center of gravity of an cell in world
         *        coordinates.
         */
        void cellCenter(const Cell &e, WorldCoord &worldCoord) const
            {
                // get the center position of the cell's reference
                // element
                Dune::GeometryType geoType = e.type();
                const LocalCoord &localCoord =
                    CellReferenceElements::general(geoType).position(0, 0);

                // get world coordinate of the cell center
                worldCoord = e.geometry().global(localCoord);
            }

        /*!
         * \brief Returns the cell map in case it is externally required.
         */
        CellMap &cellMap()
            { return *_cellMap; }
        const CellMap &cellMap() const
            { return *_cellMap; }

        /*!
         * \brief Returns the current grid's number of nodes.
         */
        int numNodes() const
            { return _nodeMap->size(); }

        /*!
         * \brief Given a node index within a cell, return the
         *        respective global index.
         */
        int nodeIndex(const Cell &e, int localNodeId) const
            { return _nodeMap->template map<GridDim>(e, localNodeId); }

        /*!
         * \brief Given a node index within a cell, return the
         *        respective global index.
         */
        int nodeIndex(const Node &v) const
            { return _nodeMap->map(v); }

        /*!
         * \brief Returns the iterator pointing to the first cell of
         *        the grid.
         */
        NodeIterator nodeBegin() {
            return _grid->template leafbegin<GridDim>();
        }

        /*!
         * \brief Returns the iterator pointing to the next to last
         *        cell of the grid.
         */
        NodeIterator nodeEnd() {
            return _grid->template leafend<GridDim>();
        }

        /*!
         * \brief Return the entity for a node given a cell and the
         *        node' local index.
         */
        const Node &node(const Cell &cell, int localVertIdx)
            { return *cell.template entity<GridDim>(localVertIdx); };

        /*!
         * \brief Given a node return its position in world coodinates.
         */
        void nodePosition(WorldCoord &worldCoord,
                            const Node &v) const
            {
                worldCoord = v.geometry()[0];
            }

        /*!
         * \brief Returns the node map in case it is externally required.
         */
        NodeMap &nodeMap()
            { return *_nodeMap; }
        const NodeMap &nodeMap() const
            { return *_nodeMap; }


        /*!
         * \brief This method must be called after something within the grid has been
         *        altered.
         *
         * Cases where this method should be called are grid
         * refinements, deformation of cells, etc.
         */
        void gridChanged()
            {
                delete _cellMap;
                delete _nodeMap;

                // create the cell and node index sets for the grid
                _cellMap = new CellMap(*_grid);
                _nodeMap = new NodeMap(*_grid, _grid->leafIndexSet());
            }

        /*!
         * \brief Set the current grid.
         *
         * This method assumes that the grid was allocated on the heap
         * and it takes ownership of it, so you may not delete it
         * externally! The grid which was previously in the Domain
         * becomes invalid after calling this method.
         */
        void setGrid(GridPointer grid)
            {
                _grid = grid;

                gridChanged();
            };

    private:
        // pointer to the grid object
        GridPointer _grid;

        // map from the grid's leafs to an index of the index set
        CellMap *_cellMap;
        // translates vertices local to a cell to global ones
        NodeMap *_nodeMap;
    };
}

#endif
