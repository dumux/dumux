/*****************************************************************************
 *   Copyright (C) 2008 by Bernd Flemisch, Andreas Lauser                    *
 *   Institute of Hydraulic Engineering                                      *
 *   University of Stuttgart, Germany                                        *
 *   email: <givenname>.<name>@iws.uni-stuttgart.de                          *
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
 * \brief Provides some convenient interfaces to elements and vertices of grids.
 */
#ifndef DUMUX_BASIC_DOMAIN_HH
#define DUMUX_BASIC_DOMAIN_HH

#include <dune/common/exceptions.hh>

#include <dune/common/fvector.hh>
#include <dune/common/fmatrix.hh>

#include <dune/grid/common/referenceelements.hh>
#include <dune/grid/common/intersectioniterator.hh>
#include <dune/grid/common/mcmgmapper.hh>

#include <dune/disc/operators/boundaryconditions.hh>

#include <iostream>

namespace Dune
{
/*!
 * \brief Provides some convenient interfaces to elements and vertices of grids.
 *
 * This class doesn't handle any user data itself, it just
 * provides a way to map form elements and vertices to their
 * respective indices.
 */
template<class GridT, class ScalarT>
class BasicDomain
{
public:
    /*!
     * \brief Provides a few default types to interface element or vert
     *        based grids.
     */
    class DomainTraits
    {
    public:
        typedef GridT       Grid;
        typedef ScalarT     Scalar;

        enum {
            dim = GridT::dimension,        //!< Dimension of the grid.
            dimWorld = GridT::dimensionworld   //!< Dimension of the world the grid is embedded into.
        };

        // coordinate stuff
        typedef typename Grid::ctype                     CoordScalar;
        typedef Dune::FieldVector<CoordScalar, dim>      LocalPosition;
        typedef Dune::FieldVector<CoordScalar, dimWorld> GlobalPosition;

        // Dune grid related stuff: iterators, elements, vertices, traits...
        typedef typename Grid::Traits::template Codim<0>                ElementTraits;
        typedef typename ElementTraits::Entity                          Element;
        typedef typename ElementTraits
        ::template Partition<All_Partition>
        ::LeafIterator                                 ElementIterator;

        typedef typename Grid::Traits::template Codim<dim>              VertexTraits;
        typedef typename VertexTraits::Entity                           Vertex;
        typedef typename VertexTraits
        ::template Partition<All_Partition>
        ::LeafIterator                                 VertexIterator;

        typedef typename Grid::template Codim<0>::LeafIntersectionIterator IntersectionIterator;

        // grid-space vector and matrix types
        typedef Dune::FieldVector<Scalar, dim>     FieldVector;
        typedef Dune::FieldMatrix<Scalar, dim,dim> FieldMatrix;
    };

private:
    // some types from the traits for convenience
    typedef typename DomainTraits::Grid                  Grid;
    typedef typename DomainTraits::Element               Element;
    typedef typename DomainTraits::ElementIterator       ElementIterator;
    typedef typename DomainTraits::Vertex                Vertex;
    typedef typename DomainTraits::VertexIterator        VertexIterator;
    typedef typename DomainTraits::LocalPosition         LocalPosition;
    typedef typename DomainTraits::GlobalPosition        GlobalPosition;

    typedef typename Grid::LeafGridView                  GridView;

    enum {
        dim = DomainTraits::dim,
        dimWorld = DomainTraits::dimWorld
    };

    // helper class for DUNE. This class tells the element mapper
    // whether data should be attached to elements of a given
    // geometry type. Here we only attach data data to any
    // geometry type as long as the entity represents a call
    // (i.e. leaf in DUNE-lingo).
    template<int dim>
    struct ElementLayout_
    { bool contains(Dune::GeometryType geoType) { return geoType.dim() == dim; } };

    // mapper from grid elements to the actual array indices
    // of the element state vector.
    typedef Dune::LeafMultipleCodimMultipleGeomTypeMapper<Grid, ElementLayout_> ElementMap;

    // the set of indices for the grid's elements
    typedef typename Grid::Traits::LeafIndexSet ElementIdxSet_;

    // mapper: one data element per vert
    template<int dim>
    struct VertexLayout_
    { bool contains (Dune::GeometryType gt) { return gt.dim() == 0; } };
    typedef Dune::MultipleCodimMultipleGeomTypeMapper<GridView, VertexLayout_> VertexMap;

public:
    BasicDomain()
    {
        elementMap_ = NULL;
        vertMap_ = NULL;
    };

    BasicDomain(Grid *grid)
    {
        elementMap_ = NULL;
        vertMap_ = NULL;
        setGrid(grid);
    };

    ~BasicDomain()
    {
        delete elementMap_;
        delete vertMap_;
    }


    /*!
     * \brief Returns the current grid.
     */
    /*        Grid &grid()
              { return const_cast<Grid&>(*grid_); }
    */


    /*!
     * \brief Returns the current grid.
     */
    const Grid &grid() const
    { return *grid_; }

    /*!
     * \brief Returns the current grid view.
     */
    const GridView gridView() const
    { return grid_->leafView(); }

    /*!
     * \brief Returns the current grid's number of elements.
     */
    int numElements() const
    { return elementMap_->size(); }

    /*!
     * \brief Given a element, return it's respective index.
     */
    int elementIdx(const Element &e) const
    { return elementMap_->map(e); }

    /*!
     * \brief Returns the iterator pointing to the first element of
     *        the grid.
     */
    ElementIterator elementBegin() const {
        return grid_->template leafbegin<0, All_Partition>();
    }

    /*!
     * \brief Returns the iterator pointing to the next to last
     *        element of the grid.
     */
    ElementIterator elementEnd() {
        return grid_->template leafend<0, All_Partition>();
    }

    /*!
     * \brief Calculate the center of gravity of an element in world
     *        coordinates.
     */
    void elementCenter(const Element &e, GlobalPosition &worldCoord) const
    {
        // get the center position of the element's reference
        // element
        Dune::GeometryType geoType = e.type();
        const LocalPosition &localCoord =
            DomainTraits::referenceElement(geoType).position(0, 0);

        // get world coordinate of the element center
        worldCoord = e.geometry().global(localCoord);
    }

    /*!
     * \brief Returns the element map in case it is externally required.
     */
    ElementMap &elementMap()
    { return *elementMap_; }
    const ElementMap &elementMap() const
    { return *elementMap_; }

    /*!
     * \brief Returns the current grid's number of verts.
     */
    int numVertices() const
    { return vertMap_->size(); }

    /*!
     * \brief Given a vert index within a element, return the
     *        respective global index.
     */
    int vertexIdx(const Element &e, int localVertexId) const
    { return vertMap_->map(e, localVertexId, dim); }

    /*!
     * \brief Given a vert index within a element, return the
     *        respective global index.
     */
    int vertexIdx(const Vertex &v) const
    { return vertMap_->map(v); }

    /*!
     * \brief Returns the iterator pointing to the first element of
     *        the grid.
     */
    VertexIterator vertexBegin() {
        return grid_->template leafbegin<dim, All_Partition>();
    }

    /*!
     * \brief Returns the iterator pointing to the next to last
     *        element of the grid.
     */
    VertexIterator vertexEnd() {
        return grid_->template leafend<dim, All_Partition>();
    }

    /*!
     * \brief Return the entity for a vert given a element and the
     *        vert' local index.
     */
    const Vertex &vertex(const Element &element, int localVertIdx)
    { return *element.template subEntity<dim>(localVertIdx); };

    /*!
     * \brief Given a vert return its position in world coodinates.
     */
    void vertexPosition(GlobalPosition &worldCoord,
                        const Vertex &v) const
    {
        worldCoord = v.geometry().corner(0);
    }

    /*!
     * \brief Returns the vert map in case it is externally required.
     */
    VertexMap &vertexMap()
    { return *vertMap_; }
    const VertexMap &vertexMap() const
    { return *vertMap_; }


    /*!
     * \brief This method must be called after something within the grid has been
     *        altered.
     *
     * Cases where this method should be called are grid
     * refinements, deformation of elements, etc.
     */
    void gridChanged()
    {
        delete elementMap_;
        delete vertMap_;

        // create the element and vert index sets for the grid
        elementMap_ = new ElementMap(*grid_);
        vertMap_ = new VertexMap(grid_->leafView());
    }

    /*!
     * \brief Set the current grid.
     *
     * This method assumes that the grid was allocated on the heap
     * and it takes ownership of it, so you may not delete it
     * externally! The grid which was previously in the Domain
     * becomes invalid after calling this method.
     */
    void setGrid(Grid *grid)
    {
        grid_ = grid;
        grid->loadBalance();

        gridChanged();
    };

private:
    // pointer to the grid object
    Grid  *grid_;

    // map from the grid's leafs to an index of the index set
    ElementMap *elementMap_;
    // translates vertices local to a element to global ones
    VertexMap *vertMap_;
};

} // namespace Dune

#endif
