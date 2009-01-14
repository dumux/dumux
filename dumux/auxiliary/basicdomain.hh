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
 * \brief Provides some convenient interfaces to elements and vertices of grids.
 */
#ifndef DUMUX_BASIC_DOMAIN_HH
#define DUMUX_BASIC_DOMAIN_HH

#include <dumux/auxiliary/apis.hh>

#include <dune/common/fvector.hh>
#include <dune/common/fmatrix.hh>

#include <dune/grid/common/referenceelements.hh>
#include <dune/grid/common/intersectioniterator.hh>
#include <dune/grid/io/file/dgfparser.hh>
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
            i = I::dim;
            i = I::dimWorld;

            // make sure the required types exist
            {typename I::Grid                       *x;x=NULL;}
            {typename I::Scalar                     *x;x=NULL;}
            {typename I::CoordScalar                *x;x=NULL;}
            {typename I::LocalPosition                 *x;x=NULL;}
            {typename I::GlobalPosition                 *x;x=NULL;}
            {typename I::Element                    *x;x=NULL;}
            {typename I::ElementIterator            *x;x=NULL;}
            {typename I::Vertex                     *x;x=NULL;}
            {typename I::VertexIterator             *x;x=NULL;}
            {typename I::ReferenceElement           *x;x=NULL;}
            {typename I::ReferenceElementContainer  *x;x=NULL;}
            {typename I::IntersectionIterator       *x;x=NULL;}
            {typename I::FieldVector                *x;x=NULL;}
            {typename I::FieldMatrix                *x;x=NULL;}

            {const typename I::ReferenceElementContainer *x = &I::referenceElement; x=NULL;}
        }
        END_API_DEF
    };

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
            typedef typename Grid::Traits::template Codim<0>                  ElementTraits;
            typedef typename ElementTraits::Entity                            Element;
            typedef typename ElementTraits::LeafIterator                      ElementIterator;

            typedef typename Grid::Traits::template Codim<dim>                VertexTraits;
            typedef typename VertexTraits::Entity                             Vertex;
            typedef typename VertexTraits::LeafIterator                       VertexIterator;

            // TODO: Dune::ReferenceElement uses virtual functions in
            //       order to support arbitary element types. It would be
            //       better to use ReferenceCubeContainer,
            //       ReferenceSimplexContainer,
            //       ReferencePrismContainer or
            //       ReferencePyramidContainer if the Grid only uses
            //       one element type. We need to investigate how to do
            //       this, but we probably need to beef up
            //       Dune::Capabilities a bit and use specialization
            //       for the DomainTraits...
            typedef Dune::ReferenceElementContainer<CoordScalar, dim>         ReferenceElementContainer;
            typedef typename ReferenceElementContainer::value_type            ReferenceElement;
            static const ReferenceElementContainer &referenceElement;

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
        typedef typename DomainTraits::ReferenceElement      ReferenceElement;

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
        typedef Dune::MultipleCodimMultipleGeomTypeMapper<Grid,ElementIdxSet_,VertexLayout_> VertexMap;

    public:
        BasicDomain()
            {
                Api::require<Api::BasicDomainTraits, DomainTraits>();

                elementMap_ = NULL;
                vertMap_ = NULL;
            };

        BasicDomain(Grid *grid)
            {
                Api::require<Api::BasicDomainTraits, DomainTraits>();

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
/*        const Grid &grid()
            { return *grid_; }
*/

        /*!
         * \brief Returns the current grid.
         */
        const Grid &grid() const
            { return *grid_; }

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
            return grid_->template leafbegin<0>();
        }

        /*!
         * \brief Returns the iterator pointing to the next to last
         *        element of the grid.
         */
        ElementIterator elementEnd() {
            return grid_->template leafend<0>();
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
        int vertIdx(const Element &e, int localVertexId) const
            { return vertMap_->template map<dim>(e, localVertexId); }

        /*!
         * \brief Given a vert index within a element, return the
         *        respective global index.
         */
        int vertIdx(const Vertex &v) const
            { return vertMap_->map(v); }

        /*!
         * \brief Returns the iterator pointing to the first element of
         *        the grid.
         */
        VertexIterator vertexBegin() {
            return grid_->template leafbegin<dim>();
        }

        /*!
         * \brief Returns the iterator pointing to the next to last
         *        element of the grid.
         */
        VertexIterator vertexEnd() {
            return grid_->template leafend<dim>();
        }

        /*!
         * \brief Return the entity for a vert given a element and the
         *        vert' local index.
         */
        const Vertex &vert(const Element &element, int localVertIdx)
            { return *element.template entity<dim>(localVertIdx); };

        /*!
         * \brief Given a vert return its position in world coodinates.
         */
        void vertPosition(GlobalPosition &worldCoord,
                          const Vertex &v) const
            {
                worldCoord = v.geometry().corner(0);
            }

        /*!
         * \brief Returns the vert map in case it is externally required.
         */
        VertexMap &vertMap()
            { return *vertMap_; }
        const VertexMap &vertMap() const
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
                vertMap_ = new VertexMap(*grid_, grid_->leafIndexSet());
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


    // this is butt-ugly, but the only way I could come up with to
    // initialize a static const member of a class
    template<class Grid, class ScalarT>
    const typename BasicDomain<Grid, ScalarT>::DomainTraits::ReferenceElementContainer &
         BasicDomain<Grid, ScalarT>::DomainTraits::referenceElement
         =  Dune::ReferenceElements<typename Grid::ctype,
                                    Grid::dimension>::general;

}

#endif
