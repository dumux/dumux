// $Id$
/*****************************************************************************
 *   Copyright (C) 2010 by Andreas Lauser                                    *
 *   Institute of Hydraulic Engineering                                      *
 *   University of Stuttgart, Germany                                        *
 *   email: <givenname>.<name>@iws.uni-stuttgart.de                          *
 *                                                                           *
 *   This program is free software: you can redistribute it and/or modify    *
 *   it under the terms of the GNU General Public License as published by    *
 *   the Free Software Foundation, either version 2 of the License, or       *
 *   (at your option) any later version.                                     *
 *                                                                           *
 *   This program is distributed in the hope that it will be useful,         *
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of          *
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the           *
 *   GNU General Public License for more details.                            *
 *                                                                           *
 *   You should have received a copy of the GNU General Public License       *
 *   along with this program.  If not, see <http://www.gnu.org/licenses/>.   *
 *****************************************************************************/
/*!
 * \file
 *
 * \brief Boundary types gathered on an element
 */
#ifndef DUMUX_BOX_ELEMENT_BOUNDARY_TYPES_HH
#define DUMUX_BOX_ELEMENT_BOUNDARY_TYPES_HH

#include "boxproperties.hh"

#include <dune/grid/common/geometry.hh>

#include <dumux/common/valgrind.hh>

namespace Dumux
{

/*!
 * \ingroup BoxModel
 *
 * \brief This class stores an array of BoundaryTypes objects
 */
template<class TypeTag>
class BoxElementBoundaryTypes : public std::vector<typename GET_PROP_TYPE(TypeTag, PTAG(BoundaryTypes)) >
{
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(BoundaryTypes)) BoundaryTypes;
    typedef std::vector<BoundaryTypes> ParentType;

    typedef typename GET_PROP_TYPE(TypeTag, PTAG(Problem)) Problem;

    typedef typename GET_PROP_TYPE(TypeTag, PTAG(GridView)) GridView;
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(FVElementGeometry)) FVElementGeometry;

    enum { dim = GridView::dimension };

    typedef typename GridView::template Codim<0>::Entity Element;
    typedef typename GridView::template Codim<dim>::Entity Vertex;
    typedef typename GridView::template Codim<dim>::EntityPointer VertexPointer;
    typedef typename GridView::IntersectionIterator IntersectionIterator;

    typedef typename GridView::ctype CoordScalar;
    typedef typename Dune::GenericReferenceElements<CoordScalar, dim> ReferenceElements;
    typedef typename Dune::GenericReferenceElement<CoordScalar, dim> ReferenceElement;

public:
    /*!
     * \brief Copy constructor.
     *
     * Copying a the boundary types of an element should be explicitly
     * requested
     */
    explicit BoxElementBoundaryTypes(const BoxElementBoundaryTypes &v)
        : ParentType(v)
    {}

    /*!
     * \brief Default constructor.
     */
    BoxElementBoundaryTypes()
    {
        hasDirichlet_ = false;
        hasNeumann_ = false;
        hasOutflow_ = false;
    }

    /*!
     * \brief Update the boundary types for all vertices of an element.
     *
     * \param problem The problem object which needs to be simulated
     * \param element The DUNE Codim<0> entity for which the boundary
     *                types should be collected
     * \param fvElemGeom The element's finite volume geometry
     */
    void update(const Problem &problem,
                const Element &element,
                const FVElementGeometry &fvElemGeom)
    {
        int numVerts = element.template count<dim>();
        this->resize(numVerts);

        hasDirichlet_ = false;
        hasNeumann_ = false;
        hasOutflow_ = false;
        for (int i = 0; i < numVerts; ++i) {
            (*this)[i].reset();

            if (!problem.model().onBoundary(element, i))
                continue;

            const VertexPointer vptr = element.template subEntity<dim>(i);
            problem.boundaryTypes((*this)[i], *vptr);

            hasDirichlet_ = hasDirichlet_ or (*this)[i].hasDirichlet();
            hasNeumann_ = hasNeumann_ or (*this)[i].hasNeumann();
            hasOutflow_ = hasOutflow_ or (*this)[i].hasOutflow();
        }
    };

    /*!
     * \brief Returns whether the element has a vertex which contains
     *        a Dirichlet value.
     */
    bool hasDirichlet() const
    { return hasDirichlet_; }

    /*!
     * \brief Returns whether the element potentially features a
     *        Neumann boundary segment.
     */
    bool hasNeumann() const
    { return hasNeumann_; }

    /*!
     * \brief Returns whether the element potentially features an
     *        outflow boundary segment.
     */
    bool hasOutflow() const
    { return hasOutflow_; }

protected:
    bool hasDirichlet_;
    bool hasNeumann_;
    bool hasOutflow_;
};

} // namespace Dumux

#endif
