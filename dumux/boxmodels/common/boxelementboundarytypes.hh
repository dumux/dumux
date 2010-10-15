// $Id$
/*****************************************************************************
 *   Copyright (C) 2010 by Andreas Lauser                                    *
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

    typedef typename GET_PROP(TypeTag, PTAG(ReferenceElements)) RefElemProp;
    typedef typename RefElemProp::Container ReferenceElements;
    typedef typename RefElemProp::ReferenceElement ReferenceElement;

    typedef typename GET_PROP_TYPE(TypeTag, PTAG(FVElementGeometry)) FVElementGeometry;

    enum { dim = GridView::dimension };

    typedef typename GridView::template Codim<0>::Entity Element;
    typedef typename GridView::template Codim<dim>::Entity Vertex;
    typedef typename GridView::template Codim<dim>::EntityPointer VertexPointer;
    typedef typename GridView::IntersectionIterator IntersectionIterator;

public:
    // copying a the boundary types of an element should be explicitly
    // requested
    explicit BoxElementBoundaryTypes(const BoxElementBoundaryTypes &v)
        : ParentType(v)
    {}

    /*!
     * \brief The constructor.
     */
    BoxElementBoundaryTypes()
    {
        hasDirichlet_ = false;
        hasNeumann_ = false;
    }

    /*!
     * \brief Update the boundary types for all vertices of an element.
     */
    void update(const Problem &problem,
                const Element &element,
                const FVElementGeometry &fvElemGeom)
    {
        int numVerts = element.template count<dim>();
        this->resize(numVerts);

        hasDirichlet_ = false;
        hasNeumann_ = false;
        for (int i = 0; i < numVerts; ++i) {
            (*this)[i].reset();

            if (!problem.model().onBoundary(element, i))
                continue;

            const VertexPointer vptr = element.template subEntity<dim>(i);
            problem.boundaryTypes((*this)[i], *vptr);

            hasDirichlet_ = hasDirichlet_ or (*this)[i].hasDirichlet();
            hasNeumann_ = hasNeumann_ or (*this)[i].hasNeumann();
        }
    };

    /*!
     * \brief Returns whether the element has a Dirichlet value.
     */
    bool hasDirichlet() const
    { return hasDirichlet_; }

    /*!
     * \brief Returns whether the element potentially has a Neumann
     *        boundary segment.
     */
    bool hasNeumann() const
    { return hasNeumann_; }
    
protected:
    bool hasDirichlet_;
    bool hasNeumann_;
};

} // namespace Dumux

#endif
