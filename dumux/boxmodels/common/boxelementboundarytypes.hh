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

#include <vector>
#include <dumux/common/boundarytypes.hh>

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
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(Scalar)) Scalar;

    typedef typename GET_PROP_TYPE(TypeTag, PTAG(GridView)) GridView;
    typedef typename GridView::template Codim<0>::Entity Element;
    typedef typename GridView::IntersectionIterator IntersectionIterator;

    typedef typename GET_PROP(TypeTag, PTAG(ReferenceElements)) RefElemProp;
    typedef typename RefElemProp::Container ReferenceElements;
    typedef typename RefElemProp::ReferenceElement ReferenceElement;

    typedef typename GET_PROP_TYPE(TypeTag, PTAG(FVElementGeometry)) FVElementGeometry;

    enum { dim = GridView::dimension };

public:
    /*!
     * \brief The constructor.
     */
    BoxElementBoundaryTypes()
    { }

    void update(const Problem &problem, 
                const Element &element, 
                const FVElementGeometry &fvElemGeom)
    {
        Dune::GeometryType      geoType = element.geometry().type();
        const ReferenceElement &refElem = ReferenceElements::general(geoType);
        
        int numVerts = element.template count<dim>();
        this->resize(numVerts);
        for (int i = 0; i < numVerts; ++i)
            (*this)[i].reset();
    
        // evaluate boundary conditions
        IntersectionIterator isIt = problem.gridView().template ibegin(element);
        const IntersectionIterator &endIt = problem.gridView().template iend(element);
        for (; isIt != endIt; ++isIt) {
            // Ignore non- boundary faces.
            if (!isIt->boundary())
                continue;
            
            // Set the boundary type for all vertices of the face
            int faceIdx = isIt->indexInInside();
            int numFaceVerts = refElem.size(faceIdx, 1, dim);
            for (int faceVertIdx = 0;
                 faceVertIdx < numFaceVerts;
                 faceVertIdx++)
            {
                int elemVertIdx = refElem.subEntity(faceIdx,
                                                    1,
                                                    faceVertIdx,
                                                    dim);
                int boundaryFaceIdx =
                    fvElemGeom.boundaryFaceIndex(faceIdx,
                                                 faceVertIdx);
                // set the boundary types
                problem.boundaryTypes((*this)[elemVertIdx],
                                      element,
                                      fvElemGeom,
                                      *isIt,
                                      elemVertIdx,
                                      boundaryFaceIdx);
                (*this)[elemVertIdx].checkWellPosed();
                Valgrind::CheckDefined((*this)[elemVertIdx]);
            }
        }
    };
};

} // namespace Dumux

#endif
