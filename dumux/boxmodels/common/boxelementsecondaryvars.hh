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
#ifndef DUMUX_BOX_ELEMENT_SECONDARY_VARS_HH
#define DUMUX_BOX_ELEMENT_SECONDARY_VARS_HH

#include "boxproperties.hh"

#include <vector>
#include <dumux/common/boundarytypes.hh>

namespace Dumux
{

/*!
 * \ingroup BoxModel
 *
 * \brief This class stores an array of SecondaryVars objects
 */
template<class TypeTag>
class BoxElementSecondaryVars : public std::vector<typename GET_PROP_TYPE(TypeTag, PTAG(SecondaryVars)) >
{
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(SecondaryVars)) SecondaryVars;
    typedef std::vector<SecondaryVars> ParentType;

    typedef typename GET_PROP_TYPE(TypeTag, PTAG(Problem)) Problem;
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(Scalar)) Scalar;
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(FVElementGeometry)) FVElementGeometry;
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(SolutionVector)) SolutionVector;
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(VertexMapper)) VertexMapper;
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(PrimaryVarVector)) PrimaryVarVector;

    typedef typename GET_PROP_TYPE(TypeTag, PTAG(GridView)) GridView;
    typedef typename GridView::template Codim<0>::Entity Element;
    typedef typename GridView::IntersectionIterator IntersectionIterator;

    enum { dim = GridView::dimension };

public:
    /*!
     * \brief The constructor.
     */
    BoxElementSecondaryVars()
    { }

    void update(const Problem &problem, 
                const Element &element, 
                const FVElementGeometry &fvElemGeom,
                bool oldSol)
    {
        const SolutionVector &globalSol = 
            oldSol?
            problem.model().prevSol():
            problem.model().curSol();
        const VertexMapper &vertexMapper = problem.vertexMapper();
        // we assert that the i-th shape function is
        // associated to the i-th vert of the element.
        int n = element.template count<dim>();
        this->resize(n);       
        for (int i = 0; i < n; i++) {
            const PrimaryVarVector &solI
                = globalSol[vertexMapper.map(element, i, dim)];
            (*this)[i].update(solI,
                              problem,
                              element,
                              fvElemGeom,
                              i,
                              oldSol);
        }
    };
};

} // namespace Dumux

#endif
