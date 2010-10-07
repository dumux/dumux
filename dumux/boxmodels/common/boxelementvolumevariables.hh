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
#ifndef DUMUX_BOX_ELEMENT_VOLUME_VARIABLES_HH
#define DUMUX_BOX_ELEMENT_VOLUME_VARIABLES_HH

#include "boxproperties.hh"


namespace Dumux
{

/*!
 * \ingroup BoxModel
 *
 * \brief This class stores an array of VolumeVariables objects
 */
template<class TypeTag>
class BoxElementVolumeVariables : public std::vector<typename GET_PROP_TYPE(TypeTag, PTAG(VolumeVariables)) >
{
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(VolumeVariables)) VolumeVariables;
    typedef std::vector<VolumeVariables> ParentType;

    typedef typename GET_PROP_TYPE(TypeTag, PTAG(Scalar)) Scalar;
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(Problem)) Problem;
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(FVElementGeometry)) FVElementGeometry;
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(SolutionVector)) SolutionVector;
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(VertexMapper)) VertexMapper;
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(PrimaryVariables)) PrimaryVariables;

    typedef typename GET_PROP_TYPE(TypeTag, PTAG(GridView)) GridView;
    typedef typename GridView::template Codim<0>::Entity Element;
    typedef typename GridView::IntersectionIterator IntersectionIterator;

    enum { dim = GridView::dimension };

public:
    /*!
     * \brief The constructor.
     */
    BoxElementVolumeVariables()
    { }
    /*!
     * \todo please doc me
     */
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
            const PrimaryVariables &solI
                = globalSol[vertexMapper.map(element, i, dim)];
            (*this)[i].update(solI,
                              problem,
                              element,
                              fvElemGeom,
                              i,
                              oldSol);
        }
    };

    //overloaded update function
    //possible to give solution vector
    template<typename SolVectorType>
    void update(const Problem &problem,
                    const Element &element,
                    const FVElementGeometry &fvElemGeom,
                    const SolVectorType& elementSolVector,
                    const int numEq)
        {
            int n = element.template count<dim>();
            this->resize(n);
            for (int vertexIdx= 0; vertexIdx < n; vertexIdx++)
            {
                PrimaryVariables solI(0);
                for (int eqnIdx=0; eqnIdx<numEq; eqnIdx++)
                {
                    solI[eqnIdx] = elementSolVector[vertexIdx+eqnIdx*n];
                }
                    (*this)[vertexIdx].update(solI,
                                      problem,
                                      element,
                                      fvElemGeom,
                                      vertexIdx,
                                      false);

            }
        };
};

} // namespace Dumux

#endif
