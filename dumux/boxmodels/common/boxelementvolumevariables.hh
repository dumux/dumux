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
 * \brief Volume variables gathered on an element
 */
#ifndef DUMUX_BOX_ELEMENT_VOLUME_VARIABLES_HH
#define DUMUX_BOX_ELEMENT_VOLUME_VARIABLES_HH

#include "boxproperties.hh"


namespace Dumux
{

/*!
 * \ingroup BoxModel
 *
 * \brief This class stores an array of VolumeVariables objects, one
 *        volume variables object for each of the element's vertices
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
    enum { numEq = GET_PROP_VALUE(TypeTag, PTAG(NumEq)) };

public:
    /*!
     * \brief The constructor.
     */
    BoxElementVolumeVariables()
    { }

    /*!
     * \brief Construct the volume variables for all of vertices of an element.
     *
     * \param problem The problem which needs to be simulated.
     * \param element The DUNE Codim<0> entity for which the volume variables ought to be calculated
     * \param fvElemGeom The finite volume geometry of the element
     * \param oldSol Tells whether the model's previous or current solution should be used.
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

    /*!
     * \brief Construct the volume variables for all of vertices of an
     *        element given a solution vector computed by PDELab.
     *
     * \tparam ElemSolVectorType The container type which stores the
     *                           primary variables of the element
     *                           using _local_ indices
     *
     * \param problem The problem which needs to be simulated.
     * \param element The DUNE Codim<0> entity for which the volume variables ought to be calculated
     * \param fvElemGeom The finite volume geometry of the element
     * \param elementSolVector The local solution for the element using PDELab ordering
     */
    template<typename ElemSolVectorType>
    void updatePDELab(const Problem &problem,
                      const Element &element,
                      const FVElementGeometry &fvElemGeom,
                      const ElemSolVectorType& elementSolVector)
    {
        int n = element.template count<dim>();
        this->resize(n);
        for (int vertexIdx = 0; vertexIdx < n; vertexIdx++)
        {
            PrimaryVariables solI(0);
            for (int eqnIdx=0; eqnIdx<numEq; eqnIdx++)
                solI[eqnIdx] = elementSolVector[vertexIdx + eqnIdx*n];
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
