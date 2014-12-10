// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*****************************************************************************
 *   See the file COPYING for full copying permissions.                      *
 *                                                                           *
 *   This program is free software: you can redistribute it and/or modify    *
 *   it under the terms of the GNU General Public License as published by    *
 *   the Free Software Foundation, either version 2 of the License, or       *
 *   (at your option) any later version.                                     *
 *                                                                           *
 *   This program is distributed in the hope that it will be useful,         *
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of          *
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the            *
 *   GNU General Public License for more details.                            *
 *                                                                           *
 *   You should have received a copy of the GNU General Public License       *
 *   along with this program.  If not, see <http://www.gnu.org/licenses/>.   *
 *****************************************************************************/
/*!
 * \file
 * \brief Volume variables gathered on an element
 */
#ifndef DUMUX_BOX_ELEMENT_VOLUME_VARIABLES_HH
#define DUMUX_BOX_ELEMENT_VOLUME_VARIABLES_HH

#include <dune/common/version.hh>

#include "boxproperties.hh"


namespace Dumux
{

/*!
 * \ingroup BoxModel
 * \brief This class stores an array of VolumeVariables objects, one
 *        volume variables object for each of the element's vertices
 */
template<class TypeTag>
class BoxElementVolumeVariables : public std::vector<typename GET_PROP_TYPE(TypeTag, VolumeVariables) >
{
    typedef typename GET_PROP_TYPE(TypeTag, Problem) Problem;
    typedef typename GET_PROP_TYPE(TypeTag, FVElementGeometry) FVElementGeometry;
    typedef typename GET_PROP_TYPE(TypeTag, SolutionVector) SolutionVector;
    typedef typename GET_PROP_TYPE(TypeTag, VertexMapper) VertexMapper;
    typedef typename GET_PROP_TYPE(TypeTag, PrimaryVariables) PrimaryVariables;

    typedef typename GET_PROP_TYPE(TypeTag, GridView) GridView;
    typedef typename GridView::template Codim<0>::Entity Element;

    enum { dim = GridView::dimension };
    enum { numEq = GET_PROP_VALUE(TypeTag, NumEq) };

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
     * \param fvGeometry The finite volume geometry of the element
     * \param oldSol Tells whether the model's previous or current solution should be used.
     */
    void update(const Problem &problem,
                const Element &element,
                const FVElementGeometry &fvGeometry,
                const bool oldSol)
    {
        const SolutionVector &globalSol =
            oldSol?
            problem.model().prevSol():
            problem.model().curSol();
        const VertexMapper &vertexMapper = problem.vertexMapper();
        // we assert that the i-th shape function is
        // associated to the i-th vertex of the element.
#if DUNE_VERSION_NEWER(DUNE_COMMON, 2, 4)
        int numVertices = element.subEntities(dim);
#else
        int numVertices = element.template count<dim>();
#endif
        this->resize(numVertices);
        for (int scvIdx = 0; scvIdx < numVertices; scvIdx++) {
#if DUNE_VERSION_NEWER(DUNE_COMMON, 2, 4)
            const PrimaryVariables &priVars
                = globalSol[vertexMapper.subIndex(element, scvIdx, dim)];
#else
            const PrimaryVariables &priVars
                = globalSol[vertexMapper.map(element, scvIdx, dim)];
#endif

            // reset evaluation point to zero
            (*this)[scvIdx].setEvalPoint(0);

            (*this)[scvIdx].update(priVars,
                              problem,
                              element,
                              fvGeometry,
                              scvIdx,
                              oldSol);
        }
    }

    /*!
     * \brief Construct the volume variables for all of vertices of an
     *        element given a solution vector computed by PDELab.
     *
     * \tparam ElementSolutionVector The container type which stores the
     *                           primary variables of the element
     *                           using _local_ indices
     *
     * \param problem The problem which needs to be simulated.
     * \param element The DUNE Codim<0> entity for which the volume variables ought to be calculated
     * \param fvGeometry The finite volume geometry of the element
     * \param elementSolVector The local solution for the element using PDELab ordering
     */
    template<typename ElementSolutionVector>
    void updatePDELab(const Problem &problem,
                      const Element &element,
                      const FVElementGeometry &fvGeometry,
                      const ElementSolutionVector& elementSolVector)
    {
#if DUNE_VERSION_NEWER(DUNE_COMMON, 2, 4)
        int numVertices = element.subEntities(dim);
#else
        int numVertices = element.template count<dim>();
#endif
        this->resize(numVertices);
        for (int scvIdx = 0; scvIdx < numVertices; scvIdx++)
        {
            PrimaryVariables priVars(0);
            for (int eqnIdx = 0; eqnIdx < numEq; eqnIdx++)
                priVars[eqnIdx] = elementSolVector[scvIdx + eqnIdx*numVertices];
            (*this)[scvIdx].update(priVars,
                                      problem,
                                      element,
                                      fvGeometry,
                                      scvIdx,
                                      false);

        }
    }
};

} // namespace Dumux

#endif
