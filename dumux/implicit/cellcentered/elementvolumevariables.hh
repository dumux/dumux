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
#ifndef DUMUX_CC_ELEMENT_VOLUME_VARIABLES_HH
#define DUMUX_CC_ELEMENT_VOLUME_VARIABLES_HH

#include "properties.hh"

namespace Dumux
{

/*!
 * \ingroup CCModel
 * \brief This class stores an array of VolumeVariables objects, one
 *        volume variables object for each of the element's vertices
 */
template<class TypeTag>
class CCElementVolumeVariables
{
    typedef typename GET_PROP_TYPE(TypeTag, Problem) Problem;
    typedef typename GET_PROP_TYPE(TypeTag, FVElementGeometry) FVElementGeometry;
    typedef typename GET_PROP_TYPE(TypeTag, SolutionVector) SolutionVector;
    typedef typename GET_PROP_TYPE(TypeTag, PrimaryVariables) PrimaryVariables;
    typedef typename GET_PROP_TYPE(TypeTag, VolumeVariables) VolumeVariables;
    typedef typename GET_PROP_TYPE(TypeTag, ElementBoundaryTypes) ElementBoundaryTypes;
    typedef typename GET_PROP_TYPE(TypeTag, BoundaryTypes) BoundaryTypes;

    typedef typename GET_PROP_TYPE(TypeTag, GridView) GridView;
    typedef typename GridView::template Codim<0>::Entity Element;

public:
    /*!
     * \brief Construct the volume variables for all of vertices of an element.
     *
     * \param problem The problem which needs to be simulated.
     * \param element The DUNE Codim<0> entity for which the volume variables ought to be calculated
     * \param fvGeometry The finite volume geometry of the element
     * \param oldSol Tells whether the model's previous or current solution should be used.
     */
    void update(Problem &problem,
                const Element &element,
                const FVElementGeometry &fvGeometry,
                bool oldSol)
    {
        oldSol_ = oldSol;
        problem_ = &problem;

        int numNeighbors = fvGeometry.numNeighbors;
        stencil_.resize(numNeighbors);

        for (int i = 0; i < numNeighbors; i++)
        {
            const Element& neighbor = fvGeometry.neighbors[i];
            stencil_[i] = problem.elementMapper().index(neighbor);
        }

        // only treat boundary if current solution is evaluated
        if (!oldSol)
        {
            // check if element intersects with the boundary
            ElementBoundaryTypes elemBCTypes;
            elemBCTypes.update(problem, element);
            if (elemBCTypes.hasDirichlet()
                || elemBCTypes.hasNeumann()
                || elemBCTypes.hasOutflow())
            {
                ghostVolVars_ = std::vector<VolumeVariables>(element.subEntities(1));
                isDirichlet_ = std::vector<bool>(element.subEntities(1), false);

                // add volume variables for the boundary faces
                for (const auto& intersection : intersections(problem.gridView(), element))
                {
                    if (!intersection.boundary())
                        continue;

                    BoundaryTypes bcTypes;
                    problem.boundaryTypes(bcTypes, intersection);

                    int fIdx = intersection.indexInInside();
                    if (bcTypes.hasDirichlet())
                    {
                        PrimaryVariables dirichletValues;
                        problem.dirichlet(dirichletValues, intersection);

                        ghostVolVars_[fIdx].update(dirichletValues,
                                                    problem,
                                                    element,
                                                    fvGeometry,
                                                    /*scvIdx=*/0,
                                                    oldSol);
                        isDirichlet_[fIdx] = true;
                    }
                }
            }
        }
    }

    const VolumeVariables& operator [](const unsigned int i) const
    {
        if (i >= stencil_.size())
        {
            assert(!oldSol_); // old boundary volVars don't exist
            if (isDirichlet_[i - stencil_.size()]) // ghost volVars only exist for Dirichlet boundaries
                return ghostVolVars_[i - stencil_.size()];
            else
                return problem_->model().volVars(stencil_[0]);
        }
        else
        {
            if(!oldSol_)
                return problem_->model().volVars(stencil_[i]);
            else
                return problem_->model().prevVolVars(stencil_[i]);
        }
    }

    VolumeVariables& operator [](const unsigned int i)
    {
        assert(!oldSol_); // old volVars should never be modified
        assert(i < stencil_.size()); // boundary volVars should not be modified
        return problem_->model().volVars(stencil_[i]);
    }

private:
    bool oldSol_ = false;
    std::vector<unsigned int> stencil_;
    std::vector<VolumeVariables> ghostVolVars_;
    std::vector<bool> isDirichlet_;
    Problem* problem_;
};

} // namespace Dumux

#endif
