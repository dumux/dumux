// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*****************************************************************************
 *   See the file COPYING for full copying permissions.                      *
 *                                                                           *
 *   This program is free software: you can redistribute it and/or modify    *
 *   it under the terms of the GNU General Public License as published by    *
 *   the Free Software Foundation, either version 3 of the License, or       *
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
 * \ingroup BoxDiscretization
 * \brief Boundary types gathered on an element
 */
#ifndef DUMUX_BOX_ELEMENT_BOUNDARY_TYPES_HH
#define DUMUX_BOX_ELEMENT_BOUNDARY_TYPES_HH

#include <cassert>
#include <vector>

namespace Dumux {

/*!
 * \ingroup BoxDiscretization
 * \brief This class stores an array of BoundaryTypes objects
 */
template<class BTypes>
class BoxElementBoundaryTypes
{
public:
    using BoundaryTypes = BTypes;

    /*!
     * \brief Update the boundary types for all vertices of an element.
     *
     * \param problem The problem object which needs to be simulated
     * \param element The DUNE Codim<0> entity for which the boundary
     *                types should be collected
     * \param fvGeometry The element's finite volume geometry
     */
    template<class Problem, class Element, class FVElementGeometry>
    void update(const Problem &problem,
                const Element &element,
                const FVElementGeometry &fvGeometry)
    {
        using GridGeometry = typename FVElementGeometry::GridGeometry;
        using GridView = typename GridGeometry::GridView;

        vertexBCTypes_.resize( element.subEntities(GridView::dimension) );

        hasDirichlet_ = false;
        hasNeumann_ = false;
        hasOutflow_ = false;

        for (const auto& scv : scvs(fvGeometry))
        {
            int scvIdxLocal = scv.localDofIndex();
            vertexBCTypes_[scvIdxLocal].reset();

            if (fvGeometry.gridGeometry().dofOnBoundary(scv.dofIndex()))
            {
                vertexBCTypes_[scvIdxLocal] = problem.boundaryTypes(element, scv);

                hasDirichlet_ = hasDirichlet_ || vertexBCTypes_[scvIdxLocal].hasDirichlet();
                hasNeumann_ = hasNeumann_ || vertexBCTypes_[scvIdxLocal].hasNeumann();
                hasOutflow_ = hasOutflow_ || vertexBCTypes_[scvIdxLocal].hasOutflow();
            }
        }
    }

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

    /*
     * \brief Access operator
     * \return BoundaryTypes
     */
    const BoundaryTypes& operator[] (std::size_t i) const
    {
        assert(i < vertexBCTypes_.size());
        return vertexBCTypes_[i];
    }

protected:
    std::vector< BoundaryTypes > vertexBCTypes_;
    bool hasDirichlet_ = false;
    bool hasNeumann_ = false;
    bool hasOutflow_ = false;
};

} // namespace Dumux

#endif
