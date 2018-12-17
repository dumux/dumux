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
 * \ingroup FacetCoupling
 * \brief Boundary types gathered on an element for the box scheme
 *        with coupling occurring across the element facets.
 */
#ifndef DUMUX_FACET_COUPLING_BOX_ELEMENT_BOUNDARY_TYPES_HH
#define DUMUX_FACET_COUPLING_BOX_ELEMENT_BOUNDARY_TYPES_HH

#include <cassert>
#include <vector>

#include <dumux/discretization/box/elementboundarytypes.hh>

namespace Dumux {

/*!
 * \ingroup FacetCoupling
 * \brief This class stores an array of BoundaryTypes objects
 *        on an element for the box scheme with coupling occuring
 *        across the element facets.
 */
template<class BTypes>
class BoxFacetCouplingElementBoundaryTypes : public BoxElementBoundaryTypes<BTypes>
{

public:
    /*!
     * \brief Update the boundary types for all vertices of an element.
     *
     * \param problem The problem object which needs to be simulated
     * \param element The DUNE Codim<0> entity for which the boundary
     *                types should be collected
     * \param fvGeometry The element's finite volume geometry
     *
     * \note We need this overload so that we can call different
     *       problem interfaces for domain boundary and interior boundaries.
     */
    template<class Problem, class Element, class FVElementGeometry>
    void update(const Problem &problem,
                const Element &element,
                const FVElementGeometry &fvGeometry)
    {
        using FVGridGeometry = typename FVElementGeometry::FVGridGeometry;
        using GridView = typename FVGridGeometry::GridView;

        this->vertexBCTypes_.resize( element.subEntities(GridView::dimension) );

        this->hasDirichlet_ = false;
        this->hasNeumann_ = false;
        this->hasOutflow_ = false;

        for (const auto& scv : scvs(fvGeometry))
        {
            int scvIdxLocal = scv.localDofIndex();
            this->vertexBCTypes_[scvIdxLocal].reset();

            // lambda to update the element boundary info
            auto updateElemBCInfo = [&] (const BTypes& bcTypes)
            {
                this->hasDirichlet_ = this->hasDirichlet_ || this->vertexBCTypes_[scvIdxLocal].hasDirichlet();
                this->hasNeumann_ = this->hasNeumann_ || this->vertexBCTypes_[scvIdxLocal].hasNeumann();
                this->hasOutflow_ = this->hasOutflow_ || this->vertexBCTypes_[scvIdxLocal].hasOutflow();
            };

            // We have to have Neumann-type BCs on nodes that touch interior boundaries
            if (fvGeometry.fvGridGeometry().dofOnInteriorBoundary(scv.dofIndex()))
            {
                this->vertexBCTypes_[scvIdxLocal].setAllNeumann();
                updateElemBCInfo(this->vertexBCTypes_[scvIdxLocal]);
            }

            // otherwise, let the problem decide
            else if (fvGeometry.fvGridGeometry().dofOnBoundary(scv.dofIndex()))
            {
                this->vertexBCTypes_[scvIdxLocal] = problem.boundaryTypes(element, scv);
                updateElemBCInfo(this->vertexBCTypes_[scvIdxLocal]);
            }
        }
    }
};

} // namespace Dumux

#endif
