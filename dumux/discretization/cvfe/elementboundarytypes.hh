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
 * \ingroup CVFEDiscretization
 * \brief Boundary types gathered on an element
 */
#ifndef DUMUX_CVFE_ELEMENT_BOUNDARY_TYPES_HH
#define DUMUX_CVFE_ELEMENT_BOUNDARY_TYPES_HH

#include <cassert>
#include <vector>

#include <dumux/discretization/method.hh>

namespace Dumux {

/*!
 * \ingroup CVFEDiscretization
 * \brief This class stores an array of BoundaryTypes objects
 */
template<class BTypes>
class CVFEElementBoundaryTypes
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
    template<class Problem, class FVElementGeometry>
    void update(const Problem& problem,
                const typename FVElementGeometry::Element& element,
                const FVElementGeometry& fvGeometry)
    {
        using GridGeometry = typename FVElementGeometry::GridGeometry;

        if constexpr (GridGeometry::discMethod == DiscretizationMethods::fcdiamond)
            if (!fvGeometry.hasBoundaryScvf())
                return;

        bcTypes_.resize(fvGeometry.numScv());

        hasDirichlet_ = false;
        hasNeumann_ = false;

        // boundary dofs always have boundary corresponding faces
        if constexpr (GridGeometry::discMethod == DiscretizationMethods::fcdiamond)
        {
            for (const auto& scvf : scvfs(fvGeometry))
            {
                if (scvf.boundary())
                {
                    const auto localIndex = fvGeometry.scv(scvf.insideScvIdx()).localDofIndex();
                    bcTypes_[localIndex] = problem.boundaryTypes(element, scvf);
                    hasDirichlet_ = hasDirichlet_ || bcTypes_[localIndex].hasDirichlet();
                    hasNeumann_ = hasNeumann_ || bcTypes_[localIndex].hasNeumann();
                }
            }
        }

        // generally, single vertices maybe on the boundary even if the face is not
        // in this case we need a different mechanism
        else
        {
            for (const auto& scv : scvs(fvGeometry))
            {
                const auto scvIdxLocal = scv.localDofIndex();
                bcTypes_[scvIdxLocal].reset();

                if (fvGeometry.gridGeometry().dofOnBoundary(scv.dofIndex()))
                {
                    bcTypes_[scvIdxLocal] = problem.boundaryTypes(element, scv);
                    hasDirichlet_ = hasDirichlet_ || bcTypes_[scvIdxLocal].hasDirichlet();
                    hasNeumann_ = hasNeumann_ || bcTypes_[scvIdxLocal].hasNeumann();
                }
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

    /*
     * \brief Access operator
     * \return BoundaryTypes
     * \note yields undefined behaviour of the scv is not on the boundary
     */
    template<class FVElementGeometry>
    const BoundaryTypes& get(const FVElementGeometry& fvGeometry, const typename FVElementGeometry::SubControlVolumeFace& scvf) const
    {
        assert(scvf.boundary());
        const auto localDofIdx = fvGeometry.scv(scvf.insideScvIdx()).localDofIndex();
        assert(localDofIdx < bcTypes_.size());
        return bcTypes_[localDofIdx];
    }

    /*
     * \brief Access operator
     * \return BoundaryTypes
     * \note yields undefined behaviour of the scv is not on the boundary
     */
    template<class FVElementGeometry>
    const BoundaryTypes& get(const FVElementGeometry&, const typename FVElementGeometry::SubControlVolume& scv) const
    {
        const auto localDofIdx = scv.localDofIndex();
        assert(localDofIdx < bcTypes_.size());
        return bcTypes_[localDofIdx];
    }

private:
    std::vector<BoundaryTypes> bcTypes_;
    bool hasDirichlet_ = false;
    bool hasNeumann_ = false;
};

} // namespace Dumux

#endif
