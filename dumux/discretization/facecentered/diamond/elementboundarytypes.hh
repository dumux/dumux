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
 * \ingroup DiamondDiscretization
 * \brief Boundary types gathered on an element
 */
#ifndef DUMUX_DISCRETIZATION_FACECENTERED_DIAMOND_ELEMENT_BOUNDARY_TYPES_HH
#define DUMUX_DISCRETIZATION_FACECENTERED_DIAMOND_ELEMENT_BOUNDARY_TYPES_HH

#include <vector>

namespace Dumux {

/*!
 * \ingroup DiamondDiscretization
 * \brief This class stores an array of BoundaryTypes objects
 */
template<class BTypes, class FVElementGeometry>
class FaceCenteredDiamondElementBoundaryTypes
{
    using SubControlVolume = typename FVElementGeometry::SubControlVolume;
    using SubControlVolumeFace = typename FVElementGeometry::SubControlVolumeFace;
public:
    using BoundaryTypes = BTypes;

    FaceCenteredDiamondElementBoundaryTypes() = default;

    /*!
     * \brief Update the boundary types for all faces of an element.
     *
     * \param problem The problem object which needs to be simulated
     * \param element The DUNE Codim<0> entity for which the boundary
     *                types should be collected
     * \param fvGeometry The element's finite volume geometry
     */
    template<class Problem>
    void update(const Problem& problem,
                const typename FVElementGeometry::Element& element,
                const FVElementGeometry& fvGeometry)
    {
        if (!fvGeometry.hasBoundaryScvf())
            return;

        bcTypes_.resize(fvGeometry.numScv());

        hasDirichlet_ = false;
        hasNeumann_ = false;

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

    /*!
     * \brief Returns whether the element has a face which contains
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
    const BoundaryTypes& get(const FVElementGeometry& fvGeometry, const SubControlVolumeFace& scvf) const
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
    const BoundaryTypes& get(const FVElementGeometry&, const SubControlVolume& scv) const
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

} // end namespace Dumux

#endif
