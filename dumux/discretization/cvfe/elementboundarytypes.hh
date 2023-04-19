// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
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
        bcTypes_.resize(fvGeometry.numScv());

        hasDirichlet_ = false;
        hasNeumann_ = false;

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
