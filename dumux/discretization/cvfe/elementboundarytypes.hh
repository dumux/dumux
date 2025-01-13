// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
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

#include <dumux/common/typetraits/localdofs.hh>
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

        for (const auto& localDof : cvLocalDofs(fvGeometry))
        {
            const auto localIdx = localDof.indexInElement();
            bcTypes_[localIdx].reset();

            if (fvGeometry.gridGeometry().dofOnBoundary(localDof.dofIndex()))
            {
                bcTypes_[localIdx] = problem.boundaryTypes(element, localDof.scv());
                hasDirichlet_ = hasDirichlet_ || bcTypes_[localIdx].hasDirichlet();
                hasNeumann_ = hasNeumann_ || bcTypes_[localIdx].hasNeumann();
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
     * \note yields undefined behaviour if the scv is not on the boundary
     */
    template<class FVElementGeometry,
             class ScvOrLocalDof>
    const BoundaryTypes& get(const FVElementGeometry&, const ScvOrLocalDof& scvOrLocalDof) const
    {
        const auto localDofIdx = scvOrLocalDof.indexInElement();
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
