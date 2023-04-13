// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup FaceCenteredStaggeredDiscretization
 * \brief Boundary types gathered on an element
 */
#ifndef DUMUX_DISCRETIZATION_FACECENTERED_STAGGERED_ELEMENT_BOUNDARY_TYPES_HH
#define DUMUX_DISCRETIZATION_FACECENTERED_STAGGERED_ELEMENT_BOUNDARY_TYPES_HH

#include <vector>

namespace Dumux {

/*!
 * \ingroup FaceCenteredStaggeredDiscretization
 * \brief This class stores an array of BoundaryTypes objects
 */
template<class BTypes>
class FaceCenteredStaggeredElementBoundaryTypes
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
        if (!fvGeometry.hasBoundaryScvf())
            return;

        bcTypes_.resize(fvGeometry.numScvf());

        hasDirichlet_ = false;
        hasNeumann_ = false;

        for (const auto& scvf : scvfs(fvGeometry))
        {
            if (scvf.boundary())
            {
                bcTypes_[scvf.localIndex()] = problem.boundaryTypes(element, scvf);
                hasDirichlet_ = hasDirichlet_ || bcTypes_[scvf.localIndex()].hasDirichlet();
                hasNeumann_ = hasNeumann_ || bcTypes_[scvf.localIndex()].hasNeumann();
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
     */
    const BoundaryTypes& operator[] (std::size_t i) const
    {
        assert(i < bcTypes_.size());
        return bcTypes_[i];
    }

protected:
    std::vector<BoundaryTypes> bcTypes_;
    bool hasDirichlet_ = false;
    bool hasNeumann_ = false;
};

} // end namespace Dumux

#endif
