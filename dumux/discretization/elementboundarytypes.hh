// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup Discretization
 * \brief Boundary types gathered on an element
 */
#ifndef DUMUX_ELEMENT_BOUNDARY_TYPES_HH
#define DUMUX_ELEMENT_BOUNDARY_TYPES_HH

#include <cassert>
#include <vector>

#include <dune/geometry/referenceelements.hh>

namespace Dumux {

/*!
 * \ingroup Discretization
 * \brief This class stores an array of BoundaryTypes objects.
 *        This class is not dependent on the used discretization method
 *        It only requires the function fvGeometry.intersectionIndex(scvf)
 */
template<class BTypes>
class ElementIntersectionBoundaryTypes
{
public:
    using BoundaryTypes = BTypes;

    /*!
     * \brief Update the boundary types for all element intersections.
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
        static constexpr int dim = FVElementGeometry::GridGeometry::GridView::dimension;
        using Scalar = typename FVElementGeometry::GridGeometry::GlobalCoordinate::value_type;
        // Resize according to the number of faces (intersections)
        bcTypes_.resize(Dune::referenceElement<Scalar, dim>(element.geometry().type()).size(dim-1));

        hasDirichlet_ = false;
        hasNeumann_ = false;

        for (const auto& intersection : intersections(fvGeometry.gridGeometry().gridView(), element))
        {
            const auto localIdx = intersection.indexInInside();
            bcTypes_[localIdx].reset();

            if (intersection.boundary() && !intersection.neighbor())
            {
                bcTypes_[localIdx] = problem.boundaryTypes(fvGeometry, intersection);
                hasDirichlet_ = hasDirichlet_ || bcTypes_[localIdx].hasDirichlet();
                hasNeumann_ = hasNeumann_ || bcTypes_[localIdx].hasNeumann();
            }
        }
    }

    /*!
     * \brief Returns whether the element has an intersection
     *        with Dirichlet conditions
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
     * \note yields undefined behaviour if the scvf is not on the boundary
     */
    template<class FVElementGeometry>
    const BoundaryTypes& get(const FVElementGeometry& fvGeometry, const typename FVElementGeometry::SubControlVolumeFace& scvf) const
    {
        assert(scvf.boundary());
        const auto localIdx = fvGeometry.intersectionIndex(scvf);
        assert(localIdx < bcTypes_.size());
        return bcTypes_[localIdx];
    }

    /*
     * \brief Access operator
     * \return BoundaryTypes
     * \note yields undefined behaviour if the intersection is not on the boundary
     */
    template<class FVElementGeometry>
    const BoundaryTypes& get(const FVElementGeometry&, const typename FVElementGeometry::GridGeometry::GridView::Intersection& intersection) const
    {
        assert(intersection.boundary());
        const auto localIdx = intersection.indexInInside();
        assert(localIdx < bcTypes_.size());
        return bcTypes_[localIdx];
    }

private:
    std::vector<BoundaryTypes> bcTypes_;
    bool hasDirichlet_ = false;
    bool hasNeumann_ = false;
};


} // namespace Dumux

#endif
