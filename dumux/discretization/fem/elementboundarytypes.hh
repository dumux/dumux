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
 * \ingroup FEMDiscretization
 * \brief Boundary types gathered on an element
 */
#ifndef DUMUX_FE_ELEMENT_BOUNDARY_TYPES_HH
#define DUMUX_FE_ELEMENT_BOUNDARY_TYPES_HH

#include <algorithm>
#include <cassert>
#include <vector>

#include <dune/common/exceptions.hh>

namespace Dumux {

/*!
 * \ingroup FEMDiscretization
 * \brief This class stores an array of BoundaryTypes objects
 */
template<class BTypes>
class FEElementBoundaryTypes
{
public:
    using BoundaryTypes = BTypes;

    /*!
     * \brief Update the boundary types for all intersections of an element.
     *
     * \param problem The problem object which needs to be simulated
     * \param element The DUNE Codim<0> entity for which the boundary
     *                types should be collected
     * \param feGeometry The finite element geometry
     */
    template<class Problem, class Element, class FEElementGeometry>
    void update(const Problem& problem,
                const Element& element,
                const FEElementGeometry& feGeometry)
    {
        const auto& gridGeometry = feGeometry.gridGeometry();
        const auto& feBasisLocalView = feGeometry.feBasisLocalView();
        const auto numLocalDofs = feBasisLocalView.size();

        localBCTypes_.resize( numLocalDofs );
        hasDirichlet_ = false;
        hasNeumann_ = false;
        hasOutflow_ = false;

        static constexpr int dim = Element::Geometry::mydimension;
        const auto& fe = feBasisLocalView.tree().finiteElement();

        for (unsigned int localDofIdx = 0; localDofIdx < numLocalDofs; localDofIdx++)
        {
            // skip all dofs not living on the boundary
            if (!gridGeometry.dofOnBoundary(feBasisLocalView.index(localDofIdx)))
                continue;

            const auto& localKey = fe.localCoefficients().localKey(localDofIdx);
            const auto subEntity = localKey.subEntity();
            const auto codim = localKey.codim();

            // Obtain user-defined boundary conditions
            if (codim == 1) localBCTypes_[localDofIdx] = getBoundaryTypes_<1, dim>(problem, element, subEntity);
            else if (codim == 2) localBCTypes_[localDofIdx] = getBoundaryTypes_<2, dim>(problem, element, subEntity);
            else if (codim == 3) localBCTypes_[localDofIdx] = getBoundaryTypes_<3, dim>(problem, element, subEntity);
        }

        const auto beginIt = localBCTypes_.begin();
        const auto endIt = localBCTypes_.end();
        hasDirichlet_ = std::any_of(beginIt, endIt, [] (const auto& bt) { return bt.hasDirichlet(); });
        hasNeumann_ = std::any_of(beginIt, endIt, [] (const auto& bt) { return bt.hasNeumann(); });
        hasOutflow_ = std::any_of(beginIt, endIt, [] (const auto& bt) { return bt.hasOutflow(); });

        // Outflow BCs are not supported
        if (hasOutflow_)
            DUNE_THROW(Dune::NotImplemented, "Outflow BCs for finite element scheme");
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
        assert(i < localBCTypes_.size());
        return localBCTypes_[i];
    }

protected:
    template<int codim, int dim, class Problem, class Element,
             std::enable_if_t<(codim <= dim), int> = 0>
    BoundaryTypes getBoundaryTypes_(const Problem& problem,
                                    const Element& element,
                                    unsigned int subEntity)
    { return problem.boundaryTypes(element.template subEntity<codim>(subEntity)); }

    template<int codim, int dim, class Problem, class Element,
             std::enable_if_t<(codim > dim), int> = 0>
    BoundaryTypes getBoundaryTypes_(const Problem& problem,
                                    const Element& element,
                                    unsigned int subEntity)
    { DUNE_THROW(Dune::InvalidStateException, "Codimension higher than grid dimension"); }

    std::vector< BoundaryTypes > localBCTypes_;
    bool hasDirichlet_ = false;
    bool hasNeumann_ = false;
    bool hasOutflow_ = false;
};

} // namespace Dumux

#endif
