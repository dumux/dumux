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
 * \brief Boundary types gathered on an element in the
 *        context of finite element schemes.
 */
#ifndef DUMUX_FE_ELEMENT_BOUNDARY_TYPES_HH
#define DUMUX_FE_ELEMENT_BOUNDARY_TYPES_HH

#include <cassert>
#include <vector>

#include <dune/common/exceptions.hh>

namespace Dumux {

/*!
 * \ingroup FEMDiscretization
 * \brief This class stores an array of BoundaryTypes objects, where each
 *        entry holds the boundary types defined for a degree of freedom.
 */
template<class BTypes>
class FEElementBoundaryTypes
{
public:
    using BoundaryTypes = BTypes;

    /*!
     * \brief Update the boundary types for all boundary degrees of freedom of an element.
     * \param problem The problem to be solved
     * \param element The element for which the boundary types are to be collected.
     * \param feGeometry A local view on the finite element grid geometry,
     *                   bound to the provided element.
     */
    template<class Problem, class Element, class FEElementGeometry>
    void update(const Problem& problem,
                const Element& element,
                const FEElementGeometry& feGeometry)
    {
        const auto& gridGeometry = feGeometry.gridGeometry();
        const auto& feBasisLocalView = feGeometry.feBasisLocalView();
        const auto numLocalDofs = feBasisLocalView.size();
        const auto& fe = feBasisLocalView.tree().finiteElement();

        bcTypes_.resize(numLocalDofs);
        hasDirichlet_ = false;
        hasNeumann_ = false;
        hasOutflow_ = false;

        for (unsigned int localDofIdx = 0; localDofIdx < numLocalDofs; localDofIdx++)
        {
            // skip all dofs not living on the boundary
            if (!gridGeometry.dofOnBoundary(feBasisLocalView.index(localDofIdx)))
                continue;

            const auto& localKey = fe.localCoefficients().localKey(localDofIdx);
            const auto subEntity = localKey.subEntity();
            const auto codim = localKey.codim();

            // Obtain user-defined boundary conditions
            bcTypes_[localDofIdx] = getBoundaryTypes_(problem, element, subEntity, codim);

            hasDirichlet_ = hasDirichlet_ || bcTypes_[localDofIdx].hasDirichlet();
            hasNeumann_ = hasNeumann_ || bcTypes_[localDofIdx].hasNeumann();
            hasOutflow_ = hasOutflow_ || bcTypes_[localDofIdx].hasOutflow();
        }

        // Outflow BCs are not supported
        if (hasOutflow_)
            DUNE_THROW(Dune::NotImplemented, "Outflow BCs for finite element scheme");
    }

    /*!
     * \brief Returns true if the element contains a dof with Dirichlet BCs assigned to it.
     */
    bool hasDirichlet() const
    { return hasDirichlet_; }

    /*!
     * \brief Returns true if the element contains a dof with Neumann BCs assigned to it.
     */
    bool hasNeumann() const
    { return hasNeumann_; }

    /*!
     * \brief Returns true if the element contains a dof with outflow BCs assigned to it.
     */
    bool hasOutflow() const
    { return hasOutflow_; }

    /*
     * \brief Access operator
     * \param i The element-local dof index
     * \return the boundary types for the i-th element-local dof
     */
    const BoundaryTypes& operator[] (std::size_t i) const
    {
        assert(i < bcTypes_.size());
        return bcTypes_[i];
    }

protected:
    template<class Problem, class Element>
    BoundaryTypes getBoundaryTypes_(const Problem& problem,
                                    const Element& element,
                                    unsigned int subEntity,
                                    unsigned int codim)
    {
        static constexpr int dim = Element::Geometry::mydimension;

        if (codim == 0)
            return getBoundaryTypes_<0, dim>(problem, element, subEntity);
        if (codim == 1)
            return getBoundaryTypes_<1, dim>(problem, element, subEntity);
        if constexpr (dim > 1)
            if (codim == 2)
                return getBoundaryTypes_<2, dim>(problem, element, subEntity);
        if constexpr (dim > 2)
        {
            assert(codim == 3);
            return getBoundaryTypes_<3, dim>(problem, element, subEntity);
        }

        DUNE_THROW(Dune::InvalidStateException, "Invalid codimension");
    }

    template<int codim, int dim, class Problem, class Element>
    BoundaryTypes getBoundaryTypes_(const Problem& problem,
                                    const Element& element,
                                    unsigned int subEntity)
    {
        static_assert(codim <= dim, "Invalid codimension");
        return problem.boundaryTypes(element, element.template subEntity<codim>(subEntity));
    }

    std::vector<BoundaryTypes> bcTypes_;
    bool hasDirichlet_ = false;
    bool hasNeumann_ = false;
    bool hasOutflow_ = false;
};

} // namespace Dumux

#endif
