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
        const auto& feBasisLocalView = feGeometry.feBasisLocalView();
        const auto numLocalDofs = feBasisLocalView.size();

        localBCTypes_.resize( numLocalDofs );
        hasDirichlet_ = false;
        hasNeumann_ = false;
        hasOutflow_ = false;

        using ctype = typename Element::Geometry::ctype;
        using RefElements = typename Dune::ReferenceElements<ctype, Element::Geometry::mydimension>;

        const auto& eg = element.geometry();
        const auto refElement = RefElements::general(eg.type());
        const auto& fe = feBasisLocalView.tree().finiteElement();

        for (const auto& is : intersections(feGeometry.gridGeometry().gridView(), element))
        {
            if (is.boundary())
            {
                const auto bcTypes = problem.boundaryTypes(element, is);
                if (bcTypes.hasOutflow())
                    DUNE_THROW(Dune::NotImplemented, "Outflow BCs for FEM");

                // find local dofs lying on this intersection
                std::vector<unsigned int> localDofs; localDofs.reserve(numLocalDofs);
                for (unsigned int localDofIdx = 0; localDofIdx < numLocalDofs; localDofIdx++)
                {
                    const auto& localKey = fe.localCoefficients().localKey(localDofIdx);
                    const auto subEntity = localKey.subEntity();
                    const auto codim = localKey.codim();

                    // skip interior dofs
                    if (codim == 0)
                        continue;

                    bool found = false;
                    // try to find this local dof (on entity with known codim) on the current intersection
                    for (int j = 0; j < refElement.size(is.indexInInside(), 1, codim); j++)
                    {
                        // If j-th sub entity is the sub entity corresponding to local dof, continue and assign BC
                        if (subEntity == refElement.subEntity(is.indexInInside(), 1, j, codim))
                        { localDofs.push_back(localDofIdx); found = true; }

                        if (found) break;
                    }
                }

                // Write obtained boundary types into the ones for the local dofs
                for (const auto localDof : localDofs)
                    copyBTypes_(bcTypes, localBCTypes_[localDof]);
            }
        }

        const auto beginIt = localBCTypes_.begin();
        const auto endIt = localBCTypes_.end();
        hasDirichlet_ = std::any_of(beginIt, endIt, [] (const auto& bt) { return bt.hasDirichlet(); });
        hasNeumann_ = std::any_of(beginIt, endIt, [] (const auto& bt) { return bt.hasNeumann(); });
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
    /*
     * \brief Writes new boundary types into the ones stored
     *        for a local dof. We do not overwrite previously
     *        set Dirichlet BCS with Neumann BCS, which could
     *        occur for vertex dofs.
     */
    void copyBTypes_(const auto& curBCTypes, auto& localBCTypes)
    {
        for (unsigned int i = 0; i < BTypes::size(); ++i)
        {
            const bool curIsNeumann = curBCTypes.isNeumann(i);

            // do not overwrite Dirichlet with Neumann
            // if (curIsNeumann && !localBCTypes.isDirichlet(i))
            if (curIsNeumann)
                localBCTypes.setNeumann(i);

            // always set Dirichlet BCs
            else
                localBCTypes.setDirichlet(i);
        }
    }

    std::vector< BoundaryTypes > localBCTypes_;
    bool hasDirichlet_ = false;
    bool hasNeumann_ = false;
    bool hasOutflow_ = false;
};

} // namespace Dumux

#endif
