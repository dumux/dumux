// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup PQ1Discretization
 * \brief Helper class providing degree of freedom information for the PQ1 Lagrange basis.
 */
#ifndef DUMUX_DISCRETIZATION_PQ1_LAGRANGE_DOF_HELPER_HH
#define DUMUX_DISCRETIZATION_PQ1_LAGRANGE_DOF_HELPER_HH

#include <dune/common/rangeutilities.hh>
#include <dune/geometry/type.hh>
#include <dune/geometry/referenceelements.hh>

#include <dumux/common/indextraits.hh>
#include <dumux/discretization/fem/fedofhelper.hh>
#include <dumux/discretization/cvfe/localdof.hh>

namespace Dumux {

/*!
 * \ingroup PQ1Discretization
 * \brief Helper class providing degree of freedom information for the PQ1 Lagrange basis.
 *        The pq1 basis uses vertex-only dofs.
 */
template<class GridView>
class PQ1LagrangeDofHelper : public Dumux::Experimental::FEDofHelper<GridView>
{
    using Scalar = typename GridView::ctype;
    using LocalIndexType = typename IndexTraits<GridView>::LocalIndex;
    using GridIndexType = typename IndexTraits<GridView>::GridIndex;

    static constexpr auto dim = GridView::dimension;
public:
    //! Number of local dofs related to an intersection with index iIdx (vertex dofs only)
    static auto numLocalDofsIntersection(Dune::GeometryType type, unsigned int iIdx)
    {
        return Dune::referenceElement<Scalar, dim>(type).size(iIdx, 1, dim);
    }

    //! Local dof index of the localDofIdx-th dof on intersection with index iIdx
    static auto localDofIndexIntersection(Dune::GeometryType type, unsigned int iIdx, unsigned int localDofIdx)
    {
        return Dune::referenceElement<Scalar, dim>(type).subEntity(iIdx, 1, localDofIdx, dim);
    }

    /*!
     * \brief Iterator range over all local dofs on a given boundary face.
     *        Overrides the less efficient default implementation.
     */
    template<class ElemDisc, class BoundaryFace>
    static auto localDofsOnBoundaryFace(const ElemDisc& elemDisc, const BoundaryFace& boundaryFace)
    {
        return Dune::transformedRangeView(
            Dune::range(numLocalDofsIntersection(elemDisc.element().type(), boundaryFace.intersectionIndex())),
            [&](const auto i) {
                auto localDofIdx = localDofIndexIntersection(elemDisc.element().type(), boundaryFace.intersectionIndex(), i);
                return CVFE::LocalDof(
                    static_cast<LocalIndexType>(localDofIdx),
                    static_cast<GridIndexType>(Dumux::Experimental::FEDofHelper<GridView>::dofIndex(
                        elemDisc.gridGeometry().dofMapper(),
                        elemDisc.element(),
                        elemDisc.feLocalCoefficients().localKey(localDofIdx))),
                    static_cast<GridIndexType>(elemDisc.elementIndex())
                );
            }
        );
    }
};

} // namespace Dumux

#endif // DUMUX_DISCRETIZATION_PQ1_LAGRANGE_DOF_HELPER_HH
