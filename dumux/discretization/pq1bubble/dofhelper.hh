// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup PQ1BubbleDiscretization
 * \brief Dof helper for the PQ1Bubble method, providing the number of dofs and their positions.
 */
#ifndef DUMUX_DISCRETIZATION_PQ1BUBBLE_DOF_HELPER_HH
#define DUMUX_DISCRETIZATION_PQ1BUBBLE_DOF_HELPER_HH

#include <cstddef>

#include <dune/common/rangeutilities.hh>
#include <dune/geometry/type.hh>
#include <dune/geometry/referenceelements.hh>

#include <dumux/common/indextraits.hh>
#include <dumux/discretization/fem/fedofhelper.hh>
#include <dumux/discretization/cvfe/localdof.hh>

namespace Dumux {

/*!
 * \ingroup PQ1BubbleDiscretization
 * \brief Dof helper for the PQ1Bubble method, providing the number of dofs and their positions.
 * \tparam GridView the grid view type
 * \tparam numCubeBubbleDofs number of bubble dofs for cube elements
 */
template <class GridView, std::size_t numCubeBubbleDofs>
class PQ1BubbleDofHelper : public Dumux::Experimental::FEDofHelper<GridView>
{
    using Scalar = typename GridView::ctype;
    using LocalIndexType = typename IndexTraits<GridView>::LocalIndex;
    using GridIndexType = typename IndexTraits<GridView>::GridIndex;

    using Element = typename GridView::template Codim<0>::Entity;

    static constexpr auto dim = GridView::dimension;

public:

    //! number of element dofs
    static std::size_t numElementDofs(Dune::GeometryType type)
    {
        const auto numVertexDofs = Dune::referenceElement<Scalar, dim>(type).size(dim);
        return numVertexDofs + (type.isCube() ? numCubeBubbleDofs : 1);
    }

    //! Number of local dofs related to an intersection with index iIdx
    static auto numLocalDofsIntersection(Dune::GeometryType type, unsigned int iIdx)
    {
        return Dune::referenceElement<Scalar, dim>(type).size(iIdx, 1, dim);
    }

    //! Local dof index related to a localDof, with index ilocalDofIdx, on an intersection with index iIdx
    static auto localDofIndexIntersection(Dune::GeometryType type, unsigned int iIdx, unsigned int ilocalDofIdx)
    {
        return Dune::referenceElement<Scalar, dim>(type).subEntity(iIdx, 1, ilocalDofIdx, dim);
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

#endif // DUMUX_DISCRETIZATION_PQ1BUBBLE_DOF_HELPER_HH
