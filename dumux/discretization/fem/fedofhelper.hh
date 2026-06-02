// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup FEDiscretization
 * \brief Default Dof helper for finite-element discretizations providing dof-related utility functions.
 *        Assuming that all dofs map directly to grid entities (vertices, edges, faces, elements).
 */
#ifndef DUMUX_DISCRETIZATION_FE_DOF_HELPER_HH
#define DUMUX_DISCRETIZATION_FE_DOF_HELPER_HH

#include <ranges>

#include <dune/geometry/type.hh>
#include <dune/geometry/referenceelements.hh>

#include <dumux/common/indextraits.hh>
#include <dumux/discretization/cvfe/localdof.hh>

namespace Dumux {

/*!
 * \ingroup FEDiscretization
 * \brief Default Dof helper for finite-element discretizations providing dof-related utility functions.
 *        Assuming that all dofs map directly to grid entities (vertices, edges, faces, elements).
 * \tparam GridView the grid view type
 */
template <class GridView>
class FEDofHelper
{
    using Scalar = typename GridView::ctype;
    using GlobalPosition = typename Dune::FieldVector<Scalar, GridView::dimensionworld>;
    using LocalIndexType = typename IndexTraits<GridView>::LocalIndex;
    using GridIndexType = typename IndexTraits<GridView>::GridIndex;

    using Element = typename GridView::template Codim<0>::Entity;

    static constexpr auto dim = GridView::dimension;

public:

    //! Returns true if the local dof with a given local key is on the intersection with index iIdx
    template<class LocalKey>
    static auto localDofOnIntersection(Dune::GeometryType type, unsigned int iIdx, const LocalKey& localKey)
    {
        const auto& refElement = Dune::referenceElement<Scalar, dim>(type);
        const auto numEntitiesIntersection = refElement.size(iIdx, 1, localKey.codim());
        for (std::size_t idx = 0; idx < numEntitiesIntersection; idx++)
            if (localKey.subEntity() == refElement.subEntity(iIdx, 1, idx, localKey.codim()))
                return true;
        return false;
    }

    /*!
     * \brief Iterator range over all local dofs on a given boundary face.
     *         Uses a filter over all local dofs via localDofOnIntersection.
     * \param elemDisc the element discretization (must be bound)
     * \param boundaryFace the boundary face
     */
    template<class ElemDisc, class BoundaryFace>
    static auto localDofsOnBoundaryFace(const ElemDisc& elemDisc, const BoundaryFace& boundaryFace)
    {
        const auto& gridDisc = [&]() -> const auto& {
            if constexpr (requires { elemDisc.gridDiscretization(); })
                return elemDisc.gridDiscretization();
            else
                return elemDisc.gridGeometry();
        }();

        // The following implementation is not the most efficient one, but it is general and does not require any assumptions on the ordering of the dofs.
        // Discretization-specific implementations can override this with a more efficient implementation,
        // e.g. by directly only iterating over vertices of the boundaryFace and the related intersection.
        return std::views::iota(std::size_t(0), elemDisc.numLocalDofs())
            | std::views::filter([&](std::size_t i) {
                return localDofOnIntersection(
                    elemDisc.element().type(),
                    boundaryFace.intersectionIndex(),
                    elemDisc.feLocalCoefficients().localKey(i));
            })
            | std::views::transform([&](std::size_t i) {
                return CVFE::LocalDof(
                    static_cast<LocalIndexType>(i),
                    static_cast<GridIndexType>(dofIndex(
                        gridDisc.dofMapper(),
                        elemDisc.element(),
                        elemDisc.feLocalCoefficients().localKey(i))),
                    static_cast<GridIndexType>(elemDisc.elementIndex())
                );
            });
    }

    //! global index of dof
    template<class DofMapper, class LocalKey>
    static auto dofIndex(const DofMapper& dofMapper, const Element& element, const LocalKey& localKey)
    {
        return dofMapper.subIndex(element, localKey.subEntity(), localKey.codim()) + localKey.index();
    }

    //! global dof position
    template<class Geometry, class LocalKey>
    static GlobalPosition dofPosition(const Geometry& geo, const LocalKey& localKey)
    {
        if (localKey.codim() == dim)
            return geo.corner(localKey.subEntity());
        else if (localKey.codim() == 0)
            return geo.center();
        else
            return geo.global(localDofPosition(geo.type(), localKey));
    }

    //! local dof position
    template<class LocalKey>
    static typename Element::Geometry::LocalCoordinate localDofPosition(Dune::GeometryType type, const LocalKey& localKey)
    {
        return Dune::referenceElement<Scalar, dim>(type).position(localKey.subEntity(), localKey.codim());
    }
};

} // end namespace Dumux

#endif
