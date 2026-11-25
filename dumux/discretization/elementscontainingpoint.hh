// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup Discretization
 * \brief The elements and the sub-control volumes of a grid discretization that contain a point
 */
#ifndef DUMUX_DISCRETIZATION_ELEMENTS_CONTAINING_POINT_HH
#define DUMUX_DISCRETIZATION_ELEMENTS_CONTAINING_POINT_HH

#include <cstddef>
#include <type_traits>
#include <vector>

#include <dune/common/reservedvector.hh>
#include <dune/grid/common/capabilities.hh>

#include <dumux/discretization/localview.hh>
#include <dumux/discretization/method.hh>
#include <dumux/geometry/geometrytraits.hh>
#include <dumux/geometry/intersectingentities.hh>
#include <dumux/geometry/intersectspointgeometry.hh>

namespace Dumux {

namespace Detail {

//! Whether the grid discretization locates points in the sub-control volumes of its elements
template<class GridDiscretization>
inline constexpr bool locatesPointsInScvs = DiscretizationMethods::isCVFE<typename GridDiscretization::DiscretizationMethod>
                                            && requires (const typename GridDiscretization::LocalView& localView) { scvs(localView); };

} // end namespace Detail

/*!
 * \ingroup Discretization
 * \brief The local indices of the sub-control volumes of the bound element that contain a point
 * \tparam maxNumScvs the largest number of sub-control volumes of an element
 */
template<std::size_t maxNumScvs, class ElementDiscretization, class GlobalPosition>
Dune::ReservedVector<std::size_t, maxNumScvs> scvsContainingPoint(const ElementDiscretization& elemDisc,
                                                                  const GlobalPosition& globalPos)
{
    Dune::ReservedVector<std::size_t, maxNumScvs> scvIndices;
    for (const auto& scv : scvs(elemDisc))
        if (intersectsPointGeometry(globalPos, elemDisc.geometry(scv)))
            scvIndices.push_back(scv.indexInElement());
    return scvIndices;
}

/*!
 * \ingroup Discretization
 * \brief The indices of the elements that contain a point
 *
 * For discretizations with sub-control volumes, an element contains a point if one of its sub-control volumes
 * does; sub-control volumes that overlap others are not needed for that. The search of the grid alone also
 * returns elements that the point lies just outside of, within its tolerance, which can be larger than the one
 * of the sub-control volumes, so a point source of such an element would have no sub-control volume.
 *
 * \return the element indices, see the function `element` of the grid discretization
 */
template<class GridDiscretization, class GlobalPosition>
std::vector<std::size_t> elementsContainingPoint(const GridDiscretization& gridDiscretization,
                                                 const GlobalPosition& globalPos)
{
    auto elementIndices = intersectingEntities(globalPos, gridDiscretization.boundingBoxTree());
    if constexpr (Detail::locatesPointsInScvs<GridDiscretization>)
    {
        auto elemDisc = localView(gridDiscretization);
        std::erase_if(elementIndices, [&](const auto eIdx)
        {
            const auto element = gridDiscretization.element(eIdx);

            // the sub-control volumes tile the element, so one of them contains a point inside of the element
            using Grid = typename GridDiscretization::GridView::Grid;
            using Geometry = typename std::decay_t<decltype(element)>::Geometry;
            if constexpr (Dune::Capabilities::isCartesian<Grid>::v || GeometryTraits<Geometry>::axisAligned)
            {
                const auto geometry = element.geometry();
                const auto& lower = geometry.corner(0);
                const auto& upper = geometry.corner(geometry.corners() - 1);
                bool inside = true;
                for (int i = 0; i < GlobalPosition::dimension; ++i)
                    inside = inside && lower[i] <= globalPos[i] && globalPos[i] <= upper[i];
                if (inside)
                    return false;
            }

            elemDisc.bindElement(element);
            for (const auto& scv : scvs(elemDisc))
            {
                if constexpr (requires { scv.isOverlapping(); })
                    if (scv.isOverlapping())
                        continue;

                if (intersectsPointGeometry(globalPos, elemDisc.geometry(scv)))
                    return false;
            }
            return true;
        });
    }
    return elementIndices;
}

} // end namespace Dumux

#endif
