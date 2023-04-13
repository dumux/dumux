// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup FaceCenteredStaggeredDiscretization
 * \copydoc Dumux::FaceCenteredStaggeredLocalIntersectionIndexMapper
 */
#ifndef DUMUX_DISCRETIZATION_LOCAL_INTERSECTION_INDEX_MAPPER_HH
#define DUMUX_DISCRETIZATION_LOCAL_INTERSECTION_INDEX_MAPPER_HH

#include <array>
#include <numeric>
#include <dune/common/float_cmp.hh>
#include <dumux/common/indextraits.hh>
#include <dumux/common/math.hh>
#include <dumux/common/parameters.hh>
#include <dumux/discretization/facecentered/staggered/normalaxis.hh>
#include <dumux/discretization/facecentered/staggered/consistentlyorientedgrid.hh>

namespace Dumux {

namespace Detail {

template<class GridView, bool consistentlyOrientedGrid>
class FaceCenteredStaggeredLocalIntersectionIndexMapper;

/*!
 * \ingroup FaceCenteredStaggeredDiscretization
 * \brief Provides a mapping of local intersection indices (indexInInside)
 *        such that the local indices always follow the order of a reference element,
 *        regardless of how the element is oriented.
 */
template<class GridView>
class FaceCenteredStaggeredLocalIntersectionIndexMapper<GridView, false>
{
    using SmallLocalIndexType = typename IndexTraits<GridView>::SmallLocalIndex;
    using Element = typename GridView::template Codim<0>::Entity;
    static constexpr auto numElementFaces = GridView::Grid::dimension * 2;
public:
    void update(const GridView& gv, const Element& element)
    {
        static const bool makeConsistentlyOriented = getParam<bool>("Grid.MakeConsistentlyOriented", true);
        if (!makeConsistentlyOriented)
        {
            std::iota(realToRefMap_.begin(), realToRefMap_.end(), 0);
            refToRealMap_ = realToRefMap_;
            return;
        }

        for (const auto& is : intersections(gv, element))
        {
            const auto& otherOuterNormal = is.centerUnitOuterNormal();
            const auto idx = normalAxis(otherOuterNormal);
            const int positveOrientation = !std::signbit(otherOuterNormal[idx]);
            const auto refIdx = idx * 2 + positveOrientation;
            const auto realIdx = is.indexInInside();
            realToRefMap_[realIdx] = refIdx;
            refToRealMap_[refIdx] = realIdx;
        }
    }

    //! Return the intersection's actual local indexInElement given a local reference index.
    SmallLocalIndexType realToRefIdx(const SmallLocalIndexType localIsIdx) const
    { return realToRefMap_[localIsIdx]; }

    //! Return the intersection's local reference indexInElement given an actual local index.
    SmallLocalIndexType refToRealIdx(const SmallLocalIndexType localIsIdx) const
    { return refToRealMap_[localIsIdx]; }

private:
    std::array<SmallLocalIndexType, numElementFaces> realToRefMap_ = {};
    std::array<SmallLocalIndexType, numElementFaces> refToRealMap_ = {};
};

/*!
 * \ingroup FaceCenteredStaggeredDiscretization
 * \brief Provides a mapping of local intersection indices (indexInInside)
 *        such that the local indices always follow the order of a reference element,
 *        regardless of how the element in oriented.
 * \note This specialization is used for grids not supporting rotated elements.
 *       No mapping needs to be done here.
 */
template<class GridView>
class FaceCenteredStaggeredLocalIntersectionIndexMapper<GridView, true>
{
    using SmallLocalIndexType = typename IndexTraits<GridView>::SmallLocalIndex;
    using Element = typename GridView::template Codim<0>::Entity;
public:

    //! Update the map for getting the corresponding local face indices in another element.
    void update(const GridView&, const Element&) {}

    //! Return the intersection's actual local indexInElement given a local reference index.
    SmallLocalIndexType realToRefIdx(const SmallLocalIndexType localIsIdx) const
    { return localIsIdx; }

    //! Return the intersection's local reference indexInElement given an actual local index.
    SmallLocalIndexType refToRealIdx(const SmallLocalIndexType localIsIdx) const
    { return localIsIdx; }
};

} // end namespace Detail

/*!
 * \ingroup FaceCenteredStaggeredDiscretization
 * \brief Provides a mapping of local intersection indices (indexInInside)
 *        such that the local indices always follow the order of a reference element,
 *        regardless of how the element is oriented. Does not do anything for grids
 *        not supporting rotated elements (such as Dune::YaspGrid).
 */
template<class GridView>
using FaceCenteredStaggeredLocalIntersectionIndexMapper =
      Detail::FaceCenteredStaggeredLocalIntersectionIndexMapper<GridView, ConsistentlyOrientedGrid<typename GridView::Grid>{}>;

} // end namespace Dumux

#endif
