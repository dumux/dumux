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
 * \ingroup CCMpfaDiscretization
 * \todo Doc me!
 */
#ifndef DUMUX_DISCRETIZATION_CCMPFA_O_SUB_CELLS_HELPER_HH
#define DUMUX_DISCRETIZATION_CCMPFA_O_SUB_CELLS_HELPER_HH

#include <array>
#include <ranges>
#include <concepts>
#include <algorithm>
#include <utility>

#include <dune/common/exceptions.hh>
#include <dune/geometry/type.hh>
#include <dune/geometry/referenceelement.hh>

#include <dumux/geometry/volume.hh>
#include <dumux/experimental/new_assembly/dumux/discretization/ccmpfa-o/subcellscorners.hh>

namespace Dumux::CCMpfaO {

#ifndef DOXYGEN
namespace Detail {

template<typename Geometry>
auto getReferenceElement(const Geometry& geo)
{
    using Dune::referenceElement;
    return referenceElement(geo);
}

template<typename Geometry>
using ReferenceElementType = std::decay_t<decltype(getReferenceElement(std::declval<const Geometry&>()))>;

} // namespace Detail
#endif // DOXYGEN


template<typename Geometry,
         typename ReferenceElement = Detail::ReferenceElementType<Geometry>>
class SubCellsHelper
{
    static constexpr int dim = Geometry::mydimension;
    static_assert(dim == 2 || dim == 3);

public:
    static constexpr int numSubFacesPerSubCell = dim;
    static constexpr int numCornersPerSubFace = dim == 2 ? 2 : 4;

    using Coordinate = typename Geometry::GlobalCoordinate;
    using SubCellBasis = std::array<Coordinate, numSubFacesPerSubCell>;
    using SubFaceCorners = std::array<Coordinate, numCornersPerSubFace>;

    explicit SubCellsHelper(const Geometry& geo)
    : geo_(geo)
    , refElement_(Detail::getReferenceElement(geo))
    {}

    SubCellsHelper(const Geometry& geo, ReferenceElement&& refElement)
    : geo_(geo)
    , refElement_(std::move(refElement))
    {}

    constexpr int numSubCells() const
    {
        return dispatchLocalCornersHelper_([&] <typename Helper> (const Helper&) {
            return Helper::numSubCells;
        });
    }

    constexpr std::ranges::range decltype(auto) facetIndices(unsigned int subCellIndex) const
    {
        return dispatchLocalCornersHelper_([&] <typename Helper> (const Helper&) {
            return Helper::facetIndices(subCellIndex);
        });
    }

    constexpr SubCellBasis makeBasis(unsigned int subCellIndex) const
    {
        return dispatchLocalCornersHelper_([&] <typename Helper> (const Helper&) {
            auto result = makePositions_(1, Helper::facetIndices(subCellIndex));
            const auto& center = geo_.center();
            std::ranges::for_each(result, [&] (auto& p) { p -= center; });
            return result;
        });
    }

    constexpr std::ranges::range decltype(auto) subFaceCornerIndices(unsigned int subCellIndex,
                                                                     unsigned int subCellSubFaceIndex) const
    {
        assert(subCellSubFaceIndex < numSubFacesPerSubCell);
        return dispatchLocalCornersHelper_([&] <typename Helper> (const Helper&) {
            return Helper::subFaceCorners(subCellIndex)[subCellSubFaceIndex];
        });
    }

    constexpr SubFaceCorners subFaceCorners(unsigned int subCellIndex,
                                            unsigned int subCellSubFaceIndex) const
    {
        assert(subCellSubFaceIndex < numSubFacesPerSubCell);
        return makePositions_(subFaceCornerIndices(subCellIndex, subCellSubFaceIndex));
    }

    constexpr auto subFaceArea(unsigned int subCellIndex,
                               unsigned int subCellSubFaceIndex) const
    {
        static_assert(dim == 2);
        if constexpr (dim == 2)
            return convexPolytopeVolume<dim-1>(
                Dune::GeometryTypes::line,
                subFaceCorners(subCellIndex, subCellSubFaceIndex)
            );
    }

private:
    template<std::size_t N>
    std::array<Coordinate, N> makePositions_(const std::array<SubEntity, N>& indices) const
    { return makePositions_(indices, std::make_index_sequence<N>{}); }

    template<std::integral I1, std::integral I2, std::size_t N>
    std::array<Coordinate, N> makePositions_(const I1& codim,
                                             const std::array<I2, N>& indices) const
    { return makePositions_(codim, indices, std::make_index_sequence<N>{}); }

    template<std::size_t N, std::size_t... Indices>
    std::array<Coordinate, N> makePositions_(const std::array<SubEntity, N>& indices,
                                             const std::index_sequence<Indices...>&) const
    { return {makePosition_(indices[Indices].codim, indices[Indices].index)...}; }

    template<std::integral I1, std::integral I2, std::size_t N, std::size_t... Indices>
    std::array<Coordinate, N> makePositions_(const I1& codim,
                                             const std::array<I2, N>& indices,
                                             const std::index_sequence<Indices...>&) const
    { return {makePosition_(codim, indices[Indices])...}; }

    template<std::integral I1, std::integral I2>
    Coordinate makePosition_(const I1& codim, const I2& index) const
    { return geo_.global(refElement_.position(index, codim)); }

    template<typename HelperFunctionCall>
    constexpr decltype(auto) dispatchLocalCornersHelper_(const HelperFunctionCall& helperCall) const
    {
        static_assert(dim == 2);
        if constexpr (dim == 2)
        {
            constexpr auto quad = Dune::GeometryTypes::quadrilateral;
            if (geo_.type() == quad)
                return helperCall(SubCellsLocalCorners<quad>{});
        }

        DUNE_THROW(Dune::NotImplemented, "Unsupported dimension & geometry type");
    }

    const Geometry& geo_;
    ReferenceElement refElement_;
};

} // end namespace Dumux::CCMpfaO

#endif
