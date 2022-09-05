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
#ifndef DUMUX_DISCRETIZATION_CCMPFA_O_TRANSMISSIBILITY_SOLVER_HH
#define DUMUX_DISCRETIZATION_CCMPFA_O_TRANSMISSIBILITY_SOLVER_HH

#include <ranges>
#include <cstdint>
#include <vector>
#include <array>
#include <concepts>
#include <algorithm>
#include <utility>
#include <optional>
#include <cmath>

#include <dune/common/fvector.hh>
#include <dune/common/dynmatrix.hh>

#include <dumux/common/math.hh>
#include <dumux/experimental/new_assembly/dumux/common/concepts.hh>
#include <dumux/experimental/new_assembly/dumux/discretization/ccmpfa-o/interactionregions.hh>
#include <dumux/experimental/new_assembly/dumux/discretization/ccmpfa-o/subcellshelper.hh>

namespace Dumux::CCMpfaO {

// TODO: These should go somewhere generic
// TODO: "contains" will be in std::ranges in cpp23
template<std::ranges::range Range, typename Value>
bool contains(const Range& range, const Value& value)
{
    return std::ranges::any_of(range, [&] (const auto& rangeValue) {
        return value == rangeValue;
    });
}

template<typename IndexType = std::size_t,
         std::ranges::range Range,
         std::invocable<std::ranges::range_reference_t<Range>> Predicate>
std::optional<IndexType> findIndex(const Range& range, const Predicate& pred)
{
    auto it = std::ranges::find_if(range, pred);
    if (it != std::ranges::end(range))
        return std::ranges::distance(std::ranges::begin(range), it);
    return {};
}
// TODO


template<typename GridGeometry>
class TransmissibilitySolver
{
    static constexpr int dim = GridGeometry::GridView::dimension;
    static constexpr int dimWorld = GridGeometry::GridView::dimensionworld;

    using ctype = typename GridGeometry::GridView::ctype;
    using Vector = Dune::FieldVector<ctype, dimWorld>;
    using LocalBasis = std::array<Vector, dim>;
    using LocalIndex = std::uint_least8_t;
    using SubCellFacetIndices = std::array<LocalIndex, dim>;

    struct Facet
    {
        LocalIndex localCellIndex;
        LocalIndex cellFacetIndex;

        friend bool operator==(const Facet& a, const Facet& b)
        {
            return a.localCellIndex == b.localCellIndex
                   && a.cellFacetIndex == b.cellFacetIndex;
        }
    };

    struct SubFace
    {
        Facet insideFacet;
        std::vector<Facet> outsideFacets;

        Vector normal;
        ctype area;
    };

    struct SubCell
    {
        const SubCellFacetIndices& facetIndices;
        std::vector<LocalIndex> localFaceIndices;
        LocalBasis gradientBasis;
    };

public:
    template<typename Scalar>
    using TransmissibilityMatrix = Dune::DynamicMatrix<Scalar>;

    explicit TransmissibilitySolver(const GridGeometry& gridGeometry,
                                    const InteractionRegion& region)
    : gridGeometry_(gridGeometry)
    , region_(region)
    { setup_(); }

    auto cellIndices() const
    {
        return region_ | std::views::transform([] (const CCMpfaO::SubCell& subCell) {
            return subCell.cellIndex;
        });
    }

    template<std::floating_point Scalar = double,
             std::invocable<std::size_t> TensorAccess,
             std::invocable<std::size_t, LocalIndex> IsDirichletFace> requires(
        std::convertible_to<bool, std::invoke_result_t<IsDirichletFace, std::size_t, LocalIndex>>)
    TransmissibilityMatrix<Scalar> solve(const TensorAccess& tensorAccess,
                                         const IsDirichletFace& isDirichletFace) const
    {
        TransmissibilityMatrix<Scalar> result(subFaces_.size(), subCells_.size());
        TransmissibilityMatrix<Scalar> A(subFaces_.size(), subFaces_.size());
        TransmissibilityMatrix<Scalar> B(subFaces_.size(), subCells_.size());
        TransmissibilityMatrix<Scalar> C(subFaces_.size(), subFaces_.size());
        TransmissibilityMatrix<Scalar> D(subFaces_.size(), subCells_.size());

        LocalIndex subFaceIndex = 0;
        for (const auto& subFace : subFaces_)
        {
            const auto insideSubCellLocalIndex = subFace.insideFacet.localCellIndex;
            const auto& insideSubCell = subCells_[insideSubCellLocalIndex];
            const auto& insideTensor = tensorAccess(region_[insideSubCellLocalIndex].cellIndex);

            for (int dir = 0; dir < dim; dir++)
            {
                const auto directionLocalFaceIndex = insideSubCell.localFaceIndices[dir];
                const auto tijk = computeTransmissibility_<Scalar>(insideSubCell, subFace, insideTensor, dir);
                A[subFaceIndex][directionLocalFaceIndex] -= tijk;
                B[subFaceIndex][insideSubCellLocalIndex] -= tijk;
                C[subFaceIndex][directionLocalFaceIndex] -= tijk;
                D[subFaceIndex][insideSubCellLocalIndex] += tijk;
            }

            for (const auto [outsideSubCellLocalIndex, cellFacetIdx] : subFace.outsideFacets)
            {
                const auto& outsideSubCell = subCells_[outsideSubCellLocalIndex];
                const auto& outsideTensor = tensorAccess(region_[outsideSubCellLocalIndex].cellIndex);

                for (int dir = 0; dir < dim; dir++)
                {
                    const auto directionLocalFaceIndex = outsideSubCell.localFaceIndices[dir];
                    const auto tijk = computeTransmissibility_<Scalar>(outsideSubCell, subFace, outsideTensor, dir);
                    A[subFaceIndex][directionLocalFaceIndex] += tijk;
                    B[subFaceIndex][outsideSubCellLocalIndex] += tijk;
                }
            }

            subFaceIndex++;
        }

        A.invert();
        C.rightmultiply(A);
        B.leftmultiply(C);
        D += B;

        std::cout << "T = \n" << D << std::endl;

        return result;
    }

private:
    template<std::integral Index>
    LocalIndex getLocalCellIndex_(const Index& cellIndex) const
    {
        const auto indices = cellIndices();
        const auto it = std::ranges::find(indices, cellIndex);
        assert(it != std::ranges::end(indices));
        return std::ranges::distance(std::ranges::begin(indices), it);
    }

    void setup_()
    {
        subFaces_.reserve(region_.size()*dim);
        subCells_.reserve(region_.size());

        LocalIndex regionSubCellIndex = 0;
        for (const auto& [cellIdx, subCellIdx] : region_)
        {
            const auto& element = gridGeometry_.element(cellIdx);
            const auto& elemGeom = element.geometry();
            const auto& helper = SubCellsHelper{elemGeom};

            makeSubCell_(
                helper,
                subCellIdx,
                makeSubCellFaces_(element, helper, subCellIdx, regionSubCellIndex)
            );
            regionSubCellIndex++;
        }

        subFaces_.shrink_to_fit();
    }

    template<typename Geometry>
    void makeSubCell_(const SubCellsHelper<Geometry>& helper,
                      unsigned int localSubCellIndex,
                      std::vector<LocalIndex>&& localFaceIndices)
    {
        auto localBasis = helper.makeBasis(localSubCellIndex);
        makeGradientBasis_(localBasis);
        subCells_.emplace_back(SubCell{
            helper.facetIndices(localSubCellIndex),
            std::move(localFaceIndices),
            std::move(localBasis)
        });
    }

    void makeGradientBasis_(LocalBasis& localBasis) const
    {
        static_assert(dim == 2);
        if constexpr (dim == 2)
        {
            rotate90Degrees_(localBasis[0]);
            rotate90Degrees_(localBasis[1], true);
            std::swap(localBasis[0], localBasis[1]);

            using std::abs;
            const auto determinant = crossProduct(localBasis[0], localBasis[1]);
            std::ranges::for_each(localBasis, [&] (auto& v) { v /= determinant; });
        }
    }

    template<std::enable_if_t<dim == 2, bool> = true>
    void rotate90Degrees_(Vector& v, bool clockwise = false) const
    {
        std::swap(v[0], v[1]);
        v[(clockwise ? 1 : 0)] *= -1.0;
    }

    template<typename Element, typename Geometry>
    std::vector<LocalIndex> makeSubCellFaces_(const Element& element,
                                              const SubCellsHelper<Geometry>& helper,
                                              unsigned int elementSubCellIndex,
                                              const LocalIndex regionSubCellIndex)
    {
        const auto& facetsToInsert = helper.facetIndices(elementSubCellIndex);
        std::vector<LocalIndex> subCellSubFaceIndices(facetsToInsert.size());

        for (const auto& is : intersections(gridGeometry_.gridView(), element))
        {
            const auto subFaceIndex = findIndex(facetsToInsert, [&] (const auto& idx) {
                return idx == is.indexInInside();
            });

            if (subFaceIndex)
            {
                const auto [inserted, localFaceIdx] = insertFace_(
                    Facet{regionSubCellIndex, makeLocalIndex_(is.indexInInside())}
                );

                if (inserted)
                {
                    subFaces_.back().normal = is.centerUnitOuterNormal();
                    subFaces_.back().area = helper.subFaceArea(elementSubCellIndex, *subFaceIndex);
                    if (is.neighbor())
                        addToCurrentNeighbors_(Facet{
                            getLocalCellIndex_(gridGeometry_.elementMapper().index(is.outside())),
                            makeLocalIndex_(is.indexInOutside())
                        });
                }

                subCellSubFaceIndices[*subFaceIndex] = localFaceIdx;
            }
        }

        return subCellSubFaceIndices;
    }

    std::pair<bool, LocalIndex> insertFace_(Facet&& facet)
    {
        const auto faceIdx = findIndex(subFaces_, [&] (const SubFace& subFace) {
            return contains(subFace.outsideFacets, facet);
        });

        if (!faceIdx)
        {
            subFaces_.emplace_back(SubFace{std::move(facet), {}, {}, {}});
            return {true, makeLocalIndex_(subFaces_.size() - 1)};
        }
        return {false, *faceIdx};
    }

    void addToCurrentNeighbors_(Facet&& facet)
    { subFaces_.back().outsideFacets.emplace_back(std::move(facet)); }

    template<std::integral I>
    LocalIndex makeLocalIndex_(const I& i) const
    { return static_cast<LocalIndex>(i); }

    template<typename Scalar, typename Tensor>
    Scalar computeTransmissibility_(const SubCell& subCell,
                                    const SubFace& subFace,
                                    const Tensor& tensor,
                                    int direction) const
    {
        return subFace.area*vtmv(
            subFace.normal, tensor, subCell.gradientBasis[direction]
        );
    }

    const GridGeometry& gridGeometry_;
    const InteractionRegion& region_;

    std::vector<SubFace> subFaces_;
    std::vector<SubCell> subCells_;
};

} // end namespace Dumux::CCMpfaO

#endif
