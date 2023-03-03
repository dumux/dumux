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
 * \ingroup Common
 * \copybrief Dumux::BSplineFunction
 */
#ifndef DUMUX_COMMON_BSPLINE_FUNCTION_HH
#define DUMUX_COMMON_BSPLINE_FUNCTION_HH
#if HAVE_DUNE_FUNCTIONS

#include <array>
#include <utility>
#include <numeric>
#include <algorithm>
#include <type_traits>

#include <dune/common/fvector.hh>
#include <dune/common/iteratorfacades.hh>
#include <dune/grid/common/mcmgmapper.hh>
#include <dune/grid/yaspgrid.hh>
#include <dune/istl/bvector.hh>

#include <dune/functions/functionspacebases/bsplinebasis.hh>
#include <dune/functions/functionspacebases/lagrangebasis.hh>

#include <dune/functions/gridfunctions/discreteglobalbasisfunction.hh>
#include <dune/functions/gridfunctions/gridviewfunction.hh>

#include <dumux/common/entitymap.hh>
#include <dumux/geometry/geometricentityset.hh>
#include <dumux/geometry/intersectingentities.hh>
#include <dumux/discretization/projection/projector.hh>

namespace Dumux {
namespace BSplineFunctionDetail {

// Exposes the intersections interface a gridview with itself
template<class GridView>
class IdentityIntersectionEntitySet
{
    using ctype = typename GridView::ctype;
    using Entity = typename GridView::template Codim<0>::Entity;

    class IntersectionEntity
    {
    public:
        IntersectionEntity(const Entity& entity) : entity_(entity) {}
        decltype(auto) geometry() const { return entity_.geometry(); }
        constexpr std::size_t numDomainNeighbors() const { return 1; }
        constexpr std::size_t numTargetNeighbors() const { return 1; }
        const Entity& domainEntity(unsigned int n = 0) const { return entity_; }
        const Entity& targetEntity(unsigned int n = 0) const { return entity_; }

    private:
        Entity entity_;
    };

    template<typename ElementIterator>
    class IntersectionEntityIterator
    : public Dune::ForwardIteratorFacade<
        IntersectionEntityIterator<ElementIterator>,
        IntersectionEntity,
        IntersectionEntity
    >
    {
    public:
        IntersectionEntityIterator(ElementIterator it) : it_(it) {}
        bool equals(const IntersectionEntityIterator& other) const { return other.it_ == it_; }
        IntersectionEntity dereference() const { return {*it_}; }
        void increment() { ++it_; }

    private:
        ElementIterator it_;
    };

public:
    IdentityIntersectionEntitySet(const GridView& gv) : gv_(gv) {}
    auto ibegin() const { return IntersectionEntityIterator{gv_.template begin<0>()}; }
    auto iend() const { return IntersectionEntityIterator{gv_.template end<0>()}; }
    std::size_t size() const { return gv_.size(0); }
    friend auto intersections(const IdentityIntersectionEntitySet& set)
    { return Dune::IteratorRange{set.ibegin(), set.iend()}; }

private:
    const GridView gv_;
};

} // namespace BSplineFunctionDetail

/*!
 * \ingroup Common
 * \brief TODO: Doc me
 */
template<typename Range, typename ctype, int dim>
class BSplineFunction
{
    template<typename T> struct IsFieldVector : public std::false_type {};
    template<typename T, int s> struct IsFieldVector<Dune::FieldVector<T, s>> : public std::true_type {};
    static_assert(IsFieldVector<Range>::value || Dune::IsNumber<Range>::value);

    template<typename T> struct RangeDimension : public std::integral_constant<int, 0> {};
    template<typename T, int s> struct RangeDimension<Dune::FieldVector<T, s>> : public std::integral_constant<int, dim> {};
    static constexpr int rangeDimension = RangeDimension<Range>::value;

    using Grid = Dune::YaspGrid<dim, Dune::EquidistantOffsetCoordinates<ctype, dim>>;
    using GridView = typename Grid::LeafGridView;
    using Domain = typename GridView::template Codim<0>::Entity::Geometry::GlobalCoordinate;
    using ElementMapper = Dune::MultipleCodimMultipleGeomTypeMapper<GridView>;
    using ElementMap = EntityMap<GridView, 0>;

    using Basis = Dune::Functions::BSplineBasis<GridView>;
    using PreBasis = typename Basis::PreBasis;

    using Coefficients = Dune::BlockVector<Range>;
    using Function = Dune::Functions::DiscreteGlobalBasisFunction<Basis, Coefficients>;
    using ProjectionRange = std::conditional_t<rangeDimension == 0, Dune::FieldVector<Range, 1>, Range>;


public:
    /*!
     * \brief TODO: Doc me
     */
    BSplineFunction(Dune::FieldVector<ctype, dim> minDomain,
                    Dune::FieldVector<ctype, dim> maxDomain,
                    std::array<std::size_t, dim> numSamples,
                    const Dune::BlockVector<Range>& samples,
                    unsigned int order = 2,
                    bool makeOpen = true)
    : min_{std::move(minDomain)}
    , max_{std::move(maxDomain)}
    , cells_{makeCellsArray_(numSamples)}
    , grid_{min_, max_, cells_}
    , mapper_{grid_.leafGridView(), Dune::mcmgElementLayout()}
    , elementMap_{grid_, mapper_}
    , basis_{
        PreBasis{
            gridView_(),
            min_,
            max_,
            convertArray_<unsigned int>(cells_),
            order,
            makeOpen
        }
    }
    , coefficients_{projectSamples_(samples)}
    , function_{Dune::Functions::makeDiscreteGlobalBasisFunction<Range>(basis_, coefficients_)}
    {}

    /*!
     * \brief TODO: Doc me
     */
    Range operator()(const Domain& x) const
    {
        const auto elementIndex = intersectingEntityCartesianGrid(x, min_, max_, cells_);
        const auto element = elementMap_[10];
        auto localF = localFunction(function_);
        localF.bind(element);
        return localF(element.geometry().local(x));
    }

private:
    GridView gridView_() const
    { return grid_.leafGridView(); }

    template<typename T, typename U>
    std::array<T, dim> convertArray_(const std::array<U, dim>& in) const
    {
        std::array<T, dim> result;
        std::transform(in.begin(), in.end(), result.begin(), [&] (const U u) { return static_cast<T>(u); });
        return result;
    }

    std::array<int, dim> makeCellsArray_(const std::array<std::size_t, dim>& numSamples) const
    {
        if (std::any_of(numSamples.begin(), numSamples.end(), [] (const auto& n) { return n < 2; }))
            DUNE_THROW(Dune::InvalidStateException, "At least 2 samples per direction required");
        auto cells = convertArray_<int>(numSamples);
        std::for_each(cells.begin(), cells.end(), [&] (int& u) { u -= 1; });
        return cells;
    }

    auto projectSamples_(const Dune::BlockVector<Range>& samples)
    {
        if (samples.size() != gridView_().size(dim))
            DUNE_THROW(
                Dune::InvalidStateException,
                "Length of the samples vector does not match the given domain dimensions"
            );

        const auto sampleBasis = Dune::Functions::BasisFactory::makeBasis(
            gridView_(),
            Dune::Functions::BasisFactory::lagrange<1, Range>()
        );
        const auto intersections = BSplineFunctionDetail::IdentityIntersectionEntitySet{gridView_()};
        const auto projector = Dumux::makeProjector(sampleBasis, basis_, intersections);
        return castBlockVector_<Range>(
            projector.project(castBlockVector_<ProjectionRange>(samples))
        );
    }

    template<typename To, typename From, std::enable_if_t<!std::is_same_v<From, To>, bool> = true>
    auto castBlockVector_(const Dune::BlockVector<From>& v) const
    {
        Dune::BlockVector<To> result(v.size());
        std::copy(v.begin(), v.end(), result.begin());
        return result;
    }

    template<typename To, typename From, std::enable_if_t<std::is_same_v<From, To>, bool> = true>
    const auto& castBlockVector_(const Dune::BlockVector<From>& v) const
    { return v; }

    Dune::FieldVector<ctype, dim> min_;
    Dune::FieldVector<ctype, dim> max_;
    std::array<int, dim> cells_;

    Grid grid_;
    ElementMapper mapper_;
    ElementMap elementMap_;

    Basis basis_;
    Coefficients coefficients_;
    Function function_;
};

} // end namespace Dumux

#endif // HAVE_DUNE_FUNCTIONS
#endif
