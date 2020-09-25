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

#ifndef DUMUX_IO_VTK_INTERSECTIONWRITER_HH
#define DUMUX_IO_VTK_INTERSECTIONWRITER_HH

#include <memory>
#include <string>

#include <dune/grid/io/file/vtk.hh>
#include <dune/grid/io/file/vtk/basicwriter.hh>
#include <dune/grid/io/file/vtk/function.hh>
#include <dune/grid/io/file/vtk/skeletonfunction.hh>
#include <dumux/common/typetraits/typetraits.hh>

namespace Dumux {

/*!
 * \brief Iterate over the GridViews boundary intersections
 * This will visit all intersections for which boundary() is true and
 * neighbor() is false.
 */
template<typename GV>
class GlobalIntersectionIterator
: public Dune::ForwardIteratorFacade<GlobalIntersectionIterator<GV>,
                                     const typename GV::Intersection,
                                     const typename GV::Intersection&,
                                     typename std::iterator_traits<typename GV::template Codim<0>::Iterator>::difference_type>
{
public:
    using DerivedType = GlobalIntersectionIterator<GV>;
    using Intersection = typename GV::Intersection;
    using Value = const Intersection;
    using ElementIterator = typename GV::template Codim<0>::Iterator;
    using IntersectionIterator = typename GV::IntersectionIterator;
    using DifferenceType = typename std::iterator_traits<ElementIterator>::difference_type;

  /*!
   * \brief Construct a GlobalIntersectionIterator
   * If end == true, construct an end iterator for the given gridview.
   * Otherwise, construct a begin iterator.
   */
    GlobalIntersectionIterator(const GV& gv, bool end = false)
    : gridView_(gv), eIt_(end ? gridView_.template end<0>() : gridView_.template begin<0>())
    {
      if (eIt_ != gridView_.template end<0>())
          iIt_.reset(new IntersectionIterator(gridView_.ibegin(*eIt_)));
    }

    const Intersection& dereference() const
    {
        if constexpr (std::is_lvalue_reference_v<decltype(*std::declval<IntersectionIterator>())>)
            return **iIt_;
        else
        {
            // Some grids only return temporary intersection objects.
            // If this is the case, story a copy of the intersection here so that
            // its address can be safely taken later on.
            intersection_ = **iIt_;
            return intersection_;
        }
    }

    bool equals(const DerivedType& other) const
    {
        if (eIt_ != other.eIt_)
            return false;

        // this is a bit tricky, since we may not compare iIt_ if we are
        // passed-the-end
        bool mePassedTheEnd = eIt_ == gridView_.template end<0>();
        bool otherPassedTheEnd = other.eIt_ == other.gridView_.template end<0>();

        // both passed-the-end => consider them equal
        if(mePassedTheEnd && otherPassedTheEnd)
            return true;

        // one passed the end => not equal
        if(mePassedTheEnd || otherPassedTheEnd)
            return false;

        // none passed-the-end => do their iIt_ iterators match?
        return *iIt_ == *other.iIt_;
    }

    void increment()
    {
        ++*iIt_;
        if (*iIt_ == gridView_.iend(*eIt_))
        {
            iIt_.reset();
            ++eIt_;
            if (eIt_ != gridView_.template end<0>())
              iIt_.reset(new IntersectionIterator(gridView_.ibegin(*eIt_)));
        }
    }

private:
    const GV gridView_;
    ElementIterator eIt_;
    std::shared_ptr<IntersectionIterator> iIt_;
    mutable Intersection intersection_;
};

template<class GridView>
class NonConformingIntersectionIteratorFactory
{
public:
    static constexpr auto dimCell = GridView::dimension-1;
    using Cell = typename GridView::Intersection;
    using CellIterator = GlobalIntersectionIterator<GridView>;
    using Corner = Dune::VTK::Corner<Cell>;
    using CornerIterator = Dune::VTK::CornerIterator<CellIterator>;
    using Point = Corner;
    using PointIterator = CornerIterator;
    using ConnectivityWriter = Dune::VTK::NonConformingConnectivityWriter<Cell>;
    using CollectiveCommunication = typename GridView::CollectiveCommunication;

    explicit NonConformingIntersectionIteratorFactory(const GridView& gv)
    : gridView_(gv) {}

    CellIterator beginCells() const
    { return CellIterator(gridView_); }

    CellIterator endCells() const
    { return CellIterator(gridView_, true); }

    CornerIterator beginCorners() const
    { return CornerIterator(beginCells(), endCells()); }

    CornerIterator endCorners() const
    { return CornerIterator(endCells()); }

    PointIterator beginPoints() const
    { return beginCorners(); }

    PointIterator endPoints() const
    { return endCorners(); }

    ConnectivityWriter makeConnectivity() const
    { return ConnectivityWriter(); }

    const CollectiveCommunication& comm() const
    { return gridView_.comm(); }

private:
    const GridView gridView_;
};


template<class GridView, class Mapper, class F>
class SkeletonFunction
{
    using Intersection = typename GridView::Intersection;
public:
    using Traits = Dune::VTK::SkeletonFunctionTraits<GridView, typename GridView::ctype>;

    SkeletonFunction(const GridView& gv, const Mapper& mapper, const F& field)
    : gv_(gv)
    , mapper_(mapper)
    , field_(field)
    , components_(1)
    {
        if constexpr (std::is_invocable_v<F, Intersection, int>)
        {
            for (const auto& element : elements(gv))
            {
                for (const auto& is : intersections(gv, element))
                {
                    if constexpr (IsIndexable<std::decay_t<decltype(field(std::declval<Intersection>(), 0))>>{})
                    {
                        if constexpr (IsIndexable<std::decay_t<decltype(field(std::declval<Intersection>(), 0)[0])>>{})
                            DUNE_THROW(Dune::InvalidStateException, "Invalid field type");
                        else
                            components_ = field(is, mapper_(is, gv_)).size();
                    }

                    return;
                }
            }
        }
        else if constexpr (IsIndexable<std::decay_t<decltype(field[0])>>{})
        {
            assert(field.size() == gv.size(1));
            if constexpr (IsIndexable<std::decay_t<decltype(field[0][0])>>{})
              DUNE_THROW(Dune::InvalidStateException, "Invalid field type");
            else
              components_ = field[0].size();
        }
    }

    //! return number of components
    unsigned dimRange() const { return components_; }

    void evaluate(const typename Traits::Cell& intersection,
                  const typename Traits::Domain& xl,
                  typename Traits::Range& result) const
    {
        assert(intersection.conforming());
        result.resize(components_);
        const auto idx = mapper_(intersection, gv_);

        auto accessEntry = [&](auto i)
        {
            if constexpr (std::is_invocable_v<F, Intersection, int>)
            {
                if constexpr (IsIndexable<std::decay_t<decltype(field_(intersection, idx))>>{})
                    return field_(intersection, idx)[i];
                else
                    return field_(intersection, idx);
            }
            else
            {
                if constexpr (IsIndexable<std::decay_t<decltype(std::declval<F>()[0])>>{})
                    return field_[idx][i];
                else
                    return field_[idx];
            }
        };

        for (int i = 0; i < components_; ++i)
            result[i] = accessEntry(i);
    }

private:
    const GridView gv_;
    const Mapper& mapper_;

    // If F is callable we store a copy of it, otherwise we assume that F is a container
    // stored somewhere else, thus we only keep a reference here
    std::conditional_t<std::is_invocable_v<F,Intersection, int>, const F, const F&> field_;

    std::size_t components_;
};

template<class GridView>
class ConformingIntersectionWriter
: public NonConformingIntersectionIteratorFactory<GridView>
, public Dune::VTK::BasicWriter<NonConformingIntersectionIteratorFactory<GridView>>
{
    using Factory = NonConformingIntersectionIteratorFactory<GridView>;
    using Base = Dune::VTK::BasicWriter<Factory>;

public:
    ConformingIntersectionWriter(const GridView& gridView, const std::string& paramGroup = "")
    : Factory(gridView), Base(static_cast<const Factory&>(*this)), gridView_(gridView)
    {
        static bool addProcessRank = getParamFromGroup<bool>(paramGroup, "Vtk.AddProcessRank");
        if (addProcessRank)
        {
            auto getRank = [rank = gridView_.comm().rank()](const auto& is, const auto idx)
            { return rank; };

            auto mapper = getStandardMapper();
            auto fun = std::make_shared<SkeletonFunction<GridView, decltype(mapper), decltype(getRank)>>(gridView_, mapper, getRank);
            addCellData(fun, "processRank");
        }
    }

    using Base::addCellData;

    static auto getStandardMapper()
    {
        return [](const auto& is, const GridView& gridView){ return gridView.indexSet().subIndex(is.inside(), is.indexInInside(), 1); };
    }

    template<class F, class Mapper = decltype(getStandardMapper())>
    auto makeSkeletonFunction(const F& f, const Mapper& mapper = getStandardMapper()) const
    {
        return std::make_shared<SkeletonFunction<GridView, decltype(mapper), decltype(f)>>(gridView_, mapper, f);
    }

    template<class Func>
    void addCellData(const std::shared_ptr<Func>& p, const std::string& name)
    {
      addCellData(std::shared_ptr<typename Base::FunctionWriter>
                  (new Dune::VTK::SkeletonFunctionWriter<Func>(p, name)));
    }

    template<class Func>
    void addCellData(Func* p, const std::string& name)
    {
        addCellData(std::shared_ptr<Func>(p), name);
    }

    template<class F>
    void addField(const F& field, const std::string& name)
    {
        auto mapper = getStandardMapper();
        addCellData(makeSkeletonFunction(field, mapper), name);
    }

    using Base::addPointData;

    template<class Func>
    void addPointData(const std::shared_ptr<Func>& p, const std::string& name)
    {
        addPointData(std::shared_ptr<typename Base::FunctionWriter>
                     (new Dune::VTK::SkeletonFunctionWriter<Func>(p, name)));
    }

    template<class Func>
    void addPointData(Func* p, const std::string& name)
    {
        addPointData(std::shared_ptr<Func>(p), name);
    }

private:
    const GridView gridView_;
};

} // namespace Dumux

#endif // DUMUX_IO_VTK_INTERSECTIONWRITER_HH
