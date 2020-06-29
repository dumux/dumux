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

#include <dune/grid/io/file/vtk/basicwriter.hh>
#include <dune/grid/io/file/vtk/function.hh>
#include <dune/grid/io/file/vtk/skeletonfunction.hh>

namespace Dumux {

  //! iterate over the GridViews boundary intersections
/**
 * This will visit all intersections for which boundary() is true and
 * neighbor() is false.
 */
template<typename GV>
class GlobalIntersectionIterator
  : public Dune::ForwardIteratorFacade
    < GlobalIntersectionIterator<GV>,
        const typename GV::Intersection,
        const typename GV::Intersection&,
        typename std::iterator_traits<typename GV::template Codim<0>::
            Iterator>::difference_type>
{
public:
  // reiterator the facades typedefs here
  typedef GlobalIntersectionIterator<GV> DerivedType;
  typedef const typename GV::Intersection Value;
  typedef Value& Reference;
  typedef typename GV::template Codim<0>::Iterator ElementIterator;
  typedef typename GV::IntersectionIterator IntersectionIterator;
  typedef typename std::iterator_traits<ElementIterator>::difference_type
  DifferenceType;

private:
  typedef Dune::ForwardIteratorFacade<DerivedType, Value, Reference, DifferenceType> Facade;

  const GV* gv;
  ElementIterator eit;
  std::shared_ptr<IntersectionIterator> iit;

  bool valid() const {
    return true;
    // TODO maybe visit intersection only once
    // // we're valid if we're passed-the-end
    // if(eit == gv->template end<0>()) return true;
    // // or if we're on a boundary
    // if((*iit)->boundary() && !(*iit)->neighbor()) return true;
    // // otherwise we're invalid
    // return false;
  }

  void basic_increment() {
    ++*iit;
    if(*iit == gv->iend(*eit)) {
      iit.reset();
      ++eit;
      if(eit != gv->template end<0>())
        iit.reset(new IntersectionIterator(gv->ibegin(*eit)));
    }
  }

public:
  Reference dereference() const {
    return **iit;
  }
  bool equals(const DerivedType& other) const {
    if(eit != other.eit) return false;

    // this is a bit tricky, since we may not compare iit if we are
    // passed-the-end
    bool mePassedTheEnd = eit == gv->template end<0>();
    bool otherPassedTheEnd = other.eit == other.gv->template end<0>();

    // both passed-the-end => consider them equal
    if(mePassedTheEnd && otherPassedTheEnd) return true;

    // one passed the end => not equal
    if(mePassedTheEnd || otherPassedTheEnd) return false;

    // none passed-the-end => do their iit iterators match?
    return *iit == *other.iit;
  }

  void increment() {
    basic_increment();
    while(!valid()) basic_increment();
  }

  //! construct a GlobalIntersectionIterator
  /**
   * If end == true, construct an end iterator for the given gridview.
   * Otherwise, construct a begin iterator.
   */
  GlobalIntersectionIterator(const GV& gv_, bool end = false)
    : gv(&gv_), eit(end ? gv->template end<0>() : gv->template begin<0>())
  {
    intersectionVisited_.resize(gv_.size(1), false);

    if(eit != gv->template end<0>())
      iit.reset(new IntersectionIterator(gv->ibegin(*eit)));

    while(!valid()) basic_increment();
  }

private:
    std::vector<bool> intersectionVisited_;
};

template<class GV>
class NonConformingIntersectionIteratorFactory
{
  const GV& gv;

public:
  static const unsigned dimCell = GV::dimension-1;

  typedef typename GV::Intersection Cell;
  typedef GlobalIntersectionIterator<GV> CellIterator;

  typedef Dune::VTK::Corner<Cell> Corner;
  typedef Dune::VTK::CornerIterator<CellIterator> CornerIterator;

  typedef Corner Point;
  typedef CornerIterator PointIterator;

  typedef Dune::VTK::NonConformingConnectivityWriter<Cell> ConnectivityWriter;
  typedef typename GV::CollectiveCommunication CollectiveCommunication;

  explicit NonConformingIntersectionIteratorFactory(const GV& gv_)
    : gv(gv_)
  { }

    CellIterator beginCells() const
    { return CellIterator(gv); }

    CellIterator endCells() const
    { return CellIterator(gv, true); }

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
    { return gv.comm(); }

};

template<class GridView>
class ConformingIntersectionWriter
: public NonConformingIntersectionIteratorFactory<GridView>
, public Dune::VTK::BasicWriter<NonConformingIntersectionIteratorFactory<GridView>>
{
    using Factory = NonConformingIntersectionIteratorFactory<GridView>;
    using Base = Dune::VTK::BasicWriter<Factory>;

public:
  ConformingIntersectionWriter(const GridView& gridView)
    : Factory(gridView), Base(static_cast<const Factory&>(*this)), gridView_(gridView)
  { }

    using Base::addCellData;

  //  /**
  //    * @brief Add a grid function (represented by container) that lives on the cells of
  //    * the grid to the visualization.
  //    *
  //    * The container has to have random access via operator[] (e.g. std::vector). The
  //    * value of the grid function for an arbitrary element
  //    * will be accessed by calling operator[] with the index (corresponding
  //    * to the index from the MGMC mapper on the grid view) of the element.
  //    * For vector valued data all components for an element are assumed to
  //    * be consecutive.
  //    *
  //    * @param v The container with the values of the grid function for each cell.
  //    * @param name A name to identify the grid function.
  //    * @param ncomps Number of components (default is 1).
  //    */
    // template<class Container>
    // void addCellData (const Container& v, const std::string &name, int ncomps = 1,
    //                   Dune::VTK::Precision prec = Dune::VTK::Precision::float32)
    // {
    //   typedef Dune::P0VTKFunction<GridView, Container> Function;
    //   for (int c=0; c<ncomps; ++c) {
    //     std::stringstream compName;
    //     compName << name;
    //     if (ncomps>1)
    //       compName << "[" << c << "]";
    //     Dune::VTKFunction* p = new Function(gridView_, v, compName.str(), ncomps, c, prec);
    //     addCellData(std::shared_ptr< const VTKFunction >(p));
    //   }
    // }


    template<class Func>
    void addCellData(const std::shared_ptr<Func>& p, const std::string& name) {
      std::cout << "adding celldata 1" << std::endl;
      addCellData(std::shared_ptr<typename Base::FunctionWriter>
                    (new Dune::VTK::SkeletonFunctionWriter<Func>(p, name)));
    }

  template<class Func>
  void addCellData(Func* p, const std::string& name) {
    std::cout << "adding celldata 2" << std::endl;
    addCellData(std::shared_ptr<Func>(p), name);
  }

  using Base::addPointData;

  template<class Func>
  void addPointData(const std::shared_ptr<Func>& p, const std::string& name) {
    addPointData(std::shared_ptr<typename Base::FunctionWriter>
                    (new Dune::VTK::SkeletonFunctionWriter<Func>(p, name)));
  }

  template<class Func>
  void addPointData(Func* p, const std::string& name) {
    addPointData(std::shared_ptr<Func>(p), name);
  }

private:
    const GridView gridView_;

};

} // namespace Dumux

#endif // DUMUX_IO_VTK_INTERSECTIONWRITER_HH
