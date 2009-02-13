// $Id$

#ifndef DUNE_ONEDINNDGRID_INDEXSETS_HH
#define DUNE_ONEDINNDGRID_INDEXSETS_HH

/** \file
    \brief The index and id sets for the OneDInNDGrid class
*/

#include <vector>

namespace Dune {

/** \todo Please doc me! */

template <class GridImp>
struct OneDInNDGridLevelIndexSetTypes
{
  //! The types
  template<int cd>
  struct Codim
  {
    template<PartitionIteratorType pitype>
    struct Partition
    {
      typedef typename GridImp::Traits::template Codim<cd>::template Partition<pitype>::LevelIterator Iterator;
    };
  };
};

/** \todo Please doc me! */

template<class GridImp>
class OneDInNDGridLevelIndexSet : public IndexSetDefaultImplementation<GridImp,OneDInNDGridLevelIndexSet<GridImp>,OneDInNDGridLevelIndexSetTypes<GridImp> >
{
  typedef IndexSetDefaultImplementation<GridImp,OneDInNDGridLevelIndexSet<GridImp>,OneDInNDGridLevelIndexSetTypes<GridImp> > Base;
  enum { dimworld = GridImp::dimensionworld };
public:

    /** \brief Constructor for a given level of a given grid
    */
    OneDInNDGridLevelIndexSet (const GridImp& grid, int level)
        : grid_(&grid), level_(level)
    {}

  //! get index of an entity
  template<int cd>
  int index (const typename GridImp::Traits::template Codim<cd>::Entity& e) const
  {
      return grid_->getRealImplementation(e).levelIndex();
  }

  //! get index of subentity of a codim 0 entity
  template<int cc>
  int subIndex (const typename GridImp::Traits::template Codim<0>::Entity& e, int i) const
  {
      return grid_->getRealImplementation(e).template subLevelIndex<cc>(i);
  }

  //! get number of entities of given type and on this level
  int size (GeometryType type) const
  {
      return grid_->size(level_,type);
  }

  //! get number of entities of given codim, type and on this level
  int size (int codim) const
  {
      return grid_->size(level_,codim);
  }

    /** \brief Deliver all geometry types used in this grid */
  const std::vector<GeometryType>& geomTypes (int codim) const
  {
    return myTypes_[codim];
  }

  //! one past the end on this level
  template<int cd, PartitionIteratorType pitype>
  typename Base::template Codim<cd>::template Partition<pitype>::Iterator begin () const
  {
    return grid_->template lbegin<cd,pitype>(level_);
  }

  //! Iterator to one past the last entity of given codim on level for partition type
  template<int cd, PartitionIteratorType pitype>
  typename Base::template Codim<cd>::template Partition<pitype>::Iterator end () const
  {
    return grid_->template lend<cd,pitype>(level_);
  }

    /** \todo Should be private */
    void update() {

        // ///////////////////////////////
        //   Init the element indices
        // ///////////////////////////////
        numElements_ = 0;
        OneDInNDEntityImp<1,dimworld>* eIt;
        for (eIt = grid_->elements[level_].begin; eIt!=NULL; eIt = eIt->succ_)
            eIt->levelIndex_ = numElements_++;

        // //////////////////////////////
        //   Init the vertex indices
        // //////////////////////////////

        numVertices_ = 0;
        OneDInNDEntityImp<0,dimworld>* vIt;
        for (vIt = grid_->vertices[level_].begin; vIt!=NULL; vIt = vIt->succ_)
            vIt->levelIndex_ = numVertices_++;

        // ///////////////////////////////////////////////
        //   Update the list of geometry types present
        // ///////////////////////////////////////////////
        if (numElements_>0) {
            myTypes_[0].resize(1);
            myTypes_[0][0] = GeometryType(1);
        } else
            myTypes_[0].resize(0);

        if (numVertices_>0) {
            myTypes_[1].resize(1);
            myTypes_[1][0] = GeometryType(0);
        } else
            myTypes_[1].resize(0);
    }

private:
  const GridImp* grid_;
  int level_;

    int numElements_;
    int numVertices_;

    /** \brief The GeometryTypes present for each codim */
  std::vector<GeometryType> myTypes_[2];
};

/** \todo Please doc me! */

template <class GridImp>
struct OneDInNDGridLeafIndexSetTypes
{
  //! The types
  template<int cd>
  struct Codim
  {
    template<PartitionIteratorType pitype>
    struct Partition
    {
      typedef typename GridImp::Traits::template Codim<cd>::template Partition<pitype>::LeafIterator Iterator;
    };
  };
};

/** \todo Please doc me! */

template<class GridImp>
class OneDInNDGridLeafIndexSet :
  public IndexSetDefaultImplementation<GridImp,OneDInNDGridLeafIndexSet<GridImp>,OneDInNDGridLeafIndexSetTypes<GridImp> >
{
  typedef IndexSetDefaultImplementation<GridImp,OneDInNDGridLeafIndexSet<GridImp>,OneDInNDGridLeafIndexSetTypes<GridImp> > Base;
public:
  //! constructor stores reference to a grid and level
  OneDInNDGridLeafIndexSet (const GridImp& g) : grid_(g)
  {
  }

  //! get index of an entity
  /*
    We use the RemoveConst to extract the Type from the mutable class,
    because the const class is not instatiated yet.
  */
  template<int cd>
  int index (const typename remove_const<GridImp>::type::Traits::template Codim<cd>::Entity& e) const
  {
      return grid_.getRealImplementation(e).leafIndex();
  }
//  template<int cd>
//  int index (const typename RemoveConst<GridImp>::Type::Traits::template Codim<cd>::Entity& e) const
//  {
//      return grid_.getRealImplementation(e).leafIndex();
//  }

  //! get index of subentity of a codim 0 entity
  /*
    We use the RemoveConst to extract the Type from the mutable class,
    because the const class is not instatiated yet.
  */
  template<int cc>
  int subIndex (const typename remove_const<GridImp>::type::Traits::template Codim<0>::Entity& e, int i) const
  {
      return grid_.getRealImplementation(e).template subLeafIndex<cc>(i);
  }
//  template<int cc>
//  int subIndex (const typename RemoveConst<GridImp>::Type::Traits::template Codim<0>::Entity& e, int i) const
//  {
//      return grid_.getRealImplementation(e).template subLeafIndex<cc>(i);
//  }

  //! get number of entities of given codim, type on the leaf level
  int size (GeometryType type) const
  {
      if (type.isVertex()) {

          return numVertices_;

      } else if (type.isLine()) {

          return numElements_;

      }

      return 0;
  }

  //! get number of entities of given codim, type on the leaf level
  int size(int codim) const
  {
    return Base::size(codim);
  }

  /** deliver all geometry types used in this grid */
  const std::vector<GeometryType>& geomTypes (int codim) const
  {
    return myTypes_[codim];
  }

  //! one past the end on this level
  template<int cd, PartitionIteratorType pitype>
  typename Base::template Codim<cd>::template Partition<pitype>::Iterator begin () const
  {
    return grid_.template leafbegin<cd,pitype>();
  }

  //! Iterator to one past the last entity of given codim on level for partition type
  template<int cd, PartitionIteratorType pitype>
  typename Base::template Codim<cd>::template Partition<pitype>::Iterator end () const
  {
    return grid_.template leafend<cd,pitype>();
  }

    /** \todo Should be private */
    void update() {

        // ///////////////////////////////
        //   Init the element indices
        // ///////////////////////////////
        numElements_ = 0;
        typename GridImp::Traits::template Codim<0>::LeafIterator eIt    = grid_.template leafbegin<0>();
        typename GridImp::Traits::template Codim<0>::LeafIterator eEndIt = grid_.template leafend<0>();

        for (; eIt!=eEndIt; ++eIt) {

            grid_.getRealImplementation(*eIt).target_->leafIndex_ = numElements_++;

        }

        // //////////////////////////////
        //   Init the vertex indices
        // //////////////////////////////

        numVertices_ = 0;

        for (int i=grid_.maxLevel(); i>=0; i--) {

            OneDInNDEntityImp<0, GridImp::dimensionworld>* vIt;
            for (vIt = grid_.vertices[i].begin; vIt!=NULL; vIt = vIt->succ_) {

                if (vIt->isLeaf())
                    vIt->leafIndex_ = numVertices_++;
                else
                    vIt->leafIndex_ = vIt->son_->leafIndex_;

            }

        }

        // ///////////////////////////////////////////////
        //   Update the list of geometry types present
        // ///////////////////////////////////////////////
        if (numElements_>0) {
            myTypes_[0].resize(1);
            myTypes_[0][0] = GeometryType(1);
        } else
            myTypes_[0].resize(0);

        if (numVertices_>0) {
            myTypes_[1].resize(1);
            myTypes_[1][0] = GeometryType(0);
        } else
            myTypes_[1].resize(0);

    }

private:

    const GridImp& grid_;

    int numElements_;
    int numVertices_;

    /** \brief The GeometryTypes present for each codim */
    std::vector<GeometryType> myTypes_[2];
};

/** \todo Please doc me! */

template<class GridImp>
class OneDInNDGridIdSet : public IdSet<GridImp,OneDInNDGridIdSet<GridImp>,unsigned int>
{
public:
  //! define the type used for persistent indices
  typedef unsigned int GlobalIdType;
  typedef unsigned int LocalIdType;

  //! constructor stores reference to a grid
  OneDInNDGridIdSet (const GridImp& g) : grid_(g) {}

  //! get id of an entity
  /*
    We use the RemoveConst to extract the Type from the mutable class,
    because the const class is not instatiated yet.
  */
  template<int cd>
  GlobalIdType id (const typename remove_const<GridImp>::type::Traits::template Codim<cd>::Entity& e) const
  {
      return grid_.getRealImplementation(e).globalId();
  }
//  template<int cd>
//  GlobalIdType id (const typename RemoveConst<GridImp>::Type::Traits::template Codim<cd>::Entity& e) const
//  {
//      return grid_.getRealImplementation(e).globalId();
//  }

  //! get id of subentity
  /*
    We use the RemoveConst to extract the Type from the mutable class,
    because the const class is not instatiated yet.
  */
  template<int cd>
  GlobalIdType subId (const typename remove_const<GridImp>::type::Traits::template Codim<0>::Entity& e, int i) const
  {
      return grid_.template getRealImplementation(e).template subId<cd>(i);
  }
//  template<int cd>
//  GlobalIdType subId (const typename RemoveConst<GridImp>::Type::Traits::template Codim<0>::Entity& e, int i) const
//  {
//      return grid_.template getRealImplementation(e).template subId<cd>(i);
//  }

    /** \todo Should be private */
    void update() {}

private:

  const GridImp& grid_;
};

}  // namespace Dune


#endif
