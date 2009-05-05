#ifndef DUNE_ONEDINNDGRID_INDEXSETS_HH
#define DUNE_ONEDINND_GRID_INDEXSETS_HH

/** \file
    \brief The index and id sets for the OneDInNDGrid class
*/

#include <vector>

namespace Dune {

template<class GridImp>
class OneDInNDGridLevelIndexSet : public IndexSet<GridImp,OneDInNDGridLevelIndexSet<GridImp> >
{

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
  template<int codim>
  unsigned int subIndex (const typename GridImp::Traits::template Codim<0>::Entity& e, int i) const
  {
      return subIndex(e,i,codim);
  }

  //! get index of subentity of a codim 0 entity
  unsigned int subIndex (const typename GridImp::Traits::template Codim<0>::Entity& e, 
                         int i,
                         unsigned int codim) const
  {
      return grid_->getRealImplementation(e).subLevelIndex(i,codim);
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

    /** \brief Return true if e is contained in the index set.
        
    Warning: this implementation takes O(n) time!  It also assumes that e belongs
    to the correct grid.
    */
  template <class EntityType>
    bool contains (const EntityType& e) const
    {
        enum { cd = EntityType::codimension };
        typedef typename GridImp::template Codim<cd>::template Partition<All_Partition>::LevelIterator IteratorType; 
        IteratorType iend = grid_->template lend<cd,All_Partition>(level_);
        for (IteratorType it = grid_->template lbegin<cd,All_Partition>(level_);
             it != iend; ++it)
            {
                if (it->level() == e.level() && this->template index<cd>(*it) == this->template index<cd>(e)) 
                    return true;
            }
        return false; 
    }

    /** \todo Should be private */
    void update() {

        // ///////////////////////////////
        //   Init the element indices
        // ///////////////////////////////
        numElements_ = 0;
        typename OneDInNDGridList<OneDInNDEntityImp<1, dimworld> >::const_iterator eIt;
        for (eIt = grid_->elements[level_].begin(); eIt!=grid_->elements[level_].end(); eIt = eIt->succ_)
            /** \todo Remove this const cast */
	  const_cast<OneDInNDEntityImp<1, dimworld>*>(eIt)->levelIndex_ = numElements_++;

        // //////////////////////////////
        //   Init the vertex indices
        // //////////////////////////////

        numVertices_ = 0;
        typename OneDInNDGridList<OneDInNDEntityImp<0,dimworld> >::const_iterator vIt;
        for (vIt = grid_->vertices[level_].begin(); vIt!=grid_->vertices[level_].end(); vIt = vIt->succ_)
            /** \todo Remove this const cast */
	  const_cast<OneDInNDEntityImp<0, dimworld>*>(vIt)->levelIndex_ = numVertices_++;

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

template<class GridImp>
class OneDInNDGridLeafIndexSet : 
  public IndexSet<GridImp,OneDInNDGridLeafIndexSet<GridImp> >
{
  //enum { dimworld = GridImp::dimensionworld };
public:
  //! constructor stores reference to a grid and level
  OneDInNDGridLeafIndexSet (const GridImp& g) : grid_(g)
  {
  }

  //! get index of an entity
  /*
    We use the remove_const to extract the Type from the mutable class,
    because the const class is not instantiated yet.
  */
  template<int cd>
  int index (const typename remove_const<GridImp>::type::Traits::template Codim<cd>::Entity& e) const 
  {
      return grid_.getRealImplementation(e).leafIndex(); 
  }

  //! get index of subentity of a codim 0 entity
  /*
    We use the remove_const to extract the Type from the mutable class,
    because the const class is not instantiated yet.
  */
  template<int codim>
  int subIndex (const typename remove_const<GridImp>::type::Traits::template Codim<0>::Entity& e, int i) const
  {
      return subIndex(e,i,codim);
  }

  //! get index of subentity of a codim 0 entity
  int subIndex (const typename remove_const<GridImp>::type::Traits::template Codim<0>::Entity& e, 
                int i,
                unsigned int codim) const
  {
      return grid_.getRealImplementation(e).subLeafIndex(i,codim);
  }

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
      if (codim==1) {

          return numVertices_;
          
      } else if (codim==0) {
          
          return numElements_;

      } 

      return 0;
  } 

  /** deliver all geometry types used in this grid */
  const std::vector<GeometryType>& geomTypes (int codim) const
  {
	return myTypes_[codim];
  }

        /** \brief Return true if e is contained in the index set.
        
    Warning: this implementation takes O(n) time!  It also assumes that e belongs
    to the correct grid.
    */
    template <class EntityType>
    bool contains (const EntityType& e) const
    {
        enum { cd = EntityType::codimension };
        typedef typename GridImp::template Codim<cd>::template Partition<All_Partition>::LeafIterator IteratorType; 
        IteratorType iend = grid_.template leafend<cd,All_Partition>();
        for (IteratorType it = grid_.template leafbegin<cd,All_Partition>();
             it != iend; ++it)
            {
                if (it->level() == e.level() && this->template index<cd>(*it) == this->template index<cd>(e)) 
                    return true;
            }
        return false; 
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

	  const OneDInNDEntityImp<0, GridImp::dimensionworld>* vIt;
            for (vIt = grid_.vertices[i].begin(); vIt!=grid_.vertices[i].end(); vIt = vIt->succ_) {
                
                /** \todo Remove the const casts */
                if (vIt->isLeaf())
		  const_cast<OneDInNDEntityImp<0, GridImp::dimensionworld>*>(vIt)->leafIndex_ = numVertices_++;
                else
		  const_cast<OneDInNDEntityImp<0, GridImp::dimensionworld>*>(vIt)->leafIndex_ = vIt->son_->leafIndex_;

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


template<class GridImp>
class OneDInNDGridIdSet : public IdSet<GridImp,OneDInNDGridIdSet<GridImp>,unsigned int>
{
public:
  //! define the type used for persistent indices
  typedef unsigned int IdType;

  //! constructor stores reference to a grid
  OneDInNDGridIdSet (const GridImp& g) : grid_(g) {}

  //! get id of an entity
  /*
    We use the remove_const to extract the Type from the mutable class,
    because the const class is not instantiated yet.
  */
  template<int cd>
  IdType id (const typename remove_const<GridImp>::type::Traits::template Codim<cd>::Entity& e) const 
  {
      return grid_.getRealImplementation(e).globalId();
  }

  //! get id of subentity
  /*
    We use the remove_const to extract the Type from the mutable class,
    because the const class is not instantiated yet.
  */
  template<int codim>
  IdType subId (const typename remove_const<GridImp>::type::Traits::template Codim<0>::Entity& e, int i) const
  {
      return subId(e,i,codim);
  }


  //! get id of subentity
  IdType subId (const typename remove_const<GridImp>::type::Traits::template Codim<0>::Entity& e, 
                int i,
                unsigned int codim) const
  {
      return grid_.getRealImplementation(e).subId(i,codim);
  }

private:

  const GridImp& grid_;
};

}  // namespace Dune


#endif
