// $Id: onedinndgridleveliterator.hh 972 2009-01-12 10:15:57Z lauser $

#ifndef DUNE_ONE_D_IN_N_D_GRID_LEVELITERATOR_HH
#define DUNE_ONE_D_IN_N_D_GRID_LEVELITERATOR_HH

/** \file
 * \brief The OneDInNDGridLevelIterator class
 */

#include <dune/common/dlist.hh>

namespace Dune {



//**********************************************************************
//
// --OneDInNDGridLevelIterator
// --LevelIterator
  /** \brief Iterator over all entities of a given codimension and level of a grid.
   * \ingroup OneDInNDGrid
   */
template<int codim, PartitionIteratorType pitype, class GridImp>
class OneDInNDGridLevelIterator :
    public OneDInNDGridEntityPointer <codim, GridImp>,
    public LevelIteratorDefaultImplementation <codim, pitype, GridImp, OneDInNDGridLevelIterator>
{
public:
    enum {dim=GridImp::dimension};
    enum { dimworld = GridImp::dimensionworld };
    friend class OneDInNDGridLevelIteratorFactory<codim>;
    template <int dimworld>
    friend class OneDInNDGrid;
    friend class OneDInNDGridEntity<codim,dim,GridImp>;
    friend class OneDInNDGridEntity<0,dim,GridImp>;

    typedef typename GridImp::template Codim<codim>::Entity Entity;

protected:

    /** \brief Constructor from a given iterator */
    OneDInNDGridLevelIterator<codim,pitype, GridImp>(OneDInNDEntityImp<dim-codim, dimworld>* it)
      : OneDInNDGridEntityPointer<codim, GridImp>(it)
    {
    }

public:

    //! prefix increment
    void increment() {
        this->virtualEntity_.setToTarget(this->virtualEntity_.target()->succ_);
    }
};

}  // namespace Dune

#endif
