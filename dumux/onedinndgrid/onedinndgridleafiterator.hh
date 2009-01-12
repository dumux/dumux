// $Id$

#ifndef DUNE_ONE_D_IN_N_D_GRID_LEAFITERATOR_HH
#define DUNE_ONE_D_IN_N_D_GRID_LEAFITERATOR_HH

#include "onedinndgridentitypointer.hh"

/** \file
 * \brief The OneDInNDGridLeafIterator class
 */

namespace Dune {

  /** \brief Iterator over all entities of a given codimension and level of a grid.
   * \ingroup OneDInNDGrid
   */
    template<int codim, PartitionIteratorType pitype, class GridImp>
class OneDInNDGridLeafIterator :
        public Dune::OneDInNDGridEntityPointer <codim,GridImp>
{
    enum {dim = GridImp::dimension};
    enum { dimworld = GridImp::dimensionworld };

    friend class OneDInNDGridEntity<codim,dim,GridImp>;

public:

    OneDInNDGridLeafIterator(const GridImp& grid) : grid_(&grid) {

        /** \todo Can a make the fullRefineLevel work somehow? */
        const int fullRefineLevel = 0;

        if (codim==0)
            this->virtualEntity_.setToTarget((OneDInNDEntityImp<1-codim, dimworld>*)grid_->elements[fullRefineLevel].begin);
        else
            this->virtualEntity_.setToTarget((OneDInNDEntityImp<1-codim, dimworld>*)grid_->vertices[fullRefineLevel].begin);

        if (!this->virtualEntity_.target()->isLeaf())
            increment();
    }

  //! Constructor
    OneDInNDGridLeafIterator()
    {
        this->virtualEntity_.setToTarget(NULL);
    }

    //! prefix increment
    void increment() {
        // Increment until you find a leaf entity
        do {
            globalIncrement();
        } while (this->virtualEntity_.target() && !this->virtualEntity_.target()->isLeaf());
    }

private:

    /** \brief This increment makes the iterator wander over all entities on all levels */
    void globalIncrement() {

        // Backup current level because it may not be accessible anymore after
        // setting the pointer to the next entity.
        const int oldLevel = this->virtualEntity_.level();

        // Increment on this level
        this->virtualEntity_.setToTarget(this->virtualEntity_.target()->succ_);

        // If beyond the end of this level set to first of next level
        if (!this->virtualEntity_.target() && oldLevel < grid_->maxLevel()) {

            if (codim==0)
                this->virtualEntity_.setToTarget((OneDInNDEntityImp<1-codim, dimworld>*)grid_->elements[oldLevel+1].begin);
            else
                this->virtualEntity_.setToTarget((OneDInNDEntityImp<1-codim, dimworld>*)grid_->vertices[oldLevel+1].begin);

        }

    }

    // /////////////////////////////////////
    //   Data members
    // /////////////////////////////////////
    const GridImp* grid_;

};

}  // namespace Dune

#endif
