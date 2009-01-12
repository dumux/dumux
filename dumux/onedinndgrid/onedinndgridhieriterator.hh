// $Id$

#ifndef DUNE_ONE_D_IN_N_D_GRID_HIERITERATOR_HH
#define DUNE_ONE_D_IN_N_D_GRID_HIERITERATOR_HH

/** \file
 * \brief The OneDInNDGridHierarchicIterator class
 */

//#include <dune/common/stack.hh>
#include <stack>

namespace Dune {

//**********************************************************************
//
// --OneDInNDGridHierarchicIterator
// --HierarchicIterator
  /** \brief Iterator over the descendants of an entity.
   * \ingroup OneDInNDGrid
  Mesh entities of codimension 0 ("elements") allow to visit all entities of
  codimension 0 obtained through nested, hierarchic refinement of the entity.
  Iteration over this set of entities is provided by the HIerarchicIterator,
  starting from a given entity.
  This is redundant but important for memory efficient implementations of unstru
  hierarchically refined meshes.
 */
template<class GridImp>
class OneDInNDGridHierarchicIterator :
        public Dune::OneDInNDGridEntityPointer <0,GridImp>,
        public HierarchicIteratorDefaultImplementation <GridImp, OneDInNDGridHierarchicIterator>
{
    enum { dim = GridImp::dimension };
    enum { dimworld = GridImp::dimensionworld };
    friend class OneDInNDGridEntity<0,dim,GridImp>;

    // Stack entry
    struct StackEntry {
        OneDInNDEntityImp<1, dimworld>* element;
        /** \todo Do we need the level ? */
        int level;
    };

public:

    typedef typename GridImp::template Codim<0>::Entity Entity;

  //! Constructor
    OneDInNDGridHierarchicIterator(int maxlevel) : OneDInNDGridEntityPointer<0,GridImp>(NULL),
                                               maxlevel_(maxlevel), elemStack()
    {}

    //! prefix increment
    void increment() {

        if (elemStack.empty())
            return;

        StackEntry old_target = elemStack.pop();

        // Traverse the tree no deeper than maxlevel
        if (old_target.level < maxlevel_) {

            // Load sons of old target onto the iterator stack
            if (!old_target.element->isLeaf()) {
                StackEntry se0;
                se0.element = old_target.element->sons_[0];
                se0.level   = old_target.level + 1;
                elemStack.push(se0);

                // Add the second son only if it is different from the first one
                // i.e. the son is not just a copy of the father
                if (old_target.element->sons_[0] != old_target.element->sons_[1]) {
                    StackEntry se1;
                    se1.element = old_target.element->sons_[1];
                    se1.level   = old_target.level + 1;
                    elemStack.push(se1);
                }
            }

        }

        this->virtualEntity_.setToTarget((elemStack.empty()) ? NULL : elemStack.top().element);
    }

private:

  //! max level to go down
  int maxlevel_;

  std::stack<StackEntry> elemStack;

};

}  // end namespace Dune

#endif
