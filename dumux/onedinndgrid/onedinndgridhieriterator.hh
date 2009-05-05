#ifndef DUNE_ONE_D_IN_N_D_GRID_HIERITERATOR_HH
#define DUNE_ONE_D_IN_N_D_GRID_HIERITERATOR_HH

/** \file
 * \brief The OneDGridHierarchicIterator class
 */

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
     Iteration over this set of entities is provided by the HierarchicIterator,
     starting from a given entity.
    */
  template<class GridImp>
  class OneDInNDGridHierarchicIterator :
    public Dune::OneDInNDGridEntityPointer <0,GridImp>
  {
    enum { dim = GridImp::dimension };
    enum { dimworld = GridImp::dimensionworld };
    friend class OneDInNDGridEntity<0,dim,GridImp>;
    
    // Stack entry
    typedef typename OneDInNDGridList<OneDInNDEntityImp<1,dimworld> >::iterator StackEntry;
    
  public:
    
    typedef typename GridImp::template Codim<0>::Entity Entity;
    
    //! Constructor
    OneDInNDGridHierarchicIterator(int maxlevel) : OneDInNDGridEntityPointer<0,GridImp>(OneDInNDGridNullIteratorFactory<1, dimworld>::null()),
                                               maxlevel_(maxlevel), elemStack()
    {}
    
    //! prefix increment
    void increment() {
      
      if (elemStack.empty())
	return;
      
      StackEntry old_target = elemStack.top();
      elemStack.pop();
      
      // Traverse the tree no deeper than maxlevel
      if (old_target->level_ < maxlevel_) {
	
	// Load sons of old target onto the iterator stack
	if (!old_target->isLeaf()) {
	  
	  elemStack.push(old_target->sons_[0]);
          
	  // Add the second son only if it is different from the first one
	  // i.e. the son is not just a copy of the father
	  if (old_target->sons_[0] != old_target->sons_[1])
	    elemStack.push(old_target->sons_[1]);
	  
	}
        
      }
      
      this->virtualEntity_.setToTarget((elemStack.empty()) 
				       ? OneDInNDGridNullIteratorFactory<1, dimworld>::null() : elemStack.top());
    }
    
  private:
    
    //! max level to go down 
    int maxlevel_;
    
    std::stack<StackEntry> elemStack;
    
  };
  
}  // end namespace Dune

#endif
