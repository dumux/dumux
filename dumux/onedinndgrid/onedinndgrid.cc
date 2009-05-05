#include "config.h"

#include "onedinndgrid.hh"

// ///////////////////////////////////////////////////////////////
//
//    OneDInNDGridLevelIteratorFactory, a class used to simulate
//    specialization of member templates
//
// ///////////////////////////////////////////////////////////////

namespace Dune {

template <int codim>
struct OneDInNDGridLevelIteratorFactory {};

template <>
struct OneDInNDGridLevelIteratorFactory<1>
{
  template <Dune::PartitionIteratorType PiType, int dimworld>
    static Dune::OneDInNDGridLevelIterator<1,PiType, const OneDInNDGrid<dimworld> > 
    lbegin(const Dune::OneDInNDGrid<dimworld>* g, int level) {

        return Dune::OneDInNDGridLevelIterator<1,PiType, const Dune::OneDInNDGrid<dimworld> >
	  (const_cast<Dune::OneDInNDEntityImp<0,dimworld>* >(g->vertices[level].begin()));
    }

};

template <>
struct OneDInNDGridLevelIteratorFactory<0>
{
  template <Dune::PartitionIteratorType PiType, int dimworld>
    static Dune::OneDInNDGridLevelIterator<0,PiType, const Dune::OneDInNDGrid<dimworld> >
    lbegin(const Dune::OneDInNDGrid<dimworld>* g, int level) {

        return Dune::OneDInNDGridLevelIterator<0,PiType, const Dune::OneDInNDGrid<dimworld> >
	  (const_cast<Dune::OneDInNDEntityImp<1,dimworld>* >(g->elements[level].begin()));
    }

};

}

template <int dimworld>
Dune::OneDInNDGrid<dimworld>::OneDInNDGrid(int numElements, const ctype& leftBoundary, const ctype& rightBoundary)
    : refinementType_(LOCAL),
      leafIndexSet_(*this),
      idSet_(*this),
      freeVertexIdCounter_(0),
      freeElementIdCounter_(0)
{
    if (numElements<1)
        DUNE_THROW(GridError, "Nonpositive number of elements requested!");

    if (leftBoundary >= rightBoundary)
        DUNE_THROW(GridError, "The left boundary coordinate has to be strictly less than the right boundary one!");
    // Init grid hierarchy
    vertices.resize(1);
    elements.resize(1);

    // Init vertex set
    for (int i=0; i<numElements+1; i++) {
        ctype newCoord = leftBoundary + i*(rightBoundary-leftBoundary) / numElements;

        OneDInNDEntityImp<0, dimworld> newVertex(0, newCoord);
        newVertex.id_ = getNextFreeId(1);
        vertices[0].push_back(newVertex);

    }

    // Init element set
    typename OneDInNDGridList<OneDInNDEntityImp<0, dimworld> >::iterator it = vertices[0].begin();
    for (int i=0; i<numElements; i++) {

      OneDInNDEntityImp<1, dimworld> newElement(0, getNextFreeId(0));
        newElement.vertex_[0] = it;
        it = it->succ_;
        newElement.vertex_[1] = it;

        elements[0].push_back(newElement);

    }

    setIndices();
}
// template<int dimworld>
// Dune::OneDInNDGrid::OneDInNDGrid(const std::vector<ctype>& coords) 
//     : refinementType_(LOCAL),
//       leafIndexSet_(*this),
//       idSet_(*this),
//       freeVertexIdCounter_(0),
//       freeElementIdCounter_(0)
// {
//     if (coords.size()<2)
//         DUNE_THROW(GridError, "You have to provide at least two coordinates!");

//     // Init grid hierarchy
//     vertices.resize(1);
//     elements.resize(1);

//     // Init vertex set
//     for (size_t i=0; i<coords.size(); i++) {
//       OneDInNDEntityImp<0, dimworld> newVertex(0, coords[i], getNextFreeId(1));
//         vertices[0].push_back(newVertex);
//     }

//     // Init element set
//     OneDInNDGridList<OneDInNDEntityImp<0, dimworld> >::iterator it = vertices[0].begin();
//     for (size_t i=0; i<coords.size()-1; i++) {

//       OneDInNDEntityImp<1, dimworld> newElement(0, getNextFreeId(0));
//         newElement.vertex_[0] = it;
//         it = it->succ_;
//         newElement.vertex_[1] = it;

//         if (newElement.vertex_[0]->pos_ >= newElement.vertex_[1]->pos_)
//             DUNE_THROW(GridError, "The coordinates have to be in ascending order!");

//         elements[0].push_back(newElement);

//     }

//     setIndices();
// }

template<int dimworld>
Dune::OneDInNDGrid<dimworld>::~OneDInNDGrid()
{
    // Delete all vertices
    for (unsigned int i=0; i<vertices.size(); i++) {

      typename OneDInNDGridList<OneDInNDEntityImp<0, dimworld> >::iterator v = vertices[i].begin();

        while (v) {
            
	  typename OneDInNDGridList<OneDInNDEntityImp<0, dimworld> >::iterator vSucc = v->succ_;
            vertices[i].erase(v);
            v = vSucc;

        }

    }

    // Delete all elements
    for (unsigned int i=0; i<elements.size(); i++) {

      typename OneDInNDGridList<OneDInNDEntityImp<1, dimworld> >::iterator e = elements[i].begin();

        while (e) {
            
	  typename OneDInNDGridList<OneDInNDEntityImp<1, dimworld> >::iterator eSucc = e->succ_;
            elements[i].erase(e);
            e = eSucc;

        }

    }

    // Delete levelIndexSets
    for (unsigned int i=0; i<levelIndexSets_.size(); i++)
      if (levelIndexSets_[i])
        delete levelIndexSets_[i];
}

template<int dimworld>
template <int codim>
typename Dune::OneDInNDGrid<dimworld>::Traits::template Codim<codim>::LevelIterator
Dune::OneDInNDGrid<dimworld>::lbegin(int level) const
{
    if (level<0 || level>maxLevel())
        DUNE_THROW(Dune::GridError, "LevelIterator in nonexisting level " << level << " requested!");

    return OneDInNDGridLevelIteratorFactory<codim>::template lbegin<All_Partition>(this, level);
}

template<int dimworld>
template <int codim>
typename Dune::OneDInNDGrid<dimworld>::Traits::template Codim<codim>::LevelIterator
Dune::OneDInNDGrid<dimworld>::lend(int level) const
{
    if (level<0 || level>maxLevel())
        DUNE_THROW(GridError, "LevelIterator in nonexisting level " << level << " requested!");
    
    return OneDInNDGridLevelIterator<codim,All_Partition, const Dune::OneDInNDGrid<dimworld> >(0);
}

template<int dimworld>
template <int codim, Dune::PartitionIteratorType PiType>
typename Dune::OneDInNDGrid<dimworld>::Traits::template Codim<codim>::template Partition<PiType>::LevelIterator
Dune::OneDInNDGrid<dimworld>::lbegin(int level) const
{
    if (level<0 || level>maxLevel())
        DUNE_THROW(Dune::GridError, "LevelIterator in nonexisting level " << level << " requested!");

    return OneDInNDGridLevelIteratorFactory<codim>::template lbegin<PiType>(this, level);
}

template<int dimworld>
template <int codim, Dune::PartitionIteratorType PiType>
typename Dune::OneDInNDGrid<dimworld>::Traits::template Codim<codim>::template Partition<PiType>::LevelIterator
Dune::OneDInNDGrid<dimworld>::lend(int level) const
{
    if (level<0 || level>maxLevel())
        DUNE_THROW(GridError, "LevelIterator in nonexisting level " << level << " requested!");
    
    return OneDInNDGridLevelIterator<codim,PiType, const Dune::OneDInNDGrid<dimworld> >(static_cast<Dune::OneDInNDEntityImp<dim-codim,dimworld>*>(0));
}

template<int dimworld>
template <int codim>
typename Dune::OneDInNDGrid<dimworld>::Traits::template Codim<codim>::LeafIterator
Dune::OneDInNDGrid<dimworld>::leafbegin() const
{
    return OneDInNDGridLeafIterator<codim,All_Partition,const OneDInNDGrid<dimworld> >(*this);
}

template<int dimworld>
template <int codim>
typename Dune::OneDInNDGrid<dimworld>::Traits::template Codim<codim>::LeafIterator
Dune::OneDInNDGrid<dimworld>::leafend() const
{
    return OneDInNDGridLeafIterator<codim,All_Partition, const Dune::OneDInNDGrid<dimworld> >();
}

template<int dimworld>
template <int codim, Dune::PartitionIteratorType PiType>
typename Dune::OneDInNDGrid<dimworld>::Traits::template Codim<codim>::template Partition<PiType>::LeafIterator
Dune::OneDInNDGrid<dimworld>::leafbegin() const
{
    return OneDInNDGridLeafIterator<codim,PiType, const Dune::OneDInNDGrid<dimworld> >(*this);
}

template<int dimworld>
template <int codim, Dune::PartitionIteratorType PiType>
typename Dune::OneDInNDGrid<dimworld>::Traits::template Codim<codim>::template Partition<PiType>::LeafIterator
Dune::OneDInNDGrid<dimworld>::leafend() const
{
    return OneDInNDGridLeafIterator<codim,PiType, const Dune::OneDInNDGrid<dimworld> >();
}

template<int dimworld>
typename Dune::OneDInNDGridList<Dune::OneDInNDEntityImp<0, dimworld> >::iterator
Dune::OneDInNDGrid<dimworld>::getLeftUpperVertex(const OneDInNDEntityImp<1, dimworld>* eIt)
{
  typename OneDInNDGridList<OneDInNDEntityImp<1, dimworld> >::iterator l = eIt->pred_;

    if (!l)
        return 0;

    // return NULL if there is no geometrical left neighbor
    if (l->vertex_[1] != eIt->vertex_[0])
        return 0;

    // return NULL if that neighbor doesn't have sons
    if (l->isLeaf())
        return 0;

    // return the right vertex of the right son
    return l->sons_[1]->vertex_[1];

}

template<int dimworld>
typename Dune::OneDInNDGridList<Dune::OneDInNDEntityImp<0, dimworld> >::iterator
Dune::OneDInNDGrid<dimworld>::getRightUpperVertex(const OneDInNDEntityImp<1, dimworld>* eIt)
{
  typename OneDInNDGridList<OneDInNDEntityImp<1, dimworld> >::iterator r = eIt->succ_;

    if (!r)
        return 0;

    // return NULL if there is no geometrical right neighbor
    if (r->vertex_[0]!=eIt->vertex_[1])
        return 0;

    // return NULL if that neighbor doesn't have sons
    if (r->isLeaf())
        return 0;

    // return the left vertex of the left son
    return r->sons_[0]->vertex_[0];

}

template<int dimworld>
typename Dune::OneDInNDGridList<Dune::OneDInNDEntityImp<1, dimworld> >::iterator 
Dune::OneDInNDGrid<dimworld>::getLeftNeighborWithSon(typename OneDInNDGridList<OneDInNDEntityImp<1, dimworld> >::iterator eIt)
{
  typename OneDInNDGridList<OneDInNDEntityImp<1, dimworld> >::iterator l = eIt;

    do {
        l = l->pred_;
    } while (l && l->isLeaf());

    return l;
}

template<int dimworld>
bool Dune::OneDInNDGrid<dimworld>::adapt()
{
  typename OneDInNDGridList<OneDInNDEntityImp<1, dimworld> >::iterator eIt;
    
    // for the return value:  true if the grid was changed
    bool changedGrid = false;

    // remove elements that have been marked for coarsening.
    // If one son of an element is marked for coarsening, and the other one is not,
    // then the element is not removed.
    for (int i=1; i<=maxLevel(); i++) {

        for (eIt = elements[i].begin(); eIt!=elements[i].end(); ) {

	  typename OneDInNDGridList<OneDInNDEntityImp<1, dimworld> >::iterator leftElementToBeDeleted  = eIt;
	  typename OneDInNDGridList<OneDInNDEntityImp<1, dimworld> >::iterator rightElementToBeDeleted = eIt->succ_;

            assert(eIt->succ_);
            typename OneDInNDGridList<OneDInNDEntityImp<1, dimworld> >::iterator nextElement = eIt->succ_->succ_;

            if (leftElementToBeDeleted->markState_ == OneDInNDEntityImp<1, dimworld>::COARSEN && leftElementToBeDeleted->isLeaf()
                && rightElementToBeDeleted->markState_ == OneDInNDEntityImp<1, dimworld>::COARSEN && rightElementToBeDeleted->isLeaf()) {

                assert(rightElementToBeDeleted->isLeaf());

                // Is the left vertex obsolete?
                if (leftElementToBeDeleted->pred_==NULL 
                    || leftElementToBeDeleted->pred_->vertex_[1] != leftElementToBeDeleted->vertex_[0]){

                    // If the left vertex has a father remove the reference to this vertex at this father
                    assert(leftElementToBeDeleted->father_->vertex_[0]->son_ == leftElementToBeDeleted->vertex_[0]);
                    leftElementToBeDeleted->father_->vertex_[0]->son_ = NULL;

                    vertices[i].erase(leftElementToBeDeleted->vertex_[0]);
                }

                // Is the right vertex obsolete?
                if (rightElementToBeDeleted->succ_==NULL 
                    || rightElementToBeDeleted->succ_->vertex_[0] != rightElementToBeDeleted->vertex_[1]) {

                    // If the left vertex has a father remove the reference to this vertex at this father
                    assert(rightElementToBeDeleted->father_->vertex_[1]->son_ == rightElementToBeDeleted->vertex_[1]);
                    rightElementToBeDeleted->father_->vertex_[1]->son_ = NULL;

                    vertices[i].erase(rightElementToBeDeleted->vertex_[1]);
                }

                // Delete vertex between left and right element to be deleted
                assert(leftElementToBeDeleted->vertex_[1] == rightElementToBeDeleted->vertex_[0]);
                vertices[i].erase(leftElementToBeDeleted->vertex_[1]);

                // Remove references from the father element
                assert(rightElementToBeDeleted->father_->sons_[1] == rightElementToBeDeleted);
                leftElementToBeDeleted->father_->sons_[0]  = NULL;
                rightElementToBeDeleted->father_->sons_[1] = NULL;

                // Paranoia: make sure the father is not marked for refinement
                rightElementToBeDeleted->father_->markState_ = OneDInNDEntityImp<1, dimworld>::DO_NOTHING;

                // Actually delete elements
                elements[i].erase(leftElementToBeDeleted);
                elements[i].erase(rightElementToBeDeleted);

                // The grid has been changed
                changedGrid = true;
            }

            // increment pointer
            eIt = nextElement;

        }

    }
    
    // /////////////////////////////////////////////////////////////////////////
    //  Check if one of the elements on the toplevel is marked for refinement
    //  In that case add another level
    // /////////////////////////////////////////////////////////////////////////
    bool toplevelRefinement = false;
    for (eIt = elements[maxLevel()].begin(); eIt!=elements[maxLevel()].end(); eIt=eIt->succ_) 
      if (eIt->markState_ == OneDInNDEntityImp<1, dimworld>::REFINE) {
            toplevelRefinement = true;
            break;
        }

    if (toplevelRefinement) {
      OneDInNDGridList<OneDInNDEntityImp<0, dimworld> > newVertices;
      OneDInNDGridList<OneDInNDEntityImp<1, dimworld> > newElements;
        vertices.push_back(newVertices);
        elements.push_back(newElements);
    }

    // //////////////////////////////
    // refine all marked elements
    // //////////////////////////////
    int oldMaxlevel = (toplevelRefinement) ? maxLevel()-1 : maxLevel();
    for (int i=0; i<=oldMaxlevel; i++) {

        for (eIt = elements[i].begin(); eIt!=elements[i].end(); eIt = eIt->succ_) {

	  if (eIt->markState_ == OneDInNDEntityImp<1, dimworld>::REFINE
                && eIt->isLeaf()) {

                // Does the left vertex exist on the next-higher level?
                // If no create it
	    typename OneDInNDGridList<OneDInNDEntityImp<0, dimworld> >::iterator leftUpperVertex = getLeftUpperVertex(eIt);

                if (leftUpperVertex==NULL)
		  leftUpperVertex = new OneDInNDEntityImp<0, dimworld>(i+1, 
                                                           eIt->vertex_[0]->pos_, 
                                                           eIt->vertex_[0]->id_);

                eIt->vertex_[0]->son_ = leftUpperVertex;

                // Does the right vertex exist on the next-higher level?
                // If no create it
                typename OneDInNDGridList<OneDInNDEntityImp<0, dimworld> >::iterator rightUpperVertex = getRightUpperVertex(eIt);
                
                if (rightUpperVertex==NULL) 
		  rightUpperVertex = new OneDInNDEntityImp<0, dimworld>(i+1, 
                                                            eIt->vertex_[1]->pos_, 
                                                            eIt->vertex_[1]->id_);

                eIt->vertex_[1]->son_ = rightUpperVertex;

                // Create center vertex
                ctype p = 0.5*(eIt->vertex_[0]->pos_[0] + eIt->vertex_[1]->pos_[0]);

                OneDInNDEntityImp<0, dimworld> centerVertex(i+1, p, getNextFreeId(1));

                // //////////////////////////////////////
                // Insert new vertices into vertex list
                // //////////////////////////////////////

                typename OneDInNDGridList<OneDInNDEntityImp<1, dimworld> >::iterator leftNeighbor = getLeftNeighborWithSon(eIt);

                if (leftNeighbor!=NULL) {

                    // leftNeighbor exists
                    if ( leftNeighbor->sons_[1]->vertex_[1] != leftUpperVertex)
                        vertices[i+1].insert(leftNeighbor->sons_[1]->vertex_[1]->succ_, leftUpperVertex);

                } else {
                    // leftNeighbor does not exist
                    vertices[i+1].insert(vertices[i+1].begin(), leftUpperVertex);

                }

                typename OneDInNDGridList<OneDInNDEntityImp<0, dimworld> >::iterator centerVertexIterator = vertices[i+1].insert(leftUpperVertex->succ_, centerVertex);

                // Check if rightUpperVertex is already in the list
                typename OneDInNDGridList<OneDInNDEntityImp<0, dimworld> >::iterator succOfCenter = centerVertexIterator->succ_;

                if (succOfCenter==NULL || succOfCenter != rightUpperVertex)
                    vertices[i+1].insert(centerVertexIterator->succ_, rightUpperVertex);

                // ///////////////////////
                // Create new elements
                // ///////////////////////
                OneDInNDEntityImp<1, dimworld> newElement0(i+1, getNextFreeId(0));
                newElement0.vertex_[0] = leftUpperVertex;
                newElement0.vertex_[1] = centerVertexIterator;
                newElement0.father_ = eIt;
                newElement0.isNew_ = true;

                OneDInNDEntityImp<1, dimworld> newElement1(i+1, getNextFreeId(0));
                newElement1.vertex_[0] = centerVertexIterator;
                newElement1.vertex_[1] = rightUpperVertex;
                newElement1.father_ = eIt;
                newElement1.isNew_ = true;

                // Insert new elements into element list
                if (leftNeighbor!=NULL)
                    // leftNeighbor exists
                    eIt->sons_[0] = elements[i+1].insert(leftNeighbor->sons_[1]->succ_, newElement0);
                else 
                    // leftNeighbor does not exist
                    eIt->sons_[0] = elements[i+1].insert(elements[i+1].begin(), newElement0);

                eIt->sons_[1] = elements[i+1].insert(eIt->sons_[0]->succ_, newElement1);

                // The grid has been modified
                changedGrid = true;

            }

        }

    }

    // delete uppermost level if it doesn't contain elements anymore
    if (elements[maxLevel()].size()==0) {
        assert(vertices[maxLevel()].size()==0);
        elements.pop_back();
        vertices.pop_back();
    }


    // If the refinement mode is 'COPY', fill the empty slots in the grid
    // by copying elements
    if (refinementType_ == COPY) {
        
        for (int i=0; i<maxLevel(); i++) {
            
	  typename OneDInNDGridList<OneDInNDEntityImp<1, dimworld> >::iterator eIt;
            for (eIt = elements[i].begin(); eIt!=elements[i].end(); eIt = eIt->succ_) {
                
                if (eIt->isLeaf()) {
                    
                    // Does the left vertex exist on the next-higher level?
                    // If no create it
		  typename OneDInNDGridList<OneDInNDEntityImp<0, dimworld> >::iterator leftUpperVertex = getLeftUpperVertex(eIt);
                    
                    if (leftUpperVertex==NULL)
		      leftUpperVertex = new OneDInNDEntityImp<0, dimworld>(i+1, eIt->vertex_[0]->pos_, eIt->vertex_[0]->id_);
 
                    eIt->vertex_[0]->son_ = leftUpperVertex;
                   
                    // Does the right vertex exist on the next-higher level?
                    // If no create it
                    typename OneDInNDGridList<OneDInNDEntityImp<0, dimworld> >::iterator rightUpperVertex = getRightUpperVertex(eIt);
                    
                    if (rightUpperVertex==NULL) 
		      rightUpperVertex = new OneDInNDEntityImp<0, dimworld>(i+1, eIt->vertex_[1]->pos_, eIt->vertex_[1]->id_);

                    eIt->vertex_[1]->son_ = rightUpperVertex;

                    // //////////////////////////////////////
                    // Insert new vertices into vertex list
                    // //////////////////////////////////////
                    
                    typename OneDInNDGridList<OneDInNDEntityImp<1, dimworld> >::iterator leftNeighbor = getLeftNeighborWithSon(eIt);
                    
                    if (leftNeighbor!=NULL) {
                        
                        // leftNeighbor exists
                        if ( leftNeighbor->sons_[1]->vertex_[1] != leftUpperVertex)
                            vertices[i+1].insert(leftNeighbor->sons_[1]->vertex_[1]->succ_, leftUpperVertex);
                        
                    } else {
                        // leftNeighbor does not exist
                        vertices[i+1].insert(vertices[i+1].begin(), leftUpperVertex);

                    }

                    // Check if rightUpperVertex is already in the list
                    typename OneDInNDGridList<OneDInNDEntityImp<0, dimworld> >::iterator succOfLeft = leftUpperVertex->succ_;

                    if (succOfLeft==NULL || succOfLeft != rightUpperVertex)
                        vertices[i+1].insert(leftUpperVertex->succ_, rightUpperVertex);
                    
                    // /////////////////////////
                    //   Create new element
                    // /////////////////////////
                    OneDInNDEntityImp<1, dimworld> newElement(i+1, eIt->id_);
                    newElement.vertex_[0] = leftUpperVertex;
                    newElement.vertex_[1] = rightUpperVertex;
                    newElement.father_ = eIt;
                    newElement.isNew_ = true;

                    // Insert new elements into element list
                    typename OneDInNDGridList<OneDInNDEntityImp<1, dimworld> >::iterator newElementIterator;
                    if (leftNeighbor!=NULL)
                        // leftNeighbor exists
                        newElementIterator = elements[i+1].insert(leftNeighbor->sons_[1]->succ_, newElement);
                    else 
                        // leftNeighbor does not exist
                        newElementIterator = elements[i+1].insert(elements[i+1].begin(), newElement);
                    
                    // Mark the new element as the sons of the refined element
                    eIt->sons_[0] = eIt->sons_[1] = newElementIterator;
                    
                }
                
            }
            
        }

    }

    // ////////////////////////////////////
    //   renumber vertices and elements
    // ////////////////////////////////////
    setIndices();

    return changedGrid;
}

template<int dimworld>
bool Dune::OneDInNDGrid<dimworld>::preAdapt()
{
    typename Traits::template Codim<0>::LeafIterator eIt    = leafbegin<0>();
    typename Traits::template Codim<0>::LeafIterator eEndIt = leafend<0>();

    for (; eIt!=eEndIt; ++eIt)
      if (getRealImplementation(*eIt).target_->markState_ != OneDInNDEntityImp<1, dimworld>::DO_NOTHING)
            return true;

    return false;
}

template<int dimworld>
void Dune::OneDInNDGrid<dimworld>::postAdapt()
{
    for (int i=0; i<=maxLevel(); i++) {
      typename OneDInNDGridList<OneDInNDEntityImp<1, dimworld> >::iterator eIt;
        for (eIt = elements[i].begin(); eIt!=elements[i].end(); eIt = eIt->succ_)
	  eIt->markState_ = OneDInNDEntityImp<1, dimworld>::DO_NOTHING;

    }
        
}

template<int dimworld>
void Dune::OneDInNDGrid<dimworld>::setIndices()
{
    // Add space for new LevelIndexSets if the grid hierarchy got higher
    // They are not created until they are actually requested
    for (int i=levelIndexSets_.size(); i<maxLevel()+1; i++)
        levelIndexSets_.push_back( (OneDInNDGridLevelIndexSet< const OneDInNDGrid<dimworld> > *)0 );

    // Delete old LevelIndexSets if the grid hierarchy got lower
    int excess = levelIndexSets_.size() - (maxLevel() + 1);
    for (int i=0; i<excess; i++) {
        if (levelIndexSets_.back())
            delete(levelIndexSets_.back());
        levelIndexSets_.pop_back();
    }

    for (int i=0; i<=maxLevel(); i++)
      if (levelIndexSets_[i])
        levelIndexSets_[i]->update();

    leafIndexSet_.update();

    // IdSets don't need updating
}

template<int dimworld>
void Dune::OneDInNDGrid<dimworld>::globalRefine(int refCount)
{
  for (int i=0; i<refCount; i++) {

      // mark all entities for grid refinement
      typename Traits::template Codim<0>::LeafIterator iIt    = leafbegin<0>();
      typename Traits::template Codim<0>::LeafIterator iEndIt = leafend<0>();
      
      for (; iIt!=iEndIt; ++iIt)
          mark(1, *iIt);
	  
      this->preAdapt();
      adapt();
      this->postAdapt();
  }
}

template<int dimworld>
bool Dune::OneDInNDGrid<dimworld>::mark(int refCount,
                          const typename Traits::template Codim<0>::Entity & e )
{
    if (refCount < 0) {

        if (getRealImplementation(e).target_->level_ == 0)
            return false;
        else {
	  getRealImplementation(e).target_->markState_ = OneDInNDEntityImp<1, dimworld>::COARSEN;
            return true;
        }

    } else if (refCount > 0)
      getRealImplementation(e).target_->markState_ = OneDInNDEntityImp<1,dimworld>::REFINE;
    else
      getRealImplementation(e).target_->markState_ = OneDInNDEntityImp<1,dimworld>::DO_NOTHING;

    return true;
}

template<int dimworld>
int Dune::OneDInNDGrid<dimworld>::getMark(const typename Traits::template Codim<0>::Entity & e ) const 
{
  if(getRealImplementation(e).target_->markState_ == OneDInNDEntityImp<1, dimworld>::COARSEN) 
      return -1; 
  else if(getRealImplementation(e).target_->markState_ == OneDInNDEntityImp<1, dimworld>::REFINE) 
      return 1; 
  return 0;  
}



// /////////////////////////////////////////////////////////////////////////
//   Explicitly instantiate the OnedGrid member templates.
//   gcc-4.0 wants these instantiations after the method implementations
// /////////////////////////////////////////////////////////////////////////
template Dune::OneDInNDGrid<2>::Codim<0>::LevelIterator Dune::OneDInNDGrid<2>::lbegin<0>(int level) const;
template Dune::OneDInNDGrid<2>::Codim<1>::LevelIterator Dune::OneDInNDGrid<2>::lbegin<1>(int level) const;

template Dune::OneDInNDGrid<2>::Codim<0>::LevelIterator Dune::OneDInNDGrid<2>::lend<0>(int level) const;
template Dune::OneDInNDGrid<2>::Codim<1>::LevelIterator Dune::OneDInNDGrid<2>::lend<1>(int level) const;

template Dune::OneDInNDGrid<2>::Codim<0>::Partition<Dune::Interior_Partition>::LevelIterator 
Dune::OneDInNDGrid<2>::lbegin<0,Dune::Interior_Partition>(int level) const;
template Dune::OneDInNDGrid<2>::Codim<0>::Partition<Dune::InteriorBorder_Partition>::LevelIterator 
Dune::OneDInNDGrid<2>::lbegin<0,Dune::InteriorBorder_Partition>(int level) const;
template Dune::OneDInNDGrid<2>::Codim<0>::Partition<Dune::Overlap_Partition>::LevelIterator
 Dune::OneDInNDGrid<2>::lbegin<0,Dune::Overlap_Partition>(int level) const;
template Dune::OneDInNDGrid<2>::Codim<0>::Partition<Dune::OverlapFront_Partition>::LevelIterator
Dune::OneDInNDGrid<2>::lbegin<0,Dune::OverlapFront_Partition>(int level) const;
template Dune::OneDInNDGrid<2>::Codim<0>::Partition<Dune::All_Partition>::LevelIterator
Dune::OneDInNDGrid<2>::lbegin<0,Dune::All_Partition>(int level) const;
template Dune::OneDInNDGrid<2>::Codim<0>::Partition<Dune::Ghost_Partition>::LevelIterator
Dune::OneDInNDGrid<2>::lbegin<0,Dune::Ghost_Partition>(int level) const;

template Dune::OneDInNDGrid<2>::Codim<0>::Partition<Dune::Interior_Partition>::LevelIterator 
Dune::OneDInNDGrid<2>::lend<0,Dune::Interior_Partition>(int level) const;
template Dune::OneDInNDGrid<2>::Codim<0>::Partition<Dune::InteriorBorder_Partition>::LevelIterator 
Dune::OneDInNDGrid<2>::lend<0,Dune::InteriorBorder_Partition>(int level) const;
template Dune::OneDInNDGrid<2>::Codim<0>::Partition<Dune::Overlap_Partition>::LevelIterator
 Dune::OneDInNDGrid<2>::lend<0,Dune::Overlap_Partition>(int level) const;
template Dune::OneDInNDGrid<2>::Codim<0>::Partition<Dune::OverlapFront_Partition>::LevelIterator
Dune::OneDInNDGrid<2>::lend<0,Dune::OverlapFront_Partition>(int level) const;
template Dune::OneDInNDGrid<2>::Codim<0>::Partition<Dune::All_Partition>::LevelIterator
Dune::OneDInNDGrid<2>::lend<0,Dune::All_Partition>(int level) const;
template Dune::OneDInNDGrid<2>::Codim<0>::Partition<Dune::Ghost_Partition>::LevelIterator
Dune::OneDInNDGrid<2>::lend<0,Dune::Ghost_Partition>(int level) const;

template Dune::OneDInNDGrid<2>::Codim<1>::Partition<Dune::Interior_Partition>::LevelIterator 
Dune::OneDInNDGrid<2>::lbegin<1,Dune::Interior_Partition>(int level) const;
template Dune::OneDInNDGrid<2>::Codim<1>::Partition<Dune::InteriorBorder_Partition>::LevelIterator 
Dune::OneDInNDGrid<2>::lbegin<1,Dune::InteriorBorder_Partition>(int level) const;
template Dune::OneDInNDGrid<2>::Codim<1>::Partition<Dune::Overlap_Partition>::LevelIterator
 Dune::OneDInNDGrid<2>::lbegin<1,Dune::Overlap_Partition>(int level) const;
template Dune::OneDInNDGrid<2>::Codim<1>::Partition<Dune::OverlapFront_Partition>::LevelIterator
Dune::OneDInNDGrid<2>::lbegin<1,Dune::OverlapFront_Partition>(int level) const;
template Dune::OneDInNDGrid<2>::Codim<1>::Partition<Dune::All_Partition>::LevelIterator
Dune::OneDInNDGrid<2>::lbegin<1,Dune::All_Partition>(int level) const;
template Dune::OneDInNDGrid<2>::Codim<1>::Partition<Dune::Ghost_Partition>::LevelIterator
Dune::OneDInNDGrid<2>::lbegin<1,Dune::Ghost_Partition>(int level) const;

template Dune::OneDInNDGrid<2>::Codim<1>::Partition<Dune::Interior_Partition>::LevelIterator 
Dune::OneDInNDGrid<2>::lend<1,Dune::Interior_Partition>(int level) const;
template Dune::OneDInNDGrid<2>::Codim<1>::Partition<Dune::InteriorBorder_Partition>::LevelIterator 
Dune::OneDInNDGrid<2>::lend<1,Dune::InteriorBorder_Partition>(int level) const;
template Dune::OneDInNDGrid<2>::Codim<1>::Partition<Dune::Overlap_Partition>::LevelIterator
 Dune::OneDInNDGrid<2>::lend<1,Dune::Overlap_Partition>(int level) const;
template Dune::OneDInNDGrid<2>::Codim<1>::Partition<Dune::OverlapFront_Partition>::LevelIterator
Dune::OneDInNDGrid<2>::lend<1,Dune::OverlapFront_Partition>(int level) const;
template Dune::OneDInNDGrid<2>::Codim<1>::Partition<Dune::All_Partition>::LevelIterator
Dune::OneDInNDGrid<2>::lend<1,Dune::All_Partition>(int level) const;
template Dune::OneDInNDGrid<2>::Codim<1>::Partition<Dune::Ghost_Partition>::LevelIterator
Dune::OneDInNDGrid<2>::lend<1,Dune::Ghost_Partition>(int level) const;

template Dune::OneDInNDGrid<2>::Codim<0>::LeafIterator Dune::OneDInNDGrid<2>::leafbegin<0>() const;
template Dune::OneDInNDGrid<2>::Codim<1>::LeafIterator Dune::OneDInNDGrid<2>::leafbegin<1>() const;

template Dune::OneDInNDGrid<2>::Codim<0>::LeafIterator Dune::OneDInNDGrid<2>::leafend<0>() const;
template Dune::OneDInNDGrid<2>::Codim<1>::LeafIterator Dune::OneDInNDGrid<2>::leafend<1>() const;

template Dune::OneDInNDGrid<2>::Codim<0>::Partition<Dune::Interior_Partition>::LeafIterator 
Dune::OneDInNDGrid<2>::leafbegin<0,Dune::Interior_Partition>() const;
template Dune::OneDInNDGrid<2>::Codim<0>::Partition<Dune::InteriorBorder_Partition>::LeafIterator 
Dune::OneDInNDGrid<2>::leafbegin<0,Dune::InteriorBorder_Partition>() const;
template Dune::OneDInNDGrid<2>::Codim<0>::Partition<Dune::Overlap_Partition>::LeafIterator
 Dune::OneDInNDGrid<2>::leafbegin<0,Dune::Overlap_Partition>() const;
template Dune::OneDInNDGrid<2>::Codim<0>::Partition<Dune::OverlapFront_Partition>::LeafIterator
Dune::OneDInNDGrid<2>::leafbegin<0,Dune::OverlapFront_Partition>() const;
template Dune::OneDInNDGrid<2>::Codim<0>::Partition<Dune::All_Partition>::LeafIterator
Dune::OneDInNDGrid<2>::leafbegin<0,Dune::All_Partition>() const;
template Dune::OneDInNDGrid<2>::Codim<0>::Partition<Dune::Ghost_Partition>::LeafIterator
Dune::OneDInNDGrid<2>::leafbegin<0,Dune::Ghost_Partition>() const;

template Dune::OneDInNDGrid<2>::Codim<0>::Partition<Dune::Interior_Partition>::LeafIterator 
Dune::OneDInNDGrid<2>::leafend<0,Dune::Interior_Partition>() const;
template Dune::OneDInNDGrid<2>::Codim<0>::Partition<Dune::InteriorBorder_Partition>::LeafIterator 
Dune::OneDInNDGrid<2>::leafend<0,Dune::InteriorBorder_Partition>() const;
template Dune::OneDInNDGrid<2>::Codim<0>::Partition<Dune::Overlap_Partition>::LeafIterator
 Dune::OneDInNDGrid<2>::leafend<0,Dune::Overlap_Partition>() const;
template Dune::OneDInNDGrid<2>::Codim<0>::Partition<Dune::OverlapFront_Partition>::LeafIterator
Dune::OneDInNDGrid<2>::leafend<0,Dune::OverlapFront_Partition>() const;
template Dune::OneDInNDGrid<2>::Codim<0>::Partition<Dune::All_Partition>::LeafIterator
Dune::OneDInNDGrid<2>::leafend<0,Dune::All_Partition>() const;
template Dune::OneDInNDGrid<2>::Codim<0>::Partition<Dune::Ghost_Partition>::LeafIterator
Dune::OneDInNDGrid<2>::leafend<0,Dune::Ghost_Partition>() const;

template Dune::OneDInNDGrid<2>::Codim<1>::Partition<Dune::Interior_Partition>::LeafIterator 
Dune::OneDInNDGrid<2>::leafbegin<1,Dune::Interior_Partition>() const;
template Dune::OneDInNDGrid<2>::Codim<1>::Partition<Dune::InteriorBorder_Partition>::LeafIterator 
Dune::OneDInNDGrid<2>::leafbegin<1,Dune::InteriorBorder_Partition>() const;
template Dune::OneDInNDGrid<2>::Codim<1>::Partition<Dune::Overlap_Partition>::LeafIterator
 Dune::OneDInNDGrid<2>::leafbegin<1,Dune::Overlap_Partition>() const;
template Dune::OneDInNDGrid<2>::Codim<1>::Partition<Dune::OverlapFront_Partition>::LeafIterator
Dune::OneDInNDGrid<2>::leafbegin<1,Dune::OverlapFront_Partition>() const;
template Dune::OneDInNDGrid<2>::Codim<1>::Partition<Dune::All_Partition>::LeafIterator
Dune::OneDInNDGrid<2>::leafbegin<1,Dune::All_Partition>() const;
template Dune::OneDInNDGrid<2>::Codim<1>::Partition<Dune::Ghost_Partition>::LeafIterator
Dune::OneDInNDGrid<2>::leafbegin<1,Dune::Ghost_Partition>() const;

template Dune::OneDInNDGrid<2>::Codim<1>::Partition<Dune::Interior_Partition>::LeafIterator 
Dune::OneDInNDGrid<2>::leafend<1,Dune::Interior_Partition>() const;
template Dune::OneDInNDGrid<2>::Codim<1>::Partition<Dune::InteriorBorder_Partition>::LeafIterator 
Dune::OneDInNDGrid<2>::leafend<1,Dune::InteriorBorder_Partition>() const;
template Dune::OneDInNDGrid<2>::Codim<1>::Partition<Dune::Overlap_Partition>::LeafIterator
 Dune::OneDInNDGrid<2>::leafend<1,Dune::Overlap_Partition>() const;
template Dune::OneDInNDGrid<2>::Codim<1>::Partition<Dune::OverlapFront_Partition>::LeafIterator
Dune::OneDInNDGrid<2>::leafend<1,Dune::OverlapFront_Partition>() const;
template Dune::OneDInNDGrid<2>::Codim<1>::Partition<Dune::All_Partition>::LeafIterator
Dune::OneDInNDGrid<2>::leafend<1,Dune::All_Partition>() const;
template Dune::OneDInNDGrid<2>::Codim<1>::Partition<Dune::Ghost_Partition>::LeafIterator
Dune::OneDInNDGrid<2>::leafend<1,Dune::Ghost_Partition>() const;













template Dune::OneDInNDGrid<3>::Codim<0>::LevelIterator Dune::OneDInNDGrid<3>::lbegin<0>(int level) const;
template Dune::OneDInNDGrid<3>::Codim<1>::LevelIterator Dune::OneDInNDGrid<3>::lbegin<1>(int level) const;

template Dune::OneDInNDGrid<3>::Codim<0>::LevelIterator Dune::OneDInNDGrid<3>::lend<0>(int level) const;
template Dune::OneDInNDGrid<3>::Codim<1>::LevelIterator Dune::OneDInNDGrid<3>::lend<1>(int level) const;

template Dune::OneDInNDGrid<3>::Codim<0>::Partition<Dune::Interior_Partition>::LevelIterator 
Dune::OneDInNDGrid<3>::lbegin<0,Dune::Interior_Partition>(int level) const;
template Dune::OneDInNDGrid<3>::Codim<0>::Partition<Dune::InteriorBorder_Partition>::LevelIterator 
Dune::OneDInNDGrid<3>::lbegin<0,Dune::InteriorBorder_Partition>(int level) const;
template Dune::OneDInNDGrid<3>::Codim<0>::Partition<Dune::Overlap_Partition>::LevelIterator
 Dune::OneDInNDGrid<3>::lbegin<0,Dune::Overlap_Partition>(int level) const;
template Dune::OneDInNDGrid<3>::Codim<0>::Partition<Dune::OverlapFront_Partition>::LevelIterator
Dune::OneDInNDGrid<3>::lbegin<0,Dune::OverlapFront_Partition>(int level) const;
template Dune::OneDInNDGrid<3>::Codim<0>::Partition<Dune::All_Partition>::LevelIterator
Dune::OneDInNDGrid<3>::lbegin<0,Dune::All_Partition>(int level) const;
template Dune::OneDInNDGrid<3>::Codim<0>::Partition<Dune::Ghost_Partition>::LevelIterator
Dune::OneDInNDGrid<3>::lbegin<0,Dune::Ghost_Partition>(int level) const;

template Dune::OneDInNDGrid<3>::Codim<0>::Partition<Dune::Interior_Partition>::LevelIterator 
Dune::OneDInNDGrid<3>::lend<0,Dune::Interior_Partition>(int level) const;
template Dune::OneDInNDGrid<3>::Codim<0>::Partition<Dune::InteriorBorder_Partition>::LevelIterator 
Dune::OneDInNDGrid<3>::lend<0,Dune::InteriorBorder_Partition>(int level) const;
template Dune::OneDInNDGrid<3>::Codim<0>::Partition<Dune::Overlap_Partition>::LevelIterator
 Dune::OneDInNDGrid<3>::lend<0,Dune::Overlap_Partition>(int level) const;
template Dune::OneDInNDGrid<3>::Codim<0>::Partition<Dune::OverlapFront_Partition>::LevelIterator
Dune::OneDInNDGrid<3>::lend<0,Dune::OverlapFront_Partition>(int level) const;
template Dune::OneDInNDGrid<3>::Codim<0>::Partition<Dune::All_Partition>::LevelIterator
Dune::OneDInNDGrid<3>::lend<0,Dune::All_Partition>(int level) const;
template Dune::OneDInNDGrid<3>::Codim<0>::Partition<Dune::Ghost_Partition>::LevelIterator
Dune::OneDInNDGrid<3>::lend<0,Dune::Ghost_Partition>(int level) const;

template Dune::OneDInNDGrid<3>::Codim<1>::Partition<Dune::Interior_Partition>::LevelIterator 
Dune::OneDInNDGrid<3>::lbegin<1,Dune::Interior_Partition>(int level) const;
template Dune::OneDInNDGrid<3>::Codim<1>::Partition<Dune::InteriorBorder_Partition>::LevelIterator 
Dune::OneDInNDGrid<3>::lbegin<1,Dune::InteriorBorder_Partition>(int level) const;
template Dune::OneDInNDGrid<3>::Codim<1>::Partition<Dune::Overlap_Partition>::LevelIterator
 Dune::OneDInNDGrid<3>::lbegin<1,Dune::Overlap_Partition>(int level) const;
template Dune::OneDInNDGrid<3>::Codim<1>::Partition<Dune::OverlapFront_Partition>::LevelIterator
Dune::OneDInNDGrid<3>::lbegin<1,Dune::OverlapFront_Partition>(int level) const;
template Dune::OneDInNDGrid<3>::Codim<1>::Partition<Dune::All_Partition>::LevelIterator
Dune::OneDInNDGrid<3>::lbegin<1,Dune::All_Partition>(int level) const;
template Dune::OneDInNDGrid<3>::Codim<1>::Partition<Dune::Ghost_Partition>::LevelIterator
Dune::OneDInNDGrid<3>::lbegin<1,Dune::Ghost_Partition>(int level) const;

template Dune::OneDInNDGrid<3>::Codim<1>::Partition<Dune::Interior_Partition>::LevelIterator 
Dune::OneDInNDGrid<3>::lend<1,Dune::Interior_Partition>(int level) const;
template Dune::OneDInNDGrid<3>::Codim<1>::Partition<Dune::InteriorBorder_Partition>::LevelIterator 
Dune::OneDInNDGrid<3>::lend<1,Dune::InteriorBorder_Partition>(int level) const;
template Dune::OneDInNDGrid<3>::Codim<1>::Partition<Dune::Overlap_Partition>::LevelIterator
 Dune::OneDInNDGrid<3>::lend<1,Dune::Overlap_Partition>(int level) const;
template Dune::OneDInNDGrid<3>::Codim<1>::Partition<Dune::OverlapFront_Partition>::LevelIterator
Dune::OneDInNDGrid<3>::lend<1,Dune::OverlapFront_Partition>(int level) const;
template Dune::OneDInNDGrid<3>::Codim<1>::Partition<Dune::All_Partition>::LevelIterator
Dune::OneDInNDGrid<3>::lend<1,Dune::All_Partition>(int level) const;
template Dune::OneDInNDGrid<3>::Codim<1>::Partition<Dune::Ghost_Partition>::LevelIterator
Dune::OneDInNDGrid<3>::lend<1,Dune::Ghost_Partition>(int level) const;

template Dune::OneDInNDGrid<3>::Codim<0>::LeafIterator Dune::OneDInNDGrid<3>::leafbegin<0>() const;
template Dune::OneDInNDGrid<3>::Codim<1>::LeafIterator Dune::OneDInNDGrid<3>::leafbegin<1>() const;

template Dune::OneDInNDGrid<3>::Codim<0>::LeafIterator Dune::OneDInNDGrid<3>::leafend<0>() const;
template Dune::OneDInNDGrid<3>::Codim<1>::LeafIterator Dune::OneDInNDGrid<3>::leafend<1>() const;

template Dune::OneDInNDGrid<3>::Codim<0>::Partition<Dune::Interior_Partition>::LeafIterator 
Dune::OneDInNDGrid<3>::leafbegin<0,Dune::Interior_Partition>() const;
template Dune::OneDInNDGrid<3>::Codim<0>::Partition<Dune::InteriorBorder_Partition>::LeafIterator 
Dune::OneDInNDGrid<3>::leafbegin<0,Dune::InteriorBorder_Partition>() const;
template Dune::OneDInNDGrid<3>::Codim<0>::Partition<Dune::Overlap_Partition>::LeafIterator
 Dune::OneDInNDGrid<3>::leafbegin<0,Dune::Overlap_Partition>() const;
template Dune::OneDInNDGrid<3>::Codim<0>::Partition<Dune::OverlapFront_Partition>::LeafIterator
Dune::OneDInNDGrid<3>::leafbegin<0,Dune::OverlapFront_Partition>() const;
template Dune::OneDInNDGrid<3>::Codim<0>::Partition<Dune::All_Partition>::LeafIterator
Dune::OneDInNDGrid<3>::leafbegin<0,Dune::All_Partition>() const;
template Dune::OneDInNDGrid<3>::Codim<0>::Partition<Dune::Ghost_Partition>::LeafIterator
Dune::OneDInNDGrid<3>::leafbegin<0,Dune::Ghost_Partition>() const;

template Dune::OneDInNDGrid<3>::Codim<0>::Partition<Dune::Interior_Partition>::LeafIterator 
Dune::OneDInNDGrid<3>::leafend<0,Dune::Interior_Partition>() const;
template Dune::OneDInNDGrid<3>::Codim<0>::Partition<Dune::InteriorBorder_Partition>::LeafIterator 
Dune::OneDInNDGrid<3>::leafend<0,Dune::InteriorBorder_Partition>() const;
template Dune::OneDInNDGrid<3>::Codim<0>::Partition<Dune::Overlap_Partition>::LeafIterator
 Dune::OneDInNDGrid<3>::leafend<0,Dune::Overlap_Partition>() const;
template Dune::OneDInNDGrid<3>::Codim<0>::Partition<Dune::OverlapFront_Partition>::LeafIterator
Dune::OneDInNDGrid<3>::leafend<0,Dune::OverlapFront_Partition>() const;
template Dune::OneDInNDGrid<3>::Codim<0>::Partition<Dune::All_Partition>::LeafIterator
Dune::OneDInNDGrid<3>::leafend<0,Dune::All_Partition>() const;
template Dune::OneDInNDGrid<3>::Codim<0>::Partition<Dune::Ghost_Partition>::LeafIterator
Dune::OneDInNDGrid<3>::leafend<0,Dune::Ghost_Partition>() const;

template Dune::OneDInNDGrid<3>::Codim<1>::Partition<Dune::Interior_Partition>::LeafIterator 
Dune::OneDInNDGrid<3>::leafbegin<1,Dune::Interior_Partition>() const;
template Dune::OneDInNDGrid<3>::Codim<1>::Partition<Dune::InteriorBorder_Partition>::LeafIterator 
Dune::OneDInNDGrid<3>::leafbegin<1,Dune::InteriorBorder_Partition>() const;
template Dune::OneDInNDGrid<3>::Codim<1>::Partition<Dune::Overlap_Partition>::LeafIterator
 Dune::OneDInNDGrid<3>::leafbegin<1,Dune::Overlap_Partition>() const;
template Dune::OneDInNDGrid<3>::Codim<1>::Partition<Dune::OverlapFront_Partition>::LeafIterator
Dune::OneDInNDGrid<3>::leafbegin<1,Dune::OverlapFront_Partition>() const;
template Dune::OneDInNDGrid<3>::Codim<1>::Partition<Dune::All_Partition>::LeafIterator
Dune::OneDInNDGrid<3>::leafbegin<1,Dune::All_Partition>() const;
template Dune::OneDInNDGrid<3>::Codim<1>::Partition<Dune::Ghost_Partition>::LeafIterator
Dune::OneDInNDGrid<3>::leafbegin<1,Dune::Ghost_Partition>() const;

template Dune::OneDInNDGrid<3>::Codim<1>::Partition<Dune::Interior_Partition>::LeafIterator 
Dune::OneDInNDGrid<3>::leafend<1,Dune::Interior_Partition>() const;
template Dune::OneDInNDGrid<3>::Codim<1>::Partition<Dune::InteriorBorder_Partition>::LeafIterator 
Dune::OneDInNDGrid<3>::leafend<1,Dune::InteriorBorder_Partition>() const;
template Dune::OneDInNDGrid<3>::Codim<1>::Partition<Dune::Overlap_Partition>::LeafIterator
 Dune::OneDInNDGrid<3>::leafend<1,Dune::Overlap_Partition>() const;
template Dune::OneDInNDGrid<3>::Codim<1>::Partition<Dune::OverlapFront_Partition>::LeafIterator
Dune::OneDInNDGrid<3>::leafend<1,Dune::OverlapFront_Partition>() const;
template Dune::OneDInNDGrid<3>::Codim<1>::Partition<Dune::All_Partition>::LeafIterator
Dune::OneDInNDGrid<3>::leafend<1,Dune::All_Partition>() const;
template Dune::OneDInNDGrid<3>::Codim<1>::Partition<Dune::Ghost_Partition>::LeafIterator
Dune::OneDInNDGrid<3>::leafend<1,Dune::Ghost_Partition>() const;
