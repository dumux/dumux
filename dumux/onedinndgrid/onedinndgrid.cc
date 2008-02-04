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
    static Dune::OneDInNDGridLevelIterator<1,PiType, const Dune::OneDInNDGrid<dimworld> > 
    lbegin(const Dune::OneDInNDGrid<dimworld>* g, int level) {

    return Dune::OneDInNDGridLevelIterator<1,PiType, const Dune::OneDInNDGrid<dimworld> >(g->vertices[level].begin);
    }

};

template <>
struct OneDInNDGridLevelIteratorFactory<0>
{
    template <Dune::PartitionIteratorType PiType, int dimworld>
    static Dune::OneDInNDGridLevelIterator<0,PiType, const Dune::OneDInNDGrid<dimworld> > 
    lbegin(const Dune::OneDInNDGrid<dimworld> * g, int level) {

    return Dune::OneDInNDGridLevelIterator<0,PiType, const Dune::OneDInNDGrid<dimworld> >(g->elements[level].begin);
    }

};

}


//Dune::OneDInNDGrid::OneDInNDGrid(int numElements, double leftBoundary, double rightBoundary)
//    : refinementType_(LOCAL),
//      leafIndexSet_(*this),
//      idSet_(*this),
//      freeVertexIdCounter_(0),
//      freeElementIdCounter_(0)
//{
//    if (numElements<1)
//        DUNE_THROW(GridError, "Nonpositive number of elements requested!");
//
//    if (leftBoundary >= rightBoundary)
//        DUNE_THROW(GridError, "The left boundary coordinate has to be strictly less than the right boundary one!");
//    // Init grid hierarchy
//    vertices.resize(1);
//    elements.resize(1);
//
//    // Init vertex set
//    for (int i=0; i<numElements+1; i++) {
//        double newCoord = leftBoundary + i*(rightBoundary-leftBoundary) / numElements;
//
//        OneDInNDEntityImp<0, dimworld>* newVertex = new OneDInNDEntityImp<0, dimworld>(0, newCoord);
//        newVertex->id_ = getNextFreeId(1);
//        vertices[0].insert_after(vertices[0].rbegin, newVertex);
//
//    }
//
//    // Init element set
//    OneDInNDEntityImp<0, dimworld>* it = vertices[0].begin;
//    for (int i=0; i<numElements; i++) {
//
//        OneDInNDEntityImp<1, dimworld>* newElement = new OneDInNDEntityImp<1, dimworld>(0, getNextFreeId(0));
//        newElement->vertex_[0] = it;
//        it = it->succ_;
//        newElement->vertex_[1] = it;
//
//        elements[0].insert_after(elements[0].rbegin, newElement);
//
//    }
//
//    setIndices();
//}

template<int dimworld>
Dune::OneDInNDGrid<dimworld>::OneDInNDGrid(int numElements, Dune::FieldVector<double, dimworld> leftBoundary, Dune::FieldVector<double, dimworld> rightBoundary)
    : refinementType_(LOCAL),
      leafIndexSet_(*this),
      idSet_(*this),
      freeVertexIdCounter_(0),
      freeElementIdCounter_(0)
{
    if (numElements<1)
        DUNE_THROW(GridError, "Nonpositive number of elements requested!");

    // Init grid hierarchy
    vertices.resize(1);
    elements.resize(1);

    // Init vertex set
    for (int i=0; i<numElements+1; i++) {
    	Dune::FieldVector<double, dimworld> newCoord;
    	for (int k = 0; k < dimworld; k++)
    		newCoord[k] = leftBoundary[k] + i*(rightBoundary[k]-leftBoundary[k]) / numElements;

        OneDInNDEntityImp<0, dimworld>* newVertex = new OneDInNDEntityImp<0, dimworld>(0, newCoord);
        newVertex->id_ = getNextFreeId(1);
        vertices[0].insert_after(vertices[0].rbegin, newVertex);

    }

    // Init element set
    OneDInNDEntityImp<0, dimworld>* it = vertices[0].begin;
    for (int i=0; i<numElements; i++) {

        OneDInNDEntityImp<1, dimworld>* newElement = new OneDInNDEntityImp<1, dimworld>(0, getNextFreeId(0));
        newElement->vertex_[0] = it;
        it = it->succ_;
        newElement->vertex_[1] = it;

        elements[0].insert_after(elements[0].rbegin, newElement);

    }

    setIndices();
}

//Dune::OneDInNDGrid::OneDInNDGrid(const std::vector<OneDInNDCType>& coords) 
//    : refinementType_(LOCAL),
//      leafIndexSet_(*this),
//      idSet_(*this),
//      freeVertexIdCounter_(0),
//      freeElementIdCounter_(0)
//{
//    if (coords.size()<2)
//        DUNE_THROW(GridError, "You have to provide at least two coordinates!");
//
//    // Init grid hierarchy
//    vertices.resize(1);
//    elements.resize(1);
//
//    // Init vertex set
//    for (size_t i=0; i<coords.size(); i++) {
//        OneDInNDEntityImp<0, dimworld>* newVertex = new OneDInNDEntityImp<0, dimworld>(0, coords[i], getNextFreeId(1));
//        vertices[0].insert_after(vertices[0].rbegin, newVertex);
//    }
//
//    // Init element set
//    OneDInNDEntityImp<0, dimworld>* it = vertices[0].begin;
//    for (size_t i=0; i<coords.size()-1; i++) {
//
//        OneDInNDEntityImp<1, dimworld>* newElement = new OneDInNDEntityImp<1, dimworld>(0, getNextFreeId(0));
//        newElement->vertex_[0] = it;
//        it = it->succ_;
//        newElement->vertex_[1] = it;
//
//        if (newElement->vertex_[0]->pos_ >= newElement->vertex_[1]->pos_)
//            DUNE_THROW(GridError, "The coordinates have to be in ascending order!");
//
//        elements[0].insert_after(elements[0].rbegin, newElement);
//
//    }
//
//    setIndices();
//}

template<int dimworld>
Dune::OneDInNDGrid<dimworld>::~OneDInNDGrid()
{
    // Delete all vertices
    for (unsigned int i=0; i<vertices.size(); i++) {

        OneDInNDEntityImp<0, dimworld>* v = vertices[i].begin;

        while (v) {
            
            OneDInNDEntityImp<0, dimworld>* vSucc = v->succ_;
            vertices[i].remove(v);
            delete(v);
            v = vSucc;

        }

    }

    // Delete all elements
    for (unsigned int i=0; i<elements.size(); i++) {

        OneDInNDEntityImp<1, dimworld>* e = elements[i].begin;

        while (e) {
            
            OneDInNDEntityImp<1, dimworld>* eSucc = e->succ_;
            elements[i].remove(e);
            delete(e);
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
    
    return OneDInNDGridLevelIterator<codim,PiType, const Dune::OneDInNDGrid<dimworld> >(0);
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
Dune::OneDInNDEntityImp<0, dimworld>*
Dune::OneDInNDGrid<dimworld>::getLeftUpperVertex(const OneDInNDEntityImp<1, dimworld>* eIt)
{
    OneDInNDEntityImp<1, dimworld>* l = eIt->pred_;

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
Dune::OneDInNDEntityImp<0, dimworld>*
Dune::OneDInNDGrid<dimworld>::getRightUpperVertex(const OneDInNDEntityImp<1, dimworld>* eIt)
{
    OneDInNDEntityImp<1, dimworld>* r = eIt->succ_;

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
Dune::OneDInNDEntityImp<1, dimworld>* 
Dune::OneDInNDGrid<dimworld>::getLeftNeighborWithSon(OneDInNDEntityImp<1, dimworld>* eIt)
{
    OneDInNDEntityImp<1, dimworld>* l = eIt;

    do {
        l = l->pred_;
    } while (l && l->isLeaf());

    return l;
}


template<int dimworld>
bool Dune::OneDInNDGrid<dimworld>::adapt()
{
    OneDInNDEntityImp<1, dimworld>* eIt;
    
    // for the return value:  true if the grid was changed
    bool changedGrid = false;

    // remove all elements that have been marked for coarsening
    for (int i=1; i<=maxLevel(); i++) {

        for (eIt = elements[i].begin; eIt!=NULL; ) {

            OneDInNDEntityImp<1, dimworld>* leftElementToBeDeleted  = eIt;
            OneDInNDEntityImp<1, dimworld>* rightElementToBeDeleted = eIt->succ_;

            assert(eIt->succ_);
            OneDInNDEntityImp<1, dimworld>* nextElement = eIt->succ_->succ_;

            if (leftElementToBeDeleted->markState_ == OneDInNDEntityImp<1, dimworld>::COARSEN && leftElementToBeDeleted->isLeaf()
                && rightElementToBeDeleted->markState_ == OneDInNDEntityImp<1, dimworld>::COARSEN && rightElementToBeDeleted->isLeaf()) {

                assert(rightElementToBeDeleted->isLeaf());

                // Is the left vertex obsolete?
                if (leftElementToBeDeleted->pred_==NULL 
                    || leftElementToBeDeleted->pred_->vertex_[1] != leftElementToBeDeleted->vertex_[0]){

                    // If the left vertex has a father remove the reference to this vertex at this father
                    assert(leftElementToBeDeleted->father_->vertex_[0]->son_ == leftElementToBeDeleted->vertex_[0]);
                    leftElementToBeDeleted->father_->vertex_[0]->son_ = NULL;

                    vertices[i].remove(leftElementToBeDeleted->vertex_[0]);
                    delete(leftElementToBeDeleted->vertex_[0]);
                }

                // Is the right vertex obsolete?
                if (rightElementToBeDeleted->succ_==NULL 
                    || rightElementToBeDeleted->succ_->vertex_[0] != rightElementToBeDeleted->vertex_[1]) {

                    // If the left vertex has a father remove the reference to this vertex at this father
                    assert(rightElementToBeDeleted->father_->vertex_[1]->son_ == rightElementToBeDeleted->vertex_[1]);
                    rightElementToBeDeleted->father_->vertex_[1]->son_ = NULL;

                    vertices[i].remove(rightElementToBeDeleted->vertex_[1]);
                    delete(rightElementToBeDeleted->vertex_[1]);
                }

                // Delete vertex between left and right element to be deleted
                assert(leftElementToBeDeleted->vertex_[1] == rightElementToBeDeleted->vertex_[0]);
                vertices[i].remove(leftElementToBeDeleted->vertex_[1]);
                delete(leftElementToBeDeleted->vertex_[1]);

                // Remove references from the father element
                assert(rightElementToBeDeleted->father_->sons_[1] == rightElementToBeDeleted);
                leftElementToBeDeleted->father_->sons_[0]  = NULL;
                rightElementToBeDeleted->father_->sons_[1] = NULL;

                // Paranoia: make sure the father is not marked for refinement
                rightElementToBeDeleted->father_->markState_ = OneDInNDEntityImp<1, dimworld>::NONE;

                // Actually delete elements
                elements[i].remove(leftElementToBeDeleted);
                elements[i].remove(rightElementToBeDeleted);
                delete(leftElementToBeDeleted);
                delete(rightElementToBeDeleted);

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
    for (eIt = elements[maxLevel()].begin; eIt!=NULL; eIt=eIt->succ_) 
        if (eIt->markState_ == OneDInNDEntityImp<1, dimworld>::REFINED) {
            toplevelRefinement = true;
            break;
        }

    if (toplevelRefinement) {
        List<OneDInNDEntityImp<0, dimworld> > newVertices;
        List<OneDInNDEntityImp<1, dimworld> > newElements;
        vertices.push_back(newVertices);
        elements.push_back(newElements);
    }

    // //////////////////////////////
    // refine all marked elements
    // //////////////////////////////
    int oldMaxlevel = (toplevelRefinement) ? maxLevel()-1 : maxLevel();
    for (int i=0; i<=oldMaxlevel; i++) {

        for (eIt = elements[i].begin; eIt!=NULL; eIt = eIt->succ_) {

            if (eIt->markState_ == OneDInNDEntityImp<1, dimworld>::REFINED
                && eIt->isLeaf()) {

                // Does the left vertex exist on the next-higher level?
                // If no create it
                OneDInNDEntityImp<0, dimworld>* leftUpperVertex = getLeftUpperVertex(eIt);

                if (leftUpperVertex==NULL)
                    leftUpperVertex = new OneDInNDEntityImp<0, dimworld>(i+1, 
                                                           eIt->vertex_[0]->pos_, 
                                                           eIt->vertex_[0]->id_);

                eIt->vertex_[0]->son_ = leftUpperVertex;

                // Does the right vertex exist on the next-higher level?
                // If no create it
                OneDInNDEntityImp<0, dimworld>* rightUpperVertex = getRightUpperVertex(eIt);
                
                if (rightUpperVertex==NULL) 
                    rightUpperVertex = new OneDInNDEntityImp<0, dimworld>(i+1, 
                                                            eIt->vertex_[1]->pos_, 
                                                            eIt->vertex_[1]->id_);

                eIt->vertex_[1]->son_ = rightUpperVertex;

                // Create center vertex
                FieldVector<double, dimworld> newPos = eIt->vertex_[0]->pos_;
                newPos += eIt->vertex_[1]->pos_;
                newPos *= 0.5;
                //double p = 0.5*(eIt->vertex_[0]->pos_ + eIt->vertex_[1]->pos_);

                OneDInNDEntityImp<0, dimworld>* centerVertex = new OneDInNDEntityImp<0, dimworld>(i+1, newPos, getNextFreeId(1));

                // //////////////////////////////////////
                // Insert new vertices into vertex list
                // //////////////////////////////////////

                OneDInNDEntityImp<1, dimworld>* leftNeighbor = getLeftNeighborWithSon(eIt);

                if (leftNeighbor!=NULL) {

                    // leftNeighbor exists
                    if ( leftNeighbor->sons_[1]->vertex_[1] != leftUpperVertex)
                        vertices[i+1].insert_after(leftNeighbor->sons_[1]->vertex_[1], leftUpperVertex);

                } else {
                    // leftNeighbor does not exist
                    vertices[i+1].insert_before(vertices[i+1].begin, leftUpperVertex);

                }

                vertices[i+1].insert_after(leftUpperVertex, centerVertex);

                // Check if rightUpperVertex is already in the list
                OneDInNDEntityImp<0, dimworld>* succOfCenter = centerVertex->succ_;

                if (succOfCenter==NULL || succOfCenter != rightUpperVertex)
                    vertices[i+1].insert_after(centerVertex, rightUpperVertex);

                // ///////////////////////
                // Create new elements
                // ///////////////////////
                OneDInNDEntityImp<1, dimworld>* newElement0 = new OneDInNDEntityImp<1, dimworld>(i+1, getNextFreeId(0));
                newElement0->vertex_[0] = leftUpperVertex;
                newElement0->vertex_[1] = centerVertex;
                newElement0->father_ = eIt;
                newElement0->adaptationState_ = OneDInNDEntityImp<1, dimworld>::REFINED;

                OneDInNDEntityImp<1, dimworld>* newElement1 = new OneDInNDEntityImp<1, dimworld>(i+1, getNextFreeId(0));
                newElement1->vertex_[0] = centerVertex;
                newElement1->vertex_[1] = rightUpperVertex;
                newElement1->father_ = eIt;
                newElement1->adaptationState_ = OneDInNDEntityImp<1, dimworld>::REFINED;

                // Insert new elements into element list
                if (leftNeighbor!=NULL)
                    // leftNeighbor exists
                    elements[i+1].insert_after(leftNeighbor->sons_[1], newElement0);
                else 
                    // leftNeighbor does not exist
                    elements[i+1].insert_before(elements[i+1].begin, newElement0);

                elements[i+1].insert_after(newElement0, newElement1);

                // Mark the two new elements as the sons of the refined element
                eIt->sons_[0] = newElement0;
                eIt->sons_[1] = newElement1;

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
            
            OneDInNDEntityImp<1, dimworld>* eIt;
            for (eIt = elements[i].begin; eIt!=NULL; eIt = eIt->succ_) {
                
                if (eIt->isLeaf()) {
                    
                    // Does the left vertex exist on the next-higher level?
                    // If no create it
                    OneDInNDEntityImp<0, dimworld>* leftUpperVertex = getLeftUpperVertex(eIt);
                    
                    if (leftUpperVertex==NULL)
                        leftUpperVertex = new OneDInNDEntityImp<0, dimworld>(i+1, eIt->vertex_[0]->pos_, eIt->vertex_[0]->id_);
 
                    eIt->vertex_[0]->son_ = leftUpperVertex;
                   
                    // Does the right vertex exist on the next-higher level?
                    // If no create it
                    OneDInNDEntityImp<0, dimworld>* rightUpperVertex = getRightUpperVertex(eIt);
                    
                    if (rightUpperVertex==NULL) 
                        rightUpperVertex = new OneDInNDEntityImp<0, dimworld>(i+1, eIt->vertex_[1]->pos_, eIt->vertex_[1]->id_);

                    eIt->vertex_[1]->son_ = rightUpperVertex;

                    // //////////////////////////////////////
                    // Insert new vertices into vertex list
                    // //////////////////////////////////////
                    
                    OneDInNDEntityImp<1, dimworld>* leftNeighbor = getLeftNeighborWithSon(eIt);
                    
                    if (leftNeighbor!=NULL) {
                        
                        // leftNeighbor exists
                        if ( leftNeighbor->sons_[1]->vertex_[1] != leftUpperVertex)
                            vertices[i+1].insert_after(leftNeighbor->sons_[1]->vertex_[1], leftUpperVertex);
                        
                    } else {
                        // leftNeighbor does not exist
                        vertices[i+1].insert_before(vertices[i+1].begin, leftUpperVertex);

                    }

                    // Check if rightUpperVertex is already in the list
                    OneDInNDEntityImp<0, dimworld>* succOfLeft = leftUpperVertex->succ_;

                    if (succOfLeft==NULL || succOfLeft != rightUpperVertex)
                        vertices[i+1].insert_after(leftUpperVertex, rightUpperVertex);
                    
                    // /////////////////////////
                    //   Create new element
                    // /////////////////////////
                    OneDInNDEntityImp<1, dimworld>* newElement = new OneDInNDEntityImp<1, dimworld>(i+1, eIt->id_);
                    newElement->vertex_[0] = leftUpperVertex;
                    newElement->vertex_[1] = rightUpperVertex;
                    newElement->father_ = eIt;
                    newElement->adaptationState_ = OneDInNDEntityImp<1, dimworld>::REFINED;

                    // Insert new elements into element list
                    if (leftNeighbor!=NULL)
                        // leftNeighbor exists
                        elements[i+1].insert_after(leftNeighbor->sons_[1], newElement);
                    else 
                        // leftNeighbor does not exist
                        elements[i+1].insert_before(elements[i+1].begin, newElement);
                    
                    // Mark the new element as the sons of the refined element
                    eIt->sons_[0] = eIt->sons_[1] = newElement;
                    
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
        if (getRealImplementation(*eIt).target_->markState_ != OneDInNDEntityImp<1, dimworld>::NONE)
            return true;

    return false;
}

template<int dimworld>
void Dune::OneDInNDGrid<dimworld>::postAdapt()
{
    for (int i=0; i<=maxLevel(); i++) {
        OneDInNDEntityImp<1, dimworld>* eIt;
        for (eIt = elements[i].begin; eIt!=NULL; eIt = eIt->succ_)
            eIt->markState_ = OneDInNDEntityImp<1, dimworld>::NONE;

    }
        
}

template<int dimworld>
void Dune::OneDInNDGrid<dimworld>::setIndices()
{
    // Add space for new LevelIndexSets if the grid hierarchy got higher
    // They are not created until they are actually requested
    for (int i=levelIndexSets_.size(); i<maxLevel()+1; i++)
        levelIndexSets_.push_back(0);

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

    idSet_.update();

}

template<int dimworld>
void Dune::OneDInNDGrid<dimworld>::globalRefine(int refCount)
{
  for (int i=0; i<refCount; i++) {

      // mark all entities for grid refinement
	  typename Traits::template Codim<0>::LeafIterator iIt    = leafbegin<0>();
	  typename Traits::template Codim<0>::LeafIterator iEndIt = leafend<0>();
      
      for (; iIt!=iEndIt; ++iIt)
          mark(1, iIt);
	  
      this->preAdapt();
      adapt();
      this->postAdapt();
  }
}

template<int dimworld>
bool Dune::OneDInNDGrid<dimworld>::mark(int refCount,
                          const typename Traits::template Codim<0>::EntityPointer & e )
{
    if (refCount < 0) {

        if (getRealImplementation(*e).target_->level_ == 0)
            return false;
        else {
            getRealImplementation(*e).target_->markState_ = OneDInNDEntityImp<1, dimworld>::COARSEN;
            return true;
        }

    } else if (refCount > 0)
        getRealImplementation(*e).target_->markState_ = OneDInNDEntityImp<1, dimworld>::REFINED;
    else
        getRealImplementation(*e).target_->markState_ = OneDInNDEntityImp<1, dimworld>::NONE;

    return true;
}

template<int dimworld>
int Dune::OneDInNDGrid<dimworld>::getMark(const typename Traits::template Codim<0>::EntityPointer & ep ) const 
{
  if(getRealImplementation(*ep).target_->markState_ == OneDInNDEntityImp<1, dimworld>::COARSEN) 
      return -1; 
  else if(getRealImplementation(*ep).target_->markState_ == OneDInNDEntityImp<1, dimworld>::REFINED) 
      return 1; 
  return 0;  
}



// /////////////////////////////////////////////////////////////////////////
//   Explicitly instantiate the OnedGrid member templates.
//   gcc-4.0 wants these instantiations after the method implementations
// /////////////////////////////////////////////////////////////////////////

template void Dune::OneDInNDGrid<2>::globalRefine(int refCount);
template Dune::OneDInNDGrid<2>::OneDInNDGrid(int numElements, Dune::FieldVector<double, 2> leftBoundary, Dune::FieldVector<double, 2> rightBoundary);
template Dune::OneDInNDGrid<2>::~OneDInNDGrid();

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









template void Dune::OneDInNDGrid<3>::globalRefine(int refCount);
template Dune::OneDInNDGrid<3>::OneDInNDGrid(int numElements, Dune::FieldVector<double, 3> leftBoundary, Dune::FieldVector<double, 3> rightBoundary);
template Dune::OneDInNDGrid<3>::~OneDInNDGrid();

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

