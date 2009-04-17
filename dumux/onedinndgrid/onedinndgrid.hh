// $Id: onedinndgrid.hh 1227 2009-02-13 12:01:13Z anneli $

#ifndef DUNE_ONE_D_IN_N_D_GRID_HH
#define DUNE_ONE_D_IN_N_D_GRID_HH

#include <vector>

#include <dune/common/misc.hh>
#include <dune/grid/common/capabilities.hh>
#include <dune/common/collectivecommunication.hh>
#include <dune/grid/common/grid.hh>


/** \file
 * \brief The OneDInNDGrid class
 */

namespace Dune
{
/** \todo Please doc me! */

// forward declarations
    template<int codim, int dim, class GridImp> class OneDInNDGridEntity;
    template<int codim, class GridImp> class OneDInNDGridEntityPointer;
    template<int codim, PartitionIteratorType pitype, class GridImp> class OneDInNDGridLevelIterator;

    template<int mydim, int coorddim, class GridImp>            class OneDInNDGridGeometry;
    template<class GridImp>            class OneDInNDGridHierarchicIterator;
    template<class GridImp, bool LeafIterator> class OneDInNDGridIntersectionIterator;

    template <int dimworld>
    class OneDInNDGrid;

    template<int codim>                        class OneDInNDGridLevelIteratorFactory;

}  // namespace Dune

#include "onedinndgridentity.hh"
#include "onedinndgridentitypointer.hh"
#include "onedinndgridgeometry.hh"
#include "onedinndintersectionit.hh"
#include "onedinndgridleveliterator.hh"
#include "onedinndgridleafiterator.hh"
#include "onedinndgridhieriterator.hh"
#include "onedinndgridindexsets.hh"

namespace Dune {

    // A simple double linked list
    template<class T>
    class List
    {

    public:

        List() : numelements(0), begin(0), rbegin(0) {}

        int size() const {return numelements;}

        T* insert_after (T* i, T* t) {

            // Teste Eingabe
            if (i==0 && begin!=0)
                DUNE_THROW(DoubleLinkedListError, "invalid iterator for insert_after");

            // einfuegen
            if (begin==0) {
    // einfuegen in leere Liste
    begin = t;
                rbegin = t;
            }
            else
                {
                    // nach Element i.p einsetzen
                    t->pred_ = i;
                    t->succ_ = i->succ_;
                    i->succ_ = t;

                    if (t->succ_!=0)
                        t->succ_->pred_ = t;

                    // tail neu ?
                    if (rbegin==i)
                        rbegin = t;
                }

            // Groesse und Rueckgabeiterator
            numelements = numelements+1;

            return t;
        }

        T* insert_before (T* i, T* t) {

            // Teste Eingabe
            if (i==0 && begin!=0)
                DUNE_THROW(DoubleLinkedListError,
                           "invalid iterator for insert_before");

            // einfuegen
            if (begin==0)
                {
                    // einfuegen in leere Liste
                    begin=t;
                    rbegin=t;
                }
            else
                {
                    // vor Element i.p einsetzen
                    t->succ_ = i;
                    t->pred_ = i->pred_;
                    i->pred_ = t;

                    if (t->pred_!=0)
                        t->pred_->succ_ = t;
                    // head neu ?
                    if (begin==i)
                        begin = t;
                }

            // Groesse und Rueckgabeiterator
            numelements = numelements+1;
            return t;
        }

        void remove (T* i)
        {
            // Teste Eingabe
            if (i==0)
                return;

            // Ausfaedeln
            if (i->succ_!=0)
                i->succ_->pred_ = i->pred_;
            if (i->pred_!=0)
                i->pred_->succ_ = i->succ_;

            // head & tail
            if (begin==i)
                begin=i->succ_;
            if (rbegin==i)
                rbegin = i->pred_;

            // Groesse
            numelements = numelements-1;
        }


        int numelements;

        T* begin;
        T* rbegin;

    };

template<int dim, int dimw>
struct OneDInNDGridFamily
{
    typedef GridTraits< dim,
                         dimw,
                         Dune::OneDInNDGrid<dimw>,
                         OneDInNDGridGeometry,
                         OneDInNDGridEntity,
                         OneDInNDGridEntityPointer,
                         OneDInNDGridLevelIterator,
                         OneDInNDGridLeafIntersectionIterator, // leaf  intersection iter
                         OneDInNDGridLevelIntersectionIterator, // level intersection iter
                         OneDInNDGridLeafIntersectionIterator, // leaf  intersection iter
                         OneDInNDGridLevelIntersectionIterator, // level intersection iter
                         OneDInNDGridHierarchicIterator,
                         OneDInNDGridLeafIterator,
                         OneDInNDGridLevelIndexSet<const OneDInNDGrid<dimw> >,
                         OneDInNDGridLevelIndexSetTypes<const OneDInNDGrid<dimw> >,
                         OneDInNDGridLeafIndexSet<const OneDInNDGrid<dimw> >,
                         OneDInNDGridLeafIndexSetTypes<const OneDInNDGrid<dimw> >,
                         OneDInNDGridIdSet<const OneDInNDGrid<dimw> >,
                         unsigned int,
                         OneDInNDGridIdSet<const OneDInNDGrid<dimw> >,
                         unsigned int,
                         CollectiveCommunication<Dune::OneDInNDGrid<dimw> >
                       >
  Traits;
};

//**********************************************************************
//
// --OneDInNDGrid
//
//**********************************************************************

/**
 \brief [<em> provides Grid </em>]
 Onedimensional adaptive grid
 \ingroup GridImplementations

 This implementation of the grid interface provides one-dimensional
 grids only.  No matter what the values of dim and dimworld may be,
 you'll always get a 1D-grid in a 1D-world.  Unlike SGrid, however,
 which can also be instantiated in 1D, the OneDInNDGrid is nonuniform
 and provides local mesh refinement and coarsening.
 */
template<int dimworld>
class OneDInNDGrid : public GridDefaultImplementation <1, dimworld, double,OneDInNDGridFamily<1, dimworld> >
{
    // Grid and world dimension are hardwired in this grid
    enum {dim = 1};

    typedef OneDInNDGrid<dimworld> ThisClass;

    friend class OneDInNDGridLevelIteratorFactory <0>;
    friend class OneDInNDGridLevelIteratorFactory <1>;
    friend class OneDInNDGridEntity <0,dim,ThisClass>;
    friend class OneDInNDGridEntity <dim,dim,ThisClass>;
    friend class OneDInNDGridHierarchicIterator<ThisClass>;
    friend class OneDInNDGridLeafIntersectionIterator<ThisClass>;
    friend class OneDInNDGridLevelIntersectionIterator<ThisClass>;

    friend class OneDInNDGridLevelIndexSet<const ThisClass>;
    friend class OneDInNDGridLeafIndexSet<const ThisClass>;
    friend class OneDInNDGridIdSet<const ThisClass>;

    template <int codim_, PartitionIteratorType PiType_, class GridImp_>
    friend class OneDInNDGridLeafIterator;

    template<int codim_, int dim_, class GridImp_, template<int,int,class> class EntityImp_>
    friend class Entity;

    /** \brief The type used by to store coordinates */
    typedef double OneDInNDCType;

    // **********************************************************
    // The Interface Methods
    // **********************************************************

public:
    /** \brief GridFamily of OneDInNDGrid */
    typedef OneDInNDGridFamily<dim,dimworld> GridFamily;

    /** \brief Provides the standard grid types */
    typedef typename GridFamily::Traits Traits;
    //typedef typename OneDInNDGridFamily<dim,dimworld>::Traits Traits;

    /** \brief Constructor with an explicit set of coordinates */
//    OneDInNDGrid(const std::vector<OneDInNDCType>& coords);

    /** \brief Constructor for a uniform grid */
//    OneDInNDGrid(int numElements, double leftBoundary, double rightBoundary);

    /** \brief Constructor for a uniform grid */
    OneDInNDGrid(int numElements, FieldVector<double, dimworld> leftBoundary, FieldVector<double, dimworld> rightBoundary);

    //! Destructor
    ~OneDInNDGrid();

    /** \brief Return maximum level defined in this grid.

    Levels are numbered 0 ... maxlevel with 0 the coarsest level.
    */
    int maxLevel() const {return vertices.size()-1;}

  //! Iterator to first entity of given codim on level
  template<int codim>
  typename Traits::template Codim<codim>::LevelIterator lbegin (int level) const;

  //! one past the end on this level
  template<int codim>
  typename Traits::template Codim<codim>::LevelIterator lend (int level) const;

    //! Iterator to first entity of given codim on level
    template<int codim, PartitionIteratorType PiType>
    typename Traits::template Codim<codim>::template Partition<PiType>::LevelIterator lbegin (int level) const;

    //! one past the end on this level
    template<int codim, PartitionIteratorType PiType>
    typename Traits::template Codim<codim>::template Partition<PiType>::LevelIterator lend (int level) const;

  //! Iterator to first entity of given codim on leaf level
  template<int codim>
  typename Traits::template Codim<codim>::LeafIterator leafbegin () const;

  //! one past the end on leaf level
  template<int codim>
  typename Traits::template Codim<codim>::LeafIterator leafend () const;

        //! Iterator to first entity of given codim on level
    template<int codim, PartitionIteratorType PiType>
    typename Traits::template Codim<codim>::template Partition<PiType>::LeafIterator leafbegin() const;

    //! one past the end on this level
    template<int codim, PartitionIteratorType PiType>
    typename Traits::template Codim<codim>::template Partition<PiType>::LeafIterator leafend() const;

    /** \brief Number of grid entities per level and codim
     */
    int size (int level, int codim) const {
        if (codim<0 || codim>1)
            DUNE_THROW(GridError, "There are no codim " << codim << " entities in a OneDInNDGrid!");

        if (codim==0)
            return elements[level].size();

        return vertices[level].size();
    }



  //! number of leaf entities per codim in this process
  int size (int codim) const
  {
      return leafIndexSet().size(codim);
  }

  //! number of entities per level and geometry type in this process
  int size (int level, GeometryType type) const
  {
      // There is only one type for each codim
      return size(level,1-type.dim());
  }

  //! number of leaf entities per geometry type in this process
  int size (GeometryType type) const
  {
      return leafIndexSet().size(type);
  }

    /** \brief The processor overlap for parallel computing.  Always zero because
        this is a strictly sequential grid */
    int overlapSize(int codim) const {
        return 0;
    }

    /** \brief The processor ghost overlap for parallel computing.  Always zero because
        this is a strictly sequential grid */
    int ghostSize(int codim) const {
        return 0;
    }

    /** \brief The processor overlap for parallel computing.  Always zero because
        this is a strictly sequential grid */
    int overlapSize(int level, int codim) const {
        return 0;
    }

    /** \brief The processor ghost overlap for parallel computing.  Always zero because
        this is a strictly sequential grid */
    int ghostSize(int level, int codim) const {
        return 0;
    }

    /** \brief Get the set of global ids */
    const typename Traits::GlobalIdSet& globalIdSet() const
    {
        return idSet_;
    }

    /** \brief Get the set of local ids */
    const typename Traits::LocalIdSet& localIdSet() const
    {
        return idSet_;
    }

    /** \brief Get an index set for the given level */
    const typename Traits::LevelIndexSet& levelIndexSet(int level) const
    {
        if (! levelIndexSets_[level]) {
            levelIndexSets_[level] =
                new OneDInNDGridLevelIndexSet<const ThisClass>(*this, level);
            levelIndexSets_[level]->update();
        }

        return * levelIndexSets_[level];
    }

    /** \brief Get an index set for the leaf level */
    const typename Traits::LeafIndexSet& leafIndexSet() const
    {
        return leafIndexSet_;
    }


    /** \brief Mark entity for refinement
     *
     * \param refCount if >0 mark for refinement, if <0 mark for coarsening
     * \param e EntityPointer to the entity you want to mark
     *
     * \return True, if marking was successfull
     */
    bool mark(int refCount, const typename Traits::template Codim<0>::EntityPointer& e );

    /** \brief return current adaptation marker of given entity

        \param e Entity to the entity you want to mark

        \return int current adaptation marker of entity pointer e
    */
    int getMark(const typename Traits::template Codim<0>::EntityPointer & e ) const;

    //! Does nothing except return true if some element has been marked for refinement
    bool preAdapt();

    //! Triggers the grid refinement process
    bool adapt();

    /** \brief Adaptation post-processing: Reset all adaptation state flags */
    void postAdapt();

    /** \brief grid identification */
    std::string name () const { return "OneDInNDGrid"; }

    // **********************************************************
    // End of Interface Methods
    // **********************************************************

       /** \brief The different forms of grid refinement supported by OneDInNDGrid */
    enum RefinementType {
        /** \brief New level consists only of the refined elements */
        LOCAL,
        /** \brief New level consists of the refined elements and the unrefined ones, too */
        COPY};

   /** \brief Sets the type of grid refinement */
    void setRefinementType(RefinementType type) {
        refinementType_ = type;
    }

    /** \brief Does one uniform refinement step
     *
     * \param refCount I don't know what this is good for.  It doesn't
     *        actually do anything.
     */
    void globalRefine(int refCount);

  // dummy parallel functions

  template<class DataHandle>
  void communicate (DataHandle& data, InterfaceType iftype, CommunicationDirection dir, int level) const
  {
  }

  template<class DataHandle>
  void communicate (DataHandle& data, InterfaceType iftype, CommunicationDirection dir) const
  {
  }

  const CollectiveCommunication<ThisClass>& comm () const
  {
  return ccobj;
  }


private:

  CollectiveCommunication<ThisClass> ccobj;

    /** \brief Update all indices and ids */
    void setIndices();

    unsigned int getNextFreeId(int codim) {
        return (codim==0) ? freeElementIdCounter_++ : freeVertexIdCounter_++;
    }

   //! The type of grid refinement currently in use
    RefinementType refinementType_;

    OneDInNDEntityImp<0, dimworld>* getLeftUpperVertex(const OneDInNDEntityImp<1, dimworld>* eIt);

    OneDInNDEntityImp<0, dimworld>* getRightUpperVertex(const OneDInNDEntityImp<1, dimworld>* eIt);

    /** \brief Returns an iterator the the first element on the left of
        the input element which has sons.
    */
    OneDInNDEntityImp<1, dimworld>* getLeftNeighborWithSon(OneDInNDEntityImp<1, dimworld>* eIt);

    // The vertices of the grid hierarchy
    std::vector<List<OneDInNDEntityImp<0, dimworld> > > vertices;

    // The elements of the grid hierarchy
    std::vector<List<OneDInNDEntityImp<1, dimworld> > > elements;

    // Our set of level indices
    mutable std::vector<OneDInNDGridLevelIndexSet<const ThisClass>* > levelIndexSets_;

    OneDInNDGridLeafIndexSet<const ThisClass> leafIndexSet_;

    OneDInNDGridIdSet<const ThisClass> idSet_;

    unsigned int freeVertexIdCounter_;

    unsigned int freeElementIdCounter_;

}; // end Class OneDInNDGrid

namespace Capabilities
{
  /** \struct hasBackupRestoreFacilities
  \ingroup OneDInNDGrid
  */

  /** \struct IsUnstructured
  \ingroup OneDInNDGrid
  */

  /** \brief OneDInNDGrid has entities for all codimension
  \ingroup OneDInNDGrid
  */
  template<int dimworld, int cdim>
  struct hasEntity< OneDInNDGrid<dimworld>, cdim >
  {
    static const bool v = true;
  };

  /** \brief OneDInNDGrid is not parallel
  \ingroup OneDInNDGrid
  */
  template<int dimworld>
  struct isParallel< OneDInNDGrid<dimworld> >
  {
    static const bool v = false;
  };

  /** \brief OneDInNDGrid is levelwise conforming
  \ingroup OneDInNDGrid
  */
  template<int dimworld>
  struct isLevelwiseConforming< OneDInNDGrid<dimworld> >
  {
    static const bool v = true;
  };

  /** \brief OneDInNDGrid is leafwise conforming
  \ingroup OneDInNDGrid
  */
  template<int dimworld>
  struct isLeafwiseConforming< OneDInNDGrid<dimworld> >
  {
    static const bool v = true;
  };

  /** \brief OneDInNDGrid does not support hanging nodes
  \ingroup OneDInNDGrid
  */
  template<int dimworld>
  struct hasHangingNodes< OneDInNDGrid<dimworld> >
  {
    static const bool v = false;
  };

}

} // namespace Dune

#endif
