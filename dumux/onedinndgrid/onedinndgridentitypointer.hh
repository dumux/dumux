#ifndef DUNE_ONEDINNDGRID_ENTITY_POINTER_HH
#define DUNE_ONEDINNDGRID_ENTITY_POINTER_HH


namespace Dune {

/*! Acts as a pointer to an  entities of a given codimension.
 */
template<int codim, class GridImp>
class OneDInNDGridEntityPointer
{
  enum { dim = GridImp::dimension };
  enum { dimworld = GridImp::dimensionworld };
    template <class GridImp_>
    friend class OneDInNDGridLevelIntersectionIterator;
    template <class GridImp_>
    friend class OneDInNDGridLeafIntersectionIterator;
    friend class OneDInNDGridEntity<0,dim,GridImp>;
    
public:
  typedef typename GridImp::template Codim<codim>::Entity Entity;
  typedef OneDInNDGridEntityPointer<codim,GridImp> Base;

    /** \brief The type of the class itself
        Do we really need this?
    */
    typedef OneDInNDGridEntityPointer<codim,GridImp> EntityPointerImp;

    //! codimension of entity pointer 
  enum { codimension = codim };

  //! equality
  bool equals(const OneDInNDGridEntityPointer<codim,GridImp>& other) const {
      return other.virtualEntity_.target() == virtualEntity_.target();
  }

  //! dereferencing
  Entity& dereference() const {return virtualEntity_;}
  
  //! ask for level of entity
  int level () const {return virtualEntity_.level();}

    OneDInNDGridEntityPointer() {}

  /** \brief Constructor from a given entity  */
  OneDInNDGridEntityPointer(const OneDInNDGridEntity<codim, dim, GridImp> & entity) 
  {
      virtualEntity_.setToTarget(entity.target_);
  }

  //! empty method since internal entity is not a pointer 
  void compactify () {}

protected:

    /** \brief Constructor from a given iterator */
  OneDInNDGridEntityPointer(OneDInNDEntityImp<dim-codim, dimworld>* it) {
        virtualEntity_.setToTarget(it);
    };

protected:

  mutable OneDInNDEntityWrapper<codim,GridImp::dimension,GridImp,GridImp::dimensionworld> virtualEntity_;

};


} // end namespace Dune

#endif
