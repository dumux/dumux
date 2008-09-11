// $Id$ 

#ifndef DUNE_ONE_D_IN_N_D_GEOMETRY_HH
#define DUNE_ONE_D_IN_N_D_GEOMETRY_HH

/** \file
 * \brief The OneDInNDGridElement class and its specializations
 */

namespace Dune {

    template<int mydim, class GridImp>
    class OneDInNDMakeableGeometry : public Geometry<mydim, GridImp::dimensionworld, GridImp, OneDInNDGridGeometry>
    {
    	enum {dimworld = GridImp::dimensionworld};
    public:

        OneDInNDMakeableGeometry() :
            Geometry<mydim, dimworld, GridImp, OneDInNDGridGeometry>(OneDInNDGridGeometry<mydim, dimworld, GridImp>())
        {};

        void setToTarget(OneDInNDEntityImp<mydim,dimworld>* target) {
            this->realGeometry.target_ = target;
        }

        void setPosition(double p) {
        	std::cout << "Here in setPosition with p = " << p << std::endl;
            this->realGeometry.storeCoordsLocally_ = true;
            this->realGeometry.pos_[0] = p;
        }
    };

    template<class GridImp>
    class OneDInNDMakeableGeometry<1,GridImp> : public Geometry<1, GridImp::dimensionworld, GridImp, OneDInNDGridGeometry>
    {
    	enum { dimworld = GridImp::dimensionworld };
    public:

        OneDInNDMakeableGeometry() :
            Geometry<1, dimworld, GridImp, OneDInNDGridGeometry>(OneDInNDGridGeometry<1, dimworld, GridImp>())
        {};

        void setToTarget(OneDInNDEntityImp<1, dimworld>* target) {
            this->realGeometry.target_ = target;
        }

        void setPositions(double p1, double p2) {
        	std::cout << "Here in setPositions with p1 = " << p1 << ", p2 = " << p2 << std::endl;
            this->realGeometry.storeCoordsLocally_ = true;
            this->realGeometry.pos_[0][0] = p1;
            this->realGeometry.pos_[1][0] = p2;
        }
    };

    // forward declaration
    template <int codim, int dim, class GridImp>
    class OneDInNDGridEntity;


template<class GridImp, int dimworld>  
class OneDInNDGridGeometry <0, dimworld, GridImp> : 
        public GeometryDefaultImplementation <0, dimworld, GridImp,OneDInNDGridGeometry>
{ 
    template <int codim_, int dim_, class GridImp_>
    friend class OneDInNDGridEntity;
    template <int mydim_, int coorddim_, class GridImp_>
    friend class OneDInNDGridGeometry;

public:

    OneDInNDGridGeometry() : storeCoordsLocally_(false) {}

    //! return the element type identifier (vertex)
    GeometryType type () const {return GeometryType(0);}

    //! return the number of corners of this element (==1)
    int corners () const {return 1;}

    //! access to coordinates of corners. Index is the number of the corner 
    const FieldVector<typename GridImp::ctype, dimworld>& operator[] (int i) const {
        return (storeCoordsLocally_) ? pos_ : target_->pos_;
    }

    /** \brief Maps a local coordinate within reference element to 
     * global coordinate in element  */
    FieldVector<typename GridImp::ctype, dimworld> global (const FieldVector<typename GridImp::ctype, 0>& local) const {
        return (storeCoordsLocally_) ? pos_ : target_->pos_;
    }

    /** \brief Maps a global coordinate within the element to a 
     * local coordinate in its reference element */
    FieldVector<typename GridImp::ctype, 0> local (const FieldVector<typename GridImp::ctype, dimworld>& global) const {
        FieldVector<typename GridImp::ctype, 0> l;
        return l;
    }

    /** \brief Returns true if the point is in the current element

    This method really doesn't make much sense for a zero-dimensional
    object.  It always returns 'true'.
    */
    bool checkInside(const FieldVector<typename GridImp::ctype, 0> &local) const {
        return true;
    }


    /** \brief !!!

    This method really doesn't make much sense for a zero-dimensional
    object.  It always returns '1'.
    */
    typename GridImp::ctype integrationElement (const FieldVector<typename GridImp::ctype, 0>& local) const {
        return 1;
    }

    //! The Jacobian matrix of the mapping from the reference element to this element
    const FieldMatrix<typename GridImp::ctype,0,0>& jacobianInverseTransposed (const FieldVector<typename GridImp::ctype, 0>& local) const {
        return jacInverse_;
    }

    //private:
    bool storeCoordsLocally_;

    // Stores the element corner positions if it is returned as geometryInFather
    FieldVector<typename GridImp::ctype, dimworld> pos_;

    OneDInNDEntityImp<0, dimworld>* target_;
    
    FieldMatrix<typename GridImp::ctype,0,0> jacInverse_;
};

//**********************************************************************
//
// --OneDInNDGridGeometry
  /** \brief Defines the geometry part of a mesh entity. 
   * \ingroup OneDInNDGrid
*/
template<int mydim, int coorddim, class GridImp>  
class OneDInNDGridGeometry : 
public GeometryDefaultImplementation <mydim, coorddim, GridImp, OneDInNDGridGeometry>
{ 
	enum { dimworld = GridImp::dimensionworld };

	template <int codim_, int dim_, class GridImp_>
    friend class OneDInNDGridEntity;

    template <int dimworld>
    friend class OneDInNDGrid;

    template <int cc_, int dim_, class GridImp_>
    friend class OneDInNDGridSubEntityFactory;

    template <class GridImp_>
    friend class OneDInNDGridLevelIntersectionIterator;
    template <class GridImp_>
    friend class OneDInNDGridLeafIntersectionIterator;

public:

    OneDInNDGridGeometry() : storeCoordsLocally_(false) {}

    /** \brief Return the element type identifier 
     *
     * OneDInNDGrid obviously supports only lines
     */
    GeometryType type () const {return GeometryType(1);}

    //! return the number of corners of this element. Corners are numbered 0...n-1
    int corners () const {return 2;}

    //! access to coordinates of corners. Index is the number of the corner 
    const FieldVector<typename GridImp::ctype, coorddim>& operator[](int i) const {
        assert(i==0 || i==1);
        return (storeCoordsLocally_) ? pos_[i] : target_->vertex_[i]->pos_;
    }

    /** \brief Maps a local coordinate within reference element to 
     * global coordinate in element  */
    FieldVector<typename GridImp::ctype, coorddim> global (const FieldVector<typename GridImp::ctype, mydim>& local) const {
        FieldVector<typename GridImp::ctype, coorddim> g;
        for (int k = 0; k < coorddim; k++)
        	g[k] = (storeCoordsLocally_)
            ? pos_[0][k] * (1-local[0]) + pos_[1][k] * local[0]
            : target_->vertex_[0]->pos_[k] * (1-local[0]) + target_->vertex_[1]->pos_[k] * local[0];
        return g;
    }

    /** \brief Maps a global coordinate within the element to a 
     * local coordinate in its reference element */
    FieldVector<typename GridImp::ctype, mydim> local (const FieldVector<typename GridImp::ctype, coorddim>& global) const {
        FieldVector<typename GridImp::ctype, mydim> l;
        if (storeCoordsLocally_) {
            l[0] = (global[0] - pos_[0][0]) / (pos_[1][0] - pos_[0][0]);
        } else {
            const double& v0 = target_->vertex_[0]->pos_[0];
            const double& v1 = target_->vertex_[1]->pos_[0];
            l[0] = (global[0] - v0) / (v1 - v0);
        }
        return l;
    }
    
    //! Returns true if the point is in the current element
    bool checkInside(const FieldVector<typename GridImp::ctype, coorddim> &global) const {
        return (storeCoordsLocally_)
            ? pos_[0][0] <= global[0] && global[0] <= pos_[1][0]
            : target_->vertex_[0]->pos_[0] <= global[0] && global[0] <= target_->vertex_[1]->pos_[0];
    }

    /** ???
   */
    typename GridImp::ctype integrationElement (const FieldVector<typename GridImp::ctype, mydim>& local) const {
        return (storeCoordsLocally_)
            ? pos_[1][0] - pos_[0][0]
            : target_->vertex_[1]->pos_[0] - target_->vertex_[0]->pos_[0];
    }

    //! The Jacobian matrix of the mapping from the reference element to this element
    const FieldMatrix<typename GridImp::ctype,mydim,mydim>& jacobianInverseTransposed (const FieldVector<typename GridImp::ctype, mydim>& local) const {
        if (storeCoordsLocally_)
            jacInverse_[0][0] = 1 / (pos_[1][0] - pos_[0][0]);
        else
            jacInverse_[0][0] = 1 / (target_->vertex_[1]->pos_[0] - target_->vertex_[0]->pos_[0]);

        return jacInverse_;
    }


    //private:
    OneDInNDEntityImp<1, dimworld>* target_;

    bool storeCoordsLocally_;

    // Stores the element corner positions if it is returned as geometryInFather
    FieldVector<typename GridImp::ctype,coorddim> pos_[2];

    //! The jacobian inverse
    // \todo Should be of different dimension
    mutable FieldMatrix<typename GridImp::ctype,mydim,mydim> jacInverse_;

};

}  // namespace Dune

#endif
