// $Id$
/*****************************************************************************
 *   Copyright (C) <YEARS> by <ADD_AUTHOR_HERE>                              *
 *   Institute of Hydraulic Engineering                                      *
 *   University of Stuttgart, Germany                                        *
 *   email: <givenname>.<name>@iws.uni-stuttgart.de                          *
 *                                                                           *
 *   This program is free software; you can redistribute it and/or modify    *
 *   it under the terms of the GNU General Public License as published by    *
 *   the Free Software Foundation; either version 2 of the License, or       *
 *   (at your option) any later version, as long as this copyright notice    *
 *   is included in its original form.                                       *
 *                                                                           *
 *   This program is distributed WITHOUT ANY WARRANTY.                       *
 *****************************************************************************/

#ifndef DUNE_CRSHAPEFUNCTIONS_HH
#define DUNE_CRSHAPEFUNCTIONS_HH

#include<iostream>
#include<dune/common/fvector.hh>
#include<dune/common/exceptions.hh>
#include<dune/common/misc.hh>
#include<dune/common/geometrytype.hh>
#include"shapefunctions.hh"

#include"cr/cubeshapefunctions.hh"
#include"cr/simplexshapefunctions.hh"

/**
 * @file
 * @brief  define CR type shape functions
 * @author Peter Bastian
 */
namespace Dune
{
/** @addtogroup DISC_Shapefnkt
 *
 * @{
 */

/***********************************************************
 * The interface for CR shape functions.
 ***********************************************************/

/** \brief A scalar (N=1) ShapeFunction extended by a method providing a position
    \param C Type used for coordinates in the reference element
    \param T Type used for the shape function values
    \param d Dimension of the element
*/
template<typename C, typename T, int d>
class CRShapeFunction : public ShapeFunction<C,T,d,1>
{
public:
    // compile time sizes
    enum { dim=d };
    enum { comps=1 };
    typedef C CoordType;
    typedef T ResultType;

    //! interpolation point associated with shape function
    virtual const FieldVector<CoordType,dim>& position () const = 0;
};



/** \brief A scalar (N=1) ShapeFunctionSet that returns a CRShapeFunction
    \param C Type used for coordinates in the reference element
    \param T Type used for the shape function values
    \param d Dimension of the element
*/
template<typename C, typename T, int d>
class CRShapeFunctionSet : public ShapeFunctionSet<C,T,d,1>
{
public:
    // compile time sizes
    enum { dim=d };
    enum { comps=1 };

    // exported types
    typedef C CoordType;
    typedef T ResultType;
    typedef CRShapeFunction<CoordType,ResultType,dim> value_type;

    //! random access to i'th CRShapeFunction
    virtual const value_type& operator[] (int i) const = 0;
};


/***********************************************************
 * Wrappers
 *
 ***********************************************************/

/*! wrap inlinable implementation into class that is derived
 * from abstract base class and implement the functions with
 * the given implementation
 */
template<typename Imp>
class CRShapeFunctionWrapper :
        public CRShapeFunction<typename Imp::CoordType,typename Imp::ResultType,Imp::dim>,
        private Imp
{
public:

    // compile time sizes
    enum { dim=Imp::dim };
    enum { comps=1 };

    // exported types
    typedef typename Imp::CoordType CoordType;
    typedef typename Imp::ResultType ResultType;
    typedef Imp ImplementationType;

    //! assignment from implementation type (this class has no data)
    CRShapeFunctionWrapper& operator= (const Imp& imp)
    {
        Imp::operator=(imp);
        return *this;
    }

    //! evaluate component comp at point x
    virtual ResultType evaluateFunction (int comp, const FieldVector<CoordType,dim>& x) const
    {
        return Imp::evaluateFunction(comp,x);
    }

    //! evaluate derivative of component comp in direction dir at point x
    virtual ResultType evaluateDerivative (int comp, int dir, const FieldVector<CoordType,dim>& x) const
    {
        return Imp::evaluateDerivative(comp,dir,x);
    }

    //! consecutive number of associated dof within element
    virtual int localindex (int comp) const
    {
        return Imp::localindex(comp);
    }

    //! codim of associated dof
    virtual int codim () const
    {
        return Imp::codim();
    }

    //! entity (of codim) of associated dof
    virtual int entity () const
    {
        return Imp::entity();
    }

    //! consecutive number of dof within entity
    virtual int entityindex () const
    {
        return Imp::entityindex();
    }

    //! interpolation point associated with shape function
    virtual const FieldVector<CoordType,dim>& position () const
    {
        return Imp::position();
    }
};




/*! wrap inlinable implementation into class that is derived
 * from abstract base class and implement the functions with
 * the given implementation
 */
template<typename Imp>
class CRShapeFunctionSetWrapper :
        public CRShapeFunctionSet<typename Imp::CoordType,typename Imp::ResultType,Imp::dim>,
        private Imp
{
public:

    // compile time sizes
    enum { dim=Imp::dim };
    enum { comps=1 };     // must be available at compile time

    // exported types
    typedef typename Imp::CoordType CoordType;
    typedef typename Imp::ResultType ResultType;
    typedef Imp ImplementationType;
    typedef CRShapeFunction<CoordType,ResultType,dim> value_type; // Note: Imp::value_type references
    // must be convertible to references to
    // this type !
    //! return total number of shape functions
    virtual int size () const
    {
        return Imp::size();
    }

    //! total number of shape functions associated with entity in codim
    virtual int size (int entity, int codim) const
    {
        return Imp::size(entity,codim);
    }

    //! random access to i'th ShapeFunction
    virtual const value_type& operator[] (int i) const
    {
        return Imp::operator[](i);
    }

    //! return order
    virtual int order () const
    {
        return Imp::order();
    }

    //! return type of element
    virtual GeometryType type () const
    {
        return Imp::type();
    }
};

/***********************************************************
 * The general container for CR shape functions.
 * All containers are accessible in a singleton.
 ***********************************************************/

/** \brief This are CR shape functions
    \param C Type used for coordinates in the reference element
    \param T Type used for the shape function values
    \param d Dimension of the element
*/
template<typename C, typename T, int d>
class CRShapeFunctionSetContainer : public ShapeFunctionSetContainer<C,T,d,1,Power_m_p<3,d>::power >
{
public:
    // compile time sizes
    enum { dim=d };
    enum { comps=1 };
    enum { maxsize=Power_m_p<3,dim>::power };

    // exported types
    typedef C CoordType;
    typedef T ResultType;

    //! type of objects in the container
    typedef CRShapeFunctionSet<C,T,d> value_type;

    const value_type& operator() (GeometryType type, int order) const
    {
        if ( type.isCube() )
        {
            if (order==1) return wrappedcube;
            DUNE_THROW(RangeError, "order not available for cubes");
        }

        if ( type.isSimplex() )
        {
            if (order==1) return wrappedsimplex;
            DUNE_THROW(RangeError, "order not available for simplex");
        }

        if ( type.isPyramid() )
        {
            DUNE_THROW(RangeError, "No pyramid for this dimension");
        }

        if ( type.isPrism() )
        {
            DUNE_THROW(RangeError, "No prism for this dimension ");
        }

        DUNE_THROW(RangeError, "type or order not available");
    }

private:
    // the cubes
    typedef CRShapeFunctionWrapper<CRCubeShapeFunction<C,T,d> > WrappedCubeShapeFunction;
    typedef CRCubeShapeFunctionSet<C,T,d,WrappedCubeShapeFunction> CRCubeWrappedShapeFunctionSet;
    typedef CRShapeFunctionSetWrapper<CRCubeWrappedShapeFunctionSet> WrappedCRCubeShapeFunctionSet;
    WrappedCRCubeShapeFunctionSet wrappedcube;

    // the simplices
    typedef CRShapeFunctionWrapper<CRSimplexShapeFunction<C,T,d> > WrappedSimplexShapeFunction;
    typedef CRSimplexShapeFunctionSet<C,T,d,WrappedSimplexShapeFunction> CRSimplexWrappedShapeFunctionSet;
    typedef CRShapeFunctionSetWrapper<CRSimplexWrappedShapeFunctionSet> WrappedCRSimplexShapeFunctionSet;
    WrappedCRSimplexShapeFunctionSet wrappedsimplex;
};

//! This are CR shape functions for any element type and order (in the future ... )
template<typename C, typename T>
class CRShapeFunctionSetContainer<C,T,3> : public ShapeFunctionSetContainer<C,T,3,1,Power_m_p<3,3>::power >
{
public:
    // compile time sizes
    enum { dim=3 };
    enum { comps=1 };
    enum { maxsize=Power_m_p<3,dim>::power };

    // exported types
    typedef C CoordType;
    typedef T ResultType;

    //! type of objects in the container
    typedef CRShapeFunctionSet<C,T,3> value_type;

    const value_type& operator() (GeometryType type, int order) const
    {
        if ( type.isCube() )
        {
            if (order==1) return wrappedcube;
            DUNE_THROW(RangeError, "order not available for cubes");
        }

        if ( type.isSimplex() )
        {
            if (order==1) return wrappedsimplex;
            DUNE_THROW(RangeError, "order not available for simplex");
        }

        DUNE_THROW(RangeError, "type or order not available");
    }

private:
    // the cubes
    typedef CRShapeFunctionWrapper<CRCubeShapeFunction<C,T,dim> > WrappedCubeShapeFunction;
    typedef CRCubeShapeFunctionSet<C,T,dim,WrappedCubeShapeFunction> CRCubeWrappedShapeFunctionSet;
    typedef CRShapeFunctionSetWrapper<CRCubeWrappedShapeFunctionSet> WrappedCRCubeShapeFunctionSet;
    WrappedCRCubeShapeFunctionSet wrappedcube;

    // the simplices
    typedef CRShapeFunctionWrapper<CRSimplexShapeFunction<C,T,dim> > WrappedSimplexShapeFunction;
    typedef CRSimplexShapeFunctionSet<C,T,dim,WrappedSimplexShapeFunction> CRSimplexWrappedShapeFunctionSet;
    typedef CRShapeFunctionSetWrapper<CRSimplexWrappedShapeFunctionSet> WrappedCRSimplexShapeFunctionSet;
    WrappedCRSimplexShapeFunctionSet wrappedsimplex;
};


// singleton holding several reference element containers
template<typename C, typename T, int d>
struct CRShapeFunctions {
    static CRCubeShapeFunctionSetContainer<C,T,d> cube;
    static CRSimplexShapeFunctionSetContainer<C,T,d> simplex;
    static CRShapeFunctionSetContainer<C,T,d> general;
};

/** @} */
}
#endif
