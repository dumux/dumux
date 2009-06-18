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

#ifndef DUNE_CRSIMPLEXSHAPEFUNCTIONS_HH
#define DUNE_CRSIMPLEXSHAPEFUNCTIONS_HH

#include<iostream>
#include<dune/common/fvector.hh>
#include<dune/common/exceptions.hh>
#include<dune/common/misc.hh>
#include<dune/common/geometrytype.hh>
#include<dune/grid/common/referenceelements.hh>

/**
 * @file
 * @brief  define Lagrange type shape functions
 * @author Peter Bastian
 */
namespace Dune
{
/** @addtogroup DISC_Shapefnkt
 *
 * @{
 */
/**
 * @brief Shape functions for Lagrange finite elements
 *
 */

/***********************************************************
 * CR shape functions for the simplex of any dimension
 * The basic classes implement the interface of a
 * LagrangeShapeFunction(Set) but are NOT derived from
 * it. They can be used in simple cases (fixed element type and
 * order) to avoid virtual function calls. The generally usable
 * classes are then wrapped up with the Wrappers to supply
 * the virtual functions of the abstract interface.
 ***********************************************************/

/*!
 * A class for piecewise bilinear shape functions in the simplex.
 * This class implements the interface of ShapeFunction but
 * is NOT derived from it to avoid making the functions virtual at this
 * stage.
 *
 * Let \f$i=(i_{dim-1},...,i_1,i_0)\f$ be the binary representation of the shape
 * function number. Then the corresponding shape function can be written as
 *
 *    \f[phi_i (x) = \prod_{j=0}^{dim-1} [ 1-i_j + x_j*(2*i_j-1) ]\f]
 *
 * and its derivative is
 *
 *    \f[d/dx_k phi_i (x) = (2*i_k-1) * \prod_{j!=k} [ 1-i_j + x_j*(2*i_j-1) ]\f]
 *
 * The coefficients \f$a_{ij} = 1-i_j\f$ and \f$b_{ij} = 2*i_j-1\f$ are precomputed.
 */
template<typename C, typename T, int d>
class CRSimplexShapeFunction
{
public:

    // compile time sizes
    enum { dim=d };    // maps from R^d
    enum { comps=1 };      // to R^1

    enum { m=1+dim }; // total number of basis functions

    // export types
    typedef C CoordType;
    typedef T ResultType;
    typedef CRSimplexShapeFunction ImplementationType;

    /*! \brief make a shape function object
     *
     *  \todo check the local ordering for 2D and 3D
     */
    CRSimplexShapeFunction (int i) // make it the i'th shape function
    {
        number = i;
        if (dim == 2) {
            switch (i) {
            case 2:
                pos[0] = 0.5; pos[1] = 0;
                a[0] = 1; a[1] = 0; a[2] = -2;
                break;
            case 0:
                pos[0] = 0.5; pos[1] = 0.5;
                a[0] = -1; a[1] = 2; a[2] = 2;
                break;
            case 1:
                pos[0] = 0; pos[1] = 0.5;
                a[0] = 1; a[1] = -2; a[2] = 0;
                break;
            }
        }
        else if (dim == 3) {
            switch (i) {
            case 0:
                pos[0] = 0.5; pos[1] = 0.5; pos[2] = 0;
                a[0] = 1; a[1] = 0; a[2] = 0; a[3] = -2;
                break;
            case 1:
                pos[0] = 0; pos[1] = 0.5; pos[2] = 0.5;
                a[0] = 1; a[1] = -2; a[2] = 0; a[3] = 0;
                break;
            case 2:
                pos[0] = 0.5; pos[1] = 0; pos[2] = 0.5;
                a[0] = 1; a[1] = 0; a[2] = -2; a[3] = 0;
                break;
            case 3:
                pos[0] = 0.5; pos[1] = 0.5; pos[2] = 0.5;
                a[0] = -2; a[1] = 2; a[2] = 2; a[3] = 2;
                break;
            }
        }
    }

    //! must be defaultconstructible
    CRSimplexShapeFunction ()
    {}

    //! evaluate shape function in local coordinates
    ResultType evaluateFunction (int comp, const FieldVector<CoordType,d>& x) const
    {
        ResultType phi = a[0];
        for (int i = 1; i < dim+1; i++)
            phi += a[i]*x[i-1];
        return phi;
    }

    //! evaluate gradient in local coordinates
    ResultType evaluateDerivative (int comp, int dir, const FieldVector<CoordType,d>& x) const
    {
        return(a[dir+1]);
    }

    //! consecutive number of associated dof within element
    int localindex (int comp) const
    {
        return number;
    }

    //! codim of associated dof
    int codim () const
    {
        return 1;
    }

    //! entity (of codim) of associated dof
    int entity () const
    {
        return number;
    }

    //! consecutive number of dof within entity
    int entityindex () const
    {
        return 0;
    }

    //! interpolation point associated with shape function
    const FieldVector<CoordType,dim>& position () const
    {
        return pos;
    }

private:
    int number;
    ResultType a[dim+1]; // store coefficients for this shape function
    FieldVector<CoordType,d> pos;
};


/*! CRSimplexShapeFunctionSet implements the interface of
 * LagrangeShapeFunctionSet but is NOT derived from it.
 * S is either CRSimplexShapeFunction<..> or LagrangeShapeFunctionWrapper<CRSimplexShapeFunction<..> >
 */
template<typename C, typename T, int d, typename S>
class CRSimplexShapeFunctionSet
{
public:

    // compile time sizes
    enum { dim=d };    // maps from R^d
    enum { comps=1 };      // to R^1

    enum { m=1+dim }; // total number of basis functions

    // export types
    typedef C CoordType;
    typedef T ResultType;
    typedef S value_type;
    typedef typename S::ImplementationType Imp; // Imp is either S or derived from S

    //! make a shape function object
    CRSimplexShapeFunctionSet ()
    {
        for (int i=0; i<m; i++)
            sf[i] = Imp(i); // assignment of derived class objects defined in wrapper
    }

    //! return total number of shape functions
    int size () const
    {
        return m;
    }

    //! total number of shape functions associated with entity in codim
    int size (int entity, int codim) const
    {
        if (codim==1) return 1; else return 0;
    }

    //! random access to shape functions
    const value_type& operator[] (int i) const
    {
        return sf[i]; // ok derived class reference goes for base class reference
    }

    //! return order
    int order () const
    {
        return 1;
    }

    //! return type of element
    GeometryType type () const
    {
        static GeometryType simplex(GeometryType::simplex, dim);
        return simplex;
    }

private:
    S sf[m];
};



//! This are CR shape functions in the simplex without virtual functions
template<typename C, typename T, int d>
class CRSimplexShapeFunctionSetContainer
{
public:
    // compile time sizes
    enum { dim=d };
    enum { comps=1 };
    enum { maxsize=dim+1 };

    // exported types
    typedef C CoordType;
    typedef T ResultType;
    typedef CRSimplexShapeFunctionSet<C,T,d,CRSimplexShapeFunction<C,T,d> > value_type;

    const value_type& operator() (GeometryType type, int order) const
    {
        if (type.isSimplex())
            return simplex;
        DUNE_THROW(NotImplemented, "type not implemented yet");
    }
private:
    value_type simplex;
};



/** @} */
}
#endif
