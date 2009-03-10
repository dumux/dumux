// $Id$

#ifndef DUNE_SHALLOWPROBLEMBASE_HH
#define DUNE_SHALLOWPROBLEMBASE_HH

#include<iostream>
#include<iomanip>

#include<dune/common/fvector.hh>
#include<dune/common/fmatrix.hh>
#include<dune/common/exceptions.hh>
#include<dune/grid/common/grid.hh>
#include<dune/grid/common/referenceelements.hh>
#include<dune/disc/operators/boundaryconditions.hh>
#include<dumux/shallowwater/solidsurfacebase.hh>
#include<dumux/shallowwater/shallowvariableclass.hh>

/**
 * @file
 * @brief  Base class for defining an instance of the transport problem
 * @author Bernd Flemisch
 */

namespace Dune
{
/**  \ingroup transport
 *   \defgroup transportProblems Problems
 */
/*! \ingroup transportProblems
 *   *- Template parameters are:
 *   *- Grid  a DUNE grid type
 *	- DT    type used for return values
 */

template<class Grid, class Scalar, class VC> class ShallowProblemBase
{
    enum
    {   dim=Grid::dimension, m=1, blocksize=2*Grid::dimension};
    typedef typename Grid::Traits::template Codim<0>::Entity Element;
    typedef FieldVector<Scalar,dim> LocalPosition;
    typedef FieldVector<Scalar,dim> GlobalPosition;
    typedef Dune::FieldVector<Scalar, dim> VelType;
    typedef Dune::FieldVector<Scalar,dim+1> SystemType;

public:
    typedef Dune::SolidSurfaceBase<Grid,Scalar> Surface;

    // return type of boundary condition at the given global coordinate

    virtual BoundaryConditions::Flags bctype(const GlobalPosition& globalPos,
            const Element& element, const LocalPosition& localPos) const = 0;

    //! evaluate Dirichlet boundary condition at given position
     
    virtual Scalar dirichlet(const GlobalPosition& globalPos,
            const Element& element, const LocalPosition& localPos) const = 0;

    // evaluate Neumann boundary condition at given position
  
    virtual SystemType neumann(const GlobalPosition& globalPos,
            const Element& element, const LocalPosition& localPos) const = 0;

    // evaluate initial boundary condition at given position
     
    virtual Scalar setInitWDepth(const GlobalPosition& globalPos,
            const Element& element, const LocalPosition& localPos) const = 0;

    virtual VelType setInitVel(const GlobalPosition& globalPos,
            const Element& element, const LocalPosition& localPos) const = 0;

    //Definition of sources

    virtual Scalar setSource(const GlobalPosition& globalPos,
            const Element& element, const LocalPosition& localPos) const = 0;

    ShallowProblemBase(VC& variableobject, Surface& surfaceobject) :
        variables(variableobject), surface(surfaceobject)
    {
    }

    //! always define virtual destructor in abstract base class
    virtual ~ShallowProblemBase()
    {
    }

    VC& variables;
    Surface& surface;

};
}
#endif
