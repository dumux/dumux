// $Id:$

#ifndef DUNE_SHALLOWWATER_HH
#define DUNE_SHALLOWWATER_HH

#include "dumux/shallowwater/shallowproblembase.hh"

/**
 * @file
 * @brief  Base class for defining an instance of a numerical diffusion model
 * @author Bernd Flemisch
 * \defgroup transport Transport
 */

namespace Dune
{
//! \ingroup transport
//! Base class for defining an instance of a numerical transport model.

template<class Grid, class Scalar, class VC> class ShallowWater
{
public:

    enum
    {   dim = Grid::dimension};
    typedef Dune::BlockVector<Dune::FieldVector<Scalar,dim+1> > SolutionType;
    typedef Dune::BlockVector<Dune::FieldVector<Scalar,dim+1> >
            RepresentationType;
    const Grid& grid;

    typename Dune::ShallowProblemBase<Grid,Scalar,VC>& problem;

    //updates the solution vector
    virtual int
            update(const Scalar t, Scalar& dt, RepresentationType& updateVec) = 0;

    void initial()
    {
        initialize();
        return;
    }

    //initializes the variables 
    virtual void initialize() = 0;

    virtual void postProcessUpdate(Scalar t, Scalar dt)
    {
        return;
    }

    //! return reference to globalSolution vector
    virtual RepresentationType& operator*()
    {
        return problem.variables.globalSolution;
    }

    virtual void vtkout(const char* name, int k) const
    {
        problem.variables.vtkout(name, k);
        return;
    }

    //! always define virtual destructor in abstract base class
    virtual ~ShallowWater()
    {
    }

    /*! @brief constructor
     *  @param g a DUNE grid object
     *  @param prob an object of class TransportProblem or derived
     */
    ShallowWater(const Grid& grid, ShallowProblemBase<Grid, Scalar, VC>& problem) :
        grid(grid), problem(problem)
    {
    }

};

}
#endif
