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
/*! An interface for defining a numerical transport model for the
 *  solution of equations of the form

 - Grid      a DUNE grid type
 - RT        type used for return values
 - RepresentationType   type of the vector holding the saturation values
 - VelType   type of the vector holding the velocity values

 */
template<class Grid, class Scalar, class VC> class ShallowWater
{
public:
	
	enum {dim = Grid::dimension};
	typedef Dune::BlockVector<Dune::FieldVector<Scalar,dim+1> > SolutionType;
	typedef Dune::BlockVector<Dune::FieldVector<Scalar,dim+1> > RepresentationType;
	const Grid& grid;
	
	typename Dune::ShallowProblemBase<Grid,Scalar,VC>& problem;
	
	
	
	virtual int update(const Scalar t, Scalar& dt, RepresentationType& updateVec,
			Scalar& CLFFac) = 0;

	void initial()
	{
		initialize();
		return;
	}

	//! \brief Sets the initial solution 
	virtual void initialize() = 0 ;

	//! return const reference to globalSolution vector
	//virtual const RepresentationType& operator*() const
	//{
		//return problem.variables.globalSolution;
	//}

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
	ShallowWater(const Grid& grid, ShallowProblemBase<Grid, Scalar, VC>& problem):
		grid(grid), problem(problem)
	{
	}

		
};

}
#endif
