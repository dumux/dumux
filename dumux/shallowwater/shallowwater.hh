// $Id: transport.hh 566 2008-09-11 11:38:31Z bernd $ 

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
 *  \f$S_t - \text{div}\, (f_\text{w}(S) \boldsymbol{v}_\text{total}) = 0\f$,
 * \f$S = g\f$ on \f$\Gamma_1\f$, and \f$S(t = 0) = S_0\f$. Here,
 * \f$S\f$ denotes the wetting phase saturation,
 * \f$\boldsymbol{v}_\text{total}\f$ the total velocity,
 * and \f$f_\text{w}\f$ the wetting phase fractional flow function.

 - Grid      a DUNE grid type
 - RT        type used for return values
 - RepresentationType   type of the vector holding the saturation values
 - VelType   type of the vector holding the velocity values

 */
template<class G, class DT, class VC> class ShallowWater
{
public:
	
	enum {dim = G::dimension};
	typedef Dune::BlockVector<Dune::FieldVector<DT,dim+1> > SolutionType;
	typedef Dune::BlockVector<Dune::FieldVector<DT,dim+1> > RepresentationType;
	const G& grid;
	
	typename Dune::ShallowProblemBase<G,DT,VC>& problem;
	
	
	
	virtual int update(const DT t, DT& dt, RepresentationType& updateVec,
			DT& CLFFac) = 0;

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
	ShallowWater(const G& g, ShallowProblemBase<G, DT, VC>& problem):
		grid(g), problem(problem)
	{
	}

		
};

}
#endif
