// $Id$ 

#ifndef EXPLICITEULERSTEP_HH
#define EXPLICITEULERSTEP_HH

#include <dumux/timedisc/timestep.hh>

namespace Dune
{
	//! \brief For models with postupdate routine
	/** Class Model must contain a method postupdate(double t, double dt)
	 * with t the time at begin of timestep and dt the timestep length.
	 */
	template<class G, class Model>
	class ExplicitEulerStep : public TimeStep<G, Model>
	{
	public:
			/** Performs a simple Euler just as the RungeKuttaStep with stages == 1.
			 *  After updating the primary variable vector of the model, the postupdate routine
			 *  of the model is called.
			 */
			void execute(Model& model, double t, double& dt,
							double maxDt, double tEnd, double cFLFactor)
			{
				  // allocate temporary vectors for the updates
				  typedef typename Model::RepresentationType RepresentationType;
				  RepresentationType k1 = *model;

				  // obtain the first update and the time step size
				  model.update(t, dt, k1);

				  // scale dt with safety factor
				  dt *= cFLFactor;
				  dt = std::min( dt, maxDt );
				  dt = std::min( dt, tEnd - t);

			      // explicit Euler: Sat <- Sat + dt*N(Sat)
				  *model += (k1 *= dt);

				  model.postupdate(t, dt);
			}

	};

}//end namespace

#endif /*EXPLICITEULERSTEP_HH*/
