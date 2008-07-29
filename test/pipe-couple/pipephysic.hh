/*
 * Calculates Lambda for a given fluid and hydrodynamic situation 
 * 
*/
template<class Fluid>
class Lambda {
public:
  Fluid fluid;

/*! \brief Implements the Lambda function for pipe roughness
 *   if Re < 2000 laminar flow otherwise turbulent flow 
 *   otherwise turbulent flow
 * 	 for the turbulent flow Zigrang and Sylvester eqn is used
 *  \param 
 *  \return
 */
  double calculate (double pressure, double velocity, double T, double roughness, double diameter) 
  { 
	double viscosity = fluid.viscosity(T,pressure); 
	double Re = fabs(velocity) * diameter / viscosity;
	double lambda;
	if (Re < 2000) // laminar
	{
		if (velocity != 0.0)
		 lambda = 64/Re;
		else
		 lambda = 0.0;
		
		return lambda;
	}
	else // turbulent
	{
		if (velocity != 0.0)
		 lambda = pow( -2*log( 2*(roughness/diameter)/3.7 - 5.02/Re*log( 2* (roughness/diameter)/3.7+13/Re ) ) ,-2 );
		else
		 lambda = 0.0;
		
		return lambda;
	}
  }
  
  Lambda()
  {
    ;
  }
};
