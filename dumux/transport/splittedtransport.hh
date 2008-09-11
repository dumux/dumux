// $Id$ 

#ifndef DUNE_SPLITTEDTRANSPORT_HH
#define DUNE_SPLITTEDTRANSPORT_HH

#include "dumux/transport/transportproblem.hh"

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
  template<class G, class RT, class VC, class HyperbolicType, class ParabolicType>
  class SplittedTransport {
  public:
	  typedef BlockVector< FieldVector<RT,1> > RepresentationType;
	  typedef typename HyperbolicType::RepresentationType HyperbolicRepresentationType;
	  typedef typename ParabolicType::RepresentationType ParabolicRepresentationType;
	    
    const G& grid;
    RepresentationType sat; //!< vector of saturation values
    HyperbolicType& hyperbolicPart; //!< for the hyperbolic part
    ParabolicType& parabolicPart; //!< for the parabolic part
    
    virtual void transferHyperbolicToParabolic(const HyperbolicRepresentationType& hyperSat, 
					       ParabolicRepresentationType& paraSat) = 0;

    virtual void transferParabolicToHyperbolic(const ParabolicRepresentationType& paraSat, 
					       HyperbolicRepresentationType& hyperSat) = 0;

    virtual void transferParabolicToRepresentationType(const ParabolicRepresentationType& paraSat, 
					       RepresentationType& rSat) = 0;

    virtual void transferRepresentationTypeToParabolic(const RepresentationType& rSat, 
					       ParabolicRepresentationType& paraSat) = 0;


	  
	//! \brief Calculate the update vector.
	/*!
	 *  \param[in]  t         time 
	 *  \param[out] dt        time step size
	 *  
	 *  Calculate the update vector, i.e., the discretization 
	 *  of \f$\text{div}\, (f_\text{w}(S) \boldsymbol{v}_t)\f$.
	 */
    virtual int update(RT t, RT& dt, RepresentationType& updateVec, double cFLFactor)
    {
      // parabolic part:
      transferRepresentationTypeToParabolic(sat, parabolicPart.transproblem.variables.saturation);
      parabolicPart.update(t, dt, updateParabolic, cFLFactor);
      ParabolicRepresentationType helpVector1(updateParabolic);
      helpVector1 *= dt*cFLFactor;
      parabolicPart.transproblem.variables.saturation += helpVector1;
      transferParabolicToHyperbolic(parabolicPart.transproblem.variables.saturation, hyperbolicPart.transproblem.variables.saturation);

      // hyperbolic part:
      double dummyDT;
      hyperbolicPart.update(t, dummyDT, updateHyperbolic, cFLFactor);
      //printvector(std::cout, updateHyperbolic, "updateHyperbolic", "row", 200, 1);    
      //updateHyperbolic = 0;
      
      // combine:
      ParabolicRepresentationType helpVector2(parabolicPart.transproblem.variables.saturation.size());
      transferHyperbolicToParabolic(updateHyperbolic, helpVector2);
      updateParabolic += helpVector2; 
      transferParabolicToRepresentationType(updateParabolic, updateVec);

      return (0);
    } 
		
    //! \brief Sets the initial solution \f$S_0\f$.
    virtual void initial() 
    {
      hyperbolicPart.initial();
      parabolicPart.initial();
      transferParabolicToRepresentationType(parabolicPart.transproblem.variables.saturation, sat);
      return;
    }
		
    //! generate vtk output
    virtual void vtkout (const char* name, int k) const = 0;
	
    //! return const reference to saturation vector
    const RepresentationType& operator* () const
    {
      return sat;
    }

    //! return reference to saturation vector
    RepresentationType& operator* ()
    {
      return sat;
    }

    //! always define virtual destructor in abstract base class
    virtual ~SplittedTransport () {}
    
    /*! @brief constructor
     *  @param g a DUNE grid object
     *  @param prob an object of class TransportProblem or derived
     */
    SplittedTransport(const G& g, HyperbolicType& hyper, ParabolicType& para) 
      : grid(g), hyperbolicPart(hyper), parabolicPart(para), level_(g.maxLevel())
    { 
      updateHyperbolic.resize(hyperbolicPart.transproblem.variables.saturation.size());
      updateParabolic.resize(parabolicPart.transproblem.variables.saturation.size());
    }
	
    //! returns the level on which the transport eqution is solved.
    int hyperbolicLevel() const
    {
      return hyperbolicPart.level();
    }
    
    //! returns the level on which the transport eqution is solved.
    int parabolicLevel() const
    {
      return parabolicPart.level();
    }
    
	
  protected:
    const int level_;
    HyperbolicRepresentationType updateHyperbolic;
    ParabolicRepresentationType updateParabolic;
  };

}
#endif
