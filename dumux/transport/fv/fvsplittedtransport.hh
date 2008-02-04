#ifndef DUNE_FVSPLITTEDTRANSPORT_HH
#define DUNE_FVSPLITTEDTRANSPORT_HH

#include "dumux/transport/splittedtransport.hh"
#include "dumux/transport/fv/fvtransport.hh"

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
  template<class G, class RT>
  class FVSplittedTransport : public SplittedTransport< 
                                          G, RT, 
                                          BlockVector< 
                                            Dune::FieldVector<RT,1> >,
							BlockVector< Dune::FieldVector<Dune::FieldVector<RT, G::dimension>, 2*G::dimension> >, 
							FVTransport<G, RT>, FVTransport<G, RT> > 
  {
  public:
    typedef typename FVTransport<G, RT>::RepresentationType HyperbolicRepresentationType;
    typedef typename FVTransport<G, RT>::RepresentationType ParabolicRepresentationType;
    typedef typename FVTransport<G, RT>::RepresentationType RepresentationType;
    typedef FVTransport<G, RT>  HyperbolicType;
    typedef FVTransport<G, RT>  ParabolicType;
   

    virtual void transferHyperbolicToParabolic(const HyperbolicRepresentationType& hyperSat, 
					       ParabolicRepresentationType& paraSat)
    {
      paraSat = hyperSat;
      return;
    }

    virtual void transferParabolicToHyperbolic(const ParabolicRepresentationType& paraSat, 
					       HyperbolicRepresentationType& hyperSat)
    {
      hyperSat = paraSat;
      return;
    }

    virtual void transferParabolicToRepresentationType(const ParabolicRepresentationType& paraSat, 
					       RepresentationType& rSat) 
    {
      rSat = paraSat;
      return;
    }

    virtual void transferRepresentationTypeToParabolic(const RepresentationType& rSat, 
					       ParabolicRepresentationType& paraSat) 
    {
      paraSat = rSat;
      return;
    }

    //! generate vtk output
    virtual void vtkout (const char* name, int k) const 
    {
      const typename G::template Codim<0>::LevelIndexSet& iset(this->grid.levelIndexSet(this->parabolicLevel()) );
      Dune::VTKWriter<G, typename G::template Codim<0>::LevelIndexSet> 
	vtkwriter(this->grid, 
		  iset );
      char fname[128];	
      sprintf(fname,"%s-%05d",name,k);
      vtkwriter.addCellData(this->sat,"saturation");
      vtkwriter.write(fname,Dune::VTKOptions::ascii);		
    }

    /*! @brief constructor
     *  @param g a DUNE grid object
     *  @param prob an object of class TransportProblem or derived
     */
    FVSplittedTransport(const G& g, HyperbolicType& hyper, ParabolicType& para) 
      : SplittedTransport< G, RT, BlockVector< Dune::FieldVector<RT,1> >,
					BlockVector< Dune::FieldVector<Dune::FieldVector<RT, G::dimension>, 2*G::dimension> >, 
					FVTransport<G, RT>, FVTransport<G, RT> >(g, hyper, para)
    { 
      this->sat.resize(this->parabolicPart.sat.size());
    }
  };

}
#endif
