#ifndef DUNE_IMPESMS_HH
#define DUNE_IMPESMS_HH

#include "dumux/fractionalflow/impes/impes.hh"
#include "dumux/transport/fv/diffusivepart.hh"

/**
 * @file
 * @brief  upscaled IMPES scheme
 * @author Bernd Flemisch, Jochen Fritz
 */

namespace Dune
{
  template<class G, class Diffusion, class Transport>
  class IMPESMS : public IMPES<G, Diffusion, Transport> 
  {
      enum{dim = G::dimension, dimworld = G::dimensionworld};
	  typedef typename Diffusion::RepresentationType PressType;
	  typedef typename Diffusion::NumberType RT;

  public:
	typedef typename G::ctype ct;
	typedef typename Transport::RepresentationType RepresentationType;

	virtual void totalVelocity(const RT t=0) 
	{
	  if ( this->level() == this->diffusion.level() )
		  this->diffusion.totalVelocity(this->problem.velocity,t);
	  else
		  this->diffusion.totalVelocity(this->problem.velocity,t,this->level());
	}

	
	virtual void initial()
	{
		double t = 0;
		Transport::initial();

		this->pressure(t);		
		totalVelocity(t);
	}
	
		
	int update(const RT t, RT& dt, RepresentationType& updateVec, RT cFLFactor = 1)
	{
		  typedef typename G::template Codim<0>::LevelIterator ElementLevelIterator;
		  typedef typename G::template Codim<0>::HierarchicIterator HierarchicIterator;
		  const typename G::Traits::LevelIndexSet& isetC ( this->grid.levelIndexSet(this->level()) );
		  const typename G::Traits::LevelIndexSet& isetF ( this->grid.levelIndexSet(this->diffusion.level()) );
		  
//		  RepresentationType saturation(isetF.size(0));
		  // entity pointer type
		  typedef typename G::template Codim<0>::EntityPointer ElementEntityPointer;
		  int maxlevel = this->grid.maxLevel();
		
		  int pressSize = (*(this->diffusion)).size();
		  int satSize = saturation.size();
		  PressType pressOldIter(*(this->diffusion));
		  PressType pressHelp(pressSize);  
		  RepresentationType satOldIter(saturation);
		  RepresentationType satHelp(satSize);
		  RepresentationType satDiff(satSize);
		  RepresentationType updateOldIter(satSize);
		  RepresentationType updateHelp(satSize);
		  RepresentationType updateDiff(satSize);

		    bool converg = false;
		    int iter = 0;
		    int iterTot = 0;
		    updateOldIter = 0;
		    while (!converg)
		    {
		    	iter++;
		    	iterTot++;
		    	if (!(this->diffusion).problem.materialLaw.isLinear()) 
		    	{ // update pressure 
		    		pressure(t);
		    		totalVelocity(t);
		    	}
		    	Transport::update(t, dt, updateVec);
		    	if (this->iterFlag)
		    	{   // only needed if iteration has to be done
		    		*(this->diffusion) *= this->omega;
		    		pressHelp = pressOldIter;
		    		pressHelp *= (1-this->omega);
		    		*(this->diffusion) += pressHelp;
		    		updateHelp = updateVec;
		    		saturation = this->sat;
		    		saturation += (updateHelp *= dt*cFLFactor);
		    		saturation *= this->omega;
		    		satHelp = satOldIter;
		    		satHelp *= (1-this->omega);
		    		saturation += satHelp;
		    		updateDiff = updateVec;
		    		updateDiff -= updateOldIter;
		    		satOldIter = saturation;
		    		pressOldIter = *(this->diffusion);
		    		updateOldIter = updateVec;
		    	}
		    // break criteria for iteration loop
		    	if ( this->iterFlag==2 && dt*updateDiff.two_norm()/(saturation).two_norm() <= this->maxDefect ) 
		    		converg = true;
		    	else if ( this->iterFlag==2 && iter > this->nIter ) 
		    	{
		    		std::cout << "Nonlinear loop in IMPES.update exceeded nIter = " 
		    		<< this->nIter << " iterations." << std::endl;
		    	return 1;
		    	}
		    	else if ( this->iterFlag==1 && iter > this->nIter ) 
		    		converg = true;
		    	else if ( this->iterFlag==0 ) 
		    		converg = true;
		    }
		    // outputs
		    if (this->iterFlag==2) 
		    	std::cout << "Iteration steps: "<< iterTot << std::endl;
		    std::cout.setf (std::ios::scientific, std::ios::floatfield);
  
		    return 0;
		  
	}
			
	virtual void vtkout (const char* name, int k) const 
	{
	  char fname[128];	

	  Transport::vtkout (name,k);
	  
//	  Dune::VTKWriter<G, typename G::template Codim<0>::LevelIndexSet> vtkWriterSat(this->grid, this->grid.levelIndexSet(this->level()));
//	  sprintf(fname,"%s-%05d",name,k);
//	  vtkWriterSat.addCellData(this->sat,"saturation");
//	  vtkWriterSat.write(fname,Dune::VTKOptions::ascii);		

	  Dune::VTKWriter<G, typename G::template Codim<0>::LevelIndexSet> vtkWriterPress(this->grid, this->grid.levelIndexSet(this->diffusion.level()));
	  sprintf(fname,"%s-press.%05d",name,k);
	  vtkWriterPress.addCellData(*(this->diffusion),"total pressure p~");
	  vtkWriterPress.write(fname,Dune::VTKOptions::ascii);		
	}
	
	//! Construct an IMPES object.
	IMPESMS (Diffusion& diff, Transport& trans, int flag = 1, int nIt = 2, double maxDef = 1e-5, 
			double om = 1)
	: IMPES<G, Diffusion, Transport>(diff, trans, flag, nIt, maxDef, om)
	{  
		saturation.resize(trans.grid.size(0));
	}
	  
  private:
	  RepresentationType saturation;
  };
}
#endif
