#ifndef DUNE_IMPES_HH
#define DUNE_IMPES_HH

#include <dune/grid/io/file/vtk/vtkwriter.hh>
#include "dumux/fractionalflow/fractionalflow.hh"

/**
 * @file
 * @brief  IMPES scheme
 * @author Bernd Flemisch
 */

namespace Dune
{
 /**
  * \ingroup fracflow
  * @brief IMplicit Pressure Explicit Saturation (IMPES) scheme for the solution of
  * coupled diffusion/transport problems
  */

  template<class G, class Diffusion, class Transport>
  class IMPES : public FractionalFlow<G, Diffusion, Transport> {
	  typedef typename Diffusion::RepresentationType PressType;
	  typedef typename Diffusion::VelType VelType;
	  typedef typename Diffusion::NumberType RT;

  public:
	typedef typename Transport::RepresentationType RepresentationType;

	virtual void totalVelocity(const RT t=0)
	{
		this->diffusion.totalVelocity(this->problem.velocity, t);
	}

	virtual void initial()
	{
		double t = 0;
		Transport::initial();
		this->pressure(t);
		totalVelocity(t);

	}
	
	
	virtual int update(const RT t, RT& dt, RepresentationType& updateVec, RT cFLFactor = 1) 
	{
		int pressSize = (*(this->diffusion)).size();
		int satSize = this->sat.size();
		PressType pressOldIter(*(this->diffusion));
		PressType pressHelp(pressSize);  
		RepresentationType saturation(this->sat);
		RepresentationType satOldIter(this->sat);
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
            if (!(this->diffusion).problem.materialLaw.isLinear()) { // update pressure 
            	pressure(t);
            	totalVelocity(t);
            }
            Transport::update(t, dt, updateVec);
            if (iterFlag)
            {   // only needed if iteration has to be done
                *(this->diffusion) *= omega;
                pressHelp = pressOldIter;
                pressHelp *= (1-omega);
                *(this->diffusion) += pressHelp;
                updateHelp = updateVec;
                saturation = this->sat;
                saturation += (updateHelp *= dt*cFLFactor);
                saturation *= omega;
                satHelp = satOldIter;
                satHelp *= (1-omega);
                saturation += satHelp;
                updateDiff = updateVec;
                updateDiff -= updateOldIter;
                satOldIter = saturation;
                pressOldIter = *(this->diffusion);
                updateOldIter = updateVec;
            }
            // break criteria for iteration loop
            if ( iterFlag==2 && dt*updateDiff.two_norm()/(saturation).two_norm() <= maxDefect ) 
            	converg = true;
            else if ( iterFlag==2 && iter > nIter ) {
            	std::cout << "Nonlinear loop in IMPES.update exceeded nIter = " 
            			<< nIter << " iterations." << std::endl;
            	return 1;
            }
            else if ( iterFlag==1 && iter > nIter ) 
            	converg = true;
            else if ( iterFlag==0 ) 
            	converg = true;
        }
        // outputs
        if (iterFlag==2) 
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
	IMPES (Diffusion& diff, Transport& trans, int flag = 0, int nIt = 2, double maxDef = 1e-5, 
			double om = 1)
	: FractionalFlow<G, Diffusion, Transport>(diff, trans), 
	  iterFlag(flag), nIter(nIt), maxDefect(maxDef), omega(om)
	{  }
	
  protected: 
	  const int iterFlag;
	  const int nIter;
	  const double maxDefect;
	  const double omega;
	  
  };
}
#endif
