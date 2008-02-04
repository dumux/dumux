#ifndef DUNE_FVTRANSPORT2P2C_HH
#define DUNE_FVTRANSPORT2P2C_HH

#include <dune/common/helpertemplates.hh>
#include <dune/common/typetraits.hh>
#include <dune/grid/common/mcmgmapper.hh>
#include <dune/grid/io/file/vtk/vtkwriter.hh>
#include <dune/grid/utility/intersectiongetter.hh>
#include <dune/istl/bvector.hh>
#include <dumux/transport/transport.hh>
#include <dumux/transport/fv/numericalflux.hh>
#include <dumux/transport/fv/diffusivepart.hh>
#include <dumux/transport/transportproblem2p2c.hh>

namespace Dune
{
  //! \ingroup transport
  //! The finite volume model for the solution of the transport equation
  template<class G, class RT>
  class FVTransport2p2c 
  : public Transport< G, RT, BlockVector< Dune::FieldVector<RT,1> >,
  					   BlockVector< Dune::FieldVector<Dune::FieldVector<RT, G::dimension>, 2*G::dimension> > > 
  {
	  template<int dim>
	  struct ElementLayout
	  {
		  bool contains (Dune::GeometryType gt)
	      {
			  return gt.dim() == dim;
	      }
	  }; 
	  
	  enum{dim = G::dimension};	
	  enum{dimworld = G::dimensionworld};
	  typedef BlockVector< Dune::FieldVector<RT,1> > PressType;
	  typedef BlockVector< FieldVector<FieldVector<RT, G::dimension>, 2*G::dimension> > VelType;
	  typedef typename G::Traits::template Codim<0>::Entity Entity;
	  typedef typename G::Traits::LevelIndexSet IS;
	  typedef typename IS::template Codim<0>::template Partition<All_Partition>::Iterator Iterator;
	  typedef Dune::MultipleCodimMultipleGeomTypeMapper<G,IS,ElementLayout> EM;
	  typedef typename G::template Codim<0>::EntityPointer EntityPointer;
	  typedef typename IntersectionIteratorGetter<G,LevelTag>::IntersectionIterator IntersectionIterator;
	  typedef typename G::ctype ct; 
	  typedef BlockVector< Dune::FieldVector<RT,dim> > SlopeType;

  public:
	  typedef BlockVector< Dune::FieldVector<RT,1> > RepresentationType;
	  /*! 
	   *  \param t time 
	   *  \param dt time step size to estimate 
	   *  \param update vector to be filled with the update 
	   *
	   *  This method calculates the update vector, i.e., the FV discretization 
	   *  of \f$\text{div}\, (f_\text{w} \boldsymbol{v}_t)\f$. The total velocity 
	   *  \f$\boldsymbol{v}_t)\f$ has to be given by the block vector \a velocity, 
	   *  containing values at each cell face. The fractional flow function \f$f_\text{w}\f$ 
	   *  is realized by the numerical flux function \a numericalFlux such that the \f$i\f$-th entry 
	   *  of the vector update is obtained by calculating 
	   *  \f[ \sum_{j \in \mathcal{N}(i)} v_{ij}g(S_i, S_j) - v_{ji}g(S_j, S_i), \f]
	   *  where \f$\mathcal{N}(i)\f$ denotes the set of neighbors of the cell \f$i\f$ and 
	   *  \f$g\f$ stands for \a numericalFlux. The normal velocities \f$v_{ij}\f$ and \f$v_{ji}\f$ 
	   *  are given by 
	   *  \f{align*} v_{ij} = \int_{\Gamma_{ij}} \max(\boldsymbol{v}_\text{t}{\cdot}\boldsymbol{n}_{ij}, \, 0), \qquad
	   *  v_{ji} = \int_{\Gamma_{ij}} \min(\boldsymbol{v}_\text{t}{\cdot}\boldsymbol{n}_{ij}, \, 0), \f}
	   *  where \f$\boldsymbol{n}_{ij}\f$ denotes the unit normal vector from cell \f$i\f$ towards cell \f$j\f$. 
	   *
	   *  Additionally to the \a update vector, the recommended time step size \a dt is calculated 
	   *  employing the usual CFL condition. 
	   */
	   int update(const RT t, RT& dt, RepresentationType& updateVec);
		
	   void initial(); 
	   
	   void initialguess();
	   
	   void flashCalculation(double C1, double p, double temp, double poro, double& sat, double& Cw1, double& Cn1, double& Cw2);
	   
	   void satFlash(double sat, double p, double temp, double poro, double& C1, double& Cw1, double& Cn1, double& Cw2);
	   
	   void postupdate(double t, double dt);
	   
	   RepresentationType& operator*()
	   {
		   return totalConcentration;
	   }
	
	
	   /*! @brief prints the saturation and the concentrations to a VTK file
	    * 
	    *  The file name is "<name>-<k>.vtu" where k is an integer number.
	    *  @param name specifies the name of the VTK file
	    *  @param k specifies a number
	    */	
	   void vtkout (const char* name, int k) const 
	   {
			const typename G::template Codim<0>::LevelIndexSet& iset(this->grid.levelIndexSet(this->level()) );
			Dune::VTKWriter<G, typename G::template Codim<0>::LevelIndexSet> 
				vtkwriter(this->grid, 
						iset );
			char fname[128];	
			sprintf(fname,"%s-%05d",name,k);
			vtkwriter.addCellData(this->sat,"saturation");
			vtkwriter.addCellData(totalConcentration, "total concentration");
			vtkwriter.addCellData(wet_c1, "concentration in wetting phase");
			vtkwriter.addCellData(nonwet_c1, "concentration in non-wetting phase");
			vtkwriter.addCellData(wet_c2, "concentration2 in wetting phase");
			vtkwriter.write(fname,Dune::VTKOptions::ascii);		
	}
	
		/*! @brief constructor
		 * 
		 * @param g a DUNE grid object
		 * @param prob an object of class TransportProblem or derived
		 * @param lev the grid level on which the Transport equation is to be solved.
		 * @param diffPart an object of class DiffusivePart or derived. This determines the diffusive flux incorporated in the transport.
		 * @param rec flag to switch on linear reconstruction (second order TVD)
		 * @param amax alphamax parameter for slope limiter in TVD
		 * @param numFl an object of class Numerical Flux or derived
		 */
		FVTransport2p2c(G& g, TransportProblem2p2c<G, RT, VelType>& prob, int lev = 0,
				   	  DiffusivePart<G, RT>& diffPart = *(new DiffusivePart<G, RT>),
				   	  bool rec = false, double amax = 0.8, const NumericalFlux<RT>& numFl = *(new Upwind<RT>))
		: Transport<G, RT, RepresentationType, VelType>(g, prob,lev), 
		  elementmapper(g, g.levelIndexSet(lev)), indexset(g.levelIndexSet(lev)), reconstruct(rec), 
		  numFlux(numFl), diffusivePart(diffPart), alphamax(amax), problem(prob)
		{
			this->sat.resize(elementmapper.size());
			problem.velocity.resize(elementmapper.size());
			totalConcentration.resize(elementmapper.size());
			wet_c1.resize(elementmapper.size());
			wet_c2.resize(elementmapper.size());
			nonwet_c1.resize(elementmapper.size());
		}

  private:
	  
	  	void CalculateSlopes(SlopeType& slope, RT t);
	  	
	  	bool line1 (FieldVector<RT,dim>& pos)
	  	{
	  		FieldVector<RT,dim> ll(0);
	  		if (pos[0] < ll[0] + 1e-8) return true;
	  		else return false;
	  	}
	  	bool line2 (FieldVector<RT,dim>& pos)
	  	{
	  		FieldVector<RT,dim> ur(300);
	  		if (pos[0] > ur[0] - 1e-8) return true;
	  		else return false;
	  	}
	  
  private:
	  	const IS& indexset;
	  	bool reconstruct;
	  	const NumericalFlux<RT>& numFlux;
	  	const DiffusivePart<G, RT>& diffusivePart;
	  	double alphamax;
  public:
	  EM elementmapper;
	  RepresentationType totalConcentration;
	   RepresentationType wet_c1, nonwet_c1; //wetting phase and nonwetting phase concentrations
	   RepresentationType wet_c2;
		TransportProblem2p2c<G,RT,VelType>& problem;
	    double flux1;
	    double flux2;
  };
  
  
  
  template<class G, class RT>
	int FVTransport2p2c<G,RT>::update(const RT t, RT& dt, RepresentationType& updateVec) 
	{
		// initialize dt very large
		dt = 1E100;
		
		flux1 = 0;
		flux2 = 0;

		// set update vector to zero
		updateVec = 0;

		// printvector(std::cout,saturation,"saturation","row",200,1);
		SlopeType slope(elementmapper.size());
		CalculateSlopes(slope,t);
		// printvector(std::cout,slope,"slope","row",200,1);

		// compute update vector 
		Iterator eendit = indexset.template end<0,All_Partition>();
		for (Iterator it = indexset.template begin<0,All_Partition>(); it != eendit; ++it)
		{
			// cell geometry type
			Dune::GeometryType gt = it->geometry().type();

			// cell center in reference element
			const Dune::FieldVector<ct,dim>& 
				local = Dune::ReferenceElements<ct,dim>::general(gt).position(0,0);

			// cell volume, assume linear map here
			double volume = it->geometry().integrationElement(local)
				*Dune::ReferenceElements<ct,dim>::general(gt).volume();

			// cell index
			int indexi = elementmapper.map(*it);
			// 	  std::cout << "i = " << indexi << ":" << std::endl;

			// for time step calculation
			double sumfactor = 0;
			double sumfactor2 = 0;
			double sumDiff = 0;
			double sumDiff2 = 0;

			// run through all intersections with neighbors and boundary
			IntersectionIterator endit = IntersectionIteratorGetter<G,LevelTag>::end(*it);
			for (IntersectionIterator is = IntersectionIteratorGetter<G,LevelTag>::begin(*it); 
				  is!=endit; ++is)
			{
				// local number of facet 
				int numberInSelf = is.numberInSelf();

				// get geometry type of face
				Dune::GeometryType gtf = is.intersectionSelfLocal().type();
				  
				// center in face's reference element
				const Dune::FieldVector<ct,dim-1>& 
					facelocal = Dune::ReferenceElements<ct,dim-1>::general(gtf).position(0,0);

				// center of face inside volume reference element
				const Dune::FieldVector<ct,dim>& 
					facelocalDim = Dune::ReferenceElements<ct,dim>::general(gtf).position(is.numberInSelf(),1);
					    
				// get normal vector scaled with volume
				Dune::FieldVector<ct,dimworld> integrationOuterNormal 
					= is.integrationOuterNormal(facelocal);
				integrationOuterNormal 
					*= Dune::ReferenceElements<ct,dim-1>::general(gtf).volume();

				// compute factor occuring in flux formula	
				double velocityIJ = std::max(this->problem.vTotal(*it, numberInSelf)*integrationOuterNormal/(volume), 0.0);
				double factor, diffFactor, factorC1;

				// handle interior face
				if (is.neighbor())
				{
					// access neighbor
					EntityPointer outside = is.outside();
					int indexj = elementmapper.map(*outside);
					  
					// compute flux from one side only
					// this should become easier with the new IntersectionIterator functionality!
					if ( it->level()>=outside->level() )
					{
						// compute factor in neighbor
						Dune::GeometryType nbgt = outside->geometry().type();
						const Dune::FieldVector<ct,dim>& 
							nblocal = Dune::ReferenceElements<ct,dim>::general(nbgt).position(0,0);
					      
						double nbvolume = outside->geometry().integrationElement(nblocal)
							*Dune::ReferenceElements<ct,dim>::general(nbgt).volume();
					    double velocityJI = std::max(-(this->problem.vTotal(*it, numberInSelf)*integrationOuterNormal/nbvolume), 0.0);
					      
					    // cell center in global coordinates
					    Dune::FieldVector<ct,dimworld> global = it->geometry().global(local);
					      
					    // neighbor cell center in global coordinates
					    Dune::FieldVector<ct,dimworld> nbglobal = outside->geometry().global(nblocal);
					      
					    // distance vector between barycenters
					    Dune::FieldVector<ct,dimworld> distVec = global - nbglobal;
					      
					    // compute distance between cell centers
					    double dist = distVec.two_norm();
					      
					    // get saturation and concentration value at cell center
					    double satI = this->sat[indexi];
					    double cw1_I  = wet_c1[indexi];
					    double cn1_I = nonwet_c1[indexi];

					    // get saturation and concentrtion value at neighbor cell center
					    double satJ = this->sat[indexj];
					    double cw1_J = wet_c1[indexj];
					    double cn1_J = nonwet_c1[indexj];
						
					    // calculate the cocentration gradients 
					    Dune::FieldVector<ct,dim> cn1Gradient = distVec;		
					    cn1Gradient *= (cn1_J - cn1_I)/(dist*dist);
					    Dune::FieldVector<ct,dim> cw1Gradient = distVec;		
					    cw1Gradient *= (cw1_J - cw1_I)/(dist*dist);
					    
					    // the arithmetic average 
					    double satAvg = 0.5*(satI + satJ);
					    
					    // get the diffusive part
					    double n1_diffPart = this->diffusivePart(*it, numberInSelf, satAvg, cn1Gradient, t)*integrationOuterNormal;
					    double w1_diffPart = this->diffusivePart(*it, numberInSelf, satAvg, cw1Gradient, t)*integrationOuterNormal;

					    // CAREFUL: works only for axisymmetric grids 
					    if (reconstruct) {
					    	for (int k = 0; k < dim; k++) 
					    		if (fabs(distVec[k]) > 0.5*dist) 
					    		{
					    			satI -= fabs(distVec[k])/distVec[k]*0.5*dist*slope[indexi][k];
					    			satJ += fabs(distVec[k])/distVec[k]*0.5*dist*slope[indexj][k];
					    		}
					    }
					    
					    double fwI = this->problem.materialLaw.fractionalW(satI);
					    double fwJ = this->problem.materialLaw.fractionalW(satJ);
					    double fnI = this->problem.materialLaw.fractionalN(1.0-satI);
					    double fnJ = this->problem.materialLaw.fractionalN(1.0-satJ);
					    					    
					    // for timestep control
					    factor = velocityJI -velocityIJ;
					    
					  //  diffFactor =diffPart / volume; TODO include diffusion into timestep control
					    factorC1 = (n1_diffPart + w1_diffPart) / volume
					    		+ velocityJI * cw1_J */*fwJ*/ numFlux(satJ, satI, fwJ, fwI)
					    		- velocityIJ * cw1_I */*fwI*/ numFlux(satI, satJ, fwI, fwJ)
					    		+ velocityJI * cn1_J */*fnJ*/ numFlux(1.0-satJ, 1.0-satI, fnJ, fnI)
					    		- velocityIJ * cn1_I */*fnI;*/ numFlux(1.0-satI, 1.0-satJ, fnI, fnJ); 
					    
					    Dune::FieldVector<ct,dimworld> faceglobal = is.intersectionGlobal().global(facelocal);
					    if (line1(faceglobal)) flux1 += factorC1 * volume;
					    if (line2(faceglobal)) flux2 += factorC1 * volume;
					    
					}
				}
				  
				// handle boundary face
				if (is.boundary()) 
				{
					// cell center in global coordinates
					Dune::FieldVector<ct,dimworld> global = it->geometry().global(local);
				    
					Dune::FieldVector<ct,dim> faceglobal = is.intersectionGlobal().global(facelocal);
					
					// distance vector between barycenters
					Dune::FieldVector<ct,dimworld> distVec = global - faceglobal;
					
					// get saturation value at cell center
					double satI = this->sat[indexi];

					double velocityJI = std::max(-(this->problem.vTotal(*it, numberInSelf)*integrationOuterNormal/volume), 0.0);
					
					//get boundary conditions
					double satBound, c1Bound, cw1Bound, cn1Bound, dummy;
					BoundaryConditions2p2c::Flags bctype = this->problem.cbctype(faceglobal, *it, facelocalDim);
					if (bctype == BoundaryConditions2p2c::saturation)
					{
						satBound = this->problem.g(faceglobal, *it, facelocalDim);
						satFlash(satBound, this->problem.pressBC(faceglobal, *it, facelocalDim), 283.15, this->problem.porosity(), c1Bound, cw1Bound, cn1Bound, dummy);
					}
					if (bctype == BoundaryConditions2p2c::concentration)
					{
						c1Bound = this->problem.g(faceglobal, *it, facelocalDim);
						flashCalculation(c1Bound, this->problem.pressBC(faceglobal, *it, facelocalDim), 283.15, this->problem.porosity(), satBound, cw1Bound, cn1Bound, dummy);
					}
				      
					// compute distance between cell centers
					double dist = distVec.two_norm();
					
				    // calculate the saturation gradient 
				    Dune::FieldVector<ct,dim> satGradient = distVec;		
				    satGradient *= (satBound - satI)/(dist*dist);
				    
				    // the arithmetic average 
				    double satAvg = 0.5*(satI + satBound);
				    
				    // get the diffusive part
				    double diffPart = 0;// diffusivePart(*it, numberInSelf, satAvg, satGradient, t)*integrationOuterNormal;

					// CAREFUL: works only for axisymmetric grids 
					if (reconstruct) {
						for (int k = 0; k < dim; k++) 
							if (fabs(distVec[k]) > 0.5*dist) 
							{
								satI -= fabs(distVec[k])/distVec[k]*dist*slope[indexi][k];
							}
					}
					
			    	double cw1_I = wet_c1[indexi];
			    	double cn1_I = nonwet_c1[indexi];
			    	
			    	double fwI     = this->problem.materialLaw.fractionalW(satI);
			    	double fwBound = this->problem.materialLaw.fractionalW(satBound);
			    	double fnI     = this->problem.materialLaw.fractionalN(1.-satI);
			    	double fnBound = this->problem.materialLaw.fractionalN(1.-satBound);
			    	
			    	factor = velocityJI - velocityIJ;
			    	
					diffFactor = diffPart / volume;
			    	factorC1 = diffFactor
							+ velocityJI * cw1Bound * numFlux(satBound, satI, fwBound, fwI)
							- velocityIJ * cw1_I * numFlux(satI, satBound, fwI, fwBound)
							+ velocityJI * cn1Bound * numFlux(1.0-satBound, 1.0-satI, fnBound, fnI)
							- velocityIJ * cn1_I * numFlux(1.0-satI, 1.0-satBound, fnI, fnBound); 
			    	
				    if (line1(faceglobal)) flux1 += factorC1*volume;
				    if (line2(faceglobal)) flux2 += factorC1*volume;
				}

				// add to update vector 
				updateVec[indexi] += factorC1;

				// for time step calculation
				if (factor>=0) 
					sumfactor += factor;
				else 
					sumfactor2 += (-factor);
				if (diffFactor>=0)
					sumDiff += diffFactor;
				else
					sumDiff += (-diffFactor);
			} // end all intersections    
			// compute dt restriction
			sumfactor = std::max(sumfactor,sumfactor2);
			sumDiff = std::max(sumDiff,sumDiff2);
			sumfactor = std::max(sumfactor,100*sumDiff);
			dt = std::min(dt,1.0/sumfactor);     
			
		} // end grid traversal                 
		dt = dt*this->problem.porosity();
		updateVec /= this->problem.porosity();
		
		return 0;
	}
  
  template<class G, class RT>
  void FVTransport2p2c<G,RT>::postupdate(double t, double dt)
  {    
	  	// iterate through leaf grid an evaluate c0 at cell center
	    Iterator eendit = indexset.template end<0,All_Partition>();
	    for (Iterator it = indexset.template begin<0,All_Partition>(); it != eendit; ++it)
		{
	    	// get geometry information of cell
  			Dune::GeometryType gt = it->geometry().type();
  			const Dune::FieldVector<ct,dim>& 
  				local = Dune::ReferenceElements<ct,dim>::general(gt).position(0,0);
  			Dune::FieldVector<ct,dimworld> global = it->geometry().global(local);
  			
	    	int indexi = elementmapper.map(*it);
			flashCalculation(totalConcentration[indexi], this->problem.press(global, *it, local), 283.15, this->problem.porosity(), this->sat[indexi][0], wet_c1[indexi][0], nonwet_c1[indexi][0], wet_c2[indexi][0]);
		}
  }
  
  template<class G, class RT>
  void FVTransport2p2c<G,RT>::initial() 
  	{ 
  		// iterate through leaf grid an evaluate c0 at cell center
  	    Iterator eendit = indexset.template end<0,All_Partition>();
  	    for (Iterator it = indexset.template begin<0,All_Partition>(); it != eendit; ++it)
  		{
  	    	int indexi = elementmapper.map(*it);
  	    	
  			// get geometry information of cell
  			Dune::GeometryType gt = it->geometry().type();
  			const Dune::FieldVector<ct,dim>& 
  				local = Dune::ReferenceElements<ct,dim>::general(gt).position(0,0);
  			Dune::FieldVector<ct,dimworld> global = it->geometry().global(local);
  			
  			double sat_0, C1_0;
  			// get type of initial condition
  			Dune::BoundaryConditions2p2c::Flags ictype = this->problem.ictype(global, *it, local);
  			// saturation initial condition
  			if (ictype == Dune::BoundaryConditions2p2c::saturation)
  			{
  				sat_0 = this->problem.S0(global, *it, local);
  				satFlash(sat_0, this->problem.press(global, *it, local), 283.15, this->problem.porosity(), C1_0, wet_c1[indexi][0], nonwet_c1[indexi][0], wet_c2[indexi][0]);
  			}
  			// concentration initial condition
  			if (ictype == Dune::BoundaryConditions2p2c::concentration)
  			{
  				C1_0 = this->problem.C1_0(global, *it, local);
  				flashCalculation(C1_0, this->problem.press(global, *it, local), 283.15, this->problem.porosity(), sat_0, wet_c1[indexi][0], nonwet_c1[indexi][0], wet_c2[indexi][0]);
  			}
  			// initialize cell concentration
  			totalConcentration[indexi] = C1_0;
  			this->sat[indexi][0] = sat_0;
  		}

  		return;
  	}
  
  template<class G, class RT>
  void FVTransport2p2c<G,RT>::initialguess()
  {
	// iterate through leaf grid an evaluate c0 at cell center
	Iterator eendit = indexset.template end<0,All_Partition>();
	for (Iterator it = indexset.template begin<0,All_Partition>(); it != eendit; ++it)
	{
		int indexi = elementmapper.map(*it);
		
		// get geometry information of cell
		GeometryType gt = it->geometry().type();
		const FieldVector<ct,dim>& 
			local = ReferenceElements<ct,dim>::general(gt).position(0,0);
		FieldVector<ct,dimworld> global = it->geometry().global(local);
		
		double sat_0, C1_0;
		// get type of initial condition
		BoundaryConditions2p2c::Flags ictype = this->problem.ictype(global, *it, local);
		// saturation initial condition
		if (ictype == BoundaryConditions2p2c::saturation)
			sat_0 = this->problem.S0(global, *it, local);
		// concentration initial condition
		if (ictype == BoundaryConditions2p2c::concentration)
		{
			C1_0 = this->problem.C1_0(global, *it, local);
			sat_0 = C1_0 / this->problem.materialLaw.wettingPhase.density(283.15, 1e6);
		}
		// initialize cell concentration
		totalConcentration[indexi] = C1_0;
		this->sat[indexi][0] = sat_0;
	}	
	return;
  }
  
  template<class G, class RT>
  void FVTransport2p2c<G,RT>::CalculateSlopes(SlopeType& slope, RT t)
  	  {
  	    Iterator endit = this->grid.template lend<0>(this->level());   
  	    for (Iterator it = this->grid.template lbegin<0>(this->level()); it!=endit; ++it)
  	  	{
  	  	  // get some cell properties
  	  	  Dune::GeometryType gt = it->geometry().type();
  	  	  const Dune::FieldVector<ct,dim>& 
  	  		local = Dune::ReferenceElements<ct,dim>::general(gt).position(0,0);
  	  	  Dune::FieldVector<ct,dimworld>  global = it->geometry().global(local);
  	  	  int indexi = elementmapper.map(*it);

  	  	  // vector containing the distances to the neighboring cells
  	  	  Dune::FieldVector<double, 2*dim> dist;

  	  	  // vector containing the saturations of the neighboring cells
  	  	  Dune::FieldVector<double, 2*dim> saturation;

  	  	  // location[k], k = 0,...,2dim-1, contains the local index w.r.t. IntersectionIterator, 
  	  	  // i.e. the numberInSelf, which is east, west, north, south, top, bottom to the cell center 
  	  	  Dune::FieldVector<int, 2*dim> location;

  	  	  // run through all intersections with neighbors and boundary
  	  	  IntersectionIterator isend = IntersectionIteratorGetter<G,LevelTag>::end(*it);
  	  	  for (IntersectionIterator is = IntersectionIteratorGetter<G,LevelTag>::begin(*it); is!=isend; ++is)
  	  		{
  	  		  // local number of facet 
  	  		  int numberInSelf = is.numberInSelf();

  	  		  // handle interior face
  	  		  if (is.neighbor())
  	  		    {
  	  		      // access neighbor
  	  		      EntityPointer outside = is.outside();
  	  		      int indexj = elementmapper.map(*outside);
  	  		      
  	  		      // get saturation value 
  	  		      saturation[numberInSelf] = this->sat[indexj];

  	  		      // compute factor in neighbor
  	  		      Dune::GeometryType nbgt = outside->geometry().type();
  	  		      const Dune::FieldVector<ct,dim>& 
  	  			  nblocal = Dune::ReferenceElements<ct,dim>::general(nbgt).position(0,0);
  	  		      
  	  		      // neighboring cell center in global coordinates
  	  		      Dune::FieldVector<ct,dimworld> 
  	  			  nbglobal = outside->geometry().global(nblocal);

  	  		      // distance vector between barycenters
  	  		      Dune::FieldVector<ct,dimworld> distVec = global - nbglobal;
  	  		      
  	  		      // compute distance between cell centers
  	  		      dist[numberInSelf] = distVec.two_norm();

  	  		      // CAREFUL: works only for axiparallel grids 
  	  		      for (int k = 0; k < dim; k++) 
  	  			  if (nbglobal[k] - global[k] > 0.5*dist[numberInSelf]) 
  	  			  {
  	  			    location[2*k] = numberInSelf;
  	  			  }
  	  			  else if (nbglobal[k] - global[k] < -0.5*dist[numberInSelf]) 
  	  			  {
  	  			    location[2*k + 1] = numberInSelf;
  	  			  }
  	  		    }
  	  		  
  	  		  // handle boundary face
  	  		  if (is.boundary()) 
  	  		  {
  	  		      // get geometry type of face
  	  		      Dune::GeometryType gtf = is.intersectionSelfLocal().type();
  	  		      
  	  		      // center in face's reference element
  	  		      const Dune::FieldVector<ct,dim-1>& 
  	  			  facelocal = Dune::ReferenceElements<ct,dim-1>::general(gtf).position(0,0);

  	  			  // center of face inside volume reference element
  	  			  const Dune::FieldVector<ct,dim>& 
  	  				facelocalDim = Dune::ReferenceElements<ct,dim>::general(gtf).position(is.numberInSelf(),1);
  	  			  
  	  		      // center of face in global coordinates
  	  		      Dune::FieldVector<ct,dimworld> 
  	  			  faceglobal = is.intersectionGlobal().global(facelocal);
  	  			
  	  		      // get saturation value 
  	  		      bool dir;
  	  		      
  	  		      saturation[numberInSelf] = this->problem.g(faceglobal, *it, facelocalDim);

  	  		      // distance vector between barycenters
  	  		      Dune::FieldVector<ct,dimworld> distVec = global - faceglobal;
  	  		      
  	  		      // compute distance 
  	  		      dist[numberInSelf] = distVec.two_norm();
  	  			
  	  		      // CAREFUL: works only for axiparallel grids 
  	  		      for (int k = 0; k < dim; k++) 
  	  			  if (faceglobal[k] - global[k] > 0.5*dist[numberInSelf]) 
  	  			  {
  	  			    location[2*k] = numberInSelf;
  	  			  }
  	  			  else if (faceglobal[k] - global[k] < -0.5*dist[numberInSelf]) 
  	  			  {
  	  			    location[2*k + 1] = numberInSelf;
  	  			  }
  	  		    }
  	  		} // end all intersections 

  	  	    for (int k = 0; k < dim; k++) {
  	  	    double slopeIK = (saturation[location[2*k]] -  saturation[location[2*k + 1]])
  	  	                      /(dist[location[2*k]] + dist[location[2*k + 1]]);

  	  	    double alphaIK = 1.0;
  	  	    if (fabs(slopeIK) > 1e-8*dist[location[2*k]]) 
  	  	    {
  	  	      double satI = this->sat[indexi];
  	  	      double alphaRight = fabs(2.0/(dist[location[2*k]]*slopeIK)*(saturation[location[2*k]] - satI));
  	  	      double alphaLeft = fabs(2.0/(dist[location[2*k + 1]]*slopeIK)*(satI - saturation[location[2*k + 1]]));
  	  	      alphaIK = std::min(std::min(std::max(alphaRight, 0.0), std::max(alphaLeft, 0.0)), alphamax);
  	  	    }
  	  	    
  	  	    double gabagabahey = alphaIK*slopeIK;
  	  	    slope[indexi][k] = alphaIK*slopeIK;
  	  	    
  	  	  }
  	  	} // end grid traversal 
  	    
  	    return;
  	  }
  
  template<class G, class RT>
  void FVTransport2p2c<G,RT>::flashCalculation(double C1, double p, double temp, double poro, double& sat, double& Cw1, double& Cn1, double& Cw2)
  {
	  double K1 = this->problem.materialLaw.wettingPhase.vaporPressure(temp) / p;
      double K2 = 1. / (p * this->problem.henry(temp));
      double xw1 = (1. - K2) / (K1 -K2);
      double xn1 = xw1 * K1;
      double Xw1 = xw1 * this->problem.materialLaw.wettingPhase.molarMass() 
                  / ( xw1 * this->problem.materialLaw.wettingPhase.molarMass() + (1.-xw1) * this->problem.materialLaw.nonwettingPhase.molarMass() );
      double Xn1 = xn1 * this->problem.materialLaw.nonwettingPhase.molarMass() 
                  / ( xn1 * this->problem.materialLaw.nonwettingPhase.molarMass() + (1.-xn1) * this->problem.materialLaw.wettingPhase.molarMass() );
      K1 = Xn1 / Xw1;
      K2 = (1.-Xn1) / (1.-Xw1);
      double rho_w = this->problem.materialLaw.wettingPhase.density(temp, p);
      double rho_n = this->problem.materialLaw.nonwettingPhase.density(temp, p);
      sat = (C1/poro * (K1 - K2) - K1 * rho_n * (1. - K2)) / ((K2 - 1.) * (K1 * rho_n - rho_w));
      
      if (sat < 0.) 
      {
        sat = 0.;
        Cw1 = 0.;
        Cn1 = C1 ;
		  Cw2 = 0;
      }
      else if (sat > 1.)
      {
        sat = 1.;
        Cw1 = C1 ;
        Cn1 = 0.;
		  Cw2 = rho_w - C1;
      }
      else
      {
        Cw1 = rho_w * C1 
        	/ (rho_w * sat + rho_n * (1.-sat) * K1) / poro;
        Cn1 = rho_n * K1 * C1 
            / (rho_w * sat + rho_n * (1.-sat) * K1) / poro;
		  Cw2 = rho_w - Cw1; 
      }  
  }
  
  template <class G, class RT>
  void FVTransport2p2c<G,RT>::satFlash(double sat, double p, double temp, double poro, double& C1, double& Cw1, double& Cn1, double& Cw2)
  {
	  if (sat <= 0 || sat >= 1)
		  DUNE_THROW(RangeError,
		  	"FVTransport2p2c :: saturation initial and boundary conditions may not equal zero or one!");
	  double K1 = this->problem.materialLaw.wettingPhase.vaporPressure(temp) / p;
      double K2 = 1. / (p * this->problem.henry(temp));
      double xw1 = (1. - K2) / (K1 -K2);
      double xn1 = xw1 * K1;
      double Xw1 = xw1 * this->problem.materialLaw.wettingPhase.molarMass() 
                  / ( xw1 * this->problem.materialLaw.wettingPhase.molarMass() + (1.-xw1) * this->problem.materialLaw.nonwettingPhase.molarMass() );
      double Xn1 = xn1 * this->problem.materialLaw.nonwettingPhase.molarMass() 
                  / ( xn1 * this->problem.materialLaw.nonwettingPhase.molarMass() + (1.-xn1) * this->problem.materialLaw.wettingPhase.molarMass() );
      K1 = Xn1 / Xw1;
      K2 = (1.-Xn1) / (1.-Xw1);
      double rho_w = this->problem.materialLaw.wettingPhase.density(temp, p);
      double rho_n = this->problem.materialLaw.nonwettingPhase.density(temp, p);
      
      Cw1 = rho_w * Xw1;
      Cn1 = rho_n * Xn1;
      Cw2 = rho_w - Cw1;
      C1  = poro* (sat * Cw1 + (1-sat) * Cn1);
  }

}
#endif
