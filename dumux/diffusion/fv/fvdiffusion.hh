#ifndef DUNE_FVDIFFUSION_HH
#define DUNE_FVDIFFUSION_HH

#include <dune/common/helpertemplates.hh>
#include <dune/common/typetraits.hh>
#include <dune/grid/common/mcmgmapper.hh>
#include <dune/grid/io/file/vtk/vtkwriter.hh>
#include <dune/grid/utility/intersectiongetter.hh>
#include <dune/istl/operators.hh>
#include <dune/istl/solvers.hh>
#include <dune/istl/preconditioners.hh>
#include "dumux/diffusion/diffusion.hh"
#include "dumux/pardiso/pardiso.hh"
#include "dumux/diffusion/problems/uniformproblem.hh"
#include "dumux/transport/problems/simpleproblem.hh"

/**
 * @file
 * @brief  Finite Volume Diffusion Model
 * @author Bernd Flemisch, Jochen Fritz
 */

namespace Dune
{
  //! \ingroup diffusion
  //! Finite Volume Diffusion Model
  /*! Provides a Finite Volume implementation for the evaluation 
   * of equations of the form 
   * \f$ - \text{div}\, (\lambda K \text{grad}\, p ) = q, \f$, 
   * \f$p = g\f$ on \f$\Gamma_1\f$, and 
   * \f$\lambda K \text{grad}\, p \cdot \mathbf{n} = J\f$ 
   * on \f$\Gamma_2\f$. Here, 
   * \f$p\f$ denotes the pressure, \f$K\f$ the absolute permeability, 
   * and \f$\lambda\f$ the total mobility, possibly depending on the 
   * saturation, \f$q\f$ the source term.
	Template parameters are:

	- G         a DUNE grid type
	- RT        type used for return values 
   */
  template<class G, class RT>
  class FVDiffusion 
  : public Diffusion< G, RT, BlockVector< FieldVector<RT,1> >,
  					   BlockVector< FieldVector<FieldVector<RT, G::dimension>, 2*G::dimension> > > 
  {
	  template<int dim>
	  struct ElementLayout
	  {
		  bool contains (GeometryType gt)
	      {
			  return gt.dim() == dim;
	      }
	  }; 
	  
	  enum{dim = G::dimension};	
	  enum{dimworld = G::dimensionworld};
	  
	  typedef typename G::Traits::template Codim<0>::Entity Entity;
	  typedef typename G::Traits::LevelIndexSet IS;
	  typedef typename IS::template Codim<0>::template Partition<All_Partition>::Iterator Iterator;
	  typedef typename G::template Codim<0>::HierarchicIterator HierarchicIterator;
	  typedef MultipleCodimMultipleGeomTypeMapper<G,IS,ElementLayout> EM;
	  typedef typename G::template Codim<0>::EntityPointer EntityPointer;
	  typedef typename IntersectionIteratorGetter<G,LevelTag>::IntersectionIterator IntersectionIterator;
	  typedef typename G::ctype ct; 
	  typedef FieldMatrix<double,1,1> MB;
	  typedef BCRSMatrix<MB> MatrixType;
	  typedef FieldVector<double, 1> VB;
	  typedef BlockVector<VB> Vector;
	  typedef FieldVector<double,dim> R2;
	  typedef BlockVector<R2> R3;
	  typedef BlockVector<R3> LocVelType;
	  
  public:
	typedef BlockVector< FieldVector<FieldVector<RT, G::dimension>, 2*G::dimension> > VelType;
	typedef BlockVector< FieldVector<RT,1> > RepresentationType;

	void assemble(const RT t); 

	void solve(); 

	void pressure(const RT t=0)
	{
		assemble(t);
		solve();
		return;
	}

	void totalVelocity(VelType& velocity, const RT t) const;
	
    void totalVelocity(VelType& velocity, const RT t, const int level) const;
		
	void vtkout (const char* name, int k) const 
	{
		int level(this->level());
		VTKWriter<G, typename G::template Codim<0>::LevelIndexSet> 
			vtkwriter(this->grid, this->grid.levelIndexSet( level ));
		char fname[128];	
		sprintf(fname,"%s-%05d",name,k);
		vtkwriter.addCellData(this->press,"total pressure p~");
		vtkwriter.write(fname,VTKOptions::ascii);		
	}
	
	void initializeMatrix();
	
	FVDiffusion(G& g, 
			    DiffusionProblem<G, RT>& prob, 
				TransportProblem<G, RT, VelType>& satprob = *(new SimpleProblem<G, RT>), 
				int lev = -1)
	            : Diffusion<G, RT, RepresentationType, VelType>(g, prob, lev == -1 ? g.maxLevel() : lev), 
	              satProblem(satprob), 
	              elementmapper(g, g.levelIndexSet(this->level())), 
	              indexset(g.levelIndexSet(this->level())), 
	              A(g.size(this->level(), 0), g.size(this->level(), 0), (2*dim+1)*g.size(this->level(), 0), BCRSMatrix<MB>::random), 
	              f(g.size(this->level(), 0)) , solverName_("BiCGSTAB"), preconditionerName_("SeqILU0")
	{
		this->press.resize(g.size(this->level(), 0));
		this->press = 0;
		initializeMatrix();
	}
	
	FVDiffusion(G& g, 
			    DiffusionProblem<G, RT>& prob, 
			    std::string solverName, 
			    std::string preconditionerName, 
				TransportProblem<G, RT, VelType>& satprob = *(new SimpleProblem<G, RT>), 
				int lev = -1)
	            : Diffusion<G, RT, RepresentationType, VelType>(g, prob, lev == -1 ? g.maxLevel() : lev), 
	              satProblem(satprob), 
	              elementmapper(g, g.levelIndexSet(this->level())), 
	              indexset(g.levelIndexSet(this->level())), 
	              A(g.size(this->level(), 0), g.size(this->level(), 0), (2*dim+1)*g.size(this->level(), 0), BCRSMatrix<MB>::random), 
	              f(g.size(this->level(), 0)), solverName_(solverName), preconditionerName_(preconditionerName)
	{
		this->press.resize(g.size(this->level(), 0));
		this->press = 0;
		initializeMatrix();
	}
	
  private:
	  TransportProblem<G, RT, VelType>& satProblem; //!< problem data
	  EM elementmapper;
	  const IS& indexset;
	  MatrixType A;
	  RepresentationType f;
	  std::string solverName_;
	  std::string preconditionerName_;
  };

  
  
  template<class G, class RT>
  void FVDiffusion<G, RT>::initializeMatrix()
  {
	    // determine matrix row sizes 
	    Iterator eendit = indexset.template end<0,All_Partition>();
	    for (Iterator it = indexset.template begin<0,All_Partition>(); it != eendit; ++it)
	      {
			// cell index
			int indexi = elementmapper.map(*it);
	
			// initialize row size
			int rowSize = 1;
	
			// run through all intersections with neighbors 
			IntersectionIterator endit = IntersectionIteratorGetter<G,LevelTag>::end(*it);
			for (IntersectionIterator is = IntersectionIteratorGetter<G,LevelTag>::begin(*it); 
				  is!=endit; ++is)
			    if (is.neighbor()) 
			      rowSize++;
			A.setrowsize(indexi, rowSize);
	      }
	    A.endrowsizes();

	    // determine position of matrix entries 
	    for (Iterator it = indexset.template begin<0,All_Partition>(); it != eendit; ++it)
	      {
			// cell index
			int indexi = elementmapper.map(*it);
	
			// add diagonal index
			A.addindex(indexi, indexi);
	
			// run through all intersections with neighbors 
			IntersectionIterator endit = IntersectionIteratorGetter<G,LevelTag>::end(*it);
			for (IntersectionIterator is = IntersectionIteratorGetter<G,LevelTag>::begin(*it); 
			  	  is!=endit; ++is)
			    if (is.neighbor()) 
			      {
					// access neighbor
					EntityPointer outside = is.outside();
					int indexj = elementmapper.map(*outside);
		
					// add off diagonal index
					A.addindex(indexi, indexj);
			      }
	      }
	    A.endindices();		

	    return;
  }
  
  template<class G, class RT>
  void FVDiffusion<G, RT>::assemble(const RT t=0)
	{
	    // initialization: set matrix A to zero	   
        A = 0;

        // find out whether gravity effects are relevant
        bool hasGravity = false;
        const FieldVector<ct,dim>& gravity = this->problem.gravity();
        for (int k = 0; k < dim; k++)
        	if (gravity[k] != 0) 
        		hasGravity = true;
        
        Iterator eendit = indexset.template end<0,All_Partition>();
	    for (Iterator it = indexset.template begin<0,All_Partition>(); it != eendit; ++it)
	      {	
			// cell geometry type
			GeometryType gt = it->geometry().type();
			
			// cell center in reference element
			const FieldVector<ct,dim>& 
			  local = ReferenceElements<ct,dim>::general(gt).position(0,0);
			
			// get global coordinate of cell center
			FieldVector<ct,dim> global = it->geometry().global(local);
			
			// cell index
			int indexi = elementmapper.map(*it);
			

			// cell volume, assume linear map here
			double volume = it->geometry().integrationElement(local)
				*ReferenceElements<ct,dim>::general(gt).volume();

			// set right side to zero
			f[indexi] = volume*this->problem.q(global,*it,local);
	
			// get absolute permeability 
			FieldMatrix<ct,dim,dim> Ki(this->problem.K(global,*it,local));
	
			//compute total mobility
			double lambdaI, fractionalWI;
			double sati = this->problem.sat(global, *it, local);
			lambdaI = this->problem.materialLaw.mobTotal(sati);
			if (hasGravity) 
				fractionalWI = this->problem.materialLaw.fractionalW(sati);
			

			IntersectionIterator endit = IntersectionIteratorGetter<G,LevelTag>::end(*it);
			for (IntersectionIterator is = IntersectionIteratorGetter<G,LevelTag>::begin(*it); 
				  is!=endit; ++is)
			  {
		    
			    // get geometry type of face
			    GeometryType gtf = is.intersectionSelfLocal().type();
			    
			    // center in face's reference element
			    const FieldVector<ct,dim-1>& 
			      facelocal = ReferenceElements<ct,dim-1>::general(gtf).position(0,0);
			    
				// center of face inside volume reference element
			    const FieldVector<ct,dim>& 
			      facelocalDim = ReferenceElements<ct,dim>::general(gtf).position(is.numberInSelf(),1);
				    
			    // get normal vector 
			    FieldVector<ct,dimworld> unitOuterNormal 
			      = is.unitOuterNormal(facelocal);
			    
			    // get normal vector scaled with volume
			    FieldVector<ct,dimworld> integrationOuterNormal 
			      = is.integrationOuterNormal(facelocal);
			    integrationOuterNormal 
			      *= ReferenceElements<ct,dim-1>::general(gtf).volume();
	
			    // get face volume 
			    // double faceVol = ReferenceElements<ct,dim-1>::general(gtf).volume();
			    double faceVol = is.intersectionGlobal().volume();
			    
				// compute directed permeability vector Ki.n
				FieldVector<ct,dim> Kni(0);
				Ki.umv(unitOuterNormal, Kni);

				// handle interior face
			    if (is.neighbor()) 
			    {
					// access neighbor
					EntityPointer outside = is.outside();
					int indexj = elementmapper.map(*outside);
					
					// compute factor in neighbor
					GeometryType nbgt = outside->geometry().type();
					const FieldVector<ct,dim>& 
					  nblocal = ReferenceElements<ct,dim>::general(nbgt).position(0,0);
		
					// neighbor cell center in global coordinates
					FieldVector<ct,dimworld> 
					  nbglobal = outside->geometry().global(nblocal);
		
					// distance vector between barycenters
					FieldVector<ct,dimworld> 
					  distVec = global - nbglobal;
		
					// compute distance between cell centers
					double dist = distVec.two_norm();
		
					// get absolute permeability 
					FieldMatrix<ct,dim,dim> Kj(this->problem.K(nbglobal, *outside, nblocal));
						
					// compute vectorized permeabilities
	                FieldVector<ct,dim> Knj(0);
	                Kj.umv(unitOuterNormal, Knj);
	                // compute permeability normal to intersection and take harmonic mean
//					FieldVector<ct,dim> K(0);
//					for (int l = 0; l < dim; l++) {
//						double factor = (Kni[l] + Knj[l]);
//						if (factor)
//							K[l] = 2*Kni[l]*Knj[l]/factor;
//					}
 	                double K_n_i = Kni * unitOuterNormal;
 	                double K_n_j = Knj * unitOuterNormal;
 	                double Kn    = 2 * K_n_i * K_n_j / (K_n_i + K_n_j);
 	                // compute permeability tangential to intersection and take arithmetic mean
 	                FieldVector<ct,dim> uON = unitOuterNormal;
 	                FieldVector<ct,dim> K_t_i = Kni - (uON *= K_n_i);
 	                uON = unitOuterNormal;
 	                FieldVector<ct,dim> K_t_j = Knj - (uON *= K_n_j);
 	                FieldVector<ct,dim> Kt = (K_t_i += K_t_j);
 	                Kt *= 0.5;
 	                // Build vectorized averaged permeability
 	                uON = unitOuterNormal;
 	                FieldVector<ct,dim> K = (Kt += (uON *=Kn));
	
					//compute total mobility
					double lambdaJ, fractionalWJ;
					double satj = this->problem.sat(nbglobal, *outside, nblocal);
					lambdaJ = this->problem.materialLaw.mobTotal(satj);
					if (hasGravity) 
						fractionalWJ = this->problem.materialLaw.fractionalW(satj);
					
					// compute averaged total mobility
					// CAREFUL: Harmonic weightig can generate zero matrix entries, 
			                // use arithmetic weighting instead: 
					double fractionalW;
					double lambda = 0.5*(lambdaI + lambdaJ);
					if (hasGravity) 
						fractionalW = 0.5*(fractionalWI + fractionalWJ);
						
					// update diagonal entry 
					double entry = fabs(lambda*faceVol*(K*distVec)/(dist*dist));
					A[indexi][indexi] += entry;
		
					// set off-diagonal entry 
					A[indexi][indexj] = -entry;
		
					if (hasGravity) {
						double factor = fractionalW*(this->problem.materialLaw.wettingPhase.density()) 
							+ (1 - fractionalW)*(this->problem.materialLaw.nonwettingPhase.density());
						f[indexi] += factor*lambda*faceVol*(K*gravity);
					}
					
					if (this->problem.capillary) {
					  // calculate saturation gradient
					  FieldVector<ct,dim> satGradient = distVec;		
					  satGradient *= (satj - sati)/(dist*dist);
		
					  // arithmetic average of the permeability
					  K = ((Kni + Knj) *= 0.5);
					  
					  // capillary pressure w.r.t. saturation 
					  double pCI = this->problem.materialLaw.pC(sati);
					  double pCJ = this->problem.materialLaw.pC(satj);
		
					  // mobility of the nonwetting phase 
					  double lambdaN = 0.5*(this->problem.materialLaw.mobN(1 - sati) 
							  				+ this->problem.materialLaw.mobN(1 - satj));
		
		
					  // calculate capillary pressure gradient
					  FieldVector<ct,dim> pCGradient = distVec;		
					  pCGradient *= -(pCJ - pCI)/(dist*dist);
		
					  f[indexi] += lambdaN*faceVol*(K*pCGradient);
					}
			    }
			    // boundary face 
			    else 
			      { 
					// center of face in global coordinates
					FieldVector<ct,dimworld> 
					  faceglobal = is.intersectionGlobal().global(facelocal);
					  
					// compute total mobility
					double fractionalW = 1.;
					double lambda = lambdaI;
						if (hasGravity) fractionalW = fractionalWI;
		
					//get boundary condition for boundary face center
					BoundaryConditions::Flags bctype = this->problem.bctype(faceglobal, *it, facelocalDim);
					if (bctype == BoundaryConditions::dirichlet) 
					{ 
						FieldVector<ct,dimworld> distVec(global - faceglobal);
						double dist = distVec.two_norm();
						A[indexi][indexi] -= lambda*faceVol*(Kni*distVec)/(dist*dist);
						double g = this->problem.g(faceglobal, *it, facelocalDim);
						f[indexi] -= lambda*faceVol*g*(Kni*distVec)/(dist*dist);
		
						if (hasGravity) {
							double factor = fractionalW*(this->problem.materialLaw.wettingPhase.density()) 
								+ (1 - fractionalW)*(this->problem.materialLaw.nonwettingPhase.density());
							f[indexi] += factor*lambda*faceVol*(Kni*gravity);
						}

					} 
					else
					{
						double J = this->problem.J(faceglobal, *it, facelocalDim);
						f[indexi] += faceVol*J;
					}
		
					if (this->problem.capillary) {
					  double satj = satProblem.g(faceglobal, *it, facelocalDim);
					  
					  // distance vector between barycenters
					  FieldVector<ct,dimworld> 
					    distVec = global - faceglobal;
					  
					  // compute distance between cell centers
					  double dist = distVec.two_norm();
					  
					  // calculate saturation gradient
					  FieldVector<ct,dim> satGradient = distVec;		
					  satGradient *= (satj - sati)/(dist*dist);
					  
					  // capillary pressure w.r.t. saturation 
					  double pCI = this->problem.materialLaw.pC(sati);
					  double pCJ = this->problem.materialLaw.pC(satj);
					  
					  // mobility of the nonwetting phase 
					  double lambdaN = 0.5*(this->problem.materialLaw.mobN(1 - sati)
							  				+ this->problem.materialLaw.mobN(1 - satj));
					  
					  
					  // calculate capillary pressure gradient
					  FieldVector<ct,dim> pCGradient = distVec;		
					  pCGradient *= -(pCJ - pCI)/(dist*dist);
					  
					  f[indexi] += lambdaN*faceVol*(Kni*pCGradient);
					}
			      }
			  } // end all intersections         
			
	      } // end grid traversal 
	    return;
	}
	
	
  template<class G, class RT>
  void FVDiffusion<G, RT>::solve()
  {
	  MatrixAdapter<MatrixType,Vector,Vector> op(A); 
	  InverseOperatorResult r;
	  
	  if (preconditionerName_ == "SeqILU0") {
	      SeqILU0<MatrixType,Vector,Vector> preconditioner(A, 1.0);
	      if (solverName_ == "CG") {
	    	  CGSolver<Vector> solver(op, preconditioner, 1E-14, 10000, 1);
	    	  solver.apply(this->press, f, r);
	      }
	      else if (solverName_ == "BiCGSTAB") {
	    	  BiCGSTABSolver<Vector> solver(op, preconditioner, 1E-14, 10000, 1);
	    	  solver.apply(this->press, f, r);
	      }
	      else 
			  DUNE_THROW(NotImplemented, "FVDiffusion :: solve : combination " << preconditionerName_ 
					  << " and " << solverName_ << ".");
	  }
	  else if (preconditionerName_ == "SeqPardiso") {
	      SeqPardiso<MatrixType,Vector,Vector> preconditioner(A);
	      if (solverName_ == "Loop") {
	    	  LoopSolver<Vector> solver(op, preconditioner, 1E-14, 10000, 1);
	    	  solver.apply(this->press, f, r);
	      }
	      else 
	    	  DUNE_THROW(NotImplemented, "FVDiffusion :: solve : combination " << preconditionerName_ 
	    			  << " and " << solverName_ << ".");
	  }
	  else 
		  DUNE_THROW(NotImplemented, "FVDiffusion :: solve : preconditioner " << preconditionerName_ << ".");

	  return;
  }
	
  template<class G, class RT>
  void FVDiffusion<G, RT>::totalVelocity(VelType& velocity, const RT t=0) const
	{
      // find out whether gravity effects are relevant
      bool hasGravity = false;
      const FieldVector<ct,dim>& gravity = this->problem.gravity();
      for (int k = 0; k < dim; k++)
      	if (gravity[k] != 0) 
      		hasGravity = true;
      
	    Iterator eendit = indexset.template end<0,All_Partition>();
	    for (Iterator it = indexset.template begin<0,All_Partition>(); it != eendit; ++it)
	    {
	    	// cell geometry type
	    	GeometryType gt = it->geometry().type();
		      
	    	// cell center in reference element
	    	const FieldVector<ct,dim>& 
				local = ReferenceElements<ct,dim>::general(gt).position(0,0);
		      
			// cell center in global coordinates
			FieldVector<ct,dimworld> 
				global = it->geometry().global(local);
		      
		    // cell index
			int indexi = elementmapper.map(*it);
		    
			// get pressure and permeability in element
			double pressi = this->press[indexi];
		      
			// get absolute permeability 
			FieldMatrix<ct,dim,dim> Ki(this->problem.K(global,*it,local));
	
			//compute total mobility
			double lambdaI, fractionalWI;
			double sati = this->problem.sat(global, *it, local);
			lambdaI = this->problem.materialLaw.mobTotal(sati);
			if (hasGravity) 
				fractionalWI = this->problem.materialLaw.fractionalW(sati);
			
			double faceVol[2*dim];

			// run through all intersections with neighbors and boundary
			IntersectionIterator endit = IntersectionIteratorGetter<G,LevelTag>::end(*it);
			for (IntersectionIterator is = IntersectionIteratorGetter<G,LevelTag>::begin(*it); 
				  is!=endit; ++is)
			{
			  // get geometry type of face
			  GeometryType gtf = is.intersectionSelfLocal().type();

			  //Geometry dg = is.intersectionSelfLocal();
			  // local number of facet 
			  int numberInSelf = is.numberInSelf();

			  faceVol[numberInSelf] = is.intersectionGlobal().volume();
			    
			  // center in face's reference element
			  const FieldVector<ct,dim-1>& 
			    facelocal = ReferenceElements<ct,dim-1>::general(gtf).position(0,0);
			  
			  // center of face inside volume reference element
			  const FieldVector<ct,dim>& 
			  	facelocalDim = ReferenceElements<ct,dim>::general(gtf).position(numberInSelf,1);
				    
			  // get normal vector
			  FieldVector<ct,dimworld> unitOuterNormal 
			    = is.unitOuterNormal(facelocal);
			  
			  // center of face in global coordinates
			  FieldVector<ct,dimworld> 
			    faceglobal = is.intersectionGlobal().global(facelocal);
			  
			  // handle interior face
			  if (is.neighbor()) 
			    {
			      // access neighbor
			      EntityPointer outside = is.outside();
			      int indexj = elementmapper.map(*outside);
			      
			      // get neighbor pressure and permeability
			      double pressj = this->press[indexj];
			      
			      // compute factor in neighbor
			      GeometryType nbgt = outside->geometry().type();
			      const FieldVector<ct,dim>& 
			      	nblocal = ReferenceElements<ct,dim>::general(nbgt).position(0,0);
			      
			      // neighbor cell center in global coordinates
			      FieldVector<ct,dimworld> 
			      	nbglobal = outside->geometry().global(nblocal);
			      
			      // distance vector between barycenters
			      FieldVector<ct,dimworld> distVec = global - nbglobal;
			      
			      // compute distance between cell centers
			      double dist = distVec.two_norm();
			      
			      // get absolute permeability 
			      FieldMatrix<ct,dim,dim> Kj(this->problem.K(nbglobal, *outside, nblocal));

			      // compute vectorized permeabilities
                  FieldVector<ct,dim> Kni(0);
                  FieldVector<ct,dim> Knj(0);
                  Ki.umv(unitOuterNormal, Kni);
                  Kj.umv(unitOuterNormal, Knj);
                  // compute permeability normal to intersection and take harmonic mean
                  double K_n_i = Kni * unitOuterNormal;
                  double K_n_j = Knj * unitOuterNormal;
                  double Kn    = 2 * K_n_i * K_n_j / (K_n_i + K_n_j);
                  // compute permeability tangential to intersection and take arithmetic mean
                  FieldVector<ct,dim> uON = unitOuterNormal;
                  FieldVector<ct,dim> K_t_i = Kni - (uON *= K_n_i);
                  uON = unitOuterNormal;
                  FieldVector<ct,dim> K_t_j = Knj - (uON *= K_n_j);
                  FieldVector<ct,dim> Kt = (K_t_i += K_t_j);
                  Kt *= 0.5;
                  // Build vectorized averaged permeability
                  uON = unitOuterNormal;
                  FieldVector<ct,dim> K = (Kt += (uON *=Kn));

			      //compute total mobility
			      double lambdaJ, fractionalWJ;
			      double satj = this->problem.sat(nbglobal, *outside, nblocal);
			    	  lambdaJ = this->problem.materialLaw.mobTotal(satj);
			    	  if (hasGravity) 	
			    		  fractionalWJ = this->problem.materialLaw.fractionalW(satj);
			      
			      // compute averaged total mobility
			      // CAREFUL: Harmonic weightig can generate zero matrix entries, 
			      // use arithmetic weighting instead: 
			      double lambda = 1;
			      double fractionalW;
			      lambda = 0.5*(lambdaI + lambdaJ); 
			      if (hasGravity) 
			    	  fractionalW = 0.5*(fractionalWI + fractionalWJ);
			      
			      FieldVector<ct,dimworld> vTotal(K);
			      vTotal *= lambda*(pressi - pressj)/dist;
			      if (hasGravity) {
			    	  Ki += Kj;
			    	  Ki *= 0.5;
			    	  FieldVector<ct,dimworld> gEffect(0);
			    	  Ki.umv(gravity, gEffect);
			    	  double factor = fractionalW*(this->problem.materialLaw.wettingPhase.density()) 
							+ (1 - fractionalW)*(this->problem.materialLaw.nonwettingPhase.density());
			    	  gEffect *= lambda*factor;
			    	  vTotal -= gEffect;
			      }
			      velocity[indexi][numberInSelf] = vTotal;		
			    }
			  // boundary face 
			  else 
			    { 
			      //get boundary condition for boundary face center
				  BoundaryConditions::Flags bctype = this->problem.bctype(faceglobal, *it, facelocalDim);
				  if (bctype == BoundaryConditions::dirichlet) {
				      // distance vector between barycenters
				      FieldVector<ct,dimworld> distVec = global - faceglobal;
				      
					  double dist = distVec.two_norm();
					  distVec /= dist;
					  
				      // compute directed permeability vector Ki.n
					  FieldVector<ct,dim> Kni(0);
					  Ki.umv(distVec, Kni);

					  // compute averaged total mobility
					  double lambda = 1.;
					  double fractionalW = 1.;
						  lambda = lambdaI;
						  if (hasGravity) fractionalW = fractionalWI;
						  
					  double g = this->problem.g(faceglobal, *it, facelocalDim);
					  
					  FieldVector<ct,dim> vTotal(Kni);
					  vTotal *= lambda*(g-pressi)/dist;
				      if (hasGravity) {
				    	  FieldVector<ct,dimworld> gEffect(0);
				    	  Ki.umv(gravity, gEffect);
				    	  double factor = fractionalW*(this->problem.materialLaw.wettingPhase.density()) 
								+ (1 - fractionalW)*(this->problem.materialLaw.nonwettingPhase.density());
				    	  gEffect *= lambda*factor;
				    	  vTotal -= gEffect;
				      }
					  velocity[indexi][numberInSelf] = vTotal;
			      } 
			      else
			      {		  
			    	  double J = this->problem.J(faceglobal, *it, facelocalDim);
			    	  FieldVector<ct,dimworld> unitOuterNormal 
			    	  	= is.unitOuterNormal(facelocal);
			    	  velocity[indexi][numberInSelf] = unitOuterNormal; 
			    	  velocity[indexi][numberInSelf] *= -J; 
			      }
			      
			    }
			  
			} // end all intersections  
			if (dim == 2) 
			{
				double sum = (fabs(velocity[indexi][0][0]*faceVol[0]) 
						+ fabs(velocity[indexi][1][0]*faceVol[1]) 
						+ fabs(velocity[indexi][2][1]*faceVol[2]) 
						+ fabs(velocity[indexi][3][1]*faceVol[3]));
			  double diff = fabs(velocity[indexi][0][0]*faceVol[0] 
					     - velocity[indexi][1][0]*faceVol[1]
					     + velocity[indexi][2][1]*faceVol[2] 
					     - velocity[indexi][3][1]*faceVol[3])
			              /sum;
				if (diff > 1e-6 && sum > 1e-9) 
				{
					std::cout << "NOT conservative!!! diff = " << diff << ", indexi = " << indexi << std::endl;
					std::cout << velocity[indexi][0][0]*faceVol[0] << ", " 
						  << velocity[indexi][1][0]*faceVol[1] << ", " 
						  << velocity[indexi][2][1]*faceVol[2] << ", " 
						  << velocity[indexi][3][1]*faceVol[3] <<  std::endl;
				}
			}
	    } // end grid traversal          
		  
	    return;
	}


  //! Returns the velocities on the edges of the grid ofthe specified level.
  /**CAUTION: works for axisymetric grids only!
   * CAUTION: specified level must be lower than level of the class!
   * \param velocity the return vector
   * \param t simulation time
   * \param lev the grid level on which the velocity is to be calculated. Must be lower than the level of the class!
  */
  template<class G, class RT>
  void FVDiffusion<G, RT>::totalVelocity(VelType& velocity, const RT t, const int lev) const
  {
    if (lev >= this->level())
      DUNE_THROW(NotImplemented,
	"FVDiffusion :: totalVelocity (VelType&, SatType&, RT, int) \n , lev >= this->level(), use (VelType&, SatType&, RT) instead");
    typedef typename G::template Codim<0>::HierarchicIterator HierarchicIterator;
    const IS& coarseIset = this->grid.levelIndexSet(0);
    velocity.resize(coarseIset.size(0));
    
    const int blocksize = 2 * dim;
    const int nFine = (int)pow(2, ((this->level()-lev)*(dim-1)) ); //fine element edges per coarse element edge
    int hitcount[blocksize];
    typedef FieldVector<ct,dim> R2;          // velocity for fine scale elemnt edge
    typedef BlockVector<R2> R3;              // vector holding velocities of all fine edges situated on one coarse edge
    BlockVector<R3> fineOnCoarse(blocksize); // vector for all faces of the coarse scale element
    R3 fineVelocity(blocksize);                    // vector holding all velocities of a fine-scale-element
    
    for (int i=0; i<blocksize ; i++)
    {
      fineOnCoarse[i].resize(nFine);
    } 
    
    Iterator eendit = this->grid.template lend<0>(lev);
    for (Iterator cit = this->grid.template lbegin<0>(lev); cit != eendit; ++cit)
    {
      for (int i=0; i<blocksize ; i++)  hitcount[i] = 0 ;
      int coarseIndex = coarseIset.index(*cit);
    	                 
      HierarchicIterator endit = cit-> hend(this->level());
      for (HierarchicIterator it = cit->hbegin(this->level()); it != endit; ++it)
      {   
        if (indexset.contains(*it)) 
        {
          //get some cell properties
          GeometryType gt = it->geometry().type();
          const FieldVector<ct,dim>& 
              local = ReferenceElements<ct,dim>::general(gt).position(0,0);
          FieldVector<ct,dimworld> global = it->geometry().global(local);   //global coordinates of cell center
          int indexi = this->indexset.index(*it);                           // index of fine-scale cell

          // get pressure and permeability and total mobility in fine-scale element
          double pressi = this->press[indexi];
          FieldMatrix<ct,dim,dim> Ki =  this->problem.K(global,*it,local);
          double sati = this->problem.sat(global, *it, local);
          double lambdaI = this->problem.materialLaw.mobTotal(sati);

	      double faceVol[2*dim];
          // run through all intersections with neighbors and boundary
          IntersectionIterator isend = it->ilevelend();
          for (IntersectionIterator is = it->ilevelbegin(); is!=isend; ++is)
          {
            // get some face properties
            GeometryType gtf = is.intersectionSelfLocal().type();
            int numberInSelf = is.numberInSelf();
	        faceVol[numberInSelf] = is.intersectionGlobal().volume();       // volume of face
            const FieldVector<ct,dim-1>&  facelocal       
              = ReferenceElements<ct,dim-1>::general(gtf).position(0,0);
            FieldVector<ct,dim> unitOuterNormal 
              = is.unitOuterNormal(facelocal);                              // normal vector of unit length
            FieldVector<ct,dim> faceglobal      
              = is.intersectionGlobal().global(facelocal);                  // global coordinate of face center
             // handle interior face
            if (is.neighbor()) 
            {
              // neighbor's properties
              EntityPointer outside = is.outside();
              int indexj = indexset.index(*outside);                        // neigbor's fine-scale cell index   
              GeometryType nbgt = outside->geometry().type();
              const FieldVector<ct,dim>& 
                nblocal = ReferenceElements<ct,dim>::general(nbgt).position(0,0);
              FieldVector<ct,dimworld> 
                nbglobal = outside->geometry().global(nblocal);             // global coordinate of neighbor's cell center

              // get neighbor pressure and permeability
              double pressj = this->press[indexj];
              FieldMatrix<ct,dim,dim> Kj =   this->problem.K(nbglobal,*outside,nblocal);

              // compute vectorized permeabilities
              FieldVector<ct,dim> Kni(0);
              FieldVector<ct,dim> Knj(0);
              Ki.umv(unitOuterNormal, Kni);
              Kj.umv(unitOuterNormal, Knj);
              // compute permeability normal to intersection and take harmonic mean
              double K_n_i = Kni * unitOuterNormal;
              double K_n_j = Knj * unitOuterNormal;
              double Kn    = 2 * K_n_i * K_n_j / (K_n_i + K_n_j);
              // compute permeability tangential to intersection and take arithmetic mean
              FieldVector<ct,dim> uON = unitOuterNormal;
              FieldVector<ct,dim> K_t_i = Kni - (uON *= K_n_i);
              uON = unitOuterNormal;
              FieldVector<ct,dim> K_t_j = Knj - (uON *= K_n_j);
              FieldVector<ct,dim> Kt = (K_t_i += K_t_j);
              Kt *= 0.5;
              // Build vectorized averaged permeability
              uON = unitOuterNormal;
              FieldVector<ct,dim> K = (Kt += (uON *=Kn));
             
              // distance between cell centers
              FieldVector<ct,dimworld> 
                  distVec = global - nbglobal;
              double dist = distVec.two_norm();
              // get averaged total mobility
              double satj = this->problem.sat(nbglobal,*outside,nblocal);
              double lambdaJ = this->problem.materialLaw.mobTotal(satj);
              double mob = 0.5 * (lambdaI + lambdaJ);
              // compute total velocity with Darcy's Law
              fineVelocity[numberInSelf] = K;		
              fineVelocity[numberInSelf] *= mob * (pressi-pressj) / dist;	
            }
            // boundary face 
            else 
            { 
              const FieldVector<ct,dim>& 
                facelocalDim = ReferenceElements<ct,dim>::general(gtf).position(numberInSelf,1);
              //get boundary condition for boundary face center
              BoundaryConditions::Flags bctype = this->problem.bctype(faceglobal, *it, facelocalDim);
              if (bctype == BoundaryConditions::dirichlet) 
              {
                // distance vector between barycenters
                FieldVector<ct,dimworld> distVec = global - faceglobal;
                double dist = distVec.two_norm();
                distVec /= dist;

                // compute directed permeability vector Ki.n
                FieldVector<ct,dim> Kni(0);
                Ki.umv(distVec, Kni);

                // compute averaged total mobility
                double lambda = lambdaI;
               
                double g = this->problem.g(faceglobal, *it, facelocalDim);
            					  
                FieldVector<ct,dim> vTotal(Kni);
                vTotal *= lambda*(g-pressi)/dist;
                fineVelocity[numberInSelf] = vTotal;
              } 
              else
              {		  
                double J = this->problem.J(faceglobal, *it, facelocalDim);
	    	    FieldVector<ct,dimworld> unitOuterNormal 
	    	  	  = is.unitOuterNormal(facelocal);
	    	    fineVelocity[numberInSelf] = unitOuterNormal; 
	    	    fineVelocity[numberInSelf] *= -J; 
              }			      
            }
            
            // check if fine scale face *is lies on any face *cis of the coarse scale element
            IntersectionIterator cisend = cit->ilevelend();
            for (IntersectionIterator cis = cit->ilevelbegin(); cis != cisend;  ++cis)
            {
              // checking if fine scale element face *is lies on coarse-scale element face *cis
              GeometryType cgtf = cis.intersectionSelfLocal().type();
              FieldVector<ct,dim-1> cfacelocal = ReferenceElements<ct,dim-1>::general(cgtf).position(0,0);
              bool isInCis;
              if ( cis.neighbor() ) 
              {
                /* if face *is lies in face *cis, *is must be defined both in the inside and outside
                   element of *cis! for the inside element this is true anyway, because of the hierarchic 
                   Iteration. */
                isInCis = cis.outside()->geometry().checkInside( cis.outside()->geometry().local(faceglobal));
              }
              else if (cis.boundary() && is.boundary())
              {
                /* if face *is lies in face *cis and
                 *cis lies on a boundary, *is must lie on a boundry too!
                 to be sure, that it is the same boundary, check if the normals are the same. */
                isInCis = ( cis.unitOuterNormal(cfacelocal) * is.unitOuterNormal(facelocal) ) > 0.999; 
              }
              else isInCis = false;
              if ( isInCis )
              {
                int coarseNumberInSelf = cis.numberInSelf();
                int hits =  hitcount[coarseNumberInSelf];
                fineOnCoarse[coarseNumberInSelf][hits] = fineVelocity[numberInSelf];
                // DEBUG--->
                FieldVector<ct,dim> gabagabahey = fineVelocity[numberInSelf];
                // <---DEBUG
                hitcount[coarseNumberInSelf] += 1;
              }
            }
          
          } // end intersection traversal
          
          // check for conservativity
	      if (dim == 2) 
	      {
	        double diff = fabs(fineVelocity[0][0]*faceVol[0] 
				 - fineVelocity[1][0]*faceVol[1]
				 + fineVelocity[2][1]*faceVol[2] 
				 - fineVelocity[3][1]*faceVol[3])
		          /(fabs(fineVelocity[0][0]*faceVol[0]) 
			    + fabs(fineVelocity[1][0]*faceVol[1]) 
			    + fabs(fineVelocity[2][1]*faceVol[2]) 
			    + fabs(fineVelocity[3][1]*faceVol[3]));
	        if (diff > 1e-6) 
		    {
		      std::cout << "NOT conservative!!! diff = " << diff << ", indexi = " << indexi << std::endl;
		      std::cout << fineVelocity[0][0]*faceVol[0] << ", " 
			    << fineVelocity[1][0]*faceVol[1] << ", " 
			    << fineVelocity[2][1]*faceVol[2] << ", " 
			    << fineVelocity[3][1]*faceVol[3] <<  std::endl;
		    }
	      }
        } // end the if( it2->isLeaf) environment
      }// end hierarchic iteration

      // evaluate mean velocities and standard deviations at coarse element edges 
      // and write them to the velocity struct. 
      for (int i = 0; i<blocksize; i++)
      {
        FieldVector<ct,dim> mean(0);
        // the mean
        for (int j = 0; j<nFine; j++)
        {
          mean += fineOnCoarse[i][j];
        }
        mean/= nFine;

        // write to velocity struct
        velocity[coarseIndex][i] = mean;
      }  
      
    } // end grid traversal          
  
    return;
  }// end method totalVelocity
  
}
#endif
