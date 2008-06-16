#ifndef DUNE_FVBRINKMAN_HH
#define DUNE_FVBRINKMAN_HH

#include <dune/common/helpertemplates.hh>
#include <dune/common/typetraits.hh>
#include <dune/grid/common/mcmgmapper.hh>
#include <dune/grid/io/file/vtk/vtkwriter.hh>
#include <dune/grid/utility/intersectiongetter.hh>
#include <dune/istl/operators.hh>
#include <dune/istl/solvers.hh>
#include <dune/istl/preconditioners.hh>
#include "dumux/brinkman/brinkman.hh"

/**
 * @file
 * @brief  Finite Volume Brinkman Model
 * @author Bernd Flemisch, Jochen Fritz
 */

namespace Dune
{
  //! \ingroup diffusion
  //! Finite Volume Brinkman Model
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
  class FVBrinkman 
  : public Brinkman< G, RT, BlockVector< FieldVector<RT,1> >, BlockVector< FieldVector<RT,G::dimension> > >
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
	  typedef BCRSMatrix<MB> PressureMatrixType;
	  typedef FieldMatrix<double,dim,dim> MBV;
	  typedef BCRSMatrix<MBV> VelocityMatrixType;
	  typedef FieldVector<double, 1> VB;
	  typedef BlockVector<VB> Vector;
	  
  public:
	typedef BlockVector< FieldVector<FieldVector<RT, G::dimension>, 2*G::dimension> > VelType;
	typedef BlockVector< FieldVector<RT,1> > RepresentationType;

	void updateVelocityRHS(); 

	void solveVelocitySytem(); 

	void computeVelocity()
	{
		updateVelocityRHS();
		solveVelocitySytem();
		return;
	}

	void vtkout (const char* name, int k) const 
	{
		VTKWriter<G> vtkwriter(this->grid);
		char fname[128];	
		sprintf(fname,"%s-%05d",name,k);
		vtkwriter.addCellData(this->pressure,"total pressure p~");
		vtkwriter.write(fname,VTKOptions::ascii);		
	}
	
	void initializeMatrices();
	
	void assembleMatrices();
	
	FVBrinkman(G& g, BrinkmanProblem<G, RT>& prob)
	            : Brinkman<G, RT, RepresentationType, VelType>(g, prob), 
	              elementmapper(g, g.leafIndexSet()), 
	              indexset(g.leafIndexSet()), 
	              AV(g.size(0), g.size(0), (2*dim+1)*g.size(0), BCRSMatrix<MB>::random), 
	              AP(g.size(0), g.size(0), (2*dim+1)*g.size(0), BCRSMatrix<MB>::random), 
	              f(g.size(0))
	{
		this->pressure.resize(g.size(0));
		this->pressure = 0;
		this->pressureCorrection.resize(g.size(0));
		this->pressureCorrection = 0;
		this->velocity.resize(g.size(0));
		this->velocity = 0;
		this->velocityCorrection.resize(g.size(0));
		this->velocityCorrection = 0;
		initializeMatrices();
		assembleMatrices();
	}
	
  private:
	  EM elementmapper;
	  const IS& indexset;
	  VelocityMatrixType AV;
	  PressureMatrixType AP;
	  RepresentationType f;
  };

  
  
  template<class G, class RT>
  void FVBrinkman<G, RT>::initializeMatrices()
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
			    if (is->neighbor()) 
			      rowSize++;
			AV.setrowsize(indexi, rowSize);
			AP.setrowsize(indexi, rowSize);
	      }
	    AV.endrowsizes();
	    AP.endrowsizes();

	    // determine position of matrix entries 
	    for (Iterator it = indexset.template begin<0,All_Partition>(); it != eendit; ++it)
	      {
			// cell index
			int indexi = elementmapper.map(*it);
	
			// add diagonal index
			AV.addindex(indexi, indexi);
			AP.addindex(indexi, indexi);
	
			// run through all intersections with neighbors 
			IntersectionIterator endit = IntersectionIteratorGetter<G,LevelTag>::end(*it);
			for (IntersectionIterator is = IntersectionIteratorGetter<G,LevelTag>::begin(*it); 
			  	  is!=endit; ++is)
			    if (is->neighbor()) 
			      {
					// access neighbor
					EntityPointer outside = is->outside();
					int indexj = elementmapper.map(*outside);
		
					// add off diagonal index
					AV.addindex(indexi, indexj);
					AP.addindex(indexi, indexj);
			      }
	      }
	    AV.endindices();		
	    AP.endindices();		

	    return;
  }
  
  template<class G, class RT>
  void FVBrinkman<G, RT>::assembleMatrices()
  {
        // initialization: set matrix A to zero	   
        AV = 0;
        AP = 0;

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
	    
	    // cell volume 
	    double volume = it->geometry().integrationElement(local)
	      *ReferenceElements<ct,dim>::general(gt).volume();
	    
	    // get absolute permeability 
	    FieldMatrix<ct,dim,dim> Kinv(this->problem.K(global,*it,local));
	    
	    // get effective viscosity
	    RT muEffI = this->problem.muEff(global,*it,local);
	    
	    // get viscosity
	    RT mu = this->problem.mu(global,*it,local);
	    
	    Kinv *= volume*mu;
	    
	    AV[indexi][indexi] = Kinv;
	    
	    IntersectionIterator endit = IntersectionIteratorGetter<G,LevelTag>::end(*it);
	    for (IntersectionIterator is = IntersectionIteratorGetter<G,LevelTag>::begin(*it); 
		 is!=endit; ++is)
	      {
		
		// get geometry type of face
		GeometryType gtf = is->intersectionSelfLocal().type();
		
		// center in face's reference element
		const FieldVector<ct,dim-1>& 
		  facelocal = ReferenceElements<ct,dim-1>::general(gtf).position(0,0);
		
		// center of face inside volume reference element
		const FieldVector<ct,dim>& 
		  facelocalDim = ReferenceElements<ct,dim>::general(gtf).position(is->numberInSelf(),1);
		
		// get normal vector 
		FieldVector<ct,dimworld> unitOuterNormal 
		  = is->unitOuterNormal(facelocal);
		
		// get normal vector scaled with volume
		FieldVector<ct,dimworld> integrationOuterNormal 
		  = is->integrationOuterNormal(facelocal);
		
		// get face volume 
		double faceVol = is->intersectionGlobal().volume();
		
		// handle interior face
		if (is->neighbor()) 
		  {
		    // access neighbor
		    EntityPointer outside = is->outside();
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
		    
		    // get the effective viscosity 
		    RT muEffJ = this->problem.muEff(nbglobal, *outside, nblocal);
		    
		    // average the effective viscosity 
		    RT muEff = 2.0*muEffI*muEffJ/(muEffI + muEffJ);
		    
		    FieldMatrix<RT,dim,dim> gradUn(0);
		    for (int k = 0; k < dim; k++) {
		    	gradUn[k] = distVec; 
		    	gradUn[k] /= dist*dist; 
		    	for (int l = 0; l < dim; l++)
		    		gradUn[k][l] *= unitOuterNormal[l];
		    }
		    
		    // update diagonal entry 
		    AV[indexi][indexi] += gradUn;
		    
		    // set off-diagonal entry 
		    AV[indexi][indexj] = -gradUn;
		  }
		// boundary face 
		else 
		  { 
			// center of face in global coordinates
		    FieldVector<ct,dimworld> 
		      faceglobal = is->intersectionGlobal().global(facelocal);
		    
		    //get boundary condition for boundary face center
		    BoundaryConditions::Flags bctype = this->problem.bctype(faceglobal, *it, facelocalDim);
		    if (bctype == BoundaryConditions::dirichlet) 
		      { 
		    	FieldVector<ct,dimworld> distVec(global - faceglobal);
		    	double dist = distVec.two_norm();
		    	
			    FieldMatrix<RT,dim,dim> gradUn(0);
			    for (int k = 0; k < dim; k++) {
			    	gradUn[k] = distVec; 
			    	gradUn[k] /= dist*dist; 
			    	for (int l = 0; l < dim; l++)
			    		gradUn[k][l] *= unitOuterNormal[l];
			    }
			    
			    // update diagonal entry 
			    AV[indexi][indexi] += gradUn;
			    
			    // set off-diagonal entry 
			    AV[indexi][indexj] = -gradUn;

		    	A[indexi][indexi] -= lambda*faceVol*(Kni*distVec)/(dist*dist);
		    	double g = this->problem.g(faceglobal, *it, facelocalDim);
		    	f[indexi] -= lambda*faceVol*g*(Kni*distVec)/(dist*dist);
			
		      } 
		    else
		      {
		    	;
		      }
		    
		  }
	      } // end all intersections         
	    
	  } // end grid traversal 
	return;
	}
	
	
  template<class G, class RT>
  void FVBrinkman<G, RT>::solve()
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
			  DUNE_THROW(NotImplemented, "FVBrinkman :: solve : combination " << preconditionerName_ 
					  << " and " << solverName_ << ".");
	  }
	  else if (preconditionerName_ == "SeqPardiso") {
	      SeqPardiso<MatrixType,Vector,Vector> preconditioner(A);
	      if (solverName_ == "Loop") {
	    	  LoopSolver<Vector> solver(op, preconditioner, 1E-14, 10000, 1);
	    	  solver.apply(this->press, f, r);
	      }
	      else 
	    	  DUNE_THROW(NotImplemented, "FVBrinkman :: solve : combination " << preconditionerName_ 
	    			  << " and " << solverName_ << ".");
	  }
	  else 
		  DUNE_THROW(NotImplemented, "FVBrinkman :: solve : preconditioner " << preconditionerName_ << ".");

	  return;
  }
	
  
}
#endif
