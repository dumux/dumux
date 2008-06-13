#ifndef DUNE_PARALLELBOXDIFFUSION_HH
#define DUNE_PARALLELBOXDIFFUSION_HH

#include <dune/disc/shapefunctions/lagrangeshapefunctions.hh>
#include <dune/disc/functions/p1function.hh>
//#include <dune/disc/operators/p1operator.hh>
#include <dune/istl/io.hh>
#include <dune/common/timer.hh>
#include <dune/istl/bvector.hh>
#include <dune/istl/vbvector.hh>
#include <dune/istl/bcrsmatrix.hh>
#include <dune/istl/io.hh>
#include <dune/istl/gsetc.hh>
#include <dune/istl/ilu.hh>
#include <dune/istl/operators.hh>
#include <dune/istl/solvers.hh>
#include <dune/istl/preconditioners.hh>
#include <dune/istl/scalarproducts.hh>
#include <dune/istl/paamg/amg.hh>
#include <dune/istl/owneroverlapcopy.hh>
#include "dumux/nonlinear/nonlinearmodel.hh"
#include "dumux/fvgeometry/fvelementgeometry.hh"
#include "dumux/nonlinear/newtonmethod.hh"
#include "parallelboxdiffusionjacobian.hh"
#include "diffusionparameters.hh"
#include "schwarzcommunication.hh"
#include "schwarzpreconditioner.hh"
#include "dumux/pardiso/pardiso.hh"

namespace Dune
{
  template<class G, class RT, class ProblemType, class LocalJacobian, 
            class FunctionType, class OperatorAssembler>
  class ParallelBoxDiffusion 
  : public NonlinearModel<G, RT, ProblemType, LocalJacobian, FunctionType, OperatorAssembler> 
  {
  public:	
	typedef NonlinearModel<G, RT, ProblemType, LocalJacobian, 
	                          FunctionType, OperatorAssembler> NonlinearModel;
	
	ParallelBoxDiffusion(const G& g, ProblemType& prob)
	: NonlinearModel(g, prob), uOldTimeStep(g)
	{ }
	
	ParallelBoxDiffusion(const G& g, ProblemType& prob, int level)
	: NonlinearModel(g, prob, level), uOldTimeStep(g, level)
	{ 	}
	
	virtual void initial() = 0;
	
	virtual void update(double& dt) = 0;
	
	virtual void solve() = 0;
	
	FunctionType uOldTimeStep;
  };


  
  
  
  template<class G, class RT, int m=1>
  class LeafP1ParallelBoxDiffusion : public ParallelBoxDiffusion<G, RT, DiffusionParameters<G, RT>, ParallelBoxDiffusionJacobian<G, RT>, 
                                        LeafP1Function<G, RT, m>, LeafP1OperatorAssembler<G, RT, m> >
  {
  public:
	  // define the function type:
	  typedef LeafP1Function<G, RT, m> FunctionType;

	  // define the operator assembler type:
	  typedef LeafP1OperatorAssembler<G, RT, m> OperatorAssembler;

	  typedef ParallelBoxDiffusion<G, RT, DiffusionParameters<G, RT>, ParallelBoxDiffusionJacobian<G, RT>, 
	                          FunctionType, OperatorAssembler> ParallelBoxDiffusion;
	  
	  typedef LeafP1ParallelBoxDiffusion<G, RT, m> ThisType;

	  typedef ParallelBoxDiffusionJacobian<G, RT> LocalJacobian;
	  
	  // mapper: one data element per vertex
	  template<int dim>
	  struct P1Layout
	  {
		  bool contains (Dune::GeometryType gt)
		  {
			  return gt.dim() == 0;
		  }
	  }; 

	  typedef typename G::Traits::LeafIndexSet IS;
	  typedef MultipleCodimMultipleGeomTypeMapper<G,IS,P1Layout> VertexMapper;
	  typedef typename IntersectionIteratorGetter<G,LeafTag>::IntersectionIterator IntersectionIterator;
		typedef typename ThisType::FunctionType::RepresentationType VectorType;
		typedef typename ThisType::OperatorAssembler::RepresentationType MatrixType;
		typedef MatrixAdapter<MatrixType,VectorType,VectorType> Operator; 

	  LeafP1ParallelBoxDiffusion (const G& g, DiffusionParameters<G, RT>& prob) 
	  : ParallelBoxDiffusion(g, prob), grid(g), vertexmapper(g, g.leafIndexSet()), 
	    size((*(this->u)).size())
	  { }
	  
	  virtual void initial() 
	  {
		  typedef typename G::Traits::template Codim<0>::Entity Entity;
		  typedef typename G::ctype DT;
		  typedef typename IS::template Codim<0>::template Partition<All_Partition>::Iterator Iterator;
		  enum{dim = G::dimension};
		  enum{dimworld = G::dimensionworld};
		  
		  const IS& indexset(grid.leafIndexSet());
		  std::cout << "initializing solution." << std::endl;
		  // iterate through leaf grid an evaluate c0 at cell center
		  Iterator eendit = indexset.template end<0, All_Partition>();
		  for (Iterator it = indexset.template begin<0, All_Partition>(); it != eendit; ++it)
		  {
			  // get geometry type
			  Dune::GeometryType gt = it->geometry().type();

			  // get entity 
			  const Entity& entity = *it;

			  const typename Dune::LagrangeShapeFunctionSetContainer<DT,RT,dim>::value_type& 
		      	sfs=Dune::LagrangeShapeFunctions<DT,RT,dim>::general(gt, 1);
		      int size = sfs.size();

		      for (int i = 0; i < size; i++) {
	    		  int globalId = vertexmapper.template map<dim>(entity, sfs[i].entity());
	    	  
	    		  // initialize cell concentration
	    		  (*(this->u))[globalId] = 0;
		      }
		  }

		  // set Dirichlet boundary conditions
		  for (Iterator it = indexset.template begin<0, All_Partition>(); it != eendit; ++it)
		  {
			  // get geometry type
			  Dune::GeometryType gt = it->geometry().type();

			  // get entity 
			  const Entity& entity = *it;

			  const typename Dune::LagrangeShapeFunctionSetContainer<DT,RT,dim>::value_type& 
		      	sfs=Dune::LagrangeShapeFunctions<DT,RT,dim>::general(gt, 1);
		      int size = sfs.size();

		      // set type of boundary conditions 
		      this->localJacobian.fvGeom.update(entity);
		      this->localJacobian.template assembleBC<LeafTag>(entity);

// 			  for (int i = 0; i < size; i++) 
// 			    std::cout << "bc[" << i << "] = " << this->localJacobian.bc(i) << std::endl;
			  
			  IntersectionIterator endit = IntersectionIteratorGetter<G,LeafTag>::end(entity);
			  for (IntersectionIterator is = IntersectionIteratorGetter<G,LeafTag>::begin(entity); 
			       is!=endit; ++is)
				  if (is->boundary())
				  {
				    for (int i = 0; i < size; i++) 
					  // handle subentities of this face
					  for (int j = 0; j < ReferenceElements<DT,dim>::general(gt).size(is->numberInSelf(), 1, sfs[i].codim()); j++)
						if (sfs[i].entity() == ReferenceElements<DT,dim>::general(gt).subEntity(is->numberInSelf(), 1, j, sfs[i].codim()))
						{
							if (this->localJacobian.bc(i)[0] == BoundaryConditions::dirichlet) 
							{
								// get cell center in reference element
								Dune::FieldVector<DT,dim> local = sfs[i].position();

								// get global coordinate of cell center
								Dune::FieldVector<DT,dimworld> global = it->geometry().global(local);

								int globalId = vertexmapper.template map<dim>(entity, sfs[i].entity());
		    	  
								FieldVector<BoundaryConditions::Flags, m> bctype = this->problem.bctype(global, entity, is, local);
// 								std::cout << "global = " << global << ", id = " << globalId << std::endl;			  
								if (bctype[0] == BoundaryConditions::dirichlet) {
									(*(this->u))[globalId] = this->problem.g(global, entity, is, local);
								}
								else {
									std::cout << global << " is considered to be a Neumann node." << std::endl;
								}
							}
						}
			  }
		  }

		  *(this->uOldTimeStep) = *(this->u);
		  return;
	  }


	virtual void update(double& dt)
	{
		this->localJacobian.setDt(dt);
		this->localJacobian.setOldSolution(this->uOldTimeStep);
		NewtonMethod<G, ThisType> newtonMethod(this->grid, *this);
		newtonMethod.execute();
		dt = this->localJacobian.getDt();
		*(this->uOldTimeStep) = *(this->u);

		return;		
	}
		
	virtual void solve()
	{
		typedef typename G::Traits::GlobalIdSet::IdType GlobalIdType;
		typedef typename Dune::LeafP1Function<G,RT>::P1IndexInfoFromGrid P1IndexInfoFromGrid;

	  Dune::MatrixAdapter<MatrixType,VectorType,VectorType> op(*(this->A));
	  Dune::SeqILU0<MatrixType,VectorType,VectorType> ilu0(*(this->A),1.0);
	  //Dune::SeqSSOR<MatrixType,VectorType,VectorType> ssor(*(this->A),1,1.0);
	  
#if HAVE_MPI
	  // set up parallel solvers
	  Dune::IndexInfoFromGrid<GlobalIdType,int> indexinfo; 
	  (this->u).fillIndexInfoFromGrid(indexinfo);
	  typedef Dune::OwnerOverlapCopyCommunication<GlobalIdType,int> CommunicationType;
	  CommunicationType oocc(indexinfo,grid.comm());
	  int verbose=0;
	  if (grid.comm().rank() == 0) 
	    verbose = 1;
	  Dune::OverlappingSchwarzOperator<MatrixType,VectorType,VectorType,CommunicationType> oop(*(this->A),oocc);
	  Dune::OverlappingSchwarzScalarProduct<VectorType,CommunicationType> osp(oocc);
	  Dune::BlockPreconditioner<VectorType,VectorType,CommunicationType> parprec(ilu0,oocc);
	  Dune::CGSolver<VectorType> parcg(oop,osp,parprec,1E-12,1000,verbose);

	  // solve system
	  Dune::InverseOperatorResult r;	
	  parcg.apply(*(this->u), *(this->f), r);
#endif
	  return;				
	}

	virtual void vtkout (const char* name, int k) 
	{
		VTKWriter<G> vtkwriter(this->grid);
		vtkwriter.addVertexData(*(this->u),"pressure");
		vtkwriter.write(name, VTKOptions::ascii);
	}

    void globalDefect(FunctionType& defectGlobal)
    {   
      typedef typename G::Traits::template Codim<0>::Entity Entity;
      typedef typename G::ctype DT;
      typedef typename IS::template Codim<0>::template Partition<All_Partition>::Iterator Iterator;
      enum{dim = G::dimension};
      
      const IS& indexset(this->grid.leafIndexSet());
      (*defectGlobal)=0;
           
      // iterate through leaf grid 
      Iterator eendit = indexset.template end<0, All_Partition>();
      for (Iterator it = indexset.template begin<0, All_Partition>(); it != eendit; ++it)
	{
	  // get geometry type
	  Dune::GeometryType gt = it->geometry().type();
	  
	  const typename Dune::LagrangeShapeFunctionSetContainer<DT,RT,dim>::value_type& 
	    sfs=Dune::LagrangeShapeFunctions<DT,RT,dim>::general(gt, 1);
	  int size = sfs.size();
	  Dune::FieldVector<RT,1> defhelp[size];
	  
	  // get entity 
	  const Entity& entity = *it;
	  
	  this->localJacobian.fvGeom.update(entity);
	  
	  this->localJacobian.getLocalDefect(entity,defhelp);
	  //std::cout<<" defhelp: "<<*defhelp<<std::endl;
	  // begin loop over vertices
	  for(int i=0; i < size; i++)
	    {
	      int globalId = this->vertexmapper.template map<dim>(entity, sfs[i].entity());
	    
	      if (this->localJacobian.bc(i)[0] == BoundaryConditions::neumann)
		(*defectGlobal)[globalId] += defhelp[i];
	      else 
		(*defectGlobal)[globalId] = 0;

	    }
	}
    }
protected:
  const G& grid;
  VertexMapper vertexmapper;
  int size;
};

}
#endif
