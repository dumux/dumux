#ifndef FVCA5TEST8PROBLEM_HH
#define FVCA5TEST8PROBLEM_HH

#include "dumux/diffusion/diffusionproblem.hh"

namespace Dune
{
//! \ingroup diffusionProblems
//! example class for diffusion problems
	template<class G, class RT>
	class FVCA5Test8Problem : public DiffusionProblem<G,RT>
	{
		template<int dim>
	    struct ElementLayout
	    {
	      bool contains (Dune::GeometryType gt)
	      {
	    	  return gt.dim() == dim;
	      }
	    }; 

	  typedef typename G::ctype DT;
	  enum {n=G::dimension};
	  typedef typename G::Traits::template Codim<0>::Entity Entity;
	  typedef typename G::Traits::LevelIndexSet IS;
	  typedef MultipleCodimMultipleGeomTypeMapper<G,IS,ElementLayout> EM;
	
	public:
	  FVCA5Test8Problem(G& grid)
	    : DiffusionProblem<G,RT>(), elementmapper(grid, grid.levelIndexSet(grid.maxLevel()))
	  { 
		  permloc_[0][0] = permloc_[1][1] = 1.0; 
		  permloc_[0][1] = permloc_[1][0] = 0.0;
	  }
	
	  const FieldMatrix<DT,n,n>& K (const FieldVector<DT,n>& x, const Entity& e, 
					  const FieldVector<DT,n>& xi) 
	  {
		  return permloc_;
	  }
	
	  RT q   (const FieldVector<DT,n>& x, const Entity& e, 
					  const FieldVector<DT,n>& xi)
	  {
		  int elemId = elementmapper.map(e);
		  
		  if (elemId == 60) {
			  double volume = e.geometry().volume();
			  return (1.0/volume);
		  }
		  else 
			  return (0.0);  
	  }
	
	  typename BoundaryConditions::Flags bctype (const FieldVector<DT,n>& x, const Entity& e, 
						   const FieldVector<DT,n>& xi) const
	  {
	      return BoundaryConditions::dirichlet;
	  }
	
	  RT g (const FieldVector<DT,n>& x, const Entity& e, 
					const FieldVector<DT,n>& xi) const
	  {
		  return (0.0);
	  }
		  
		
	  RT J (const FieldVector<DT,n>& x, const Entity& e, 
					const FieldVector<DT,n>& xi) const
	  {
		return 0;
	  }

	  RT exact (const FieldVector<DT,n>& x) const
	  {
		  return 0;
	  }

	  FieldVector<RT,n> exactGrad (const FieldVector<DT,n>& x) const
	  {	
		  FieldVector<RT,n> grad(0);
		  return grad;
	  }
	
	private:
		FieldMatrix<DT,n,n> permloc_;
		EM elementmapper;
	};
}

#endif
