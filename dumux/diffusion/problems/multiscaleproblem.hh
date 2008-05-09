#ifndef MULTISCALEPROBLEM_HH
#define MULTISCALEPROBLEM_HH

#include "dumux/diffusion/diffusionproblem.hh"
#include <dune/istl/bvector.hh>

namespace Dune
{
//! \ingroup diffusionProblems
//! converter for diffusion and transport classes on different scales.
//	template <class G, class RT, class Problem>
//	class MultiscaleProblem:public DiffusionProblem<G, RT>
//	{
////	  typedef typename Problem::ReturnType RT;
////	  typedef typename Problem::GridType::ctype DT;
////	  typedef typename Problem::GridType::Traits::template Codim<0>::Entity Entity;
////	  enum {n=Problem::GridType::dimension};
//	  typedef typename G::Traits::template Codim<0>::Entity Entity;
//	  enum {n=G::dimension};
//	  typedef typename G::ctype DT;
//	public:
//	  MultiscaleProblem(Problem& p, const int diff_level, const int trans_level, TwoPhaseRelations& law = *(new LinearLaw), 
//				const bool cap = false) 
//	  : DiffusionProblem<G, RT>(law, cap), problem(p), l_diffusion(diff_level), l_transport(trans_level)
//	  { 
//		if (l_diffusion < l_transport) DUNE_THROW(NotImplemented, "In class MultiscaleProblem: level of transport class must be lower than level in diffusion class!");
//	  }
//	  
//	  RT sat (const Dune::FieldVector<DT, n>& x, const Entity& e, 
//					  const Dune::FieldVector<DT,n>& xi)
//	  {
//		 if (e.level() == l_transport) 
//			 return problem.sat(x, e, xi);
//		 typename Entity::EntityPointer f = e.father();
//		 while (f.level() > l_transport) f = f->father();
//		 return problem.sat(x,*f,xi);
//	  }
//	  
//	  const FieldMatrix<DT,n,n>& K (const FieldVector<DT,n>& x, const Entity& e, 
//	  					const FieldVector<DT,n>& xi)
//	  {
//		  return problem.K(x,e,xi);
//	  }
//	  
//	  RT q   (const FieldVector<DT,n>& x, const Entity& e, 
//	  					const FieldVector<DT,n>& xi)
//	  {
//		  return problem.q(x,e,xi);
//	  }
//	  
//	  BoundaryConditions::Flags bctype (const FieldVector<DT,n>& x, const Entity& e, 
//	  					   const FieldVector<DT,n>& xi)
//	  {
//		  return problem.bctype(x,e,xi);
//	  }
//	  
//	  RT g (const FieldVector<DT,n>& x, const Entity& e, 
//	  				  const FieldVector<DT,n>& xi)
//	  {
//	  	  return problem.g(x,e,xi);
//	  }
//	  
//	  RT J (const FieldVector<DT,n>& x, const Entity& e, 
//	  				  const FieldVector<DT,n>& xi)
//	  {
//		  return problem.J(x,e,xi);
//	  }
//	  
//	  const FieldVector<DT,n>& gravity()
//	  {
//		  return problem.gravity();
//	  }
//	  
//	private:
//		Problem& problem;
//		int l_transport;
//		int l_diffusion;
//	};

template<class G, class RT>
	class MultiscaleProblem : public DiffusionProblem<G,RT>
	{
	  typedef typename G::ctype DT;
	  enum {n=G::dimension};
	  typedef typename G::Traits::template Codim<0>::Entity Entity;
	
	public:
	  MultiscaleProblem(G& g, DiffusionProblem<G,RT>& p,const int diff_level, const int trans_level, TwoPhaseRelations& law = *(new LinearLaw), const bool cap = false)
	    : DiffusionProblem<G,RT>(p.materialLaw, p.capillary), l_diffusion(diff_level), l_transport(trans_level)
	  { 
		  problem = &p;
	  }
	
	  const Dune::FieldMatrix<DT,n,n>& K (const Dune::FieldVector<DT,n>& x, const Entity& e, 
					  const Dune::FieldVector<DT,n>& xi) 
	  {
		  return problem->K(x,e,xi);
	  }
	  
	  RT sat (const Dune::FieldVector<DT,n>& x, const Entity& e, 
					  const Dune::FieldVector<DT,n>& xi)
	  {
		  if (e.level() == l_transport) 
			  return problem->sat(x, e, xi);
		  typename Entity::EntityPointer f = e.father();
		  while (f.level() > l_transport) f = f->father();
		  return problem->sat(x,*f,xi);
	  }
	
	  RT q   (const Dune::FieldVector<DT,n>& x, const Entity& e, 
					  const Dune::FieldVector<DT,n>& xi)
	  {
		return problem->q(x,e,xi);
	  }
	
	  typename Dune::BoundaryConditions::Flags bctype (const Dune::FieldVector<DT,n>& x, const Entity& e, 
						   const Dune::FieldVector<DT,n>& xi) const
	  {
	    return problem->bctype(x,e,xi);
	  }
	
	  RT g (const Dune::FieldVector<DT,n>& x, const Entity& e, 
					const Dune::FieldVector<DT,n>& xi) const
	  {
		  return problem->g(x,e,xi);
	  }
		  
		
	  RT J (const Dune::FieldVector<DT,n>& x, const Entity& e, 
					const Dune::FieldVector<DT,n>& xi) const
	  {
		return problem->J(x,e,xi);
	  }
		  
	  
	private:
		DiffusionProblem<G,RT>* problem;
	public:
				int l_transport;
				int l_diffusion;
	};

}

#endif
