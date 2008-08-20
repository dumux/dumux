#ifndef TESTPROBLEM_2P2C_HH
#define TESTPROBLEM_2P2C_HH

#include "dumux/transport/transportproblem2p2c.hh"
//#include <dumux/material/randompermeability.hh>

namespace Dune
{

  template<class G, class RT>
  class Testproblem_2p2c 
  : public TransportProblem2p2c<G, RT> 
  {
  	template<int dim>
		struct ElementLayout
		{
			bool contains (GeometryType gt)
			{
				return gt.dim() == dim;
			}
		}; 
		  
	  typedef typename G::ctype DT;
	  enum {n=G::dimension, m=1, blocksize=2*G::dimension};
	  typedef typename G::Traits::template Codim<0>::Entity Entity;
	  typedef FieldVector<double, n> R1;
	  typedef typename G::Traits::LevelIndexSet IS;
	  typedef MultipleCodimMultipleGeomTypeMapper<G,IS,ElementLayout> EM;

  private:
	  EM elementmapper;

  public: 
		 
		BoundaryConditions2p2c::Flags cbctype (const FieldVector<DT,n>& x, const Entity& e, 
						   const FieldVector<DT,n>& xi) const
		{
			
			return BoundaryConditions2p2c::concentration;
		}
		
		BoundaryConditions2p2c::Flags ictype (const FieldVector<DT,n>& x, const Entity& e, 
							   const FieldVector<DT,n>& xi) const
		{
			return BoundaryConditions2p2c::concentration;
		}
	  
		BoundaryConditions::Flags pbctype (const Dune::FieldVector<DT,n>& x, const Entity& e, 
						   const Dune::FieldVector<DT,n>& xi) const
	  {
//	    if (x[0] > 300-1E-6) 
	      return Dune::BoundaryConditions::dirichlet;
	    // all other boundaries
	    return Dune::BoundaryConditions::neumann;
	  }
		
//		const Dune::FieldMatrix<DT,n,n>& K (const Dune::FieldVector<DT,n>& x, const Entity& e, 
//							  const Dune::FieldVector<DT,n>& xi) 
//	  {
//		  return permeability.K(e);
//	  }
		
	  RT gPress (const FieldVector<DT,n>& x, const Entity& e, const FieldVector<DT,n>& xi) const
	  {
		  return (x[0] < 1e-6) ? 1e5 : 1e5;
	  }
		
		RT gZ (const FieldVector<DT,n>& x, const Entity& e, 
			   const FieldVector<DT,n>& xi) const 
		{
			if (x[0] < 15) 
				return 0;
			else
				return 0;
		}
		
		RT gS (const FieldVector<DT,n>& x, const Entity& e, 
					   const FieldVector<DT,n>& xi) const 
		{
			if (x[0] < 15) 
				return 0;
			else
				return 0;
		}
		
		virtual FieldVector<RT,2> J (const FieldVector<DT,n>& x, const Entity& e, 
				   const FieldVector<DT,n>& xi) const
		{
			FieldVector<RT,2> J_(0);
			if (x[0]<1e-6) J_[0] = 0.;
			return J_;
    }
		
		virtual FieldVector<RT,2> q (const FieldVector<DT,n>& x, const Entity& e, 
					   const FieldVector<DT,n>& xi) const 
		{
			FieldVector<RT,2> q_(0);
			FieldVector<DT,n> center(150); //center[1] = 200;
			if ((x-center).two_norm()<8) q_[1] = 0.001;
			return q_;
		}
		  
		RT S0 (const FieldVector<DT,n>& x, const Entity& e, 
				const FieldVector<DT,n>& xi) const 
		{
			return 0.999;
		}
		
		RT Z1_0 (const FieldVector<DT,n>& x, const Entity& e, 
				const FieldVector<DT,n>& xi) const
		{
			FieldVector<DT,n> center(110); center[1] = 200;
//			if ((x-center).two_norm()<30) return 0;
			if (fabs(x[0]-center[0])<30) return 1;
			return 1;
		}
	
		Testproblem_2p2c(G& g, Dune::VariableClass2p2c<G, RT> var, liquid_gl& liq, gas_gl& gas, Matrix2p<G, RT> s, int level, TwoPhaseRelations<G, RT>& law = *(new TwoPhaseRelations<G, RT>), 
				 const bool cap = false)
		: TransportProblem2p2c<G, RT>(var, liq, gas, s, law, cap), 
		  elementmapper(g, g.levelIndexSet(level)),
		  grid(g) 
//		  permeability(g, name, create)
		{	
		}
		
//		RandomPermeability<G> permeability;
		
  	private:
  		G& grid;
  };

}
#endif
