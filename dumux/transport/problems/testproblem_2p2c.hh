#ifndef TESTPROBLEM_2P2C_HH
#define TESTPROBLEM_2P2C_HH

#include "dumux/transport/transportproblem2p2c.hh"
#include "dumux/diffusion/diffusion.hh"

namespace Dune
{

  template<class G, class RT>
  class Testproblem_2p2c 
  : public TransportProblem2p2c<G, RT, Dune::BlockVector< Dune::FieldVector<Dune::FieldVector<double, G::dimension>, 2*G::dimension> > > {
		template<int dim>
		struct ElementLayout
		{
			bool contains (Dune::GeometryType gt)
			{
				return gt.dim() == dim;
			}
		}; 
		  
	  typedef typename G::ctype DT;
	  enum {n=G::dimension, m=1, blocksize=2*G::dimension};
	  typedef typename G::Traits::template Codim<0>::Entity Entity;
	  typedef Dune::FieldVector<double, n> R1;
	  typedef Dune::BlockVector< Dune::FieldVector<R1, blocksize> > VelType;
	  typedef typename G::Traits::LevelIndexSet IS;
	  typedef Dune::MultipleCodimMultipleGeomTypeMapper<G,IS,ElementLayout> EM;

  private:
	  DT left;
	  DT right;
	  EM elementmapper;

  public:
	  
	  const Dune::FVDiffusion<G,RT>* diffusion;  
	  
	BoundaryConditions::Flags bctype (const FieldVector<DT,n>& x, const Entity& e, 
	   const FieldVector<DT,n>& xi) const
	{
		return BoundaryConditions::dirichlet;
	}
	   
	  
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

	RT g (const FieldVector<DT,n>& x, const Entity& e, 
		   const FieldVector<DT,n>& xi) const 
	{
		if (x[0] < left+1e-8) 
			return 1000;
		else
			return 1000;
	}
	  
	RT S0 (const FieldVector<DT,n>& x, const Entity& e, 
			const FieldVector<DT,n>& xi) const 
	{
		return 1.0;
	}
	
	RT C1_0 (const FieldVector<DT,n>& x, const Entity& e, 
			const FieldVector<DT,n>& xi) const 
	{
		FieldVector<DT,n> center(150); center[1] = 150;
		double radius = 30;
		if ((x-center).two_norm() < radius) return 950;
		return 1000;
	}
	  
	const FieldVector<DT,n>& vTotal (const Entity& e, const int numberInSelf)
	{
		int elemId = elementmapper.map(e);
		
		return(this->velocity[elemId][numberInSelf]);
	}
	
	RT press (const FieldVector<DT,n>& x, const Entity& e, 
			  const FieldVector<DT,n>& xi) const 
	{
		int index = elementmapper.map(e);
		return diffusion->press[index];
	}
	
	RT pressBC (const FieldVector<DT,n>& x, const Entity& e, 
			  const FieldVector<DT,n>& xi) const
	{
		return 2e6 - xi[0] * 1e6 / 300;
	}

	Testproblem_2p2c(Henry& h, const G& g, TwoPhaseRelations& law /*= *(new LinearLaw)*/, 
								const int level = 0, const bool cap = false) 
	: TransportProblem2p2c<G, RT, VelType>(h, law, cap), left((g.lowerLeft())[0]), right((g.upperRight())[0]), 
	  elementmapper(g, g.levelIndexSet(level))
	{	
		this->velocity.resize(elementmapper.size());
	}
  };

}
#endif
