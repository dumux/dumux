#ifndef DUNE_BUCKLEYLEVERETTPROBLEM_HH
#define DUNE_BUCKLEYLEVERETTPROBLEM_HH

#include "dumux/transport/transportproblem.hh"

namespace Dune
{
//! \ingroup transportProblems
//! @brief example class for a transport problem
  template<class G, class RT>
  class BuckleyLeverettProblem 
  : public TransportProblem<G, RT, Dune::BlockVector< Dune::FieldVector<Dune::FieldVector<double, G::dimension>, 2*G::dimension> > > {
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
	BoundaryConditions::Flags bctype (const FieldVector<DT,n>& x, const Entity& e, 
					   const FieldVector<DT,n>& xi) const
	{
		if (x[0] > right-1E-8 || x[0] < left+1e-8) 
			return Dune::BoundaryConditions::dirichlet;
		else
			return Dune::BoundaryConditions::neumann;
	}

	RT g (const FieldVector<DT,n>& x, const Entity& e, 
		   const FieldVector<DT,n>& xi) const 
	{
		if (x[0] < left+1e-8) 
			return 0.8;
		else
			return 0.2;
	}
	  
	RT S0 (const FieldVector<DT,n>& x, const Entity& e, 
			const FieldVector<DT,n>& xi) const 
	{
		return 0.2;
	}
	  
	const FieldVector<DT,n>& vTotal (const Entity& e, const int numberInSelf)
	{
		int elemId = elementmapper.map(e);
		
		return(this->velocity[elemId][numberInSelf]);
	}

	BuckleyLeverettProblem(const G& g, TwoPhaseRelations& law = *(new LinearLaw), 
								const int level = 0, const bool cap = false) 
	: TransportProblem<G, RT, VelType>(law, cap), left((g.lowerLeft())[0]), right((g.upperRight())[0]), 
	  elementmapper(g, g.levelIndexSet(level))
	{	
		this->velocity.resize(elementmapper.size());
	}
  };

}
#endif
