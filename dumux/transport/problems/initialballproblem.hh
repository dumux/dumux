// $Id$ 

#ifndef DUNE_INITIALBALLPROBLEM_HH
#define DUNE_INITIALBALLPROBLEM_HH

#include "dumux/transport/transportproblem.hh"

namespace Dune
{
//! \ingroup transportProblems
//! @brief example class for a transport problem
  template<class G, class RT, class VC>
  class InitialBallProblem 
  : public TransportProblem<G, RT, VC> {
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
	  DT bottom;
	  DT top;
	  EM elementmapper;

  public:
	BoundaryConditions::Flags bctype (const FieldVector<DT,n>& x, const Entity& e, 
					   const FieldVector<DT,n>& xi) const
	{
		if (x[1] > top-1E-8 || x[1] < bottom+1e-8) 
			return Dune::BoundaryConditions::dirichlet;
		else
			return Dune::BoundaryConditions::neumann;
	}

	RT g (const FieldVector<DT,n>& x, const Entity& e, 
		   const FieldVector<DT,n>& xi) const 
	{
			return 0.8;
	}
	  
	RT S0 (const FieldVector<DT,n>& x, const Entity& e, 
			const FieldVector<DT,n>& xi) const 
	{
		FieldVector<DT,n> dist(0.5*(left + right));
		dist[1] = top - 2.0;
		dist -= x;
		if (dist.two_norm() < 0.25)
			return 0.2;
		else
			return 0.8;	
	}
	  
	const FieldVector<DT,n>& vTotal (const Entity& e, const int numberInSelf)
	{
		int elemId = elementmapper.map(e);
		
		return(this->velocity[elemId][numberInSelf]);
	}

	InitialBallProblem(VC& variables, const G& g, TwoPhaseRelations& law = *(new LinearLaw), 
								const int level = 0) 
	: TransportProblem<G, RT, VC>(variables, law), left((g.lowerLeft())[0]), right((g.upperRight())[0]), 
	  bottom((g.lowerLeft())[1]), top((g.upperRight())[1]),
	  elementmapper(g, g.levelIndexSet(level))
	{	
		this->variables.velocity.resize(elementmapper.size());
	}
  };

}
#endif
