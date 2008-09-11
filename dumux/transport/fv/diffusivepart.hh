// $Id$ 

#ifndef DUNE_DIFFUSIVEPART_HH
#define DUNE_DIFFUSIVEPART_HH

//! \ingroup transport
//! \defgroup diffPart Diffusive transport
/**
 * @file
 * @brief  Base class for defining the diffusive part of an advection-diffusion equation
 * @author Bernd Flemisch
 */
namespace Dune
{
	/*!\ingroup diffPart
	 * @brief  Base class for defining the diffusive part of an advection-diffusion equation 
	 */
	template<class G, class RT>
	class DiffusivePart 
	{
	private:
		enum{dim = G::dimension};	
		typedef typename G::Traits::template Codim<0>::Entity Entity;
		typedef FieldVector<RT, dim> FieldVector;
		
	public:
		virtual FieldVector operator() (const Entity& entity, const int numberInSelf, 
						const RT satIntersection, const FieldVector& satGradient, const RT time) const
		{
			FieldVector trivial(0);
			return trivial;
		}

		virtual FieldVector operator() (const Entity& entity, const int numberInSelf, 
						const RT satIntersection, const FieldVector& satGradient, const RT time, 
						RT satI, RT satJ) const
		{
			FieldVector trivial(0);
			return trivial;
		}

	  virtual ~DiffusivePart()
	  { }
	};
}

#endif
