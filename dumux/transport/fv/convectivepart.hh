// $Id:

#ifndef DUNE_CONVECTIVEPART_HH
#define DUNE_CONVECTIVEPART_HH

namespace Dune
{
	template<class G, class RT>
	class ConvectivePart
	{
	private:
		enum{dim = G::dimension};
		typedef typename G::Traits::template Codim<0>::Entity Entity;
		typedef FieldVector<RT, dim> Vector;

	public:

		virtual double operator() (const Entity& entity, const RT satI, Vector faceGlobal) const
		{
			double trivial(0);
			return trivial;
		}

	  virtual ~ConvectivePart()
	  { }
	};
}

#endif
