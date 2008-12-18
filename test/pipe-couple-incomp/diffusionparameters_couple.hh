#ifndef DIFFUSIONPARAMETERS_HH
#define DIFFUSIONPARAMETERS_HH

//#include "dumux/material/properties.hh"

template<class G, class RT, class GlobalToPipeMapper, class VertexMapper, class VertexVectorOnLineType>
class DiffusionParameters
{
	typedef typename G::ctype DT;
	enum {n=G::dimension, m=1};
	typedef typename G::Traits::template Codim<0>::Entity Entity;
	typedef typename Dune::IntersectionIteratorGetter<G,Dune::LeafTag>::IntersectionIterator IntersectionIterator;
	typedef Dune::BlockVector<Dune::FieldVector<RT,1> > BlockVector;

public:
	double alphaEx;
	BlockVector& pipePressure;
	VertexVectorOnLineType vertexVectorOnLine;
	GlobalToPipeMapper& mapGlobalIDtoPipeID;
	VertexMapper& vertexMapper;
	Medium& medium;

	DiffusionParameters (double alpha, BlockVector& pipePress, VertexVectorOnLineType *vertexVectorOnL,
			GlobalToPipeMapper& map, VertexMapper& vm, Medium& med)
	: alphaEx(alpha), pipePressure(pipePress),
	mapGlobalIDtoPipeID(map), vertexMapper(vm), medium(med)
	{
		for (int i=0; i<n; i++)
			for (int j=0; j<n; j++)
				if (i==j)
					large[i][j] = 1e-10/medium.viscosity();
				else
					large[i][j] = 0/medium.viscosity();

		vertexVectorOnLine = *vertexVectorOnL;
	}

	const Dune::FieldMatrix<DT,n,n>& K (const Dune::FieldVector<DT,n>& x, const Entity& e,
			const Dune::FieldVector<DT,n>& xi) const
			{
		return large;
			}

	RT q   (const Dune::FieldVector<DT,n>& x, const Entity& e,
			const Dune::FieldVector<DT,n>& xi, const int& node) const
			{
		int globalId = vertexMapper.template map<n>(e, node);
		for (unsigned k = 0; k < vertexVectorOnLine.size(); k++)
		{
			if (globalId == vertexVectorOnLine[k].globalId)
			{
				int pipeId = ((mapGlobalIDtoPipeID).find(globalId))->second;
				return alphaEx * pipePressure[pipeId];
			}
		}

		return 0;
			}

	Dune::BoundaryConditions::Flags bctype (const Dune::FieldVector<DT,n>& x, const Entity& e,
			const IntersectionIterator& intersectionIt,
			const Dune::FieldVector<DT,n>& xi) const
			{
		Dune::FieldVector<Dune::BoundaryConditions::Flags, 1> values; //(Dune::BoundaryConditions::neumann);

		//		std::cout << "global coordinate "<< x << "bctype: boundaryId = " << intersectionIt.boundaryId() << std::endl;

		switch (intersectionIt.boundaryId()) {
		case 1: case 2: case 3: case 4: case 6:
			values = Dune::BoundaryConditions::neumann;
			break;
		case 5:
			values = Dune::BoundaryConditions::dirichlet;
			break;
		}

		return values;
			}

	virtual void dirichletIndex(const Dune::FieldVector<DT,n>& x, const Entity& e,
			const IntersectionIterator& intersectionIt,
			const Dune::FieldVector<DT,n>& xi, Dune::FieldVector<int,m>& dirichletIndex) const
			{
		for (int i = 0; i < m; i++)
			dirichletIndex[i]=i;
		return;
			}

	RT exact(const Dune::FieldVector<DT,n>& x) const
	{
		return (x[0]);//1.0 - x[0] + 2.0*x[1] - 3.0*x[2]);
	}

	RT g (const Dune::FieldVector<DT,n>& x, const Entity& e,
			const IntersectionIterator& intersectionIt,
			const Dune::FieldVector<DT,n>& xi) const
			{
		RT values(0);

		switch (intersectionIt.boundaryId()) {
		case 5:
			values= 1.0e+5;
			break;
		}

		return values;
			}

	Dune::FieldVector<RT,n> exactGrad(const Dune::FieldVector<DT,n>& x) const
	{
		Dune::FieldVector<DT,n> grad;

		grad[0] = -1.0;
		grad[1] = 0.0;
		grad[2] = 0.0;

		return grad;
	}

	RT J (const Dune::FieldVector<DT,n>& x, const Entity& e,
			const IntersectionIterator& intersectionIt,
			const Dune::FieldVector<DT,n>& xi) const
			{
		RT values(0);

		switch (intersectionIt.boundaryId()) {
		case 1: case 2: case 3: case 4: case 6:
			values = 0;
			break;
		}

		return values;
			}

private:
	Dune::FieldMatrix<DT,n,n> large;
};

#endif
