#ifndef DARCYPARAMETERS_HH
#define DARCYPARAMETERS_HH

template<class Grid, class Scalar>
class DarcyParameters
{
  typedef typename Grid::ctype Scalar;
  enum {dim=Grid::dimension, numEq=1};
  typedef typename Grid::Traits::template Codim<0>::Entity Element;
  typedef typename Dune::IntersectionIteratorGetter<Grid,Dune::LeafTag>::IntersectionIterator IntersectionIterator;

public:
  DarcyParameters ()
  {
      // CHANGE also in the Stokes problem!
	for (int i=0; i<dim; i++)
	  for (int j=0; j<dim; j++)
		if (i==j)
		  permeability_[i][j] = 1.0e-2;//5.88e-5;
		else
			permeability_[i][j] = 0;
  }

  const Dune::FieldMatrix<Scalar,dim,dim>& K (const Dune::FieldVector<Scalar,dim>& globalPos, const Element& element,
				  const Dune::FieldVector<Scalar,dim>& localPos) const
  {
	return permeability_;
  }

  Scalar q   (const Dune::FieldVector<Scalar,dim>& globalPos, const Element& element,
				  const Dune::FieldVector<Scalar,dim>& localPos) const
  {
	return 0;
  }

  Dune::BoundaryConditions::Flags bctype (const Dune::FieldVector<Scalar,dim>& globalPos, const Element& element,
			const IntersectionIterator& intersectionIt,
					   const Dune::FieldVector<Scalar,dim>& localPos) const
  {
 	if (globalPos[0] > 4 - 1e-6)
	  return Dune::BoundaryConditions::dirichlet;

	return Dune::BoundaryConditions::neumann;
  }

	virtual void dirichletIndex(const Dune::FieldVector<Scalar,dim>& globalPos, const Element& element,
			const IntersectionIterator& intersectionIt,
			const Dune::FieldVector<Scalar,dim>& localPos, Dune::FieldVector<int,numEq>& dirichletIndex) const
	{
		for (int i = 0; i < numEq; i++)
			dirichletIndex[i]=i;
		return;
	}

  Scalar exact(const Dune::FieldVector<Scalar,dim>& globalPos) const
  {
		return (globalPos[0]*globalPos[1]);
  }

  // Dirichlet boundary conditions
  Scalar g (const Dune::FieldVector<Scalar,dim>& globalPos, const Element& element,
			const IntersectionIterator& intersectionIt,
				const Dune::FieldVector<Scalar,dim>& localPos) const
  {
	return 0;
  }

  Dune::FieldVector<Scalar,dim> exactGrad(const Dune::FieldVector<Scalar,dim>& globalPos) const
  {
	  Dune::FieldVector<Scalar,dim> grad;

	  grad[0] = globalPos[1];
	  grad[1] = globalPos[0];
	  grad[2] = 0.0;

	  return grad;
  }

  // Neumann b.c.
  Scalar J (const Dune::FieldVector<Scalar,dim>& globalPos, const Element& element,
			const IntersectionIterator& intersectionIt,
				const Dune::FieldVector<Scalar,dim>& localPos) const
  {
	  return 0;
//	  if (globalPos[0] < 0.5 + 1e-6 && globalPos[1] > 0.4 - 1e-6)
//		  return 0;
//	  else
//	  {
//		  Dune::FieldVector<Scalar,dim> KGradU(0);
//		  large.umv(exactGrad(globalPos), KGradU);
//
//		  // ASSUMING face-wise constant normal
//		  Dune::FieldVector<Scalar, dim-1> localDimM1(0);
//		  return -(KGradU*intersectionIt->unitOuterNormal(localDimM1));
//	  }
  }

private:
  Dune::FieldMatrix<Scalar,dim,dim> permeability_;
};

#endif
