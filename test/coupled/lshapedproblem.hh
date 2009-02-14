// $Id: yxproblem.hh 733 2008-10-24 08:45:27Z bernd $

#ifndef DUNE_LSHAPEDPROBLEM_HH
#define DUNE_LSHAPEDPROBLEM_HH

#include"dumux/stokes/stokesproblem.hh"

namespace Dune {

template<class Grid, class Scalar>
class LShapedProblem : public StokesProblem<Grid, Scalar>
{
  typedef typename Grid::ctype Scalar;
  enum {dim=Grid::dimension, numEq=Grid::dimension+1};
  typedef typename Grid::Traits::template Codim<0>::Entity Element;
  typedef typename IntersectionIteratorGetter<Grid,LeafTag>::IntersectionIterator IntersectionIterator;

public:
  virtual FieldVector<Scalar,numEq> q(const FieldVector<Scalar,dim>& globalPos, const Element& element,
				const FieldVector<Scalar,dim>& localPos) const
  {
    FieldVector<Scalar,numEq> result(0);

    return result;
  }

  virtual BoundaryConditions::Flags bctype (const FieldVector<Scalar,dim>& globalPos, const Element& element,
					    const IntersectionIterator& intersectionIt,
					    const FieldVector<Scalar,dim>& localPos) const
  {
    if (globalPos[0] > 1.5 - 1e-6 && globalPos[1] < 0.5 + 1e-6)
      return BoundaryConditions::process;
   	if (globalPos[0] > 4 - 1e-6)
      return BoundaryConditions::neumann;
    else
      return BoundaryConditions::dirichlet;
  }

  virtual FieldVector<Scalar,dim> g(const FieldVector<Scalar,dim>& globalPos, const Element& element,
				const IntersectionIterator& intersectionIt,
				const FieldVector<Scalar,dim>& localPos) const
  {
    return velocity(globalPos);
  }

  // function returns the square root of the permeability (scalar) divided by alpha
  virtual Scalar beaversJosephC(const FieldVector<Scalar,dim>& globalPos, const Element& element,
		const IntersectionIterator& intersectionIt,
		const FieldVector<Scalar,dim>& localPos) const
  {
      // CHANGE also in the porous medium problem!
	  double permeability = 1.0e-2;//5.88e-5;
	  double alpha = 1.0;

	  //TODO: divide by viscosity - check?
	  return sqrt(permeability)/alpha;
  }

  //TODO: this is not called yet
  virtual Scalar mu(const FieldVector<Scalar,dim>& globalPos, const Element& element, const FieldVector<Scalar,dim>& localPos) const
  {
    return 0.01;
  }

  virtual FieldVector<Scalar,dim> velocity(const FieldVector<Scalar,dim>& globalPos) const
  {
	  FieldVector<Scalar,dim> result(0);

	  if (globalPos[0] < 1e-6)
		  result[0] = 1;//-globalPos[1];
	  if (globalPos[1] < 1e-6)
		  result[0] = 0;//-globalPos[1];
	  if (globalPos[1] > 1-1e-6)
		  result[0] = 0;//-globalPos[1];

    return result;
  }

  virtual Scalar pressure(const FieldVector<Scalar,dim>& globalPos) const
  {
    return (globalPos[0]);
  }

  virtual FieldMatrix<Scalar, dim, dim> velocityGradient(const FieldVector<Scalar,dim>& globalPos) const
  {
    FieldMatrix<Scalar, dim, dim> result(0);
//    result[0][1] = -1.0;
//    result[1][0] = -1.0;
    result[0][0] = -1.0;

    return result;
  }

  LShapedProblem()
  {}

};

}
#endif
