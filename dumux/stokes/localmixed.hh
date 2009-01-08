// $Id$

#ifndef DUNE_LOCALMIXED_HH
#define DUNE_LOCALMIXED_HH

#include<dune/grid/utility/intersectiongetter.hh>
#include"dumux/operators/localstiffnessextended.hh"
#include"dumux/stokes/stokesproblem.hh"

namespace Dune
{
template<class Grid, class Scalar, int m>
class LocalMixed : public LocalStiffness< LocalMixed<Grid,Scalar,m>, Grid, Scalar, m>
{
	typedef LocalStiffness< LocalMixed<Grid,Scalar,m>, Grid, Scalar, m> Stiffness;
	typedef typename Grid::Traits::template Codim<0>::Entity Element;
	typedef typename Element::Geometry Geometry;
	typedef typename IntersectionIteratorGetter<Grid,LeafTag>::IntersectionIterator IntersectionIterator;
	typedef StokesProblem<Grid, Scalar> Problem;
	enum {dim=Grid::dimension};

public:
	LocalMixed(Problem& problem)
	: Stiffness(), problem_(problem)
	{}

	// ASSUME a 2D uniform rectangular axiparallel grid
	template<class TypeTag>
	void assemble (const Element& element, int k = 1)
	{
		// get the number of element faces:
		int nFaces = element.template count<1>();

		// initialize everything with 0:
		for (int i = 0; i < nFaces+1; i++) {
			this->bctype[i] = BoundaryConditions::neumann;
			this->b[i] = 0;
			for (int j = 0; j < nFaces+1; j++)
				this->A[i][j] = 0;
		}

		// obtain the inverse of the Jacobian of the element transformation:
		const Geometry& geometry = element.geometry();
		GeometryType geomType = geometry.type();
		const FieldVector<Scalar,dim>& local = ReferenceElements<Scalar,dim>::general(geomType).position(0,0);
		const FieldMatrix<Scalar,dim,dim>& jacobianInv = geometry.jacobianInverseTransposed(local);

		// extract the discretization parameters:
		Scalar deltaX = 1.0/jacobianInv[0][0];
		Scalar deltaY = 1.0/jacobianInv[1][1];
		Scalar volume = deltaX*deltaY;

		// mass conservation
		// control volume: element
		// div v = 0 -> \sum_e |e| v.n = 0
		this->A[4][0] = deltaX;
		this->A[4][1] = -deltaX;
		this->A[4][2] = deltaY;
		this->A[4][3] = -deltaY;

		// momentum conservation
		// - div T + grad p
		// x-comp control volume: element shifted by deltaX/2 in x-direction
		// y-comp control volume: element shifted by deltaY/2 in y-direction
		// -> i-comp: \sum_e |e|( -T_i.n + p.n_i )
		Scalar mu = 1.0;
		this->A[0][0] = this->A[1][1] = mu*deltaX/deltaY;
		this->A[2][2] = this->A[3][3] = mu*deltaY/deltaX;
		this->A[0][1] = this->A[1][0] = -mu*deltaX/deltaY;
		this->A[2][3] = this->A[3][2] = -mu*deltaY/deltaX;
		this->A[0][4] = deltaX;
		this->A[1][4] = -deltaX;
		this->A[2][4] = deltaY;
		this->A[3][4] = -deltaY;

		// insert the source term and the boundary conditions:
		IntersectionIterator endIsIt = element.ileafend();
		for (IntersectionIterator isIt = element.ileafbegin(); isIt != endIsIt; ++isIt)
		{
			// the local face number:
			int numberInSelf = isIt->numberInSelf();

			// geometry data of the face:
			GeometryType geomType = isIt->intersectionSelfLocal().type();
			const FieldVector<Scalar,dim-1>& localDimM1 = ReferenceElements<Scalar,dim-1>::general(geomType).position(0,0);
			const FieldVector<Scalar,dim>& local = ReferenceElements<Scalar,dim>::general(geomType).position(numberInSelf, 1);
			FieldVector<Scalar,dim> global = isIt->intersectionGlobal().global(localDimM1);

			// get the source term from the problem:
			FieldVector<Scalar,dim> source = problem_.q(global, element, local);

			// set the right hand side:
			this->b[numberInSelf] = 0.5*volume*source[numberInSelf/2];

			// after this, only boundary elements should be considered:
			if (isIt->neighbor())
				continue;

			// get the boundary condition type:
			BoundaryConditions::Flags bctype = problem_.bctype(global, element, isIt, local);

			if (bctype == BoundaryConditions::dirichlet)
			{
				// get the boundary value from the problem:
				FieldVector<Scalar,dim> velocity = problem_.g(global, element, isIt, local);

				// set the right hand side and the boundary condition type:
				this->b[numberInSelf] = velocity[numberInSelf/2];
				this->bctype[numberInSelf] = BoundaryConditions::dirichlet;
			}
			else
			{
				DUNE_THROW(NotImplemented, "BoundaryConditions != Dirichlet for LocalMixed");
			}
		}

		return;
	}

private:
	Problem& problem_;
};

}
#endif
