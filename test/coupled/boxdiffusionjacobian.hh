#ifndef DUNE_BOXDIFFUSIONJACOBIAN_HH
#define DUNE_BOXDIFFUSIONJACOBIAN_HH

#include<map>
#include<iostream>
#include<iomanip>
#include<fstream>
#include<vector>
#include<sstream>

#include<dune/common/exceptions.hh>
#include<dune/grid/common/grid.hh>
#include<dune/grid/common/referenceelements.hh>
#include<dune/common/geometrytype.hh>
#include<dune/grid/common/quadraturerules.hh>
#include <dune/grid/utility/intersectiongetter.hh>

#include<dune/disc/shapefunctions/lagrangeshapefunctions.hh>
#include<dune/disc/operators/boundaryconditions.hh>

#include<dumux/operators/boxjacobian.hh>
//#include "diffusionparameters.hh"
#include "darcyparameters.hh"

namespace Dune
{
//! A class for computing local jacobian matrices
/*! A class for computing local jacobian matrix for the
    diffusion equation

        div j = q; j = -K grad u; in Omega

        u = g on Gamma1; j*n = J on Gamma2.

    Uses the box method.

    Template parameters are:

    - Grid     a DUNE grid type
    - Scalar    type used for return values
 */
template<class Grid, class Scalar, class BoxFunction = LeafP1Function<Grid, Scalar, 1> >
class BoxDiffusionJacobian
: public BoxJacobian<BoxDiffusionJacobian<Grid,Scalar,BoxFunction>,Grid,Scalar,1,BoxFunction>
{
	typedef typename Grid::Traits::template Codim<0>::Entity Element;
	typedef typename Element::Geometry Geometry;
	typedef BoxDiffusionJacobian<Grid,Scalar,BoxFunction> ThisType;
	typedef typename LocalJacobian<ThisType,Grid,Scalar,1>::VBlockType SolutionVector;
	typedef Dune::FVElementGeometry<Grid> FVElementGeometry;

public:
	enum {dim=Grid::dimension};

	//! Constructor
	BoxDiffusionJacobian (DarcyParameters<Grid,Scalar>& params,
			bool levelBoundaryAsDirichlet_, const Grid& grid,
			BoxFunction& sol,
			bool procBoundaryAsDirichlet_=true)
	: BoxJacobian<ThisType,Grid,Scalar,1,BoxFunction>(levelBoundaryAsDirichlet_, grid, sol, procBoundaryAsDirichlet_),
	problem(params)
	{
		this->analytic = false;
	}

	void clearVisited ()
	{
		return;
	}

	SolutionVector computeM (const Element& e, const SolutionVector* sol, int node, bool old = false)
	{
		SolutionVector result(0);
		return result;
	}

	SolutionVector computeA (const Element& e, const SolutionVector* sol, int face)
	{
		FieldVector<Scalar, dim> gradP(0);
		for (int k = 0; k < this->fvGeom.numVertices; k++) {
			FieldVector<Scalar,dim> grad(this->fvGeom.subContVolFace[face].grad[k]);
			grad *= sol[k];
			gradP += grad;
		}

		//TODO: implement correct viscosity, called from material law
		double viscosity(0.01);
		FieldVector<Scalar,dim> KGradP(0);
		elData.K.umv(gradP, KGradP);

		SolutionVector flux = KGradP*this->fvGeom.subContVolFace[face].normal;
		flux /= viscosity;

		return flux;
	}

	SolutionVector computeQ (const Element& e, const SolutionVector* sol, const int& node)
	{
		return problem.q(this->fvGeom.subContVol[node].global, e, this->fvGeom.subContVol[node].local);
	}

	void computeElementData (const Element& e)
	{
		// ASSUMING element-wise constant permeability, evaluate K at the cell center
		elData.K = problem.K(this->fvGeom.elementGlobal, e, this->fvGeom.elementLocal);
	};

	virtual void updateVariableData(const Element& e, const SolutionVector* sol, int i, bool old = false)
	{
		return;
	}

	void updateVariableData(const Element& e, const SolutionVector* sol, bool old = false)
	{
		return;
	}

	virtual void updateStaticData (const Element& e, const SolutionVector* sol)
	{
		return;
	}

	struct ElementData {
		FieldMatrix<Scalar,dim,dim> K;
	};

	ElementData elData;
	DarcyParameters<Grid,Scalar>& problem;
};
}
#endif
