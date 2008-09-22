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
#include<dune/disc/functions/p1function.hh>

#include<dumux/operators/boxjacobian.hh>
#include "diffusionparameters.hh"

namespace Dune
{
  //! A class for computing local jacobian matrices
  /*! A class for computing local jacobian matrix for the
	diffusion equation

	    div j = q; j = -K grad u; in Omega

		u = g on Gamma1; j*n = J on Gamma2.

	Uses the box method.

	Template parameters are:

	- G     a DUNE grid type
	- RT    type used for return values
  */
  template<class G, class RT, class BoxFunction = LeafP1Function<G, RT, 1> >
  class BoxDiffusionJacobian
    : public BoxJacobian<BoxDiffusionJacobian<G,RT,BoxFunction>,G,RT,1,BoxFunction>
  {
    typedef typename G::ctype DT;
    typedef typename G::Traits::template Codim<0>::Entity Entity;
    typedef typename Entity::Geometry Geometry;
    typedef BoxDiffusionJacobian<G,RT,BoxFunction> ThisType;
    typedef typename LocalJacobian<ThisType,G,RT,1>::VBlockType VBlockType;
 	typedef FVElementGeometry<G> FVElementGeometry;
 	enum {pIdx = 0};

  public:
    enum {n=G::dimension};
    enum {SIZE=LagrangeShapeFunctionSetContainer<DT,RT,n>::maxsize};
    struct VariableNodeData;

    //! Constructor
    BoxDiffusionJacobian (DiffusionParameters<G,RT>& params,
			      bool levelBoundaryAsDirichlet_, const G& grid,
			      BoxFunction& sol,
			      bool procBoundaryAsDirichlet_=true)
    : BoxJacobian<ThisType,G,RT,1,BoxFunction>(levelBoundaryAsDirichlet_, grid, sol, procBoundaryAsDirichlet_),
      problem(params), varNData(SIZE), oldVarNData(SIZE)
    {
      this->analytic = false;
    }

    void clearVisited ()
    {
    	return;
    }


    VBlockType computeM (const Entity& e, const VBlockType* sol, int node, const std::vector<VariableNodeData>& varData)
    {
   	 VBlockType result(0);
   	 result[0] = varData[node].density[pIdx]*elData.porosity;

   	 return result;
    };

    VBlockType computeM (const Entity& e, const VBlockType* sol, int node, bool old = false)
    {
    	if (old)
    		return computeM(e, sol, node, oldVarNData);
    	else
    		return computeM(e, sol, node, varNData);
    }


    VBlockType computeQ (const Entity& e, const VBlockType* sol, const int& node)
    {
   	 return problem.q(this->fvGeom.subContVol[node].global, e, this->fvGeom.subContVol[node].local);
    }

    VBlockType computeA (const Entity& e, const VBlockType* sol, int face)
    {
    	FieldVector<RT, n> gradP(0);
    	for (int k = 0; k < this->fvGeom.nNodes; k++) {
    		FieldVector<RT,n> grad(this->fvGeom.subContVolFace[face].grad[k]);
    		grad *= sol[k];
    		gradP += grad;
    	}

		int i = this->fvGeom.subContVolFace[face].i;
		int j = this->fvGeom.subContVolFace[face].j;
		RT densityFace = (varNData[i].density[pIdx] + varNData[j].density[pIdx]) / 2.0;

		FieldVector<RT,n> gravity = problem.gravity();

		gravity *= (-densityFace);

		gradP += gravity;

    	FieldVector<RT,n> KGradP(0);
    	elData.K.umv(gradP, KGradP);

    	KGradP *= densityFace;

    	VBlockType flux = KGradP*this->fvGeom.subContVolFace[face].normal;

    	return flux;
    }

    void computeElementData (const Entity& e)
    {
    	// ASSUMING element-wise constant permeability, evaluate K at the cell center
    	elData.K = problem.K(this->fvGeom.cellGlobal, e, this->fvGeom.cellLocal);
    	elData.porosity = problem.porosity(this->fvGeom.cellGlobal, e, this->fvGeom.cellLocal);
    };

	  // *********************************************************
	  // *																			*
	  // *	Calculation of Data at Nodes that has to be			 	*
	  // *	determined only once	(statNData)							 	*
	  // *																		 	*
	  // *********************************************************

  // analog to EvalStaticData in MUFTE
  void updateStaticData (const Entity& e, const VBlockType* sol)
  {
	  return;
  }


	  //*********************************************************
	  //*														*
	  //*	Calculation of variable Data at Nodes			 	*
	  //*	(varNData)										 	*
	  //*													 	*
	  //*********************************************************


    // analog to EvalPrimaryData in MUFTE, uses members of varNData
	void updateVariableData(const Entity& e, const VBlockType* sol,
			int i, std::vector<VariableNodeData>& varData)
    {
   		varData[i].density[pIdx] = problem.materialLaw().density(problem.Temp(), sol[i][pIdx]);
    }

	void updateVariableData(const Entity& e, const VBlockType* sol, int i, bool old = false)
	{
		if (old) {
			updateVariableData(e, sol, i, oldVarNData);
		}
		else
			updateVariableData(e, sol, i, varNData);
	}

	void updateVariableData(const Entity& e, const VBlockType* sol, bool old = false)
	{
		int size = this->fvGeom.nNodes;

		for (int i = 0; i < size; i++)
				updateVariableData(e, sol, i, old);
	}

    struct StaticNodeData
    {
    	bool visited;
    };
    // the members of the struct are defined here
    struct VariableNodeData
    {
    	VBlockType density;
    };

	struct ElementData {
		RT porosity;
    	FieldMatrix<DT,n,n> K;
   	};

   	ElementData elData;
    DiffusionParameters<G,RT>& problem;
    std::vector<StaticNodeData> statNData;
    std::vector<VariableNodeData> varNData;
    std::vector<VariableNodeData> oldVarNData;


  };
}
#endif
