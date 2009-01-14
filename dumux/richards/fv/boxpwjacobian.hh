// $Id: boxpwsnjacobian.hh 599 2008-09-18 13:57:47Z bernd $

#ifndef DUNE_BOXPWSNJACOBIAN_HH
#define DUNE_BOXPWSNJACOBIAN_HH

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
#include"dumux/operators/boxjacobian.hh"
#include"dumux/richards/richardsproblem.hh"

/**
 * @file
 * @brief  compute local jacobian matrix for conforming finite elements for diffusion equation
 * @author Peter Bastian
 */



namespace Dune
{
  /** @addtogroup DISC_Disc
   *
   * @{
   */
  /**
   * @brief compute local jacobian matrix for conforming finite elements for diffusion equation
   *
   */


  //! A class for computing local jacobian matrices
  /*! A class for computing local jacobian matrix for the
	diffusion equation

	    div j = q; j = -K grad u; in Omega

		u = g on Gamma1; j*n = J on Gamma2.

	Uses conforming finite elements with the Lagrange shape functions.
	It should work for all dimensions and element types.
	All the numbering is with respect to the reference element and the
	Lagrange shape functions

	Template parameters are:

	- Grid  a DUNE grid type
	- RT    type used for return values
  */
  template<class G, class RT, class BoxFunction = LeafP1FunctionExtended<G, RT, 1> >
  class BoxPwJacobian
    : public BoxJacobian<BoxPwJacobian<G,RT,BoxFunction>,G,RT,1,BoxFunction>
  {
    typedef typename G::ctype DT;
    typedef typename G::Traits::template Codim<0>::Entity Entity;
    typedef typename Entity::Geometry Geometry;
    typedef BoxPwJacobian<G,RT,BoxFunction> ThisType;
    typedef typename LocalJacobian<ThisType,G,RT,1>::VBlockType VBlockType;
    typedef typename LocalJacobian<ThisType,G,RT,1>::MBlockType MBlockType;
	enum {pWIdx = 0};


  public:
	    // define the number of components of your system, this is used outside
	    // to allocate the correct size of (dense) blocks with a FieldMatrix
	    enum {dim=G::dimension};
	    enum {m=1};
	    enum {SIZE=LagrangeShapeFunctionSetContainer<DT,RT,dim>::maxsize};
	    struct VariableNodeData;
	    typedef FieldMatrix<RT,dim,dim> FMatrix;
	    typedef FieldVector<RT,dim> FVector;

	     //! Constructor
    BoxPwJacobian (TwoPhaseProblem<G,RT>& params,
			      bool levelBoundaryAsDirichlet_, const G& grid,
			      BoxFunction& sol,
			      bool procBoundaryAsDirichlet_=true)
    : BoxJacobian<ThisType,G,RT,1,BoxFunction>(levelBoundaryAsDirichlet_, grid, sol, procBoundaryAsDirichlet_),
      problem(params),
      statNData(this->vertexMapper.size()), varNData(SIZE), oldVarNData(SIZE)
    {
      this->analytic = false;
    }

    void clearVisited ()
    {
    	return;
    }

    VBlockType computeM (const Entity& e, const VBlockType* sol,
    		int node, const std::vector<VariableNodeData>& varData)
    {
   	 VBlockType result;
   	 //std::cout << "rhoW = " << varData[node].density[pWIdx] << ", rhoN = " << varData[node].density[satNIdx] << std::endl;
   	 result[0] = elData.porosity*varData[node].density[pWIdx]*varData[node].dSdpC*varData[node].pC;//varData[node].density[pWIdx]*elData.porosity*varData[node].satW;//sol[node][satWIdx];

   	 return result;
    };

    VBlockType computeM (const Entity& e, const VBlockType* sol, int node, bool old = false)
    {
    	if (old)
    		return computeM(e, sol, node, oldVarNData);
    	else
    		return computeM(e, sol, node, varNData);
    }

    VBlockType computeA (const Entity& e, const VBlockType* sol, int face)
    {
		 int i = this->fvGeom.subContVolFace[face].i;
     	 int j = this->fvGeom.subContVolFace[face].j;

		  // permeability in edge direction
		  FieldVector<DT,dim> Kij(0);
		  elData.K.umv(this->fvGeom.subContVolFace[face].normal, Kij);


		  VBlockType flux;
		  for (int phase = 0; phase < m; phase++) {
	          // calculate FE gradient
	          FieldVector<RT, dim> pGrad(0);
	          for (int k = 0; k < this->fvGeom.nNodes; k++) {
	        	  FieldVector<DT,dim> grad(this->fvGeom.subContVolFace[face].grad[k]);
	        	  grad *= varNData[k].pW ;  //(phase) ? varNData[k].pN : sol[k][pWIdx];
	        	  pGrad += grad;
	          }

	          // adjust by gravity
			  FieldVector<RT, dim> gravity = problem.gravity();
			  gravity *= varNData[i].density[phase];
			  pGrad -= gravity;

			  // calculate the flux using upwind
			  RT outward = pGrad*Kij;
			  if (outward < 0)
				  flux[phase] = varNData[i].density[phase]*varNData[i].mobility[phase]*outward;
			  else
				  flux[phase] = varNData[j].density[phase]*varNData[j].mobility[phase]*outward;
		  }

		  return flux;
   };


    VBlockType computeQ (const Entity& e, const VBlockType* sol, const int& node)
    {
   	 // ASSUME problem.q already contains \rho.q
   	 return problem.q(this->fvGeom.subContVol[node].global, e, this->fvGeom.subContVol[node].local);
    };


	  // *********************************************************
	  // *														 *
	  // *	Calculation of Data at Elements			 			 *
	  // *						 								 *
	  // *														 *
	  // *********************************************************

    void computeElementData (const Entity& e)
    {
		 // ASSUMING element-wise constant permeability and porosity, evaluate at the cell center
 		 elData.K = this->problem.soil().K(this->fvGeom.cellGlobal, e, this->fvGeom.cellLocal);
 		 elData.porosity = this->problem.soil().porosity(this->fvGeom.cellGlobal, e, this->fvGeom.cellLocal);
    };


	  // *********************************************************
	  // *														 *
	  // *	Calculation of Data at Nodes that has to be			 *
	  // *	determined only once	(statNData)					 *
	  // *														 *
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


    // the members of the struct are defined here
    struct VariableNodeData
    {
    	RT satW;
    	RT pC;
    	RT pW;
    	RT dSdpC;
    	VBlockType mobility;  //Vector with the number of phases
    	VBlockType density;
    	FieldMatrix<DT,dim,dim> K;
    };

    // analog to EvalPrimaryData in MUFTE, uses members of varNData
	void updateVariableData(const Entity& e, const VBlockType* sol,
			int i, std::vector<VariableNodeData>& varData)
    {
//		this->fvGeom.update(e);
		FVector& global = this->fvGeom.subContVol[i].global;
  	   	FVector& local = this->fvGeom.subContVol[i].local;

 //		varData[i].K = problem.K(this->fvGeom.subContVol[i].global, e, this->fvGeom.subContVol[i].local);
  	   	varData[i].pW = sol[i][pWIdx];

  	   	/* pc = -pw || pc = 0 for computing Sw */
  	   	if (varData[i].pW >= 0)
  	   	{
  	   		varData[i].pC = 0.0;
  	   	}
  	   	else
  	   	{
  	   		varData[i].pC = -varData[i].pW;
  	   	}

   		varData[i].dSdpC = problem.materialLaw().dSdP(varData[i].pC, global, e, local);
   		varData[i].satW = problem.materialLaw().saturationW(varData[i].pC, global, e, local);
   		varData[i].mobility[pWIdx] = problem.materialLaw().mobW(varData[i].satW, global, e, local);
   		varData[i].density[pWIdx] = problem.wettingPhase().density(283.15, varData[i].pW);
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

    struct ElementData {
   	 RT cellVolume;
     RT porosity;
   	 RT gravity;
   	 FieldVector<RT, 4> parameters;
   	 FieldMatrix<DT,dim,dim> K;
   	 } elData;

     // parameters given in constructor
     TwoPhaseProblem<G,RT>& problem;
     std::vector<StaticNodeData> statNData;
     std::vector<VariableNodeData> varNData;
     std::vector<VariableNodeData> oldVarNData;
  };

  /** @} */
}
#endif
