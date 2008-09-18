// $Id$

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
#include"dumux/twophase/twophaseproblem.hh"

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
  template<class G, class RT, class BoxFunction = LeafP1FunctionExtended<G, RT, 2> >
  class BoxPwSnJacobian
    : public BoxJacobian<BoxPwSnJacobian<G,RT,BoxFunction>,G,RT,2,BoxFunction>
  {
    typedef typename G::ctype DT;
    typedef typename G::Traits::template Codim<0>::Entity Entity;
    typedef typename Entity::Geometry Geometry;
    typedef BoxPwSnJacobian<G,RT,BoxFunction> ThisType;
    typedef typename LocalJacobian<ThisType,G,RT,2>::VBlockType VBlockType;
    typedef typename LocalJacobian<ThisType,G,RT,2>::MBlockType MBlockType;
	enum {pWIdx = 0, satNIdx = 1};


  public:
	    // define the number of components of your system, this is used outside
	    // to allocate the correct size of (dense) blocks with a FieldMatrix
	    enum {dim=G::dimension};
	    enum {m=2};
	    enum {SIZE=LagrangeShapeFunctionSetContainer<DT,RT,dim>::maxsize};
	    struct VariableNodeData;
	    typedef FieldMatrix<RT,dim,dim> FMatrix;
	    typedef FieldVector<RT,dim> FVector;

	     //! Constructor
    BoxPwSnJacobian (TwoPhaseProblem<G,RT>& params,
			      bool levelBoundaryAsDirichlet_, const G& grid,
			      BoxFunction& sol,
			      bool procBoundaryAsDirichlet_=true)
    : BoxJacobian<ThisType,G,RT,2,BoxFunction>(levelBoundaryAsDirichlet_, grid, sol, procBoundaryAsDirichlet_),
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
   	 result[0] = varData[node].density[pWIdx]*elData.porosity*varData[node].satW;//sol[node][satNIdx];
   	 result[1] = varData[node].density[satNIdx]*elData.porosity*varData[node].satN;//*sol[node][satNIdx];

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
//		  FieldVector<DT,dim> Ki(0);
//		  varNData[i].K.umv(this->fvGeom.subContVolFace[face].normal, Ki);
//		  FieldVector<DT,dim> Kj(0);
//		  varNData[j].K.umv(this->fvGeom.subContVolFace[face].normal, Kj);
//		  FieldVector<DT,dim> Kij(0);
//		  for (int k = 0; k < dim; k++)
//			  if (fabs(Ki[k] + Kj[k]) > 0)
//				  Kij[k] = 2.0*Ki[k]*Kj[k]/(Ki[k] + Kj[k]);


		  VBlockType flux;
		  for (int phase = 0; phase < m; phase++) {
	          // calculate FE gradient
	          FieldVector<RT, dim> pGrad(0);
	          for (int k = 0; k < this->fvGeom.nNodes; k++) {
	        	  FieldVector<DT,dim> grad(this->fvGeom.subContVolFace[face].grad[k]);
	        	  grad *= (phase) ? varNData[k].pN : sol[k][pWIdx];
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
    	RT satN;
    	RT satW;
    	RT pC;
    	RT pW;
    	RT pN;
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
		varData[i].satN = sol[i][satNIdx];
		varData[i].satW = 1.0 - sol[i][satNIdx];
   		varData[i].pC = problem.materialLaw().pC(varData[i].satW, global, e, local);
		varData[i].pW = sol[i][pWIdx];
   		varData[i].pN = sol[i][pWIdx] + varData[i].pC;
   		varData[i].mobility[pWIdx] = problem.materialLaw().mobW(varData[i].satW, global, e, local);
   		varData[i].mobility[satNIdx] = problem.materialLaw().mobN(varData[i].satN, global, e, local);
   		varData[i].density[pWIdx] = problem.wettingPhase().density(283.15, varData[i].pW);
   		varData[i].density[satNIdx] = problem.nonwettingPhase().density(283.15, varData[i].pN);
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
