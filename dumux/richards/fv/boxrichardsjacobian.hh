// $Id: boxpwsnjacobian.hh 827 2008-11-24 19:49:54Z klaus $

#ifndef DUNE_RICHARDSJACOBIAN_HH
#define DUNE_RICHARDSJACOBIAN_HH

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
   * @brief compute local jacobian matrix for the box method for two-phase equation
   *
   */


  //! A class for computing local jacobian matrices
  /*! A class for computing local jacobian matrix for the
	fully coupled two-phase model with Pw and Sn as primary variables

	Uses the box scheme.
	It should work for all dimensions and element types.
	All the numbering is with respect to the reference element and the
	Lagrange shape functions

	Template parameters are:

	- Grid  a DUNE grid type
	- RT    type used for return values
  */
  template<class G, class RT, class BoxFunction = LeafP1FunctionExtended<G, RT, 2> >
  class BoxRichardsJacobian
    : public BoxJacobian<BoxRichardsJacobian<G,RT,BoxFunction>,G,RT,2,BoxFunction>
  {
    typedef typename G::ctype DT;
    typedef typename G::Traits::template Codim<0>::Entity Entity;
    typedef typename Entity::Geometry Geometry;
    typedef BoxRichardsJacobian<G,RT,BoxFunction> ThisType;
    typedef typename LocalJacobian<ThisType,G,RT,1>::VBlockType VBlockType;
    typedef typename LocalJacobian<ThisType,G,RT,1>::MBlockType MBlockType;

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
    BoxRichardsJacobian (TwoPhaseProblem<G,RT>& params,
			      bool levelBoundaryAsDirichlet_,
			      const G& grid,
			      BoxFunction& sol,
			      bool procBoundaryAsDirichlet_=true)
    : BoxJacobian<ThisType,G,RT,2,BoxFunction>(levelBoundaryAsDirichlet_, grid, sol, procBoundaryAsDirichlet_),
      problem(params), statNData(this->vertexMapper.size()), varNData(SIZE), oldVarNData(SIZE)
    {
      this->analytic = false;
    }

    // compute storage term
    RT computeM (const Entity& e, const VBlockType* sol,
    		int node, const std::vector<VariableNodeData>& varData)
    {
   	 RT result;

   	 GeometryType gt = e.geometry().type();
	 const typename LagrangeShapeFunctionSetContainer<DT,RT,dim>::value_type&
 	 sfs=LagrangeShapeFunctions<DT,RT,dim>::general(gt,1);

	 int globalIdx = this->vertexMapper.template map<dim>(e, sfs[node].entity());

	 //TODO: implement dSwdPc
   	 result = varData[node].density*sNDat[globalIdx].porosity*(-varData[node].dSwdPc)*varData[node].pressure;

   	 return result;
    };

    RT computeM (const Entity& e, const VBlockType* sol, int node, bool old = false)
    {
    	if (old)
    		return computeM(e, sol, node, oldVarNData);
    	else
    		return computeM(e, sol, node, varNData);
    }

    // compute flux term
    RT computeA (const Entity& e, const VBlockType* sol, int face)
    {
      	 int i = this->fvGeom.subContVolFace[face].i;
    	 int j = this->fvGeom.subContVolFace[face].j;

     	 // normal vector, value of the area of the scvf
    	 const FieldVector<RT,dim> normal(this->fvGeom.subContVolFace[face].normal);
    	 GeometryType gt = e.geometry().type();
    	 const typename LagrangeShapeFunctionSetContainer<DT,RT,dim>::value_type&
     	 sfs=LagrangeShapeFunctions<DT,RT,dim>::general(gt,1);

    	 // global index of the subcontrolvolume face neighbor nodes in element e
    	 int globalIdx_i = this->vertexMapper.template map<dim>(e, sfs[i].entity());
      	 int globalIdx_j = this->vertexMapper.template map<dim>(e, sfs[j].entity());

      	 // get global coordinates of nodes i,j
    	 const FieldVector<DT,dim> global_i = this->fvGeom.subContVol[i].global;
    	 const FieldVector<DT,dim> global_j = this->fvGeom.subContVol[j].global;
    	 const FieldVector<DT,dim> local_i = this->fvGeom.subContVol[i].local;
    	 const FieldVector<DT,dim> local_j = this->fvGeom.subContVol[j].local;

    	 // permeability in edge direction
    	 FieldMatrix<RT,dim,dim> Ki(0), Kj(0);
     	 Ki = this->problem.soil().K(global_i,e,local_i);
     	 Kj = this->problem.soil().K(global_j,e,local_j);
    	 const FMatrix K = harmonicMeanK(Ki, Kj);

    	 RT flux;
		 // calculate FE gradient of the pressure
		 FieldVector<RT, dim> pGrad(0);
		 for (int k = 0; k < this->fvGeom.numVertices; k++)
		 {
			  FieldVector<DT,dim> grad(this->fvGeom.subContVolFace[face].grad[k]);
			  grad *= sol[k][pWIdx]; // pressure

			  pGrad += grad;
		 }

		 // adjust by gravity
		 FieldVector<RT, dim> gravity = problem.gravity();
		 gravity *= varNData[i].densityW;
		 pGrad -= gravity; // v = grad P - rho*g

		 RT outward;
		 FieldVector<RT,dim> v_tilde(0);
		 K.mv(pGrad, v_tilde);  // v_tilde=K*gradP
		 outward = v_tilde*normal;

		 // TODO: Weighted updwind, blend function

		 // calculate the flux using fully upwind
		  if (outward < 0)
			  flux[phase] = varNData[i].densityW*varNData[i].mobilityW*outward;
		  else
			  flux[phase] = varNData[j].densityW*varNData[j].mobilityW*outward;

		  return flux;
   };


    // compute sources/sinks
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
    { };


	  // *********************************************************
	  // *														 *
	  // *	Calculation of Data at Nodes that has to be			 *
	  // *	determined only once	(statNData)					 *
	  // *														 *
	  // *********************************************************

    // analog to EvalStaticData in MUFTE
    void updateStaticData (const Entity& e, const VBlockType* sol)
    {
		// ASSUME porosity defined at nodes
		sNDat[globalIdx].porosity = problem.soil().porosity(this->fvGeom.elementGlobal, e, this->fvGeom.elementLocal);

		return;
    }

    // harmonic mean of the permeability computed directly
    virtual FMatrix harmonicMeanK (FMatrix& Ki, const FMatrix& Kj) const
    {
    	double eps = 1e-20;

    	for (int kx=0; kx<dim; kx++){
    		for (int ky=0; ky<dim; ky++){
    			if (Ki[kx][ky] != Kj[kx][ky])
    			{
    				Ki[kx][ky] = 2 / (1/(Ki[kx][ky]+eps) + (1/(Kj[kx][ky]+eps)));
    			}
    		}
    	}
   	 return Ki;
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
    	RT satN;
    	RT pressure; // pC is equal to pW, pN is zero
    	RT mobilityW;  //Vector with the number of phases
    	RT densityW;
    };

    // analog to EvalPrimaryData in MUFTE, uses members of varNData
	void updateVariableData(const Entity& e, const VBlockType* sol,
			int i, std::vector<VariableNodeData>& varData)
    {
	   	const int globalIdx = this->vertexMapper.template map<dim>(e, i);
		FVector& global = this->fvGeom.subContVol[i].global;
  	   	FVector& local = this->fvGeom.subContVol[i].local;

  	   	// TODO inversePressSat
		varData[i].pressure = sol[i][pWIdx];
		varData[i].satW = capillaryPressureInverse(varData[i].pressure);
		varData[i].satN = 1.0 - varData[i].satW;
   		varData[i].mobilityW = problem.materialLaw().mobW(varData[i].satW, global, e, local);
   		varData[i].densityW = problem.wettingPhase().density(283.15, varData[i].pW);
   		varData[i].dSwdPc = problem.materialLaw().dSwdPc(varData[i].satW, global, e, local);
   		//   		varData[i].pC = problem.materialLaw().pC(varData[i].satW, global, e, local);

   		// for output
		 (*outPressureW)[globalIdx] = varData[i].pressure;
		 (*outSaturationW)[globalIdx] = varData[i].satW;
		 (*outSaturationN)[globalIdx] = varData[i].satN;
		 (*outDensityW)[globalIdx] = varData[i].densityW;
		 (*outMobilityW)[globalIdx] = varData[i].mobilityW;
		 (*outdSwdPc)[globalIdx] = varData[i].dSwdPc;

  	   	 return;
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
		int size = this->fvGeom.numVertices;

		for (int i = 0; i < size; i++)
				updateVariableData(e, sol, i, old);
	}


    void clearVisited ()
    {
    	return;
    }

    struct StaticNodeData
    {
        RT porosity;
        FieldMatrix<DT,dim,dim> K;
    	bool visited;
    };

    struct ElementData {
   	 RT cellVolume;
   	 } elData;

     // parameters given in constructor
     TwoPhaseProblem<G,RT>& problem;
     std::vector<StaticNodeData> statNData;
     std::vector<VariableNodeData> varNData;
     std::vector<VariableNodeData> oldVarNData;

     // for output files
     BlockVector<FieldVector<RT, 1> > *outPressureW;
     BlockVector<FieldVector<RT, 1> > *outSaturationN;
     BlockVector<FieldVector<RT, 1> > *outSaturationW;
     BlockVector<FieldVector<RT, 1> > *outDensityW;
     BlockVector<FieldVector<RT, 1> > *outMobilityW;
     BlockVector<FieldVector<RT, 1> > *outdSwdPc;


  };

  /** @} */
}
#endif
