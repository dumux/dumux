#ifndef DUNE_BOXPWSNTEJACOBIAN_HH
#define DUNE_BOXPWSNTEJACOBIAN_HH

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
#include "dumux/operators/boxjacobian.hh"
#include "dumux/2penergy/2penergyproblem.hh"

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
  template<class G, class RT, class BoxFunction = LeafP1Function<G, RT, 3> >
  class BoxPwSnTeJacobian 
    : public BoxJacobian<BoxPwSnTeJacobian<G,RT,BoxFunction>,G,RT,3,BoxFunction>
  {
    typedef typename G::ctype DT;
    typedef typename G::Traits::template Codim<0>::Entity Entity;
    typedef typename Entity::Geometry Geometry;
    typedef BoxPwSnTeJacobian<G,RT,BoxFunction> ThisType;
    typedef typename LocalJacobian<ThisType,G,RT,3>::VBlockType VBlockType;
    typedef typename LocalJacobian<ThisType,G,RT,3>::MBlockType MBlockType;
 	typedef FVElementGeometry<G> FVElementGeometry;
	enum {pWIdx = 0, satNIdx = 1, teIdx = 2}; // Solution vector index
	enum {water = 0, co2 = 1, heat = 2};	// balacne equation index		
	
  public:
    // define the number of components of your system, this is used outside
    // to allocate the correct size of (dense) blocks with a FieldMatrix
    enum {n=G::dimension};
    enum {m=3};
    enum {SIZE=LagrangeShapeFunctionSetContainer<DT,RT,n>::maxsize};
    
    //! Constructor
    BoxPwSnTeJacobian (TwoPhaseHeatProblem<G,RT>& params,
			      bool levelBoundaryAsDirichlet_, const G& grid, 
			      BoxFunction& sol, 
			      bool procBoundaryAsDirichlet_=true)
    : BoxJacobian<ThisType,G,RT,m,BoxFunction>(levelBoundaryAsDirichlet_, grid, sol, procBoundaryAsDirichlet_), 
      problem(params), 
      statNData(this->vertexMapper.size())
    {
      this->analytic = false;
    }
    

    virtual void clearVisited ()
    {
    	return;
    }

    // Compute Storage Terms
    virtual VBlockType computeM (const Entity& e, const VBlockType* sol, int node)
    {
   	 VBlockType result; 
	 
   	 // storage term of water mass balance 
   	 result[water] = -varNData[node].density[pWIdx] * elData.porosity * varNData[node].satN;
   	 
   	 // storage term of CO2 mass balance
   	 result[co2] = varNData[node].density[satNIdx] * elData.porosity * varNData[node].satN;
   	 
   	 // storage term of energy equation
   	 result[heat] = elData.porosity * (varNData[node].density[pWIdx] * varNData[node].intenergy[pWIdx] * varNData[node].satW
   	             + varNData[node].density[satNIdx] * varNData[node].intenergy[satNIdx] * varNData[node].satN)
   	             + (1. - elData.porosity) * elData.soilDens * elData.soilCs * varNData[node].temp;
   	 
   	 return result;
    };
    
    // Get Source/Sink Terms
    virtual VBlockType computeQ (const Entity& e, const VBlockType* sol, const int& node)
    {
   	 // ASSUME problem.q already contains \rho.q
   	 return problem.q(this->fvGeom.subContVol[node].global, e, this->fvGeom.subContVol[node].local);
    };
    
    // Compute Flux Terms
    virtual VBlockType computeA (const Entity& e, const VBlockType* sol, int face)
    {
   	 int i = this->fvGeom.subContVolFace[face].i;
     int j = this->fvGeom.subContVolFace[face].j;

     // Data at Control Volume Face
     // Permeability still needs to be adjusted
	  // permeability in edge direction 
	  FieldVector<DT,n> Kij(0);
	  elData.K.umv(this->fvGeom.subContVolFace[face].normal, Kij);

	  
	  // Harmonic mean:
       // Heat Conductivity
     	RT lambda;
     	lambda = 2./((1./varNData[i].lambda) + (1./varNData[j].lambda));
  	  
     // Arithmetic mean:
     	// Enthalpy
     	RT enthW;
     	RT enthCO2;
     	enthW = (varNData[i].enthalpy[pWIdx] + varNData[j].enthalpy[pWIdx]) / 2.;
     	enthCO2 = (varNData[i].enthalpy[satNIdx] + varNData[j].enthalpy[satNIdx]) / 2.;
     	
     	
     // Calculate Gradients:
     // calculate FE gradient
     FieldVector<RT,n> pWGrad(0);
     FieldVector<RT,n> pCO2Grad(0);
     FieldVector<RT,n> teGrad(0);
     
     for (int comp = 0; comp < m; comp++)
     {
   	  	  FieldVector<RT, n> gravity = problem.gravity();
      	  switch (comp){
    		 case water: 
    		  	 for (int k = 0; k < this->fvGeom.nNodes; k++) {	 
    		     FieldVector<DT,n> grad(this->fvGeom.subContVolFace[face].grad[k]);
    			 grad *= varNData[k].pW;
    			 pWGrad += grad;
    			 }
    		  	 // adjust by gravity
    		  	 gravity *= varNData[i].density[comp];
    		  	 pWGrad -= gravity;
    			 break;
    		 case co2:
    			 for (int k = 0; k < this->fvGeom.nNodes; k++) {	 
    			 FieldVector<DT,n> grad(this->fvGeom.subContVolFace[face].grad[k]);
    			 grad *= varNData[k].pN;
    			 pCO2Grad += grad;
    			 }
    		  	 // adjust by gravity
    		  	 gravity *= varNData[i].density[comp];
    		  	 pCO2Grad -= gravity;
    			 break;
    		 case heat:
    			 for (int k = 0; k < this->fvGeom.nNodes; k++) {	 
    			 FieldVector<DT,n> grad(this->fvGeom.subContVolFace[face].grad[k]);
    			 grad *= varNData[k].temp;
    			 teGrad += grad;
    			 }
    			 break;
      	  }
     }
			  
	  // calculate the flux using upwind
	  RT s_w = pWGrad*Kij;
	  RT s_g = pCO2Grad*Kij;
	  RT s_h = teGrad*this->fvGeom.subContVolFace[face].normal;
	  	 s_h *= -lambda;
	  
	  VBlockType flux;
	  
	  int up_w, down_w, up_g, down_g;
	  
	  if (s_w >= 0){
		   up_w = i; down_w = j;}
	  else {
		  up_w = j; down_w = i;}
	  if (s_g >= 0){
		   up_g = i; down_g = j;}
	  else {
		  up_g = j; down_g = i;}
	  
	  // flux term of water mass balance equation
		  flux[water] = varNData[up_w].density[water] * varNData[up_w].mobility[water] * s_w;

	  // flux term of co2 mass balance equation
		  flux[co2] = varNData[up_g].density[co2] * varNData[up_g].mobility[co2] * s_g;
		  
	  // flux term of energy balance equation
		  flux[heat] = varNData[up_w].density[water] * varNData[up_w].mobility[water] * s_w * enthW
		  				+ varNData[up_g].density[co2] * varNData[up_g].mobility[co2] * s_g * enthCO2
		  				+ s_h;
		  
		  return flux;
   };
   


    
    
	  //*********************************************************
	  //*																			*
	  //*	Calculation of Data at Elements			 					*
	  //*						 													*
	  //*																		 	*
	  //*********************************************************

    virtual void computeElementData (const Entity& e)
    {
  		 // ASSUME element-wise constant parameters for the material law 
 		 elData.parameters = problem.materialLawParameters(this->fvGeom.cellGlobal, e, this->fvGeom.cellLocal);
   	 
		 // ASSUMING element-wise constant permeability, evaluate K at the cell center 
 		 elData.K = problem.K(this->fvGeom.cellGlobal, e, this->fvGeom.cellLocal);  

 		 elData.porosity = problem.porosity(this->fvGeom.cellGlobal, e, this->fvGeom.cellLocal);
 		 
 		 elData.soilDens = problem.soilParameters(this->fvGeom.cellGlobal, e, this->fvGeom.cellLocal)[0];
 		 
 		 elData.soilCs	= problem.soilParameters(this->fvGeom.cellGlobal, e, this->fvGeom.cellLocal)[1];
 		 
 		 elData.soilLDry = problem.soilParameters(this->fvGeom.cellGlobal, e, this->fvGeom.cellLocal)[2];
 		 
 		 elData.soilLSw	= problem.soilParameters(this->fvGeom.cellGlobal, e, this->fvGeom.cellLocal)[3];
 		 
 		 
    };

    
	  //*********************************************************
	  //*																			*
	  //*	Calculation of Data at Nodes that has to be			 	*
	  //*	determined only once	(statNData)							 	*
	  //*																		 	*
	  //*********************************************************

    // analog to EvalStaticData in MUFTE
    virtual void updateStaticData (const Entity& e, const VBlockType* sol)
    {
	  return;
    }


	  //*********************************************************
	  //*																			*
	  //*	Calculation of variable Data at Nodes					 	*
	  //*	(varNData)														 	*
	  //*																		 	*
	  //*********************************************************
   

    // analog to EvalPrimaryData in MUFTE, uses members of varNData
    virtual void updateVariableData (const Entity& e, const VBlockType* sol)
    {
     
   	 varNData.resize(this->fvGeom.nNodes);
   	 int size = varNData.size();

   	 for (int i = 0; i < size; i++) {
// 		this->def[i] = 0;
   		varNData[i].saturationW = 1.0 - sol[i][satNIdx];
   		// ASSUME element-wise constant parameters for the material law 
         FieldVector<RT, 4> parameters = problem.materialLawParameters(this->fvGeom.cellGlobal, e, this->fvGeom.cellLocal);

         varNData[i].temp = sol[i][teIdx];
         varNData[i].pW = sol[i][pWIdx];
         varNData[i].satN = sol[i][satNIdx];
         varNData[i].satW = 1 - sol[i][satNIdx];
         varNData[i].lambda = soil.heatConductivity(elData.soilLDry, elData.soilLSw, varNData[i].satW);
         varNData[i].pC = problem.materialLaw().pC(varNData[i].saturationW,parameters);
         varNData[i].pN = sol[i][pWIdx] + varNData[i].pC;
         varNData[i].density[pWIdx] = problem.materialLaw().wettingPhase.density(varNData[i].temp,varNData[i].pW);
         varNData[i].density[satNIdx] = problem.materialLaw().nonwettingPhase.density(varNData[i].temp,varNData[i].pN);
         varNData[i].viscosity[pWIdx] = problem.materialLaw().wettingPhase.viscosity(varNData[i].temp,sol[i][pWIdx]);
         varNData[i].viscosity[satNIdx] = problem.materialLaw().nonwettingPhase.viscosity(varNData[i].temp,varNData[i].pN,varNData[i].density[satNIdx]);
         varNData[i].mobility[pWIdx] = problem.materialLaw().mobW(varNData[i].saturationW, parameters, varNData[i].viscosity[pWIdx]);
         varNData[i].mobility[satNIdx] = problem.materialLaw().mobN(sol[i][satNIdx], parameters, varNData[i].viscosity[satNIdx]);
         varNData[i].enthalpy[pWIdx] = problem.materialLaw().wettingPhase.enthalpy(varNData[i].temp,varNData[i].pW);
         varNData[i].enthalpy[satNIdx] = problem.materialLaw().nonwettingPhase.enthalpy(varNData[i].temp,varNData[i].pN);
         varNData[i].intenergy[pWIdx] = problem.materialLaw().wettingPhase.intEnergy(varNData[i].temp,varNData[i].pW);
         varNData[i].intenergy[satNIdx] = problem.materialLaw().nonwettingPhase.intEnergy(varNData[i].temp,varNData[i].pN);
 
         // debug: std::cout  << "enthalpy " << varNData[i].enthalpy[satNIdx] << "internal energy " << varNData[i].intenergy[satNIdx] << std::endl;
   	 }   	 
    }
    
    
    struct StaticNodeData 
    {
    	bool visited;
    };
    
       // the members of the struct are defined here
    struct VariableNodeData  
    {
       RT saturationW;
       RT pC;
       RT pN;
       RT pW;
       RT temp;
       RT satN;
       RT satW;
       RT lambda;
       VBlockType mobility;  //Vector with the number of phases
       VBlockType density;
       VBlockType viscosity;
       VBlockType enthalpy;
       VBlockType intenergy;
    };
    
    struct ElementData {
   	 RT cellVolume;
     RT porosity;
   	 RT gravity;
   	 RT soilDens;
   	 RT soilCs;
   	 RT soilLDry;
   	 RT soilLSw;
   	 FieldVector<RT, 4> parameters;
   	 FieldMatrix<DT,n,n> K;
   	 } elData;
    
    // parameters given in constructor
    TwoPhaseHeatProblem<G,RT>& problem;
    Soil soil;
    std::vector<StaticNodeData> statNData;
    std::vector<VariableNodeData> varNData;
  };

  /** @} */
}
#endif
