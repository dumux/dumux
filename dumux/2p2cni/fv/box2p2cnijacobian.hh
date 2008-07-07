#ifndef DUNE_BOX2P2CNIJACOBIAN_HH
#define DUNE_BOX2P2CNIJACOBIAN_HH

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
#include<dune/grid/utility/intersectiongetter.hh>

#include<dune/disc/shapefunctions/lagrangeshapefunctions.hh>
#include<dune/disc/operators/boundaryconditions.hh>
#include<dune/disc/functions/p1function.hh>
//#include"../boxjacobian_2p2c.hh"
#include"dumux/operators/boxjacobian.hh"
#include"dumux/2p2cni/2p2cniproblem.hh"
//#include "varswitch.hh"
#include"dumux/../Stupid/Stupid/Auxilary/VtkMultiWriter.hh"

/**
 * @file
 * @brief  compute local jacobian matrix for box scheme for two-phase two-component flow equation
 * @author Bernd Flemisch, Klaus Mosthaf
 */



namespace Dune
{
  /** @addtogroup DISC_Disc
   *
   * @{
   */
  /**
   * @brief compute local jacobian matrix for the boxfile for two-phase two-component flow equation
   *
   */


  //! Derived class for computing local jacobian matrices
  /*! A class for computing local jacobian matrix for the two-phase two-component flow equation

	    div j = q; j = -K grad u; in Omega

		u = g on Gamma1; j*n = J on Gamma2.

	Uses box scheme with the Lagrange shape functions.
	It should work for all dimensions and element types.
	All the numbering is with respect to the reference element and the
	Lagrange shape functions

	Template parameters are:

	- Grid  a DUNE grid type
	- RT    type used for return values 
  */
  template<class G, class RT, class BoxFunction = LeafP1Function<G, RT, 3> >
  class Box2P2CNIJacobian 
    : public BoxJacobian<Box2P2CNIJacobian<G,RT,BoxFunction>,G,RT,3,BoxFunction>
  {
    typedef typename G::ctype DT;
    typedef typename G::Traits::template Codim<0>::Entity Entity;
    typedef typename Entity::Geometry Geometry;
    typedef Box2P2CNIJacobian<G,RT,BoxFunction> ThisType;
    typedef typename LocalJacobian<ThisType,G,RT,3>::VBlockType VBlockType;
    typedef typename LocalJacobian<ThisType,G,RT,3>::MBlockType MBlockType;
    typedef FVElementGeometry<G> FVElementGeometry;

 	enum {pWIdx = 0, satNIdx = 1, teIdx=2, numberOfComponents = 2};	// Solution vector index
	enum {wPhase = 0, nPhase = 1};					// Phase index
	enum {gasPhase = 0, waterPhase = 1, bothPhases = 2};	// Phase state
	enum {water = 0, air = 1, heat =2};				// Component index					
	
  public:
    // define the number of phases (m) and components (c) of your system, this is used outside
    // to allocate the correct size of (dense) blocks with a FieldMatrix
    enum {dim=G::dimension};
    enum {m=3, c=2};
    enum {SIZE=LagrangeShapeFunctionSetContainer<DT,RT,dim>::maxsize};

    typedef FieldMatrix<RT,dim,dim> FMatrix;
    typedef FieldVector<RT,dim> FVector;
    
    //! Constructor
    Box2P2CNIJacobian (TwoPTwoCNIProblem<G,RT>& params,
			      bool levelBoundaryAsDirichlet_, const G& grid, 
			      BoxFunction& sol, 
			      bool procBoundaryAsDirichlet_=true)
    : BoxJacobian<ThisType,G,RT,m,BoxFunction>(levelBoundaryAsDirichlet_, grid, sol, procBoundaryAsDirichlet_), 
      problem(params), 
      sNDat(this->vertexMapper.size())
  	{
      this->analytic = false;
    }

	/** @brief compute time dependent term (storage), loop over nodes / subcontrol volumes
	 *  @param e entity   
	 *  @param sol solution vector
	 *  @param node local node id
	 *  @return storage term
	 */
    virtual VBlockType computeM (const Entity& e, const VBlockType* sol, int node)
    {
    	 GeometryType gt = e.geometry().type();
    	 const typename LagrangeShapeFunctionSetContainer<DT,RT,dim>::value_type&
     	 sfs=LagrangeShapeFunctions<DT,RT,dim>::general(gt,1);
    	 
   	 int globalIdx = this->vertexMapper.template map<dim>(e, sfs[node].entity());
   	 //int globalCoord = this->fvGeom.subContVol[node].global;

   	 VBlockType result; 
   	 RT satN = vNDat[node].satN;
   	 RT satW = vNDat[node].satW;
   	    	    	                  
   	 // storage of component water
   	 result[water] = 
   		  sNDat[globalIdx].porosity*(vNDat[node].density[wPhase]*satW*vNDat[node].massfrac[water][wPhase]
   		                 +vNDat[node].density[nPhase]*satN*vNDat[node].massfrac[water][nPhase]);
   	 // storage of component air
   	 result[air] = 
   		  sNDat[globalIdx].porosity*(vNDat[node].density[nPhase]*satN*vNDat[node].massfrac[air][nPhase]
   	                    +vNDat[node].density[wPhase]*satW*vNDat[node].massfrac[air][wPhase]);

   	 // storage term of energy equation
   	 result[heat] = sNDat[globalIdx].porosity * (vNDat[node].density[wPhase] * vNDat[node].intenergy[wPhase] * satW
   	             + vNDat[node].density[nPhase] * vNDat[node].intenergy[nPhase] * satN)
   	             + (1. -  sNDat[globalIdx].porosity) * elData.soilDens * elData.soilCs * vNDat[node].temperature;
   	  // soil properties defined at the elements!!!               
   	 
   	 //std::cout << result << " " << node << std::endl;
   	 return result;
    };

    /** @brief compute diffusive/advective fluxes, loop over subcontrol volume faces
	 *  @param e entity   
	 *  @param sol solution vector
	 *  @param face face id
	 *  @return flux term
     */
    virtual VBlockType computeA (const Entity& e, const VBlockType* sol, int face)
    {
   	 int i = this->fvGeom.subContVolFace[face].i;
 	 int j = this->fvGeom.subContVolFace[face].j;
 	 
 	 // normal vector, value of the area of the scvf
	 const FieldVector<RT,dim> normal(this->fvGeom.subContVolFace[face].normal);

	 // get global coordinates of nodes i,j
	 const FieldVector<DT,dim> global_i = this->fvGeom.subContVol[i].global;
	 const FieldVector<DT,dim> global_j = this->fvGeom.subContVol[j].global;

	 FieldMatrix<RT,2,dim> pGrad(0.), xGrad(0.); 
	 FieldVector<RT,dim> teGrad(0.); 
	 FieldVector<RT,dim> temp(0.); 
//		 FieldVector<RT,dim> Kij(0);
     	 VBlockType flux(0.);
		 
		 // effective permeability in edge direction 
     	 // RT Kij = sIPDat[global_j].K_eff[face]; 
     	 
     	 // calculate harmonic mean of permeabilities of nodes i and j
		 const FMatrix K = harmonicMeanK(global_i, global_j);
		 //const FMatrix K = harmonicMeanK(e, face);		  
     	 //K.umv(normal, Kij);  // Kij=K*n

	 // Harmonic mean:
		// Heat Conductivity
     	RT lambda;
     	lambda = 2./((1./vNDat[i].lambda) + (1./vNDat[j].lambda));
  	  
     // Arithmetic mean:
     	// Enthalpy
     	RT enthW;
     	RT enthCO2;
     	enthW = (vNDat[i].enthalpy[pWIdx] + vNDat[j].enthalpy[pWIdx]) / 2.;
     	enthCO2 = (vNDat[i].enthalpy[satNIdx] + vNDat[j].enthalpy[satNIdx]) / 2.;
  

		 
		 // calculate FE gradient (grad p for each phase)
		 for (int k = 0; k < this->fvGeom.nNodes; k++) // loop over adjacent nodes
		 {	 
			 // FEGradient at node k
			 const FieldVector<DT,dim> feGrad(this->fvGeom.subContVolFace[face].grad[k]);
			 FieldVector<RT,m> pressure(0.0), massfrac(0.0);

			 pressure[wPhase] = vNDat[k].pW;
			 pressure[nPhase] = vNDat[k].pN;
      	 
		  	 // compute sum of pressure gradients for each phase
		  	 for (int phase = 0; phase < 2; phase++)
		  	 {	      		 
		  		 temp = feGrad;
		  		 temp *= pressure[phase];
		  		 pGrad[phase] += temp;
		  	 }
		  	 // Temperature gradient
		  	 temp = feGrad;
	     	 temp *= vNDat[k].temperature;
	     	 teGrad += temp;
	     	 
		  	 // for diffusion of air in wetting phase
		  	 temp = feGrad;
	     	 temp *= vNDat[k].massfrac[air][wPhase];
	     	 xGrad[wPhase] += temp;

		  	 // for diffusion of water in nonwetting phase
	     	 temp = feGrad;
	     	 temp *= vNDat[k].massfrac[water][nPhase];
	     	 xGrad[nPhase] += temp;
		 }

       // deduce gravity*density of each phase
		 FieldMatrix<RT,2,dim> contribComp(0);
		 for (int phase=0; phase<2; phase++)
		 {
			 contribComp[phase] = problem.gravity();
			 contribComp[phase] *= vNDat[i].density[phase];  
			 pGrad[phase] -= contribComp[phase]; // grad p - rho*g
		 }
	 	 
		 VBlockType outward(0);  // Darcy velocity of each phase
		 FieldVector<RT,dim> v_tilde(0);

		 // calculate the advective flux using upwind: K*n(grad p -rho*g)
		 for (int phase=0; phase<2; phase++) 
		 	{
	     	 K.umv(pGrad[phase], v_tilde);  // v_tilde=K*gradP
	     	 outward[phase] = v_tilde*normal;
		 	}
		 	outward[heat] = teGrad * normal;
	  	    outward[heat] *= lambda;  // schau mer mal
		 
		 // evaluate upwind nodes
		 int up_w, dn_w, up_n, dn_n;
		 if (outward[wPhase] <= 0) {up_w = i; dn_w = j;}
		 else {up_w = j; dn_w = i;};
		 if (outward[nPhase] <= 0) {up_n = i; dn_n = j;}
		 else {up_n = j; dn_n = i;};

		 RT alpha = 1.0;  // Upwind parameter
 		 
		 // Water conservation
 		 flux[water] =   (alpha* vNDat[up_w].density[wPhase]*vNDat[up_w].mobility[wPhase]
 				               * vNDat[up_w].massfrac[water][wPhase] 
 				    + (1-alpha)* vNDat[dn_w].density[wPhase]*vNDat[dn_w].mobility[wPhase]
 				               * vNDat[dn_w].massfrac[water][wPhase])
 				               * outward[wPhase]; 		 		
 		 flux[water] +=  (alpha* vNDat[up_n].density[nPhase]*vNDat[up_n].mobility[nPhase]
 				               * vNDat[up_n].massfrac[water][nPhase] 
 				    + (1-alpha)* vNDat[dn_n].density[nPhase]*vNDat[dn_n].mobility[nPhase]
 				               * vNDat[dn_n].massfrac[water][nPhase])
 				               * outward[nPhase];
 		 // Air conservation
 		 flux[air]   =   (alpha* vNDat[up_n].density[nPhase]*vNDat[up_n].mobility[nPhase]
 				               * vNDat[up_n].massfrac[air][nPhase] 
 				    + (1-alpha)* vNDat[dn_n].density[nPhase]*vNDat[dn_n].mobility[nPhase]
 				               * vNDat[dn_n].massfrac[air][nPhase])
 				               * outward[nPhase];
 		 flux[air]  +=   (alpha* vNDat[up_w].density[wPhase]*vNDat[up_w].mobility[wPhase]
 				               * vNDat[up_w].massfrac[air][wPhase]
 				    + (1-alpha)* vNDat[dn_w].density[wPhase]*vNDat[dn_w].mobility[wPhase]
 				               * vNDat[dn_w].massfrac[air][wPhase])
 				               * outward[wPhase];
 				               
 		 // Heat conservation
 		 	  // flux term of energy balance equation
 		 flux[heat]   =  (alpha* vNDat[up_n].density[nPhase]*vNDat[up_n].mobility[nPhase]
 				    + (1-alpha)* vNDat[dn_n].density[nPhase]*vNDat[dn_n].mobility[nPhase])
 				               * enthCO2 * outward[nPhase];
 		 flux[heat]  +=  (alpha* vNDat[up_w].density[wPhase]*vNDat[up_w].mobility[wPhase]
 				    + (1-alpha)* vNDat[dn_w].density[wPhase]*vNDat[dn_w].mobility[wPhase])
 				    		   * enthW * outward[wPhase];
 		 flux[heat]	+=	outward[heat];

		 
		 // DIFFUSION
		 VBlockType normDiffGrad;

		 RT diffusionWW, diffusionWN; // diffusion of water
		 RT diffusionAW, diffusionAN; // diffusion of air
		 VBlockType avgDensity, avgDpm;
		 avgDpm[wPhase]=1e-9; // needs to be changed !!!
		 avgDpm[nPhase]=1e-9; // water in the gasphase
		 
		 normDiffGrad[wPhase] = xGrad[wPhase]*normal;
		 normDiffGrad[nPhase] = xGrad[nPhase]*normal;

		 // calculate the arithmetic mean of densities
		 avgDensity[wPhase] = 0.5*(vNDat[i].density[wPhase] + vNDat[j].density[wPhase]);
		 avgDensity[nPhase] = 0.5*(vNDat[i].density[nPhase] + vNDat[j].density[nPhase]);

		 
		 diffusionAW = avgDpm[wPhase] * avgDensity[wPhase] * normDiffGrad[wPhase];
		 diffusionWW = - diffusionAW;
		 diffusionWN = avgDpm[nPhase] * avgDensity[nPhase] * normDiffGrad[nPhase];
		 diffusionAN = - diffusionWN;
		 
//		 std::cout << "Diffusive Flux: " << diffusionWG << std::endl; 
		 

		 // add water diffusion to flux
		 flux[water] += (diffusionWW + diffusionWN); 
//		 std::cout << "Water Flux: " << flux[water] << std::endl; 
		 // air diffusion not implemeted
		 flux[air] += (diffusionAN + diffusionAW);
//		 std::cout << "Air Flux: " << flux[air] << std::endl; 


		 return flux;
  };
    
  	/** @brief integrate sources / sinks
  	 *  @param e entity   
	 *  @param sol solution vector
	 *  @param node local node id
	 *  @return source/sink term
	 */
   	virtual VBlockType computeQ (const Entity& e, const VBlockType* sol, const int& node)
   	{
   		// ASSUME problem.q already contains \rho.q
   		return problem.q(this->fvGeom.subContVol[node].global, e, this->fvGeom.subContVol[node].local);
   	}

  	/** @brief perform variable switch
  	 *  @param global global node id   
	 *  @param sol solution vector
	 *  @param local local node id
	 */
//   	virtual void primaryVarSwitch (const Entity& e, const int global, VBlockType* sol, const int local)
   	virtual void primaryVarSwitch (const Entity& e, int global, VBlockType* sol, int local)
    {
        // neue Variable in SOL schreiben. Da Zeiger auf sol, sollte das funktionieren
        // Achtung: alte Loesung beibebehalten, falls Newton den Zeitschritt nicht schafft
   		// 1=gas, 2=water, 3=both
   		

   		int state = sNDat[global].phaseState;
        bool switched = false;
    	RT pWSat = problem.multicomp().vaporPressure(vNDat[local].temperature);
    	RT henryInv = problem.multicomp().henry(vNDat[local].temperature);
    	RT xWNmolar = problem.multicomp().xWNmolar(vNDat[local].pN, vNDat[local].temperature);
    	RT xAWmolar = problem.multicomp().xAWmolar(vNDat[local].pN, vNDat[local].temperature);
//    	double satControl = vNDat[local].satW;

    	FVector Coordinates = this->fvGeom.cellGlobal;

    	
        switch(state) 
        {
        case gasPhase :
        	RT xWNmass, pwn; 
        	xWNmass = sol[local][satNIdx];
        	xWNmolar = problem.multicomp().conversionMassToMoleFraction(xWNmass, gasPhase);
           	pwn = xWNmolar * vNDat[local].pN;
  
        	if (pwn > 1.01*pWSat && switched == false)
            {
            	// appearance of water phase
            	std::cout << "Water appears at node " << global << "  Coordinates: " << Coordinates << std::endl;
            	vNDat[local].satN = 1.0 - 1.e-6; // initialize
            	sNDat[global].phaseState = bothPhases;
//      		  	sNDat[global].switched = true;
            	sol[local][satNIdx] = vNDat[local].satN;
            	switched = true;
            }
            break;

        case waterPhase :
        	RT pbub, pNint, xAWmass;
         	xAWmass = sol[local][satNIdx];
         	xAWmolar = problem.multicomp().conversionMassToMoleFraction(xAWmass, waterPhase);
        	pbub = pWSat + xAWmolar/henryInv;
        	pNint = vNDat[local].pN;
        	
        	if (pbub > pNint && switched == false)
            {
            	// appearance of gas phase
            	std::cout << "Gas appears at node " << global << "  Coordinates: " << Coordinates << std::endl;
            	vNDat[local].satN = 1.e-6;  // initialize
            	sNDat[global].phaseState = bothPhases;
//      		  	sNDat[global].switched = true;
            	sol[local][satNIdx] = vNDat[local].satN;
            	switched = true;
            }
            break;

        case bothPhases :
      	  	if (vNDat[local].satN < 0.0  && switched == false)
      	  	{
      		  	// disappearance of gas phase
      		  	std::cout << "Gas disappears at node " << global << "  Coordinates: " << Coordinates << std::endl;
      		  	vNDat[local].satN = 0.0;  
      		  	vNDat[local].satW = 1.0;  
      		  	sNDat[global].phaseState = waterPhase;
//      		  	sNDat[global].switched = true;
      		  	vNDat[local].massfrac[air][wPhase] = 1e-6;
      		  	sol[local][satNIdx] = vNDat[local].massfrac[air][wPhase];
      		  	switched = true;
            }
      	  	else if (vNDat[local].satW < 0.0  && switched == false)
      	  	{
      	  		// disappearance of water phase
      	  		std::cout << "Water disappears at node " << global << "  Coordinates: " << Coordinates << std::endl;
      		  	vNDat[local].satN = 1.0;
      	  		vNDat[local].satW = 0.0;
      	  		sNDat[global].phaseState = gasPhase;
//      		  	sNDat[global].switched = true;
      		  	vNDat[local].massfrac[water][nPhase] = 1e-6;
      	  		sol[local][satNIdx] = vNDat[local].massfrac[water][nPhase];
            	switched = true;
      	  	}
      	  	break;
      	  	
      	  	if (switched == true) updateVariableData(e, sol);
        }
        
   	return;
    }
    
    // harmonic mean computed directly
    virtual FMatrix harmonicMeanK (const FVector global_i, const FVector global_j) const
    {
    	double eps = 1e-20;

    	FMatrix Ki, Kj;
   	 
    	Ki = this->problem.K(global_i); 
    	Kj = this->problem.K(global_j);

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
    
    
    virtual void clearVisited ()
    {
    	for (int i = 0; i < this->vertexMapper.size(); i++){
   		sNDat[i].visited = false;
//   	 	sNDat[i].switched = false;
    	}
   	 return;
   	}
         
    
	  //*********************************************************
	  //*														*
	  //*	Calculation of Data at Elements (elData) 			*
	  //*						 								*
	  //*														*
	  //*********************************************************

    virtual void computeElementData (const Entity& e)
    {

 	 
 		 elData.soilDens = problem.soilParameters(this->fvGeom.cellGlobal, e, this->fvGeom.cellLocal)[0];
 		 
 		 elData.soilCs	= problem.soilParameters(this->fvGeom.cellGlobal, e, this->fvGeom.cellLocal)[1];
 		 
 		 elData.soilLDry = problem.soilParameters(this->fvGeom.cellGlobal, e, this->fvGeom.cellLocal)[2];
 		 
 		 elData.soilLSw	= problem.soilParameters(this->fvGeom.cellGlobal, e, this->fvGeom.cellLocal)[3];
 		 
 		 
//  	 // ASSUME element-wise constant parameters for the material law 
// 		 elData.parameters = problem.materialLawParameters
// 		 (this->fvGeom.cellGlobal, e, this->fvGeom.cellLocal);
//   	 
//		 // ASSUMING element-wise constant permeability, evaluate K at the cell center 
// 		 elData.K = problem.K(this->fvGeom.cellGlobal, e, this->fvGeom.cellLocal);  
//
//		 // ASSUMING element-wise constant porosity 
// 		 elData.porosity = problem.porosity(this->fvGeom.cellGlobal, e, this->fvGeom.cellLocal);
   	 return;
    }

    
	  //*********************************************************
	  //*														*
	  //*	Calculation of Data at Nodes that has to be			*
	  //*	determined only once	(sNDat)						*
	  //*														*
	  //*********************************************************

    // analog to EvalStaticData in MUFTE
    virtual void updateStaticData (const Entity& e, VBlockType* sol)
    {
   	 // size of the sNDat vector is determined in the constructor
   	 
   	 // local to global id mapping (do not ask vertex mapper repeatedly
   	 //int localToGlobal[LagrangeShapeFunctionSetContainer<DT,RT,n>::maxsize];

   	 // get access to shape functions for P1 elements
   	 GeometryType gt = e.geometry().type();
   	 const typename LagrangeShapeFunctionSetContainer<DT,RT,dim>::value_type& 
   	 sfs=LagrangeShapeFunctions<DT,RT,dim>::general(gt,1);
   	 
   	 // get local to global id map
   	 for (int k = 0; k < sfs.size(); k++) {
  		 const int globalIdx = this->vertexMapper.template map<dim>(e, sfs[k].entity());
  		  
  		 // if nodes are not already visited
  		 if (!sNDat[globalIdx].visited) 
  		  {
  			  // ASSUME porosity defined at nodes
  			  sNDat[globalIdx].porosity = problem.porosity(this->fvGeom.cellGlobal, e, this->fvGeom.cellLocal);

 			  // global coordinates
 			  FieldVector<DT,dim> global_i = this->fvGeom.subContVol[k].global;
  			  
   			  // evaluate primary variable switch
// 			  primaryVarSwitch(e, globalIdx, sol, k);
  			  
  			  // mark elements that were already visited
  			  sNDat[globalIdx].visited = true;
  		  }
  	  }
  	  
	  return;
    }
    

	  //*********************************************************
	  //*														*
	  //*	Calculation of variable Data at Nodes				*
	  //*	(vNDat)												*
	  //*														*
	  //*********************************************************

    // for output files
//    BlockVector<FieldVector<RT, 1> > *hackyMassFracAir;
//    BlockVector<FieldVector<RT, 1> > *hackyMassFracWater;
//    BlockVector<FieldVector<RT, 1> > *hackySaturationN;

    // analog to EvalPrimaryData in MUFTE, uses members of vNDat
    virtual void updateVariableData (const Entity& e, const VBlockType* sol)
    {
   	 vNDat.resize(this->fvGeom.nNodes);
   	 int size = vNDat.size();

   	 GeometryType gt = e.geometry().type();
   	 const typename LagrangeShapeFunctionSetContainer<DT,RT,dim>::value_type& 
   	 sfs=LagrangeShapeFunctions<DT,RT,dim>::general(gt,1);
     
     const FieldVector<RT, 4> parameters = problem.materialLawParameters 
     (this->fvGeom.cellGlobal, e, this->fvGeom.cellLocal);

   	 for (int i = 0; i < size; i++) 
   	 {
   	   	 const int global = this->vertexMapper.template map<dim>(e, sfs[i].entity());
   		 int state = sNDat[global].phaseState;

   		 vNDat[i].pW = sol[i][pWIdx];
   		 if (state == bothPhases) vNDat[i].satN = sol[i][satNIdx];
   		 if (state == waterPhase) vNDat[i].satN = 0.0;
   		 if (state == gasPhase) vNDat[i].satN = 1.0;

   		 vNDat[i].satW = 1.0 - vNDat[i].satN;

   		 vNDat[i].pC = problem.materialLaw().pC(vNDat[i].satW, parameters);
   		 vNDat[i].pN = vNDat[i].pW + vNDat[i].pC;
   		 vNDat[i].temperature = sol[i][teIdx]; // in [K]

   		 // Solubilities of components in phases
   		 if (state == bothPhases){
   	   		 vNDat[i].massfrac[air][wPhase] = problem.multicomp().xAW(vNDat[i].pN, vNDat[i].temperature);
   	   		 vNDat[i].massfrac[water][nPhase] = problem.multicomp().xWN(vNDat[i].pN, vNDat[i].temperature);
   		 }
   		 if (state == waterPhase){
   	   		 vNDat[i].massfrac[water][nPhase] = 0.0;
   	   		 vNDat[i].massfrac[air][wPhase] =  sol[i][satNIdx];
   		 }
   		 if (state == gasPhase){
   	   		 vNDat[i].massfrac[water][nPhase] = sol[i][satNIdx];
   	   		 vNDat[i].massfrac[air][wPhase] = 0.0;
   		 }
   	   	 vNDat[i].massfrac[water][wPhase] = 1.0 - vNDat[i].massfrac[air][wPhase];
   	   	 vNDat[i].massfrac[air][nPhase] = 1.0 - vNDat[i].massfrac[water][nPhase];
   	   	 // for output
//   	   	 (*hackySaturationN)[global] = vNDat[i].satN;
//   	   	 (*hackyMassFracAir)[global] = vNDat[i].massfrac[air][wPhase];
//   	   	 (*hackyMassFracWater)[global] = vNDat[i].massfrac[water][nPhase];
   	   	 
   		 // Mobilities & densities
   		 vNDat[i].mobility[wPhase] = problem.materialLaw().mobW(vNDat[i].satW, parameters, vNDat[i].temperature, vNDat[i].pW);
   		 vNDat[i].mobility[nPhase] = problem.materialLaw().mobN(vNDat[i].satN, parameters, vNDat[i].temperature, vNDat[i].pN);
   		 vNDat[i].density[wPhase] = problem.materialLaw().wettingPhase.density(vNDat[i].temperature, vNDat[i].pN);
   		 vNDat[i].density[nPhase] = problem.materialLaw().nonwettingPhase.density(vNDat[i].temperature, vNDat[i].pN,
   				 vNDat[i].massfrac[air][nPhase]);
         vNDat[i].lambda = soil.heatConductivity(elData.soilLDry, elData.soilLSw, vNDat[i].satW);
         vNDat[i].enthalpy[pWIdx] = problem.materialLaw().wettingPhase.enthalpy(vNDat[i].temperature,vNDat[i].pW);
         vNDat[i].enthalpy[satNIdx] = problem.materialLaw().nonwettingPhase.enthalpy(vNDat[i].temperature,vNDat[i].pN);
         vNDat[i].intenergy[pWIdx] = problem.materialLaw().wettingPhase.intEnergy(vNDat[i].temperature,vNDat[i].pW);
         vNDat[i].intenergy[satNIdx] = problem.materialLaw().nonwettingPhase.intEnergy(vNDat[i].temperature,vNDat[i].pN);

         // CONSTANT solubility (for comparison with twophase)
//         vNDat[i].massfrac[air][wPhase] = 0.0; vNDat[i].massfrac[water][wPhase] = 1.0;
//         vNDat[i].massfrac[water][nPhase] = 0.0; vNDat[i].massfrac[air][nPhase] = 1.0;

         //std::cout << "water in gasphase: " << vNDat[i].massfrac[water][nPhase] << std::endl;
         //std::cout << "air in waterphase: " << vNDat[i].massfrac[air][wPhase] << std::endl;
   	 }   	 
   	 return;
    }

    
    struct StaticNodeData 
    {
   	 bool visited;
//   	 bool switched;
   	 int phaseState;
   	 RT cellVolume;
   	 RT porosity;
   	 FieldVector<RT, 4> parameters;
   	 FMatrix K;
    };
    
    struct VariableNodeData  
    {
   	 RT satN;
     RT satW;
     RT pW;
     RT pC;
     RT pN;
     RT temperature;
     RT lambda;
     VBlockType mobility;  //Vector with the number of phases
     VBlockType density;
     FieldMatrix<RT,c,m> massfrac;
     FieldVector<RT,2> enthalpy;
     FieldVector<RT,2> intenergy;
    };

	 struct StaticIPData
	 {
		 bool visited;
		 FMatrix K;
	 };
    
    
    struct ElementData {
     RT soilDens;
   	 RT soilCs;
   	 RT soilLDry;
   	 RT soilLSw;
   	 
//   	 RT cellVolume;
//     	 RT porosity;
//   	 RT gravity;
//   	 FieldVector<RT, 4> parameters;
//   	 FieldMatrix<RT,dim,dim> K;
   	 } elData;
   	    	 
    
    // parameters given in constructor
   	TwoPTwoCNIProblem<G,RT>& problem;
    CWaterAir multicomp;
    Soil soil;
    std::vector<StaticNodeData> sNDat;
    std::vector<StaticIPData> sIPDat;
    std::vector<VariableNodeData> vNDat;
  };  
  
}
#endif

//    // average permeability from the staticNode Vector
//    virtual FMatrix harmonicMeanK (const Entity& e, int k) const
//    {
//	 FMatrix Ki, Kj;
//	 const RT eps = 1e-20;
//	 
//    GeometryType gt = e.geometry().type();
//    const typename LagrangeShapeFunctionSetContainer<DT,RT,dim>::value_type& 
//    sfs=LagrangeShapeFunctions<DT,RT,dim>::general(gt,1);
//
//     int i = this->fvGeom.subContVolFace[k].i;
//     int j = this->fvGeom.subContVolFace[k].j;
//
//     int global_i = this->vertexMapper.template map<dim>(e, sfs[i].entity());
//     int global_j = this->vertexMapper.template map<dim>(e, sfs[j].entity());
//    
//    
//     	Ki = sNDat[global_i].K;
//     	Kj = sNDat[global_j].K;
//       	 
//     	for (int kx=0; kx<dim; kx++){
//     		for (int ky=0; ky<dim; ky++){
//     			if (Ki[kx][ky] != Kj[kx][ky])
//     			{
//     				Ki[kx][ky] = 2 / (1/(Ki[kx][ky]+eps) + (1/(Kj[kx][ky]+eps)));
//     			}
//     		}
//     	}
//     	return Ki;
//    }

