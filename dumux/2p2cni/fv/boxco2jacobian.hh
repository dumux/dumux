#ifndef DUNE_BOXCO2JACOBIAN_HH
#define DUNE_BOXCO2JACOBIAN_HH

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
#include"dumux/operators/boxjacobian.hh"
#include"dumux/2p2cni/2p2cniproblem.hh"
//#include "varswitch.hh"
#include"dumux/io/vtkmultiwriter.hh"

/**
 * @file
 * @brief  compute local jacobian matrix for box scheme for two-phase two-component flow equation
 * @author Bernd Flemisch, Klaus Mosthaf, Melanie Darcis
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
  class BoxCO2Jacobian 
    : public BoxJacobian<BoxCO2Jacobian<G,RT,BoxFunction>,G,RT,3,BoxFunction>
  {
    typedef typename G::ctype DT;
    typedef typename G::Traits::template Codim<0>::Entity Entity;
    typedef typename Entity::Geometry Geometry;
    typedef BoxCO2Jacobian<G,RT,BoxFunction> ThisType;
    typedef typename LocalJacobian<ThisType,G,RT,3>::VBlockType VBlockType;
    typedef typename LocalJacobian<ThisType,G,RT,3>::MBlockType MBlockType;

 	enum {pWIdx = 0, switchIdx = 1, teIdx=2, numberOfComponents = 2};	// Solution vector index
	enum {wPhase = 0, nPhase = 1};					// Phase index
	enum {gasPhase = 0, waterPhase = 1, bothPhases = 2};	// Phase state
	enum {water = 0, air = 1, heat = 2};				// Component index					
	
  public:
    // define the number of phases (m) and components (c) of your system, this is used outside
    // to allocate the correct size of (dense) blocks with a FieldMatrix
    enum {dim=G::dimension};
    enum {m=3, c=2};
    enum {SIZE=LagrangeShapeFunctionSetContainer<DT,RT,dim>::maxsize};
    struct VariableNodeData;

    typedef FieldMatrix<RT,dim,dim> FMatrix;
    typedef FieldVector<RT,dim> FVector;
    
    //! Constructor
    BoxCO2Jacobian (TwoPTwoCNIProblem<G,RT>& params,
			      bool levelBoundaryAsDirichlet_, const G& grid, 
			      BoxFunction& sol, 
			      bool procBoundaryAsDirichlet_=true)
    : BoxJacobian<ThisType,G,RT,m,BoxFunction>(levelBoundaryAsDirichlet_, grid, sol, procBoundaryAsDirichlet_), 
      problem(params), 
      sNDat(this->vertexMapper.size()), vNDat(SIZE), oldVNDat(SIZE)
  	{
      this->analytic = false;
    }

	/** @brief compute time dependent term (storage), loop over nodes / subcontrol volumes
	 *  @param e entity   
	 *  @param sol solution vector
	 *  @param node local node id
	 *  @return storage term
	 */
    virtual VBlockType computeM (const Entity& e, const VBlockType* sol, 
    		int node, std::vector<VariableNodeData>& varData)
    {
    	 GeometryType gt = e.geometry().type();
    	 const typename LagrangeShapeFunctionSetContainer<DT,RT,dim>::value_type&
                	 sfs=LagrangeShapeFunctions<DT,RT,dim>::general(gt,1);
    	 
   	 int globalIdx = this->vertexMapper.template map<dim>(e, sfs[node].entity());

   	 VBlockType result; 
   	 RT satN = varData[node].satN;
   	 RT satW = varData[node].satW;
   	    	    	                  
   	 // storage of component water
   	 result[water] = 
   		  sNDat[globalIdx].porosity*(varData[node].density[wPhase]*satW*varData[node].massfrac[water][wPhase]
   		                 +varData[node].density[nPhase]*satN*varData[node].massfrac[water][nPhase]);
   	 // storage of component air
   	 result[air] = 
   		  sNDat[globalIdx].porosity*(varData[node].density[nPhase]*satN*varData[node].massfrac[air][nPhase]
   	                    +varData[node].density[wPhase]*satW*varData[node].massfrac[air][wPhase]);

   	 // storage term of energy equation
   	 result[heat] = sNDat[globalIdx].porosity * (varData[node].density[wPhase] * varData[node].intenergy[wPhase] * satW
   	             + varData[node].density[nPhase] * varData[node].intenergy[nPhase] * satN)
   	             + (1. -  sNDat[globalIdx].porosity) * elData.soilDens * elData.soilCs * varData[node].temperature;
   	  // soil properties defined at the elements!!!               
   	 
   	 //std::cout << result << " " << node << std::endl;
   	 return result;
    };

    virtual VBlockType computeM (const Entity& e, const VBlockType* sol, int node, bool old = false) 
     {
     	if (old)
     		return computeM(e, sol, node, oldVNDat);
     	else 
     		return computeM(e, sol, node, vNDat);
     }

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
	 GeometryType gt = e.geometry().type();
	 const typename LagrangeShapeFunctionSetContainer<DT,RT,dim>::value_type&
 	 sfs=LagrangeShapeFunctions<DT,RT,dim>::general(gt,1);
 	 
 	 // global index of the subcontrolvolume face neighbour nodes in element e
	 int globalIdx_i = this->vertexMapper.template map<dim>(e, sfs[i].entity());
  	 int globalIdx_j = this->vertexMapper.template map<dim>(e, sfs[j].entity());

  	 // get global coordinates of nodes i,j
	 const FieldVector<DT,dim> global_i = this->fvGeom.subContVol[i].global;
	 const FieldVector<DT,dim> global_j = this->fvGeom.subContVol[j].global;

	 FieldMatrix<RT,2,dim> pGrad(0.), xGrad(0.); 
	 FieldVector<RT,dim> teGrad(0.); 
	 FieldVector<RT,dim> temp(0.); 
     VBlockType flux(0.);

		 
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
	  	 
	 // compute temperature gradient
		  	 temp = feGrad;
	     	 temp *= vNDat[k].temperature;
	     	 teGrad += temp;
	     	 
	 // compute concentration gradients 
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
	 	 

	  // effective permeability in edge direction 
      // RT Kij = sIPDat[global_j].K_eff[face]; 
      // calculate harmonic mean of permeabilities of nodes i and j
		 const FMatrix K = harmonicMeanK(global_i, global_j);
		 //const FMatrix K = harmonicMeanK(e, face);		  
     	 //K.umv(normal, Kij);  // Kij=K*n

	 // Darcy velocity of each phase
		 VBlockType outward(0); 
		 FieldVector<RT,dim> v_tilde(0);

		 // calculate the advective flux using upwind: K*n(grad p -rho*g)
		 for (int phase=0; phase<2; phase++) 
		 	{
	     	 K.umv(pGrad[phase], v_tilde);  // v_tilde=K*gradP
	     	 outward[phase] = v_tilde*normal;
		 	}
	 // Heat conduction 	 
		 // Harmonic mean:
			// Heat Conductivity
	     RT lambda;
	     lambda = 2./((1./vNDat[i].lambda) + (1./vNDat[j].lambda));
		 outward[heat] = teGrad * normal;
	  	 outward[heat] *= lambda;  
		 
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
 		    // Arithmetic mean of phase enthalpies (dissolved components neglected
 	     	RT enthW;
 	     	RT enthCO2;
 	     	enthW = (vNDat[i].enthalpy[pWIdx] + vNDat[j].enthalpy[pWIdx]) / 2.;
 	     	enthCO2 = (vNDat[i].enthalpy[switchIdx] + vNDat[j].enthalpy[switchIdx]) / 2.;
 	  
 		 flux[heat]   =  (alpha* vNDat[up_n].density[nPhase]*vNDat[up_n].mobility[nPhase]
 				    + (1-alpha)* vNDat[dn_n].density[nPhase]*vNDat[dn_n].mobility[nPhase])
 				               * enthCO2 * outward[nPhase];
 		 flux[heat]  +=  (alpha* vNDat[up_w].density[wPhase]*vNDat[up_w].mobility[wPhase]
 				    + (1-alpha)* vNDat[dn_w].density[wPhase]*vNDat[dn_w].mobility[wPhase])
 				    		   * enthW * outward[wPhase];
 		 flux[heat]	+=	outward[heat];

		 
		//Diffusive flux of components
      	// Calculate porous media diffusion coefficient
 		 	// Tortuosity
 		 	RT tauW_i, tauW_j, tauN_i, tauN_j; // tortuosity of wetting and nonwetting phase
 		 	tauW_i = pow(sNDat[globalIdx_i].porosity * vNDat[i].satW,(7/3))/
 				 (sNDat[globalIdx_i].porosity*sNDat[globalIdx_i].porosity);
 		 	tauW_j = pow(sNDat[globalIdx_j].porosity * vNDat[j].satW,(7/3))/
 		 		 (sNDat[globalIdx_j].porosity*sNDat[globalIdx_j].porosity);
 		 	tauN_i = pow(sNDat[globalIdx_i].porosity * vNDat[i].satN,(7/3))/
 				 (sNDat[globalIdx_i].porosity*sNDat[globalIdx_i].porosity);
 		 	tauN_j = pow(sNDat[globalIdx_j].porosity * vNDat[j].satN,(7/3))/
 				 (sNDat[globalIdx_j].porosity*sNDat[globalIdx_j].porosity);
 		     	 
 		 	// Get molecular diffusion coefficient
 		 	FieldVector<RT,2> D_i,D_j;
 		 	D_i = problem.D(global_i);
 		 	D_j = problem.D(global_j);
      	 
 		 	// Porous media diffusion coefficient
 		 	RT Dwg, Daw;
 		 	Dwg = (sNDat[globalIdx_i].porosity * vNDat[i].satN * tauN_i * D_i[water] +
 		 		sNDat[globalIdx_j].porosity * vNDat[j].satN * tauN_j * D_j[water])/2;
 		 	Daw = (sNDat[globalIdx_i].porosity * vNDat[i].satW * tauW_i * D_i[air] +
 				sNDat[globalIdx_j].porosity * vNDat[j].satW * tauW_j * D_j[air])/2;
		 
 			 // Calculate arithmetic mean of the densities 
 		 	 VBlockType avgDensity;
 			 avgDensity[wPhase] = 0.5*(vNDat[i].density[wPhase] + vNDat[j].density[wPhase]);
 			 avgDensity[nPhase] = 0.5*(vNDat[i].density[nPhase] + vNDat[j].density[nPhase]);
 			 
 		 // Diffusive flux
 		 VBlockType normDiffGrad;
 		 
 		  // get local to global id map
		 int state_i = sNDat[globalIdx_i].phaseState;
		 int state_j = sNDat[globalIdx_j].phaseState;
		 
		 RT diffusionAW(0.0), diffusionWW(0.0), diffusionWN(0.0), diffusionAN(0.0);
		 

		 normDiffGrad[wPhase] = xGrad[wPhase]*normal;
		 normDiffGrad[nPhase] = xGrad[nPhase]*normal;

		 if (state_i==bothPhases && state_j==bothPhases)
		 {
			 diffusionAW = Daw * avgDensity[wPhase] * normDiffGrad[wPhase];
			 diffusionWW = - diffusionAW;
			 diffusionWN = Dwg * avgDensity[nPhase] * normDiffGrad[nPhase];
			 diffusionAN = - diffusionWN;
		 }
		 else if (state_i == waterPhase || state_j == waterPhase)
		 {
			 diffusionAW = Daw * avgDensity[wPhase] * normDiffGrad[wPhase];
			 diffusionWW = - diffusionAW;
		 }
		 else if (state_i == gasPhase || state_j == gasPhase)
		 {
			 diffusionWN = Dwg * avgDensity[nPhase] * normDiffGrad[nPhase];
			 diffusionAN = - diffusionWN;		 
		 }

//		 std::cout << "Diffusive Flux: " << diffusionWG << std::endl; 
		 // Add water diffusion to flux
		 flux[water] += (diffusionWW + diffusionWN); 
//		 std::cout << "Water Flux: " << flux[water] << std::endl; 
		 // Air diffusion not implemeted
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
   	virtual void primaryVarSwitch (const Entity& e, int global, VBlockType* sol, int local)
    {
 		
   		this->fvGeom.update(e); // if switch is only called from boxco2 after each time step
   		int state = sNDat[global].phaseState;
   		bool switched = false;

        const FieldVector<RT, 4> parameters = problem.materialLawParameters 
        (this->fvGeom.cellGlobal, e, this->fvGeom.cellLocal);
        
        RT pW = sol[local][pWIdx];
        RT satW = 0.0;
        if (state == bothPhases) satW = 1.0-sol[local][switchIdx];
  		if (state == waterPhase) satW = 1.0;
  		if (state == gasPhase) satW = 0.0;

  		RT pC = problem.materialLaw().pC(satW, parameters);
  		RT pN = pW + pC;

    	FVector Coordinates = this->fvGeom.subContVol[local].global;
    	
        switch(state) 
        {
        case gasPhase :
        	RT xWNmass;
        	xWNmass = sol[local][switchIdx];
   
        	if (xWNmass > 0.001 && switched == false)
            {
            	// appearance of water phase
            	std::cout << "Water appears at node " << global << "  Coordinates: " << Coordinates << std::endl;
            	sNDat[global].phaseState = bothPhases;
            	sol[local][switchIdx] = 1.0 - 1.e-6; // initialize solution vector
            	switched = true;
            }
            break;

        case waterPhase :
        	RT xAWmax, xAWmass;
         	xAWmass = sol[local][switchIdx];
         	xAWmax = problem.multicomp().xAW(pN, sol[local][teIdx]);
             	
        	if (xAWmass > xAWmax && switched == false)
            {
            	// appearance of gas phase
            	std::cout << "Gas appears at node " << global << "  Coordinates: " << Coordinates << std::endl;
            	sNDat[global].phaseState = bothPhases;
            	sol[local][switchIdx] = 1.e-6; // initialize solution vector
            	switched = true;
            }
            break;

        case bothPhases :
        	RT satN = sol[local][switchIdx];
      	  	if (satN < 0.0  && switched == false)
      	  	{
      		  	// disappearance of gas phase
      		  	std::cout << "Gas disappears at node " << global << "  Coordinates: " << Coordinates << std::endl;
      		  	sNDat[global].phaseState = waterPhase;
      		  	sol[local][switchIdx] = 1e-6; // initialize solution vector
      		  	switched = true;
            }
      	  	else if (satW < 0.0  && switched == false)
      	  	{
      	  		// disappearance of water phase
      	  		std::cout << "Water disappears at node " << global << "  Coordinates: " << Coordinates << std::endl;
      	  		sNDat[global].phaseState = gasPhase;
      	  		sol[local][switchIdx] = 1e-6; // initialize solution vector
            	switched = true;
      	  	}
      	  	break;    	  	
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
    	}
   	 return;
   	}
    
    // updates old phase state after each time step
    virtual void updatePhaseState ()
    {
      	for (int i = 0; i < this->vertexMapper.size(); i++){
       		sNDat[i].oldPhaseState = sNDat[i].phaseState;
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
 			  primaryVarSwitch(e, globalIdx, sol, k);
  			  
  			  // mark elements that were already visited
  			  sNDat[globalIdx].visited = true;
 		  }
  	  }
  	  
	  return;
    }

    
   // for initialization of the Static Data (sets porosity)
    virtual void initiateStaticData (const Entity& e)
    {
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

//	void printVariableData()
//	{
//		for (int i = 0; i < 4; i++)
//		{
//			std::cout << "new: i = " << i << ": satN = " << vNDat[i].satN << ", satW = " << vNDat[i].satW 
//				<< ", pW = " << vNDat[i].pW << ", pC = " << vNDat[i].pC << ", pN = " << vNDat[i].pN 
//				<< ", T = " << vNDat[i].temperature << ", lambda = " << vNDat[i].lambda << std::endl; 
//			std::cout << "old: i = " << i << ": satN = " << oldVNDat[i].satN << ", satW = " << oldVNDat[i].satW 
//				<< ", pW = " << oldVNDat[i].pW << ", pC = " << oldVNDat[i].pC << ", pN = " << oldVNDat[i].pN 
//				<< ", T = " << oldVNDat[i].temperature << ", lambda = " << oldVNDat[i].lambda << std::endl; 
//		}
//	}
    

    
    struct VariableNodeData  
    {
   	 RT satN;
     RT satW;
     RT pW;
     RT pC;
     RT pN;
     RT temperature;
     RT lambda;
     FieldVector<RT,2> mobility;  //Vector with the number of phases
     FieldVector<RT,2> density;
     FieldMatrix<RT,c,2> massfrac;
     FieldVector<RT,2> enthalpy;
     FieldVector<RT,2> intenergy;
    };

    // analog to EvalPrimaryData in MUFTE, uses members of vNDat
	virtual void updateVariableData(const Entity& e, const VBlockType* sol, 
			int i, std::vector<VariableNodeData>& varData, int state) 
    {
     const FieldVector<RT, 4> parameters = problem.materialLawParameters 
     (this->fvGeom.cellGlobal, e, this->fvGeom.cellLocal);

   	   	 const int global = this->vertexMapper.template map<dim>(e, i);


   		 varData[i].pW = sol[i][pWIdx];
   		 if (state == bothPhases) varData[i].satN = sol[i][switchIdx];
   		 if (state == waterPhase) varData[i].satN = 0.0;
   		 if (state == gasPhase) varData[i].satN = 1.0;

   		 varData[i].satW = 1.0 - varData[i].satN;

   		 varData[i].pC = problem.materialLaw().pC(varData[i].satW, parameters);
   		 varData[i].pN = varData[i].pW + varData[i].pC;
   		 varData[i].temperature = sol[i][teIdx]; // in [K]

   		 // Solubilities of components in phases
   		 if (state == bothPhases){
   	   		 varData[i].massfrac[air][wPhase] = problem.multicomp().xAW(varData[i].pN, varData[i].temperature);
   	   		 varData[i].massfrac[water][nPhase] = problem.multicomp().xWN(varData[i].pN, varData[i].temperature);
   		 }
   		 if (state == waterPhase){
   	   		 varData[i].massfrac[water][nPhase] = 0.0;
   	   		 varData[i].massfrac[air][wPhase] =  sol[i][switchIdx];
   		 }
   		 if (state == gasPhase){
   	   		 varData[i].massfrac[water][nPhase] = sol[i][switchIdx];
   	   		 varData[i].massfrac[air][wPhase] = 0.0;
   		 }
   	   	 varData[i].massfrac[water][wPhase] = 1.0 - varData[i].massfrac[air][wPhase];
   	   	 varData[i].massfrac[air][nPhase] = 1.0 - varData[i].massfrac[water][nPhase];
   	   	 // for output
//   	   	 (*hackySaturationN)[global] = varData[i].satN;
//   	   	 (*hackyMassFracAir)[global] = varData[i].massfrac[air][wPhase];
//   	   	 (*hackyMassFracWater)[global] = varData[i].massfrac[water][nPhase];
   	   	 // Diffusion coefficients
   		 // Mobilities & densities

   		 varData[i].density[wPhase] = problem.materialLaw().wettingPhase.density(varData[i].temperature, varData[i].pW, varData[i].massfrac[air][wPhase]);
   		 varData[i].density[nPhase] = problem.materialLaw().nonwettingPhase.density(varData[i].temperature, varData[i].pN);
   		 varData[i].mobility[wPhase] = problem.materialLaw().mobW(varData[i].satW, parameters, varData[i].temperature, varData[i].pW);
   		 varData[i].mobility[nPhase] = problem.materialLaw().mobCO2(varData[i].satN, parameters, varData[i].temperature, varData[i].pN, varData[i].density[nPhase]);
   		 varData[i].lambda = soil.heatConductivity(elData.soilLDry, elData.soilLSw, varData[i].satW);
         varData[i].enthalpy[pWIdx] = problem.materialLaw().wettingPhase.enthalpy(varData[i].temperature,varData[i].pW);
         varData[i].enthalpy[switchIdx] = problem.materialLaw().nonwettingPhase.enthalpy(varData[i].temperature,varData[i].pN);
         varData[i].intenergy[pWIdx] = problem.materialLaw().wettingPhase.intEnergy(varData[i].temperature,varData[i].pW);
         varData[i].intenergy[switchIdx] = problem.materialLaw().nonwettingPhase.intEnergy(varData[i].temperature,varData[i].pN);
         // CONSTANT solubility (for comparison with twophase)
//         varData[i].massfrac[air][wPhase] = 0.0; varData[i].massfrac[water][wPhase] = 1.0;
//         varData[i].massfrac[water][nPhase] = 0.0; varData[i].massfrac[air][nPhase] = 1.0;

         //std::cout << "water in gasphase: " << varData[i].massfrac[water][nPhase] << std::endl;
         //std::cout << "air in waterphase: " << varData[i].massfrac[air][wPhase] << std::endl;

   		 // for output
   		 (*outPressureN)[global] = varData[i].pN;
   		 (*outCapillaryP)[global] = varData[i].pC;
  	   	 (*outSaturationW)[global] = varData[i].satW;
   	   	 (*outSaturationN)[global] = varData[i].satN;
   	   	 (*outTemperature)[global] = varData[i].temperature;
   	   	 (*outMassFracAir)[global] = varData[i].massfrac[air][wPhase];
   	   	 (*outMassFracWater)[global] = varData[i].massfrac[water][nPhase];
   	   	 (*outDensityW)[global] = varData[i].density[wPhase];
   	   	 (*outDensityN)[global] = varData[i].density[nPhase];
   	   	 (*outMobilityW)[global] = varData[i].mobility[wPhase];
   	   	 (*outMobilityN)[global] = varData[i].mobility[nPhase];
   	   	 (*outPhaseState)[global] = state;

   	   	 return;
    }

	virtual void updateVariableData(const Entity& e, const VBlockType* sol, int i, bool old = false) 
	{
		int state;
		const int global = this->vertexMapper.template map<dim>(e, i);
		if (old)
		{
	   	   	state = sNDat[global].oldPhaseState;
			updateVariableData(e, sol, i, oldVNDat, state);
		}
		else 
		{
		    state = sNDat[global].phaseState;
			updateVariableData(e, sol, i, vNDat, state);
		}
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
//   	 bool switched;
   	 int phaseState;
   	 int oldPhaseState;
   	 RT cellVolume;
   	 RT porosity;
   	 FieldVector<RT, 4> parameters;
   	 FMatrix K;
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
    std::vector<VariableNodeData> oldVNDat;

    // for output files
    BlockVector<FieldVector<RT, 1> > *outPressureN;
    BlockVector<FieldVector<RT, 1> > *outCapillaryP;
    BlockVector<FieldVector<RT, 1> > *outSaturationN;
    BlockVector<FieldVector<RT, 1> > *outSaturationW;
    BlockVector<FieldVector<RT, 1> > *outTemperature;
    BlockVector<FieldVector<RT, 1> > *outMassFracAir;
    BlockVector<FieldVector<RT, 1> > *outMassFracWater;
    BlockVector<FieldVector<RT, 1> > *outDensityW;
    BlockVector<FieldVector<RT, 1> > *outDensityN;
    BlockVector<FieldVector<RT, 1> > *outMobilityW;
    BlockVector<FieldVector<RT, 1> > *outMobilityN;
    BlockVector<FieldVector<RT, 1> > *outPhaseState;
    
  };  
  
}
#endif

//    // average permeability from the staticNode Vector
//    virtual FMatrix harmonicMeanK (const Entity& e, int k) const
//    {
//	 FMatrix Ki, Kj;
//	 const RT eps = 1e-20;
//	 
//     int i = this->fvGeom.subContVolFace[k].i;
//     int j = this->fvGeom.subContVolFace[k].j;
//
//     int global_i = this->vertexMapper.template map<dim>(e, i);
//     int global_j = this->vertexMapper.template map<dim>(e, j);
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

