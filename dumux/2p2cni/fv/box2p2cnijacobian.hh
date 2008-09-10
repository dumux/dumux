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
#include"dumux/operators/boxjacobian.hh"
#include"dumux/2p2cni/2p2cniproblem.hh"
//#include "varswitch.hh"
#include"dumux/io/vtkmultiwriter.hh"

/**
 * @file
 * @brief  compute local jacobian matrix for box scheme for two-phase two-component flow equation
 * @author Bernd Flemisch, Melanie Darcis, Klaus Mosthaf
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

 	enum {pWIdx = 0, switchIdx = 1, teIdx=2, numberOfComponents = 2};	// Solution vector index
	enum {wPhase = 0, nPhase = 1};					// Phase index
	enum {gasPhase = 0, waterPhase = 1, bothPhases = 2};	// Phase state
	enum {water = 0, air = 1, heat =2};				// Component index

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
    Box2P2CNIJacobian (TwoPTwoCNIProblem<G,RT>& params,
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
   	             + elData.heatCap * varData[node].temperature;
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

	 // get global coordinates of nodes i,j
	 const FieldVector<DT,dim> global_i = this->fvGeom.subContVol[i].global;
	 const FieldVector<DT,dim> global_j = this->fvGeom.subContVol[j].global;
	 const FieldVector<DT,dim> local_i = this->fvGeom.subContVol[i].local;
	 const FieldVector<DT,dim> local_j = this->fvGeom.subContVol[j].local;

	 FieldMatrix<RT,2,dim> pGrad(0.), xGrad(0.);
	 FMatrix Ki, Kj;
	 FieldVector<RT,dim> teGrad(0.);
	 FieldVector<RT,dim> temp(0.);
//		 FieldVector<RT,dim> Kij(0);
     	 VBlockType flux(0.);

		 // effective permeability in edge direction
     	 // RT Kij = sIPDat[global_j].K_eff[face];

     	 // calculate harmonic mean of permeabilities of nodes i and j
     	 Ki = this->problem.soil().K(global_i,e,local_i);
     	 Kj = this->problem.soil().K(global_j,e,local_j);
		 const FMatrix K = harmonicMeanK(Ki, Kj);
		 //const FMatrix K = harmonicMeanK(e, face);
     	 //K.umv(normal, Kij);  // Kij=K*n

	 // Harmonic mean:
		// Heat Conductivity
     	RT lambda;
     	lambda = 2./((1./vNDat[i].lambda) + (1./vNDat[j].lambda));

     // Arithmetic mean:
     	// Enthalpy
     	RT avgEnthW;
     	RT avgEnthN;
     	avgEnthW = (vNDat[i].enthalpy[wPhase] + vNDat[j].enthalpy[wPhase]) / 2.;
     	avgEnthN = (vNDat[i].enthalpy[nPhase] + vNDat[j].enthalpy[nPhase]) / 2.;



		 // calculate FE gradient (grad p for each phase)
		 for (int k = 0; k < this->fvGeom.nNodes; k++) // loop over adjacent nodes
		 {
			 // FEGradient at node k
			 const FieldVector<DT,dim> feGrad(this->fvGeom.subContVolFace[face].grad[k]);
			 FieldVector<RT,2> pressure(0.0), massfrac(0.0);

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
 				               * avgEnthN * outward[nPhase];
 		 flux[heat]  +=  (alpha* vNDat[up_w].density[wPhase]*vNDat[up_w].mobility[wPhase]
 				    + (1-alpha)* vNDat[dn_w].density[wPhase]*vNDat[dn_w].mobility[wPhase])
 				    		   * avgEnthW * outward[wPhase];
 		 flux[heat]	+=	outward[heat];


		 // DIFFUSION
		 VBlockType normDiffGrad;

	   	 // get local to global id map
		 int state_i = vNDat[i].phasestate;
		 int state_j = vNDat[j].phasestate;

	  	 RT diffusionWW(0.0), diffusionWN(0.0); // diffusion of water
		 RT diffusionAW(0.0), diffusionAN(0.0); // diffusion of air
		 VBlockType avgDensity, avgDpm;
		 avgDpm[wPhase]=1e-9; // needs to be changed !!!
		 avgDpm[nPhase]=1e-5; // water in the gasphase

		 normDiffGrad[wPhase] = xGrad[wPhase]*normal;
		 normDiffGrad[nPhase] = xGrad[nPhase]*normal;

		 // calculate the arithmetic mean of densities
		 avgDensity[wPhase] = 0.5*(vNDat[i].density[wPhase] + vNDat[j].density[wPhase]);
		 avgDensity[nPhase] = 0.5*(vNDat[i].density[nPhase] + vNDat[j].density[nPhase]);


		 if (state_i==2 && state_j==2)
		 {
			 diffusionAW = avgDpm[wPhase] * avgDensity[wPhase] * normDiffGrad[wPhase];
			 diffusionWW = - diffusionAW;
			 diffusionWN = avgDpm[nPhase] * avgDensity[nPhase] * normDiffGrad[nPhase];
			 diffusionAN = - diffusionWN;
		 }
		 else if ((state_i == 1 || state_j == 1) || (state_i == 1 && state_j == 1))
		 {
			 diffusionAW = avgDpm[wPhase] * avgDensity[wPhase] * normDiffGrad[wPhase];
			 diffusionWW = - diffusionAW;
		 }
		 else if ((state_i == 0 || state_j == 0) || (state_i == 0 && state_j == 0))
		 {
			 diffusionWN = avgDpm[nPhase] * avgDensity[nPhase] * normDiffGrad[nPhase];
			 diffusionAN = - diffusionWN;
		 }

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
   	virtual void primaryVarSwitch (const Entity& e, int global, VBlockType* sol, int local)
    {

        bool switched = false;
   		int state = sNDat[global].phaseState;

//        const FieldVector<RT, 4> parameters = problem.materialLawParameters
//        (this->fvGeom.cellGlobal, e, this->fvGeom.cellLocal);

        RT pW = sol[local][pWIdx];
        RT satW = 0.0;
        if (state == bothPhases) satW = 1.0-sol[local][switchIdx];
  		if (state == waterPhase) satW = 1.0;
  		if (state == gasPhase) satW = 0.0;

    	RT pC = problem.materialLaw().pC(satW, this->fvGeom.cellGlobal, e, this->fvGeom.cellLocal);
  		RT pN = pW + pC;

    	FVector Coordinates = this->fvGeom.subContVol[local].global;

        switch(state)
        {
        case gasPhase :
        	RT xWNmass, xWNmolar, pwn, pWSat;
        	xWNmass = sol[local][switchIdx];
        	xWNmolar = problem.multicomp().convertMassToMoleFraction(xWNmass, gasPhase);
           	pwn = xWNmolar * pN;
            pWSat = problem.multicomp().vaporPressure(sol[local][teIdx]);

        	if (pwn > 1.01*pWSat && switched == false)
            {
            	// appearance of water phase
            	std::cout << "Water appears at node " << global << "  Coordinates: " << Coordinates << std::endl;
            	sNDat[global].phaseState = bothPhases;
            	sol[local][switchIdx] = 1.0 - 1.e-6; // initialize solution vector
            	switched = true;
            }
            break;

        case waterPhase :
        	RT pbub, xAWmass, xAWmolar, henryInv;
           	xAWmass = sol[local][switchIdx];
         	xAWmolar = problem.multicomp().convertMassToMoleFraction(xAWmass, waterPhase);
        	henryInv = problem.multicomp().henry(sol[local][teIdx]);
            pWSat = problem.multicomp().vaporPressure(sol[local][teIdx]);
        	pbub = pWSat + xAWmolar/henryInv;

        	if (pbub > pN && switched == false)
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
        if (switched == true) updateVariableData(e, sol, local, vNDat, sNDat[global].phaseState);

   	return;
    }

    // harmonic mean computed directly
    virtual FMatrix harmonicMeanK (FMatrix& Ki, const FMatrix& Kj) const
    {
    	double eps = 1e-20;
//    	FMatrix K;

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


 		 elData.heatCap = problem.soil().heatCap(this->fvGeom.cellGlobal, e, this->fvGeom.cellLocal);


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
  			  sNDat[globalIdx].porosity = problem.soil().porosity(this->fvGeom.cellGlobal, e, this->fvGeom.cellLocal);

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
  			  sNDat[globalIdx].porosity = problem.soil().porosity(this->fvGeom.cellGlobal, e, this->fvGeom.cellLocal);

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
     int phasestate;
    };

    // analog to EvalPrimaryData in MUFTE, uses members of vNDat
	virtual void updateVariableData(const Entity& e, const VBlockType* sol,
			int i, std::vector<VariableNodeData>& varData, int state)
    {
   	   	 const int global = this->vertexMapper.template map<dim>(e, i);

   		 varData[i].pW = sol[i][pWIdx];
   		 if (state == bothPhases) varData[i].satN = sol[i][switchIdx];
   		 if (state == waterPhase) varData[i].satN = 0.0;
   		 if (state == gasPhase) varData[i].satN = 1.0;

   		 varData[i].satW = 1.0 - varData[i].satN;

   		 varData[i].pC = problem.materialLaw().pC(varData[i].satW, this->fvGeom.subContVol[i].global, e, this->fvGeom.subContVol[i].local);
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
   	   	 varData[i].phasestate = state;

   		 // Mobilities & densities
   		 varData[i].mobility[wPhase] = problem.materialLaw().mobW(varData[i].satW, this->fvGeom.subContVol[i].global, e, this->fvGeom.subContVol[i].local, varData[i].temperature, varData[i].pW);
   		 varData[i].mobility[nPhase] = problem.materialLaw().mobN(varData[i].satN, this->fvGeom.subContVol[i].global, e, this->fvGeom.subContVol[i].local, varData[i].temperature, varData[i].pN);
   		 varData[i].density[wPhase] = problem.materialLaw().wettingPhase.density(varData[i].temperature, varData[i].pN);
   		 varData[i].density[nPhase] = problem.materialLaw().nonwettingPhase.density(varData[i].temperature, varData[i].pN,
   				 varData[i].massfrac[air][nPhase]);
         varData[i].lambda = problem.soil().heatCond(this->fvGeom.subContVol[i].global, e, this->fvGeom.subContVol[i].local, varData[i].satW);
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
   	   	 (*outPhaseState)[global] = varData[i].phasestate;

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
     RT heatCap;

//   	 RT cellVolume;
//     	 RT porosity;
//   	 RT gravity;
//   	 FieldVector<RT, 4> parameters;
//   	 FieldMatrix<RT,dim,dim> K;
   	 } elData;


    // parameters given in constructor
   	TwoPTwoCNIProblem<G,RT>& problem;
    CWaterAir multicomp;
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

