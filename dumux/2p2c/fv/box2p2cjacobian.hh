#ifndef DUNE_BOX2P2CJACOBIAN_HH
#define DUNE_BOX2P2CJACOBIAN_HH

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
//#include"dumux/2p2c/boxjacobian_2p2c.hh"
#include"dumux/2p2c/2p2cproblem.hh"
//#include "varswitch.hh"

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
  template<class G, class RT, class BoxFunction = LeafP1Function<G, RT, 2> >
  class Box2P2CJacobian 
    : public BoxJacobian<Box2P2CJacobian<G,RT,BoxFunction>,G,RT,2,BoxFunction>
  {
    typedef typename G::ctype DT;
    typedef typename G::Traits::template Codim<0>::Entity Entity;
    typedef typename Entity::Geometry Geometry;
    typedef Box2P2CJacobian<G,RT,BoxFunction> ThisType;
    typedef typename LocalJacobian<ThisType,G,RT,2>::VBlockType VBlockType;
    typedef typename LocalJacobian<ThisType,G,RT,2>::MBlockType MBlockType;
    typedef FVElementGeometry<G> FVElementGeometry;

 	enum {pWIdx = 0, satNIdx = 1, numberOfComponents = 2};	// Solution vector index
	enum {wPhase = 0, nPhase = 1};									// Phase index
	enum {gasPhase = 0, waterPhase = 1, bothPhases = 2};		// Phase state
	enum {water = 0, air = 1};										// Component index					
	
  public:
    // define the number of phases (m) and components (c) of your system, this is used outside
    // to allocate the correct size of (dense) blocks with a FieldMatrix
    enum {dim=G::dimension};
    enum {m=2, c=2};
    enum {SIZE=LagrangeShapeFunctionSetContainer<DT,RT,dim>::maxsize};

    typedef FieldMatrix<RT,dim,dim> FMatrix;
    typedef FieldVector<RT,dim> FVector;
    
    //! Constructor
    Box2P2CJacobian (TwoPTwoCProblem<G,RT>& params,
			      bool levelBoundaryAsDirichlet_, const G& grid, 
			      BoxFunction& sol, 
			      bool procBoundaryAsDirichlet_=true)
    : BoxJacobian<ThisType,G,RT,2,BoxFunction>(levelBoundaryAsDirichlet_, grid, sol, procBoundaryAsDirichlet_), 
      problem(params), 
      statNData(this->vertexMapper.size()),
      statIPData(this->vertexMapper.size())
    {
      this->analytic = false;
    }


    virtual void primaryVarSwitch (const Entity& e, VBlockType* sol, int i)
    {
   	 RT Xaw, Xwg;
        //int i = this->fvGeom.cellLocal;

        // neue Variable in SOL schreiben. Da Zeiger auf sol, sollte das funktionieren
        // Achtung: alte Loesung beibebehalten, falls Newton den Zeitschritt nicht schafft

        const RT Xaw_max = 1.E-5; // Needs to be changed/specified !!!
        const RT Xwg_max = 2.E-1; // Needs to be changed/specified !!!
        
        Xaw   	 = varNData[i].massfrac[air][wPhase];
        Xwg     = varNData[i].massfrac[water][nPhase];;
        int state = statNData[i].phaseState;
        bool switched = false;
        
    	 GeometryType gt = e.geometry().type();
    	 const typename Dune::LagrangeShapeFunctionSetContainer<DT,RT,dim>::value_type& 
    	 sfs=Dune::LagrangeShapeFunctions<DT,RT,dim>::general(gt,1);

    	 // get local to global id map
  		 int globalIdx = this->vertexMapper.template map<dim>(e, sfs[i].entity());

        //std::cout << "negative non-wetting saturation!! --> Switsching Variables" << std::endl;
       
        switch(state) 
        {
        case gasPhase :

            if (varNData[i].massfrac[water][nPhase] > 1.01*Xwg_max && switched == false)
            {
            	// appearance of water phase
            	std::cout << "Water appears at node " << globalIdx << std::endl;
            	statNData[i].phaseState = bothPhases;
            	varNData[i].saturationN = 1 - 1.E-6; // initialize
            	switched = true;
            }
            break;

        case waterPhase :

            if (varNData[i].massfrac[air][wPhase] > 1.01*Xaw_max && switched == false)
            {
            	// appearance of gas phase
            	std::cout << "Gas appears at node " << globalIdx << std::endl;
            	statNData[i].phaseState = bothPhases;
            	varNData[i].saturationN = 1.E-6;  // Initialisierung
            	switched = true;
            }
            break;

        case bothPhases :

      	  	if ( varNData[i].saturationN < 0.0  && switched == false) 
      	  	{
      		  	// disappearance of gas phase
      		  	std::cout << "Gas disappears" << globalIdx << std::endl;
      		  	statNData[i].phaseState = waterPhase;
      		  	varNData[i].massfrac[air][wPhase] = 1.E-6;  // Initialisierung
            	switched = true;
            }
      	  	else if ( varNData[i].saturationW < 0.0  && switched == false) 
      	  	{
      	  		// disappearance of water phase
      	  		std::cout << "Water disappears" << globalIdx << std::endl;
      	  		statNData[i].phaseState = gasPhase;
      	  		varNData[i].massfrac[water][nPhase] = 1.E-6;  // Initialisierung
            	switched = true;
      	   }
      	  	break;
        }
        
        if (switched == true)
        {
	        // update primary variable
	        if (statNData[i].phaseState == gasPhase)  
	        {
	      	  sol[i][satNIdx] = varNData[i].massfrac[water][nPhase];
	        }
	        if (statNData[i].phaseState == waterPhase)  
	        {
	      	  sol[i][satNIdx] = varNData[i].massfrac[air][wPhase];
	        }
	        if (statNData[i].phaseState == bothPhases)  
	        {
	      	  sol[i][satNIdx] = varNData[i].saturationN;
	        }
        }

   	return;
    }
    
    // calculates the harmonic mean of K, if values are different
//    virtual double harmonicMeanK (int face, FieldVector<DT,dim> global_i, FieldVector<DT,dim> global_j)
//    {
//     	 RT auxPerm_i, auxPerm_j;
//     	 FieldVector<RT,dim> normal(this->fvGeom.subContVolFace[face].normal);
//   	 FieldVector<DT,dim> Ki, Kj;
//   	 FieldMatrix<DT,dim,dim> K;
//   	 
//   	 K = problem.K(global_i);
//   	 K.umv(normal, Ki);  // Kij=K*n
//   	 auxPerm_i = Ki*normal/(normal*normal);
//   	 
//   	 K = problem.K(global_j);
//   	 K.umv(normal, Kj);  // Kij=K*n
//   	 auxPerm_j = Ki*normal/(normal*normal);   	 
//
//     	 if (auxPerm_i != auxPerm_j){
//     		 if (auxPerm_i==0.0 || auxPerm_j==0.0)
//     			 auxPerm_i = 0.0;
//     		 else 
//     			 auxPerm_i = 2 / ((1/auxPerm_i)+(1/auxPerm_j));
//   	 }
//
//   	 return auxPerm_i;
//    }

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
    
    // or from the staticNode Vector
    virtual FMatrix harmonicMeanK (const Entity& e, int k) const
    {
	 FMatrix Ki, Kj;
	 const RT eps = 1e-20;
	 
    GeometryType gt = e.geometry().type();
    const typename Dune::LagrangeShapeFunctionSetContainer<DT,RT,dim>::value_type& 
    sfs=Dune::LagrangeShapeFunctions<DT,RT,dim>::general(gt,1);

     int i = this->fvGeom.subContVolFace[k].i;
	  int j = this->fvGeom.subContVolFace[k].j;

	  int global_i = this->vertexMapper.template map<dim>(e, sfs[i].entity());
	  int global_j = this->vertexMapper.template map<dim>(e, sfs[j].entity());
    
    

       	 Ki = statNData[global_i].K;
       	 Kj = statNData[global_j].K;
       	 
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
   	 for (int i = 0; i < this->vertexMapper.size(); i++)
   		statNData[i].visited = false;
   	 
   	 return;
    }

    void getLocalDefect(const Entity& entity,VBlockType *defhelp)
    { 
      setLocalSolution(entity);

      // set to Zero 
	  	for (int i=0; i < this->fvGeom.nNodes; i++) {
	  		this->bctype[i].assign(BoundaryConditions::neumann);
	  		this->b[i] = 0;
	  		this->def[i] = 0;
	  	}
	     
	  	this->template localDefect<LeafTag>(entity,this->u);
		  this->template assembleBC<LeafTag> (entity); 
		  
		  // add to defect 
		  for (int i=0; i < this->fvGeom.nNodes; i++) {
			  this->def[i] += this->b[i];
			  defhelp[i]=this->def[i];
      }
    }

    
    // Compute time dependent term (storage), loop over nodes / subcontrol volumes
    // ACHTUNG varNData always contains values from the NEW timestep
    virtual VBlockType computeM (const Entity& e, const VBlockType* sol, int node)
    {
    	 GeometryType gt = e.geometry().type();
    	 const typename Dune::LagrangeShapeFunctionSetContainer<DT,RT,dim>::value_type&
     	 sfs=Dune::LagrangeShapeFunctions<DT,RT,dim>::general(gt,1);
    	 
   	 int globalIdx = this->vertexMapper.template map<dim>(e, sfs[node].entity());
   	 //int globalCoord = this->fvGeom.subContVol[node].global;

   	 VBlockType result; 
   	 RT satN = sol[node][satNIdx];
   	 RT satW = 1.0 - satN;  
   	    	                  
   	 // storage of component water
   	 result[water] = 
   		 statNData[globalIdx].porosity*(varNData[node].density[wPhase]*satW*varNData[node].massfrac[water][wPhase]
   		                 +varNData[node].density[nPhase]*satN*varNData[node].massfrac[water][nPhase]);
   	 // storage of component air
   	 result[air] = 
   		 statNData[globalIdx].porosity*(varNData[node].density[nPhase]*satN*varNData[node].massfrac[air][nPhase]
   	                    +varNData[node].density[wPhase]*satW*varNData[node].massfrac[air][wPhase]);   
   	 
   	 //std::cout << result << " " << node << std::endl;
   	 return result;
    };

    // loop over subcontrol volume faces
    virtual VBlockType computeA (const Entity& e, const VBlockType* sol, int face)
    {
   	 int i = this->fvGeom.subContVolFace[face].i;
     	 int j = this->fvGeom.subContVolFace[face].j;
     	 
     	 // normal vector, value of the area of the scvf
    	 const FieldVector<RT,dim> normal(this->fvGeom.subContVolFace[face].normal);

    	 // get global coordinates of nodes i,j
    	 const FieldVector<DT,dim> global_i = this->fvGeom.subContVol[i].global;
		 const FieldVector<DT,dim> global_j = this->fvGeom.subContVol[j].global;

		 FieldMatrix<RT,m,dim> pGrad(0.), xGrad(0.); 
		 FieldVector<RT,dim> temp(0.); 
//		 FieldVector<RT,dim> Kij(0);
     	 VBlockType flux(0.);
		 
		 // effective permeability in edge direction 
     	 // RT Kij = statIPData[global_j].K_eff[face]; 
     	 
     	 // calculate harmonic mean of permeabilities of nodes i and j
		 const FMatrix K = harmonicMeanK(global_i, global_j);
		 //const FMatrix K = harmonicMeanK(e, face);		  
     	 //K.umv(normal, Kij);  // Kij=K*n
		 
		 // calculate FE gradient (grad p for each phase)
		 for (int k = 0; k < this->fvGeom.nNodes; k++) // loop over adjacent nodes
		 {	 
			 // FEGradient at node k
			 const FieldVector<DT,dim> feGrad(this->fvGeom.subContVolFace[face].grad[k]);
			 FieldVector<RT,m> pressure(0.0), massfrac(0.0);

			 pressure[wPhase] = sol[k][pWIdx];
			 pressure[nPhase] = varNData[k].pN;
      	 
      	 // compute sum of pressure gradients for each phase
      	 for (int phase = 0; phase < m; phase++)
      	 {	      		 
      		 temp = feGrad;
      		 temp *= pressure[phase];
      		 pGrad[phase] += temp;

          	 // compute sum of concentration gradient
         	 temp = feGrad;
         	 temp *= varNData[k].massfrac[air][phase];
         	 xGrad[phase] += temp;
      	 }

       }

       // deduce gravity*density of each phase
		 FieldMatrix<RT,m,dim> contribComp(0);
		 for (int phase=0; phase<m; phase++)
		 {
			 contribComp[phase] = problem.gravity();
			 contribComp[phase] *= varNData[i].density[phase];  
			 pGrad[phase] -= contribComp[phase]; // grad p - rho*g
		 }
	 	 
		 VBlockType outward(0);  // Darcy velocity of each phase
		 FieldVector<RT,dim> v_tilde(0);

		 // calculate the advective flux using upwind: K*n(grad p -rho*g)
		 for (int phase=0; phase<m; phase++) 
		 	{
			 //pGrad[phase] *= Kij;
			 //outward[phase] = Kij * pGrad[phase];// * normal;
	     	 K.umv(pGrad[phase], v_tilde);  // v_tilde=K*gradP
	     	 outward[phase] = v_tilde*normal;
		 	}
		 
		 // evaluate upwind nodes
		 int up_w, dn_w, up_n, dn_n;
		 if (outward[wPhase] <= 0) {up_w = i; dn_w = j;}
		 else {up_w = j; dn_w = i;};
		 if (outward[nPhase] <= 0) {up_n = i; dn_n = j;}
		 else {up_n = j; dn_n = i;};


		 RT alpha = 1.0;  // Upwind parameter
 		 
		 // Water conservation
 		 flux[water] =   (alpha* varNData[up_w].density[wPhase]*varNData[up_w].mobility[wPhase]
 				               * varNData[up_w].massfrac[water][wPhase] 
 				    + (1-alpha)* varNData[dn_w].density[wPhase]*varNData[dn_w].mobility[wPhase]
 				               * varNData[dn_w].massfrac[water][wPhase])
 				               * outward[wPhase]; 		 		
 		 flux[water] +=  (alpha* varNData[up_n].density[nPhase]*varNData[up_n].mobility[nPhase]
 				               * varNData[up_n].massfrac[water][nPhase] 
 				    + (1-alpha)* varNData[dn_n].density[nPhase]*varNData[dn_n].mobility[nPhase]
 				               * varNData[dn_n].massfrac[water][nPhase])
 				               * outward[nPhase];
 		 // Air conservation
 		 flux[air]   =   (alpha* varNData[up_n].density[nPhase]*varNData[up_n].mobility[nPhase]
 				               * varNData[up_n].massfrac[air][nPhase] 
 				    + (1-alpha)* varNData[dn_n].density[nPhase]*varNData[dn_n].mobility[nPhase]
 				               * varNData[dn_n].massfrac[air][nPhase])
 				               * outward[nPhase];
 		 flux[air]  +=   (alpha* varNData[up_w].density[wPhase]*varNData[up_w].mobility[wPhase]
 				               * varNData[up_w].massfrac[air][wPhase]
 				    + (1-alpha)* varNData[dn_w].density[wPhase]*varNData[dn_w].mobility[wPhase]
 				               * varNData[dn_w].massfrac[air][wPhase])
 				               * outward[wPhase];
		 
		 // DIFFUSION
		 VBlockType normDiffGrad;

		 //FieldMatrix<DT,m,c> diffusionW;
		 RT diffusionAW, diffusionWW; // Diffusion in the water phase
		 RT diffusionWN, diffusionAN; // Diffusion in the non-wetting phase
		 VBlockType avgDensity, avgDpm;
		 avgDpm[wPhase]=1e-9; // needs to be changed !!!
		 avgDpm[nPhase]=1e-9;
		 
		 normDiffGrad[wPhase] = -(xGrad[wPhase]*normal);
		 normDiffGrad[nPhase] = -(xGrad[nPhase]*normal);

		 // calculate the arithmetic mean of densities
		 avgDensity[wPhase] = 0.5*(varNData[i].density[wPhase] + varNData[j].density[wPhase]);
		 avgDensity[nPhase] = 0.5*(varNData[i].density[nPhase] + varNData[j].density[nPhase]);

		 
		 diffusionAW = avgDpm[wPhase] * avgDensity[wPhase] * normDiffGrad[wPhase];
		 diffusionWW = - diffusionAW;
		 diffusionWN = avgDpm[nPhase] * avgDensity[nPhase] * normDiffGrad[nPhase];
		 diffusionAN = - diffusionWN;
		 
//		 std::cout << "Diffusive Flux: " << diffusionWG << std::endl; 
		 

		 // add water diffusion to flux
		 flux[water] -= (diffusionWW + diffusionWN);
//		 std::cout << "Water Flux: " << flux[water] << std::endl; 
		 // air diffusion not implemeted
		 flux[air] -= (diffusionAW + diffusionAN);
//		 std::cout << "Air Flux: " << flux[air] << std::endl; 


		 return flux;
  };
    
    

   // Integrate sources / sinks
   virtual VBlockType computeQ (const Entity& e, const VBlockType* sol, const int& node)
   {
  	 // ASSUME problem.q already contains \rho.q
  	 return problem.q(this->fvGeom.subContVol[node].global, e, this->fvGeom.subContVol[node].local);
   }
       
    
	  //*********************************************************
	  //*																			*
	  //*	Calculation of Data at Elements (elData) 					*
	  //*						 													*
	  //*																		 	*
	  //*********************************************************

    virtual void computeElementData (const Entity& e)
    {
//  		 // ASSUME element-wise constant parameters for the material law 
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
	  //*																			*
	  //*	Calculation of Data at Nodes that has to be			 	*
	  //*	determined only once	(statNData)							 	*
	  //*																		 	*
	  //*********************************************************

    // analog to EvalStaticData in MUFTE
    virtual void updateStaticData (const Entity& e, const VBlockType* sol)
    {
   	 // size of the statNData vector is determined in the constructor
   	 
   	 // local to global id mapping (do not ask vertex mapper repeatedly
   	 //int localToGlobal[Dune::LagrangeShapeFunctionSetContainer<DT,RT,n>::maxsize];

   	 // get access to shape functions for P1 elements
   	 GeometryType gt = e.geometry().type();
   	 const typename Dune::LagrangeShapeFunctionSetContainer<DT,RT,dim>::value_type& 
   	 sfs=Dune::LagrangeShapeFunctions<DT,RT,dim>::general(gt,1);

   	 // get local to global id map
   	 for (int k = 0; k < sfs.size(); k++) {
  		 int globalIdx = this->vertexMapper.template map<dim>(e, sfs[k].entity());
  		  
  		 // if nodes are not already visited
  		 if (!statNData[globalIdx].visited) 
  		  {
  			  // initial phase state
  			 	statNData[globalIdx].phaseState = waterPhase;

  			  // ASSUME parameters defined at the node/subcontrol volume for the material law 
  			  statNData[globalIdx].parameters = problem.materialLawParameters
  			  (this->fvGeom.cellGlobal, e, this->fvGeom.cellLocal);

 			  // global coordinates
 			  FieldVector<DT,dim> global_i = this->fvGeom.subContVol[k].global;
  			  
  			  // ASSUME permeability defined at nodes, evaluate harmonic mean of K 
  			  //statNData[globalIdx].K = problem.K(global_i);  

  			  // ASSUME porosity defined at nodes
  			  statNData[globalIdx].porosity = problem.porosity(this->fvGeom.cellGlobal, e, this->fvGeom.cellLocal);

  			  
  			  // mark elements that were already visited
  			  statNData[globalIdx].visited = true;
  		  }
  	  }
  	  
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
   	 GeometryType gt = e.geometry().type();
   	 const typename Dune::LagrangeShapeFunctionSetContainer<DT,RT,dim>::value_type& 
     sfs=Dune::LagrangeShapeFunctions<DT,RT,dim>::general(gt,1);

   	 for (int i = 0; i < size; i++) {
   		this->def[i] = 0;

   		int globalIdx = this->vertexMapper.template map<dim>(e, sfs[i].entity());
   	 
   		varNData[i].saturationN = sol[i][satNIdx];
   		varNData[i].saturationW = 1.0 - sol[i][satNIdx];
   		varNData[i].pW = sol[i][pWIdx];
   		varNData[i].pC = problem.materialLaw().pC(varNData[i].saturationW, statNData[globalIdx].parameters);
         //std::cout << "Capillary pressure: " << varNData[i].pC << " Corresponding Sw: "<< varNData[i].saturationW << std::endl;
         varNData[i].pN = sol[i][pWIdx] + varNData[i].pC;
         //std::cout << "PressureN: " << sol[i][pWIdx] << std::endl;
         varNData[i].temperature = 283.15; // in [K]
         // Mobilities & densities
         varNData[i].mobility[wPhase] = problem.materialLaw().mobW(varNData[i].saturationW, statNData[globalIdx].parameters);
         varNData[i].mobility[nPhase] = problem.materialLaw().mobN(sol[i][satNIdx], statNData[globalIdx].parameters);
         varNData[i].density[wPhase] = problem.materialLaw().wettingPhase.density();
         varNData[i].density[nPhase] = problem.materialLaw().nonwettingPhase.density();
         // Solubilities of components in phases
         varNData[i].massfrac[air][wPhase] = problem.solu().Xaw(varNData[i].pN, varNData[i].temperature);
         varNData[i].massfrac[water][wPhase] = 1.0 - varNData[i].massfrac[air][wPhase];
         varNData[i].massfrac[water][nPhase] = problem.solu().Xwn(varNData[i].pN, varNData[i].temperature);
         varNData[i].massfrac[air][nPhase] = 1.0 - varNData[i].massfrac[water][nPhase];

         // CONSTANT solubility (for comparison with twophase)
//         varNData[i].massfrac[air][wPhase] = 0.1; varNData[i].massfrac[water][wPhase] = 0.9;
//         varNData[i].massfrac[water][nPhase] = 0.1; varNData[i].massfrac[air][nPhase] = 0.9;

         //std::cout << "water in gasphase: " << varNData[i].massfrac[water][nPhase] << std::endl;
         //std::cout << "air in waterphase: " << varNData[i].massfrac[air][wPhase] << std::endl;
   	 }   	 
   	 return;
    }
    

    
    struct StaticNodeData 
    {
   	 bool visited;
   	 int phaseState;//[numberOfComponents];
   	 RT cellVolume;
   	 RT porosity;
   	 FieldVector<RT, 4> parameters;
   	 FMatrix K;
    };
    
    struct VariableNodeData  
    {
   	 RT gravity;
   	 RT saturationN;
     RT saturationW;
     RT pW;
     RT pC;
     RT pN;
     RT temperature;
     VBlockType mobility;  //Vector with the number of phases
     VBlockType density;
     FieldMatrix<RT,c,m> massfrac;
    };

	 struct StaticIPData
	 {
		 bool visited;
		 FMatrix K;
	 };
    
    
    struct ElementData {
//   	 RT cellVolume;
//     RT porosity;
//   	 RT gravity;
//   	 FieldVector<RT, 4> parameters;
//   	 FieldMatrix<RT,dim,dim> K;
   	 } elData;
   	    	 
    
    // parameters given in constructor
    TwoPTwoCProblem<G,RT>& problem;
    Solubility solub;
    std::vector<StaticNodeData> statNData;
    std::vector<StaticIPData> statIPData;
    std::vector<VariableNodeData> varNData;
  };

//  virtual void updateStaticIPData (const Entity& e, const VBlockType* sol)
//  {
//	  // size is determined in the constructor
//	  
//	  	// get access to shape functions for P1 elements
//    	  GeometryType gt = e.geometry().type();
//    	  const typename Dune::LagrangeShapeFunctionSetContainer<DT,RT,dim>::value_type& 
//    	  sfs=Dune::LagrangeShapeFunctions<DT,RT,dim>::general(gt,1);
//
//    	  for (int k = 0; k < this->fvGeom.nEdges; ++k) // begin loop over edges / sub control volume faces
//  	  {
//  		  // get local to global id map
//    		  int globalIdx = this->vertexMapper.template map<dim>(e, sfs[k].entity());
//
//  		  // if edge is not already visited
//  		  //if (!statIPData[globalIdx].visited)
//    		  {
//    			  int i = this->fvGeom.subContVolFace[k].i;
//    			  int j = this->fvGeom.subContVolFace[k].j;
//
//    			  int globalIdx = this->vertexMapper.template map<dim>(e, sfs[j].entity());
//
//    			  // global coordinates
//    			  FieldVector<DT,dim> global_i = this->fvGeom.subContVol[i].global;
//    			  FieldVector<DT,dim> global_j = this->fvGeom.subContVol[j].global;
//    			  //int local_i = this->fvGeom.subContVol[i].local;
//    			  //int local_j = this->fvGeom.subContVol[j].local;
//    			     			 
//    			  //statIPData[globalIdx].K_eff[k] = harmonicMeanK(k, global_i,global_j); 
//    			  //statIPData[globalIdx].K = harmonicMeanK(global_i,global_j); 
//
//    			  // mark elements that were already visited
//     			  //statIPData[globalIdx].visited = true;
//    		  }
//  	  }
// 	 
// 	 //return;
//  }
  
  
}
#endif
