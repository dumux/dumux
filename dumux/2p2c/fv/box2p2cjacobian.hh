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
#include <dune/grid/utility/intersectiongetter.hh>

#include<dune/disc/shapefunctions/lagrangeshapefunctions.hh>
#include<dune/disc/operators/boundaryconditions.hh>
#include<dune/disc/functions/p1function.hh>
#include "dumux/operators/boxjacobian.hh"
#include "dumux/2p2c/2p2cproblem.hh"

//#include "dumux/2p2c/fv/varswitch.hh"

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
    
    //! Constructor
    Box2P2CJacobian (TwoPTwoCProblem<G,RT>& params,
			      bool levelBoundaryAsDirichlet_, const G& grid, 
			      BoxFunction& sol, 
			      bool procBoundaryAsDirichlet_=true)
    : BoxJacobian<ThisType,G,RT,2,BoxFunction>(levelBoundaryAsDirichlet_, grid, sol, procBoundaryAsDirichlet_), 
      problem(params), 
      statNData(this->vertexMapper.size())
    {
      this->analytic = false;
    }


    virtual void primaryVarSwitch (const VBlockType* sol, const Entity& e, int i, 
   		 FieldMatrix<RT,m,c>& massfrac, int state)
    {
        RT Sw, Xaw, Xwg, Xamax;

        Sw      = sol[i][satNIdx];
        Xamax = 0.5; // Needs to be changed !!!

//        state   = statNData[globalId].phaseState[0];
        Xaw   	 = massfrac[air][wPhase];
        Xwg     = massfrac[water][nPhase];;
     
        //std::cout << "negative non-wetting saturation!! --> Switsching Variables" << std::endl;
       
        
        switch(state) 
        {
        case 0 : /* GasPhase */

            if (sol[i][satNIdx] > 0.001)
            {
            // appearance of water phase
//          	std::cout << "Water appears" << std::endl;
//             state[i] = BothPhases;
//        		dat->co_Sw[i] = 1.E-6;
            }
            break;

        case 1 : /* WaterPhase */

            if (sol[i][satNIdx] > Xamax)
            {
            // appearance of gas phase
//          	std::cout << "Gas appears" \n x = %10.4f \n y = %10.4f \n XCO2max = %12.10f\n \n XCO2 = %12.10f\n", level, x, y, XCO2max, *var_wc);
//             state[i] = BothPhases;
//             dat->co_Sw[i] = 1.0 - 1.E-6;  /* Initialisierung */
            }
            break;

        case 2 : /* BothPhases */

      	  if ( sol[i][satNIdx] < 0.0 ) 
      	  {
      		  // disappearance of gas phase
//    			std::cout << "w,g (Level %d): gas disappears \n x = %10.4f \n y = %10.4f \n ";
//    			state[i] = WaterPhase;
//    			*var_wc = 1.E-6;  /* Initialisierung */
    		}
    		if ( 1 - sol[i][satNIdx] < 0.0 ) 
    		{
//    			/* disappearance of water phase */
//    			std::cout << "Water disappears" \n x = %10.4f \n y = %10.4f \n ", level, x, y);
//    			state[i] = GasPhase;
//    			*var_wc = 1.E-6;  /* Initialisierung */
    		}
    	break;
      }

    	/* update state */
//    	state = ...;
//
//
//       if (state[i] == gasPhase)  
//       {
//          *var_wc = dat->co_Xwg[i];
//       }
//    	if (state[i] == waterPhase)  
//    	{
//    		*var_wc = dat->co_Xaw[i];
//    	}
//    	if (state[i] == bothPhases)  
//    	{
//    		*var_wc = dat->co_Sw[i];
//    	}

        return;
    }

    
    virtual void clearVisited ()
    {
   	 for (int i = 0; i < this->vertexMapper.size(); i++)
   		 statNData[i].visited = false;
   	 
   	 return;
    }

    // Compute time dependent terms (storage)
    // ACHTUNG varNData always contains values from the NEW timestep
    virtual VBlockType computeM (const Entity& e, const VBlockType* sol, int node)
    {
   	 VBlockType result; 
   	 RT satN = sol[node][satNIdx];
   	 RT satW = 1.0 - satN;  
   	    	                  
   	 // wetting phase
   	 result[0] = 
   		 elData.porosity*(varNData[node].density[wPhase]*satW*varNData[node].massfrac[water][wPhase]
   		                 +varNData[node].density[nPhase]*satN*varNData[node].massfrac[water][nPhase]);
   	 // non-wetting phase
   	 result[1] = 
   		 elData.porosity*(varNData[node].density[nPhase]*satN*varNData[node].massfrac[air][nPhase]
   	                    +varNData[node].density[wPhase]*satW*varNData[node].massfrac[air][wPhase]);   
   	 
   	 //std::cout << result << " " << node << std::endl;
   	 return result;
    };
    
    // Compute advective and diffusive fluxes
    virtual VBlockType computeA (const Entity& e, const VBlockType* sol, int face)
    {
   	 int i = this->fvGeom.subContVolFace[face].i;
     	 int j = this->fvGeom.subContVolFace[face].j;
     	 FieldVector<RT,dim> normal(this->fvGeom.subContVolFace[face].normal);

     	 VBlockType flux;
		 FieldMatrix<RT,m,dim> pGrad(0); 
		 FieldVector<RT,dim> xGrad(0), temp(0); 
     	 
		 // permeability in edge direction 
     	 FieldVector<RT,dim> Kij(0);
		 elData.K.umv(normal, Kij);  // Kij=K*n
		 
		 // calculate FE gradient (grad p for each phase)
		 for (int k = 0; k < this->fvGeom.nNodes; k++) // loop over adjacent nodes
       {	 
      	 FieldVector<DT,dim> feGrad(this->fvGeom.subContVolFace[face].grad[k]); // FEGradient at node k
       	 FieldVector<RT,m> pressure(0), massfrac(0);

      	 pressure[wPhase] = sol[k][pWIdx];
      	 pressure[nPhase] = varNData[k].pN;
      	 
      	 massfrac[wPhase] = varNData[k].massfrac[water][nPhase]; // water in gas phase
      	 massfrac[nPhase] = varNData[k].massfrac[air][wPhase];	// air in water phase
      	 
      	 // compute sum of pressure gradients for each phase
      	 for (int phase = 0; phase < m; phase++)
      	 {	      		 
      		 temp = feGrad;
      		 temp *= pressure[phase];
      		 pGrad[phase] += temp;
      	 }

      	 // compute sum of concentration gradient
      	 // only water diffusion in the gas phase considered !!!
      	 temp = feGrad;
      	 temp *= massfrac[water];
      	 xGrad += temp;
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
		 FieldVector<RT,4> massfraction;
		 massfraction[0] = varNData[i].massfrac[water][wPhase];
		 massfraction[1] = varNData[i].massfrac[water][nPhase];
		 massfraction[2] = varNData[i].massfrac[air][wPhase];
		 massfraction[3] = varNData[i].massfrac[air][nPhase];

		 // calculate the advective flux using upwind
		 for (int phase=0; phase<m; phase++)
		 {
			 outward[phase] = pGrad[phase]*Kij;  //K*n(grad p -rho*g)  

			 if (outward[phase] <= 0)
				 temp[phase] = varNData[i].density[phase]*varNData[i].mobility[phase]*outward[phase];
			 else
				 temp[phase] = varNData[j].density[phase]*varNData[j].mobility[phase]*outward[phase];
		 }
		 
		 // DIFFUSION
		 RT normDiffGrad, diffusionWG;
		 VBlockType avgDensity, avgDpm(0);
		 avgDpm[wPhase]=1e-9; // needs to be changed !!!
		 
		 normDiffGrad = xGrad*normal;
		 avgDensity[wPhase] = 0.5*(varNData[i].density[wPhase] + varNData[j].density[wPhase]);
		 
		 diffusionWG = avgDpm[wPhase] * avgDensity[wPhase] * normDiffGrad;
//		 std::cout << "Diffusive Flux: " << diffusionWG << std::endl; 
		 
		 
		 // water conservation
		 flux[water] = /*diffusionWG +*/ massfraction[wPhase]*temp[wPhase]+massfraction[nPhase]*temp[nPhase];
		 // air conservation
		 flux[air] = massfraction[wPhase+c]*temp[wPhase]+massfraction[nPhase+c]*temp[nPhase];


		 return flux;
  };

   // Integrate sources / sinks
   virtual VBlockType computeQ (const Entity& e, const VBlockType* sol, const int& node)
   {
  	 // ASSUME problem.q already contains \rho.q
  	 return problem.q(this->fvGeom.subContVol[node].global, e, this->fvGeom.subContVol[node].local);
   };
       
    
	  //*********************************************************
	  //*																			*
	  //*	Calculation of Data at Elements (elData) 					*
	  //*						 													*
	  //*																		 	*
	  //*********************************************************

    virtual void computeElementData (const Entity& e)
    {
  		 // ASSUME element-wise constant parameters for the material law 
 		 elData.parameters = problem.materialLawParameters
 		 (this->fvGeom.cellGlobal, e, this->fvGeom.cellLocal);
   	 
		 // ASSUMING element-wise constant permeability, evaluate K at the cell center 
 		 elData.K = problem.K(this->fvGeom.cellGlobal, e, this->fvGeom.cellLocal);  

		 // ASSUMING element-wise constant porosity 
 		 elData.porosity = problem.porosity(this->fvGeom.cellGlobal, e, this->fvGeom.cellLocal);
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
  	  // local to global id mapping (do not ask vertex mapper repeatedly
  	  //int localToGlobal[Dune::LagrangeShapeFunctionSetContainer<DT,RT,n>::maxsize];

  	  // get access to shape functions for P1 elements
  	  GeometryType gt = e.geometry().type();
  	  const typename Dune::LagrangeShapeFunctionSetContainer<DT,RT,dim>::value_type& 
  	  sfs=Dune::LagrangeShapeFunctions<DT,RT,dim>::general(gt,1);

  	  // get local to global id map
  	  for (int k = 0; k < sfs.size(); k++) {
  		  int globalId = this->vertexMapper.template map<dim>(e, sfs[k].entity());
  		  
  		  // if nodes are not already visited
  		  if (!statNData[globalId].visited) 
  		  {
  			  // phase state
  			  if (1)
  			  {
  				statNData[globalId].phaseState[0] = gasPhase;
  			  }
  			  else 
  			  {
  				  
  			  }
  			  
  			  // absolute permeability
  			  //statNData[globalId].absolutePermeability = ;  // siehe oben!!

  			  
  			  // mark elements that were already visited
  			  statNData[globalId].visited = true;
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
   	 int state = bothPhases;

   	 for (int i = 0; i < size; i++) {
   		this->def[i] = 0;

   		if (sol[i][satNIdx] < 0){
   			primaryVarSwitch(sol, e, i, varNData[i].massfrac, state);
   		}

   		varNData[i].saturationW = 1.0 - sol[i][satNIdx];
   		if (varNData[i].saturationW < 0){
   			std::cout << "negative wetting saturation!!"; 
   		}
   		
         varNData[i].pC = problem.materialLaw().pC(varNData[i].saturationW, elData.parameters);
         varNData[i].pN = sol[i][pWIdx] + varNData[i].pC;
         varNData[i].temperature = 283.15; // in [K]
         // Mobilities & densities
         varNData[i].mobility[wPhase] = problem.materialLaw().mobW(varNData[i].saturationW, elData.parameters);
         varNData[i].mobility[nPhase] = problem.materialLaw().mobN(sol[i][nPhase], elData.parameters);
         varNData[i].density[wPhase] = problem.materialLaw().wettingPhase.density();
         varNData[i].density[nPhase] = problem.materialLaw().nonwettingPhase.density();
         // Solubilities of components in phases
         varNData[i].massfrac[air][wPhase] = problem.constrel().Xaw(varNData[i].pN, varNData[i].temperature);
         varNData[i].massfrac[water][wPhase] = 1.0 - varNData[i].massfrac[air][wPhase];
         varNData[i].massfrac[water][nPhase] = problem.constrel().Xwg(varNData[i].pN, varNData[i].temperature);
         varNData[i].massfrac[air][nPhase] = 1.0 - varNData[i].massfrac[water][nPhase];
         std::cout << "water in gasphase: " << varNData[i].massfrac[water][nPhase] << std::endl;
         std::cout << "air in waterphase: " << varNData[i].massfrac[air][wPhase] << std::endl;
   	 }   	 
    }
    
    
    // the members of the structs are defined here
    struct StaticNodeData 
    {
    	bool visited;
    	
    	int phaseState[numberOfComponents];
    };
    
    struct VariableNodeData  
    {
       RT saturationW;
       RT pC;
       RT pN;
       RT temperature;
       VBlockType mobility;  //Vector with the number of phases
       VBlockType density;
       FieldMatrix<RT,dim,dim> massfrac;
    };
    
    struct ElementData {
   	 RT cellVolume;
    	 RT porosity;
   	 RT gravity;
   	 FieldVector<RT, 4> parameters;
   	 FieldMatrix<RT,dim,dim> K;
   	 } elData;
    
    // parameters given in constructor
    TwoPTwoCProblem<G,RT>& problem;
    std::vector<StaticNodeData> statNData;
    std::vector<VariableNodeData> varNData;
  };

  /** @} */
}
#endif
