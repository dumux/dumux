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
#include "dumux/twophase/twophaseproblem.hh"

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
	enum {pWIdx = 0, satNIdx = 1, numberOfComponents = 2};	// Phase index
	enum {gaseous = 0, liquid = 1, solid = 2};					// Phase state
	enum {water = 0, air = 1};										// Component index					
	
  public:
    // define the number of phases (m) and components (c) of your system, this is used outside
    // to allocate the correct size of (dense) blocks with a FieldMatrix
    enum {n=G::dimension};
    enum {m=2, c=2}; 
    enum {SIZE=LagrangeShapeFunctionSetContainer<DT,RT,n>::maxsize};
    
    //! Constructor
    Box2P2CJacobian (TwoPhaseProblem<G,RT>& params,
			      bool levelBoundaryAsDirichlet_, const G& grid, 
			      BoxFunction& sol, 
			      bool procBoundaryAsDirichlet_=true)
    : BoxJacobian<ThisType,G,RT,2,BoxFunction>(levelBoundaryAsDirichlet_, grid, sol, procBoundaryAsDirichlet_), 
      problem(params), 
      statNData(this->vertexMapper.size())
    {
      this->analytic = false;
    }
    

    virtual void clearVisited ()
    {
    	for (int i = 0; i < this->vertexMapper.size(); i++)
    		statNData[i].visited = false;
    		
    	return;
    }

    // Compute time dependent terms
    // ACHTUNG varNData always contains values from the NEW timestep
    virtual VBlockType computeM (const Entity& e, const VBlockType* sol, int node)
    {
   	 VBlockType result; 
   	 
   	 RT satN = sol[node][satNIdx];
   	 RT satW = 1.0 - satN;  
   	 
   	                  
   	 // Water Phase
   	 result[0] = elData.porosity*(varNData[node].density[pWIdx]*satW*varNData[node].massfrac[water][pWIdx]
                   + varNData[node].density[satNIdx]*satN*varNData[node].massfrac[water][satNIdx]);
   	 // Gas Phase
   	 result[1] = elData.porosity*(varNData[node].density[satNIdx]*satN*varNData[node].massfrac[air][satNIdx]
   	             + varNData[node].density[pWIdx]*satW*varNData[node].massfrac[air][pWIdx]);   
   	 
   	 //std::cout << result << std::endl;
   	 return result;
    };
    
    // Compute the Fluxes
    virtual VBlockType computeA (const Entity& e, const VBlockType* sol, int face)
   {
   	 int i = this->fvGeom.subContVolFace[face].i;
     	 int j = this->fvGeom.subContVolFace[face].j;
   	 
		 // permeability in edge direction 
		 FieldVector<DT,n> Kij(0);
		 elData.K.umv(this->fvGeom.subContVolFace[face].normal, Kij);  // K*n
		 
		 VBlockType flux;
		 FieldMatrix<RT,m,n> pGrad(0); 
		 FieldVector<RT,n> temp(0);
		 // calculate FE gradient (grad p for each phase)
		 for (int k = 0; k < this->fvGeom.nNodes; k++) // loop over nodes
       {	 
      	 FieldVector<DT,n> grad(this->fvGeom.subContVolFace[face].grad[k]); // receives the FEGradient at node k
      	 FieldVector<DT,m> pressure(0);
      	 pressure[pWIdx] = sol[k][pWIdx];
      	 pressure[satNIdx] = varNData[k].pN;

      	 // compute sum of pressure gradients for each phase
      	 for (int d = 0; d < m; d++)
      	 {	      		 
      		 temp = grad;
      		 temp *= pressure[d];
      		 pGrad[d] += temp;
      	 }
       }

       // deduce gravity*density of each phase
		 FieldMatrix<RT,m,n> contribComp(0);
		 for (int d=0; d<m; d++)
		 {
			 contribComp[d] = problem.gravity();
			 contribComp[d] *= varNData[i].density[d];  
			 pGrad[d] -= contribComp[d]; // grad p -rho*g
		 }
	 	 
		 // calculate the flux using upwind
		 FieldVector<RT,m> outward(0);  // Darcy velocity of each phase
		 FieldVector<RT,4> massfraction;
		 massfraction[0] = varNData[i].massfrac[water][pWIdx];
		 massfraction[1] = varNData[i].massfrac[water][satNIdx];
		 massfraction[2] = varNData[i].massfrac[air][pWIdx];
		 massfraction[3] = varNData[i].massfrac[air][satNIdx];

		 for (int d=0; d<m; d++)
		 {
			 outward[d] = pGrad[d]*Kij;  //K*n(grad p -rho*g)  

			 if (outward[d] < 0)
			 {
		  			temp[d] = varNData[i].density[d]*varNData[i].mobility[d]*outward[d];
			 }
			 else
			 { 
		  			temp[d] = varNData[j].density[d]*varNData[j].mobility[d]*outward[d];
			 }
		 }
		 // water conservation
		 flux[water] = massfraction[pWIdx]*temp[pWIdx]+massfraction[satNIdx]*temp[satNIdx];
		 // air conservation
		 flux[air] = massfraction[pWIdx+c]*temp[pWIdx]+massfraction[satNIdx+c]*temp[satNIdx];

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
 		 elData.parameters = problem.materialLawParameters(this->fvGeom.cellGlobal, e, this->fvGeom.cellLocal);
   	 
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
  	  const typename Dune::LagrangeShapeFunctionSetContainer<DT,RT,n>::value_type& 
  	  sfs=Dune::LagrangeShapeFunctions<DT,RT,n>::general(gt,1);

  	  // get local to global id map
  	  for (int k = 0; k < sfs.size(); k++) {
  		  int globalId = this->vertexMapper.template map<n>(e, sfs[k].entity());
  		  
  		  // if nodes are not already visited
  		  if (!statNData[globalId].visited) 
  		  {
  			  // phase state
  			  if (1)
  			  {
  				statNData[globalId].phaseState[0] = gaseous;
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

   	 for (int i = 0; i < size; i++) {
   		this->def[i] = 0;
   		varNData[i].saturationW = 1.0 - sol[i][satNIdx];
         varNData[i].pC = problem.materialLaw().pC(varNData[i].saturationW, elData.parameters);
         varNData[i].pN = sol[i][pWIdx] + varNData[i].pC;
         // Mobilities & densities
         varNData[i].mobility[pWIdx] = problem.materialLaw().mobW(varNData[i].saturationW, elData.parameters);
         varNData[i].mobility[satNIdx] = problem.materialLaw().mobN(sol[i][satNIdx], elData.parameters);
         varNData[i].density[pWIdx] = problem.materialLaw().wettingPhase.density();
         varNData[i].density[satNIdx] = problem.materialLaw().nonwettingPhase.density();
         // Solubilities
         varNData[i].massfrac[water][pWIdx] = 1.0;	//problem.materialLaw().water.solubility();
         varNData[i].massfrac[air][pWIdx] = 	1.0 - varNData[i].massfrac[water][pWIdx];
         varNData[i].massfrac[water][satNIdx] = 0.0;	//problem.materialLaw().air.solubility();
         varNData[i].massfrac[air][satNIdx] = 1.0 - varNData[i].massfrac[water][satNIdx];
   	 }   	 
    }
    
    
    struct StaticNodeData 
    {
    	bool visited;
    	
    	int phaseState[numberOfComponents];
    };
    
       // the members of the struct are defined here
    struct VariableNodeData  
    {
       RT saturationW;
       RT pC;
       RT pN;
       VBlockType mobility;  //Vector with the number of phases
       VBlockType density;
       MBlockType massfrac;
    };
    
    struct ElementData {
   	 RT cellVolume;
    	 RT porosity;
   	 RT gravity;
   	 FieldVector<RT, 4> parameters;
   	 FieldMatrix<DT,n,n> K;
   	 } elData;
    
    // parameters given in constructor
    TwoPhaseProblem<G,RT>& problem;
    std::vector<StaticNodeData> statNData;
    std::vector<VariableNodeData> varNData;
  };

  /** @} */
}
#endif
