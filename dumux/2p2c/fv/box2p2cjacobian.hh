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
    // define the number of components of your system, this is used outside
    // to allocate the correct size of (dense) blocks with a FieldMatrix
    enum {n=G::dimension};
    enum {m=2};
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

    virtual VBlockType computeM (const Entity& e, const FVElementGeometry& fvGeom, const VBlockType* sol, int node)
    {
   	 VBlockType result; 
   	 
   	 
//   	 result[0] = -elData.porosity*(varNData[node].density[pWIdx]*sol[node][satNIdx]*varNData[node].massfrac[water][pWIdx]
//                   + varNData[node].density[satNIdx]*sol[node][satNIdx]*varNData[node].massfrac[water][satNIdx]);
   	 result[0] = -varNData[node].density[pWIdx]*elData.porosity*sol[node][satNIdx];
   	 result[1] = varNData[node].density[satNIdx]*elData.porosity*sol[node][satNIdx];
   	 
   	 return result;
    };
    
    virtual VBlockType computeQ (const Entity& e, const FVElementGeometry& fvGeom, const VBlockType* sol, const int& node)
    {
   	 // ASSUME problem.q already contains \rho.q
   	 return problem.q(fvGeom.subContVol[node].global, e, fvGeom.subContVol[node].local);
    };
    
    virtual VBlockType computeA (const Entity& e, const FVElementGeometry& fvGeom, const VBlockType* sol, int face)
    {
   	 int i = fvGeom.subContVolFace[face].i;
     	 int j = fvGeom.subContVolFace[face].j;
     	 
		  // permeability in edge direction 
		  FieldVector<DT,n> Kij(0);
		  elData.K.umv(fvGeom.subContVolFace[face].normal, Kij);
		  
		  VBlockType flux;
		  for (int comp = 0; comp < m; comp++) {
	          // calculate FE gradient
	          FieldVector<RT, n> pGrad(0);
	          for (int k = 0; k < fvGeom.nNodes; k++) {
	        	  FieldVector<DT,n> grad(fvGeom.subContVolFace[face].grad[k]);
	        	  grad *= (comp) ? varNData[k].pN : sol[k][pWIdx];
	        	  pGrad += grad;
	          }

	          // adjust by gravity 
			  FieldVector<RT, n> gravity = problem.gravity();
			  gravity *= varNData[i].density[comp];
			  pGrad -= gravity;
			  
			  // calculate the flux using upwind
			  RT outward = pGrad*Kij;
			  if (outward < 0)
				  flux[comp] = varNData[i].density[comp]*varNData[i].mobility[comp]*outward;
			  else 
				  flux[comp] = varNData[j].density[comp]*varNData[j].mobility[comp]*outward;
		  }
		  
		  return flux;
   };

    
    
	  //*********************************************************
	  //*																			*
	  //*	Calculation of Data at Elements			 					*
	  //*						 													*
	  //*																		 	*
	  //*********************************************************

    virtual void computeElementData (const Entity& e, const FVElementGeometry& fvGeom)
    {
  		 // ASSUME element-wise constant parameters for the material law 
 		 elData.parameters = problem.materialLawParameters(fvGeom.cellGlobal, e, fvGeom.cellLocal);
   	 
		 // ASSUMING element-wise constant permeability, evaluate K at the cell center 
 		 elData.K = problem.K(fvGeom.cellGlobal, e, fvGeom.cellLocal);  

 		 elData.porosity = problem.porosity(fvGeom.cellGlobal, e, fvGeom.cellLocal);
    };

    
	  //*********************************************************
	  //*																			*
	  //*	Calculation of Data at Nodes that has to be			 	*
	  //*	determined only once	(statNData)							 	*
	  //*																		 	*
	  //*********************************************************

    // analog to EvalStaticData in MUFTE
    virtual void updateStaticData (const Entity& e, const FVElementGeometry& fvGeom, const VBlockType* sol)
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
    virtual void updateVariableData (const Entity& e, const FVElementGeometry& fvGeom, const VBlockType* sol)
    {
   	 varNData.resize(fvGeom.nNodes);
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
         varNData[i].massfrac[water][pWIdx] = 1;	//problem.materialLaw().water.solubility();
         varNData[i].massfrac[air][pWIdx] = 	1 - varNData[i].massfrac[water][pWIdx];
         varNData[i].massfrac[satNIdx][water] = 1;	//problem.materialLaw().air.solubility();
         varNData[i].massfrac[satNIdx][air] = 1 - varNData[i].massfrac[satNIdx][water];
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
