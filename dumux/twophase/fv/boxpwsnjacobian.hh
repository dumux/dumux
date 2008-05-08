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
  class BoxPwSnJacobian 
    : public BoxJacobian<BoxPwSnJacobian<G,RT,BoxFunction>,G,RT,2,BoxFunction>
  {
    typedef typename G::ctype DT;
    typedef typename G::Traits::template Codim<0>::Entity Entity;
    typedef typename Entity::Geometry Geometry;
    typedef BoxPwSnJacobian<G,RT,BoxFunction> ThisType;
    typedef typename LocalJacobian<ThisType,G,RT,2>::VBlockType VBlockType;
    typedef typename LocalJacobian<ThisType,G,RT,2>::MBlockType MBlockType;
 	typedef FVElementGeometry<G> FVElementGeometry;
	enum {pWIdx = 0, satNIdx = 1};
	
	
  public:
    // define the number of components of your system, this is used outside
    // to allocate the correct size of (dense) blocks with a FieldMatrix
    enum {n=G::dimension};
    enum {m=2};
    enum {SIZE=LagrangeShapeFunctionSetContainer<DT,RT,n>::maxsize};
    
    //! Constructor
    BoxPwSnJacobian (TwoPhaseProblem<G,RT>& params,
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
    	return;
    }

    virtual VBlockType computeM (const Entity& e, const VBlockType* sol, int node)
    {
   	 VBlockType result; 
   	 
   	 result[0] = -varNData[node].density[pWIdx]*elData.porosity*sol[node][satNIdx];
   	 result[1] = varNData[node].density[satNIdx]*elData.porosity*sol[node][satNIdx];
   	 
   	 return result;
    };
    
    virtual VBlockType computeQ (const Entity& e, const VBlockType* sol, const int& node)
    {
   	 // ASSUME problem.q already contains \rho.q
   	 return problem.q(this->fvGeom.subContVol[node].global, e, this->fvGeom.subContVol[node].local);
    };
    
    virtual VBlockType computeA (const Entity& e, const VBlockType* sol, int face)
    {
   	 int i = this->fvGeom.subContVolFace[face].i;
     	 int j = this->fvGeom.subContVolFace[face].j;
     	 
		  // permeability in edge direction 
		  FieldVector<DT,n> Kij(0);
		  elData.K.umv(this->fvGeom.subContVolFace[face].normal, Kij);
		  
		  VBlockType flux;
		  for (int comp = 0; comp < m; comp++) {
	          // calculate FE gradient
	          FieldVector<RT, n> pGrad(0);
	          for (int k = 0; k < this->fvGeom.nNodes; k++) {
	        	  FieldVector<DT,n> grad(this->fvGeom.subContVolFace[face].grad[k]);
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

    virtual void computeElementData (const Entity& e)
    {
  		 // ASSUME element-wise constant parameters for the material law 
 		 elData.parameters = problem.materialLawParameters(this->fvGeom.cellGlobal, e, this->fvGeom.cellLocal);
   	 
		 // ASSUMING element-wise constant permeability, evaluate K at the cell center 
 		 elData.K = problem.K(this->fvGeom.cellGlobal, e, this->fvGeom.cellLocal);  

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

   		// ASSUME element-wise constant parameters for the material law 
         FieldVector<RT, 4> parameters = problem.materialLawParameters(this->fvGeom.cellGlobal, e, this->fvGeom.cellLocal);

         varNData[i].pC = problem.materialLaw().pC(varNData[i].saturationW, parameters);
         varNData[i].pN = sol[i][pWIdx] + varNData[i].pC;
         varNData[i].mobility[pWIdx] = problem.materialLaw().mobW(varNData[i].saturationW, parameters);
         varNData[i].mobility[satNIdx] = problem.materialLaw().mobN(sol[i][satNIdx], parameters);
         varNData[i].density[pWIdx] = problem.materialLaw().wettingPhase.density();
         varNData[i].density[satNIdx] = problem.materialLaw().nonwettingPhase.density();
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
       VBlockType mobility;  //Vector with the number of phases
       VBlockType density;
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
