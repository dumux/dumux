#ifndef DUNE_BOXPNSWJACOBIAN_HH
#define DUNE_BOXPNSWJACOBIAN_HH

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
#include"dumux/operators/boxjacobian.hh"
#include"dumux/twophase/twophaseproblem.hh"


namespace Dune
{

  template<class G, class RT, class BoxFunction = LeafP1Function<G, RT, 2> > class BoxPnSwJacobian 
    : public BoxJacobian<BoxPnSwJacobian<G,RT,BoxFunction>,G,RT,2,BoxFunction>
  {
    typedef typename G::ctype DT;
    typedef typename G::Traits::template Codim<0>::Entity Entity;
    typedef typename Entity::Geometry Geometry;
    typedef BoxPnSwJacobian<G,RT,BoxFunction> ThisType;
    typedef typename LocalJacobian<ThisType,G,RT,2>::VBlockType VBlockType;
    typedef typename LocalJacobian<ThisType,G,RT,2>::MBlockType MBlockType;
 	typedef FVElementGeometry<G> FVElementGeometry;
	enum {pNIdx = 0, satWIdx = 1};
	
	
  public:
    // define the number of components of your system, this is used outside
    // to allocate the correct size of (dense) blocks with a FieldMatrix
    enum {n=G::dimension};
    enum {m=2};
    enum {SIZE=LagrangeShapeFunctionSetContainer<DT,RT,n>::maxsize};
    struct VariableNodeData;
    
    //! Constructor
    BoxPnSwJacobian (TwoPhaseProblem<G,RT>& params,
			      bool levelBoundaryAsDirichlet_, const G& grid, 
			      BoxFunction& sol, 
			      bool procBoundaryAsDirichlet_=true)
    : BoxJacobian<ThisType,G,RT,2,BoxFunction>(levelBoundaryAsDirichlet_, grid, sol, procBoundaryAsDirichlet_), 
      problem(params), 
      statNData(this->vertexMapper.size()), varNData(SIZE), oldVarNData(SIZE)
    {
      this->analytic = false;
    }
    

    virtual void clearVisited ()
    {
    	return;
    }

    virtual VBlockType computeM (const Entity& e, const VBlockType* sol, 
    		int node, const std::vector<VariableNodeData>& varData)
    {
   	 VBlockType result; 
   	 
   	 result[0] = varData[node].density[pNIdx]*elData.porosity*sol[node][satWIdx];
   	 result[1] = -varData[node].density[satWIdx]*elData.porosity*sol[node][satWIdx];
   	 
   	 return result;
    };
    
    virtual VBlockType computeM (const Entity& e, const VBlockType* sol, int node, bool old = false) 
    {
    	if (old)
    		return computeM(e, sol, node, oldVarNData);
    	else 
    		return computeM(e, sol, node, varNData);
    }

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
	        	  grad *= (comp) ? sol[k][pNIdx] : varNData[k].pW;
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
   

    // the members of the struct are defined here
 struct VariableNodeData  
 {
    RT saturationN;
    RT pC;
    RT pW;
    VBlockType mobility;  //Vector with the number of phases
    VBlockType density;
 };
 
    // analog to EvalPrimaryData in MUFTE, uses members of varNData
	virtual void updateVariableData(const Entity& e, const VBlockType* sol, 
			int i, std::vector<VariableNodeData>& varData) 
    {
   		varData[i].saturationN = 1.0 - sol[i][satWIdx];

   		// ASSUME element-wise constant parameters for the material law 
         FieldVector<RT, 4> parameters = problem.materialLawParameters(this->fvGeom.cellGlobal, e, this->fvGeom.cellLocal);

         varData[i].pC = problem.materialLaw().pC(sol[i][satWIdx], parameters);
         varData[i].pW = sol[i][pNIdx] - varData[i].pC;
         varData[i].mobility[pNIdx] = problem.materialLaw().mobW(sol[i][satWIdx], parameters);
         varData[i].mobility[satWIdx] = problem.materialLaw().mobN(varData[i].saturationN, parameters);
         varData[i].density[pNIdx] = problem.materialLaw().wettingPhase.density();
         varData[i].density[satWIdx] = problem.materialLaw().nonwettingPhase.density();
    }
    
	virtual void updateVariableData(const Entity& e, const VBlockType* sol, int i, bool old = false) 
	{
		if (old)
			updateVariableData(e, sol, i, oldVarNData);
		else 
			updateVariableData(e, sol, i, varNData);
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
    std::vector<VariableNodeData> oldVarNData;
    
  };

  /** @} */
}
#endif
