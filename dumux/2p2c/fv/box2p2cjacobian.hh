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
#include "dumux/operators/localjacobian.hh"
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
  class Box2P2CLocalJacobian 
    : public LocalJacobian<Box2P2CLocalJacobian<G,RT>,G,RT,2>
  {
	  // mapper: one data element per vertex
	  template<int dim>
	  struct P1Layout
	  {
		  bool contains (Dune::GeometryType gt)
		  {
			  return gt.dim() == 0;
		  }
	  }; 

    typedef typename G::ctype DT;
    typedef typename G::Traits::template Codim<0>::Entity Entity;
    typedef typename Entity::Geometry Geometry;
    typedef typename LocalJacobian<Box2P2CLocalJacobian<G,RT>,G,RT,2>::VBlockType VBlockType;
    typedef typename LocalJacobian<Box2P2CLocalJacobian<G,RT>,G,RT,2>::MBlockType MBlockType;
 	typedef FVElementGeometry<G> FVElementGeometry;
	typedef MultipleCodimMultipleGeomTypeMapper<G, typename G::Traits::LeafIndexSet, P1Layout> VertexMapper;
	enum {pWIdx = 0, satNIdx = 1, numberOfComponents = 2};
	enum {gaseous = 0, liquid = 1, solid = 2};
	
	
  public:
    // define the number of components of your system, this is used outside
    // to allocate the correct size of (dense) blocks with a FieldMatrix
    enum {n=G::dimension};
    enum {m=2};
    enum {SIZE=LagrangeShapeFunctionSetContainer<DT,RT,n>::maxsize};
    
    //! Constructor
    Box2P2CLocalJacobian (TwoPhaseProblem<G,RT>& params,
			      bool levelBoundaryAsDirichlet_, const G& grid, 
			      BoxFunction& sol, 
			      bool procBoundaryAsDirichlet_=true)
    : problem(params),levelBoundaryAsDirichlet(levelBoundaryAsDirichlet_), 
      procBoundaryAsDirichlet(procBoundaryAsDirichlet_), 
      currentSolution(sol), oldSolution(grid), 
      vertexMapper(grid, grid.leafIndexSet()), statNData(vertexMapper.size()), dt(1)
    {
      this->analytic = false;
    }
    

    //**********************************************************
    //*																			*
    //*	Computation of the local defect								*
    //*																			*
    //**********************************************************
    
    template<class TypeTag>
    void localDefect (const Entity& e, const FVElementGeometry& fvGeom, const VBlockType* sol)
    {
      // extract some important parameters
      const Geometry& geometry = e.geometry();
      Dune::GeometryType gt = geometry.type();
      const typename Dune::LagrangeShapeFunctionSetContainer<DT,RT,n>::value_type& 
      	sfs=Dune::LagrangeShapeFunctions<DT,RT,n>::general(gt, 1);
      int size = sfs.size();
      this->setcurrentsize(size);
      
      computeElementData(e, fvGeom);
      varNData.resize(size);
      updateVariableData(e, fvGeom, sol);

      computeA();	// DUMMY
      
	  
	  // ASSUMING constant element mapping jacobian
	  FieldMatrix<DT,n,n> jacobianInverseTransposed = geometry.jacobianInverseTransposed(fvGeom.cellLocal);

	  for (int i=0; i < size; i++) // begin loop over vertices / sub control volumes
	  {
		  // implicit Euler
		  VBlockType massContrib = computeM(e, fvGeom, sol[i]);
		  VBlockType dummy = massContrib;
		  massContrib -= computeM(e, fvGeom, uold[i]);
		  massContrib *= fvGeom.subContVol[i].volume/dt;
		  this->def[i] += massContrib;
		  
		  // get source term 
		  VBlockType q = computeQ(e, fvGeom, i);
		  q *= fvGeom.subContVol[i].volume;
		  this->def[i] -= q;
	  }
	  
	  for (int face = 0; face < e.template count<n-1>(); face++) // begin loop over edges / sub control volume faces
//	  for (int i=0; i < size; i++) {
//		  for (int j=i+1; j < size; j++) // begin loop over neighboring vertices 
	  {
		  int i = fvGeom.subContVolFace[face].i;
		  int j = fvGeom.subContVolFace[face].j;
		  
			  // local coordinate of neighbor 
			  const FieldVector<DT,n> neighborLocal = sfs[j].position();

			  if (!gt.isSimplex()) {
				  // compute the local distance 
				  DT distanceLocal = (fvGeom.subContVol[i].local - neighborLocal).two_norm();
			  
				  // check whether the two vertices share a cell edge
				  if (distanceLocal > 1.01)
					  continue;
			  }
			  
			  // get global coordinate of neighbor
			  const FieldVector<DT,n> neighborGlobal = geometry.global(neighborLocal);
			  
			  // compute the edge vector
			  FieldVector<DT,n>  edgeVector = neighborGlobal - fvGeom.subContVol[i].global;
			  
			  // get distance between neighbors 
			  DT oneByDistanceGlobal = 1.0/edgeVector.two_norm(); 
			  
			  // normalize edge vector 
			  edgeVector *= oneByDistanceGlobal; 
			  
			  // get the local edge center 
			  FieldVector<DT,n> edgeLocal = fvGeom.subContVol[i].local + neighborLocal;
			  edgeLocal *= 0.5;
			  
			  // permeability in edge direction 
			  FieldVector<DT,n> Kij(0);
			  elData.K.umv(edgeVector, Kij);
			  
			  // calculate pressure differences 
//			  VBlockType pDiff;
//			  pDiff[pWIdx] = sol[j][pWIdx] - sol[i][pWIdx];
//			  pDiff[satNIdx] = pDiff[pWIdx] + varNData[j].pC - varNData[i].pC;
			  
			  VBlockType flux;
			  for (int comp = 0; comp < m; comp++) {
//				  // calculate pressure component gradient
//				  FieldVector<RT, n> pGrad(edgeVector);
//				  pGrad *= oneByDistanceGlobal*pDiff[comp];
				  
		          // calculate FE gradient
		          FieldVector<RT, n> pGrad(0);
		          for (int k = 0; k < size; k++) {
		        	  FieldVector<DT,n> grad(0),temp;
		        	  for (int l = 0; l < n; l++) 
		        		  temp[l] = sfs[k].evaluateDerivative(0, l, edgeLocal);
		        	  jacobianInverseTransposed.umv(temp, grad);
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
					  flux[comp] = varNData[i].mobility[comp]*outward;
				  else 
					  flux[comp] = varNData[j].mobility[comp]*outward;
			  }
			  
			  // get global coordinate of edge center
			  const FieldVector<DT,n> edgeGlobal = geometry.global(edgeLocal);
			  
			  // distance between cell center and edge center
			  DT distanceEdgeCell = (fvGeom.cellGlobal - edgeGlobal).two_norm();
			  
			  ////////////////////////////////////////////////////////////
			  // CAREFUL: only valid in 2D 
			  ////////////////////////////////////////////////////////////
			  // obtain integrated Flux 
			  flux *= distanceEdgeCell; 
			  
			  // add to defect 
			  //flux *= dt;
			  this->def[i] -= flux;
			  this->def[j] += flux;
//		  } // end loop over neighboring vertices 
	  } // end loop over vertices

	  // adjust by density 
	  for (int i=0; i < size; i++) {
		  this->def[i][pWIdx] *= varNData[i].density[pWIdx];
		  this->def[i][satNIdx] *= varNData[i].density[satNIdx];
	  }
	  
	  // assemble boundary conditions 
	  assembleBC<TypeTag> (e); 
	  
	  // add to defect 
	  for (int i=0; i < size; i++) {
//	    this->b[i] *= dt;
		  this->def[i] -= this->b[i];
	  }

//	  for (int i=0; i < 2*size; i++) {
//		  std::cout << "M = " << mass[i] << ", A = " << stiffness[i] << ", Q = " << neumann[i] << std::endl;
//	  }
//	  std::cout << std::endl;
		
      return;
    }
    
    void setLocalSolution (const Entity& e)
    {
        Dune::GeometryType gt = e.geometry().type();
        const typename Dune::LagrangeShapeFunctionSetContainer<DT,RT,n>::value_type& 
        	sfs=Dune::LagrangeShapeFunctions<DT,RT,n>::general(gt, 1);
        int size = sfs.size();
        this->setcurrentsize(size);
        
        for (int i = 0; i < size; i++) 
        	for (int comp = 0; comp < m; comp++) {
        		this->u[i][comp]= currentSolution.evallocal(comp, e, sfs[i].position());
        		uold[i][comp]= oldSolution.evallocal(comp, e, sfs[i].position());
        }
        
    	return;
    }
    
    void setDt (double d) 
    {
    	dt = d;
    }
  
    void setOldSolution (BoxFunction& uOld) 
    {
    	*oldSolution = *uOld;
    }
  
    virtual void clearVisited ()
    {
    	for (int i = 0; i < vertexMapper.size(); i++)
    		statNData[i].visited = false;
    		
    	return;
    }

    virtual VBlockType computeM (const Entity& e, const FVElementGeometry& fvGeom, const VBlockType& nodeSol)
    {
   	 VBlockType result; 
   	 
   	 result[0] = -elData.porosity*nodeSol[satNIdx];
   	 result[1] = elData.porosity*nodeSol[satNIdx];
   	 
   	 return result;
    };
    
    virtual VBlockType computeQ (const Entity& e, const FVElementGeometry& fvGeom, const int& nodeIdx)
    {
		 return problem.q(fvGeom.subContVol[nodeIdx].global, e, fvGeom.subContVol[nodeIdx].local);
    };
    
    virtual void computeA ()
    {
   	 
    };

    
    
	  //*********************************************************
	  //*																			*
	  //*	Calculation of Data at Elements			 					*
	  //*						 													*
	  //*																		 	*
	  //*********************************************************

    
    struct ElementData {
   	 RT cellVolume;
    	 RT porosity;
   	 RT gravity;
   	 double volumeFactor;
   	 FieldVector<RT, 4> parameters;
   	 FieldMatrix<DT,n,n> K;
   	 } elData;
    
    virtual void computeElementData (const Entity& e, const FVElementGeometry& fvGeom)
    {
  		 // ASSUME element-wise constant parameters for the material law 
 		 elData.parameters = problem.materialLawParameters(fvGeom.cellGlobal, e, fvGeom.cellLocal);
   	 

 		 // ASSUMING element-wise constant permeability, evaluate K at the cell center 
 		 elData.K = problem.K(fvGeom.cellGlobal, e, fvGeom.cellLocal);  

 		 elData.porosity = problem.porosity(fvGeom.cellGlobal, e, fvGeom.cellLocal);
 		 
 		 // ASSUMING element-wise constant porosity, evaluate at the cell center 
 		 elData.volumeFactor = elData.porosity*fvGeom.cellVolume;
    };

    
	  //*********************************************************
	  //*																			*
	  //*	Calculation of Data at Nodes that has to be			 	*
	  //*	determined only once	(statNData)							 	*
	  //*																		 	*
	  //*********************************************************

    
    struct StaticNodeData 
    {
    	bool visited;
    	
    	int phaseState[numberOfComponents];
    };
    
    
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
  		  int globalId = vertexMapper.template map<n>(e, sfs[k].entity());
  		  
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
   

    // the members of the struct are defined here
    struct VariableNodeData  
    {
       RT saturationW;
       RT pC;
       RT pN;
       VBlockType mobility;  //Vector with the number of phases
       VBlockType density;
    };
    

    // analog to EvalPrimaryData in MUFTE, uses members of varNData
    virtual void updateVariableData (const Entity& e, const FVElementGeometry& fvGeom, const VBlockType* sol)
    {
   	 int size = varNData.size();
     // cell center in reference element
//       const Dune::FieldVector<DT,n> 
//     cellLocal = Dune::ReferenceElements<DT,n>::general(gt).position(0,0);
  		
  	  // get global coordinate of cell center
//  	  const Dune::FieldVector<DT,n> cellGlobal = geometry.global(cellLocal);
   	  
        // ASSUME element-wise constant parameters for the material law 
//        FieldVector<RT, 4> parameters = problem.materialLawParameters(cellGlobal, e, cellLocal);

   	 
   	 for (int i = 0; i < size; i++) {
   		this->def[i] = 0;
   		varNData[i].saturationW = 1.0 - sol[i][satNIdx];

   		// ASSUME element-wise constant parameters for the material law 
         FieldVector<RT, 4> parameters = problem.materialLawParameters(fvGeom.cellGlobal, e, fvGeom.cellLocal);

         varNData[i].pC = problem.materialLaw().pC(varNData[i].saturationW, parameters);
         varNData[i].pN = sol[i][pWIdx] + varNData[i].pC;
         varNData[i].mobility[pWIdx] = problem.materialLaw().mobW(varNData[i].saturationW, parameters);
         varNData[i].mobility[satNIdx] = problem.materialLaw().mobN(sol[i][satNIdx], parameters);
         varNData[i].density[pWIdx] = problem.materialLaw().wettingPhase.density();
         varNData[i].density[satNIdx] = problem.materialLaw().nonwettingPhase.density();
   	 }   	 
    }
    

    template<class TypeTag>
		void assembleBC (const Entity& e)
		{
		  // extract some important parameters
		  Dune::GeometryType gt = e.geometry().type();
		  const typename Dune::LagrangeShapeFunctionSetContainer<DT,RT,n>::value_type& 
			sfs=Dune::LagrangeShapeFunctions<DT,RT,n>::general(gt,1);
		  setcurrentsize(sfs.size());

		  for (int i = 0; i < sfs.size(); i++)
			  this->b[i] = 0;
		  
		  // determine quadrature order
		  int p=2;
		  // evaluate boundary conditions via intersection iterator
		  typedef typename IntersectionIteratorGetter<G,TypeTag>::IntersectionIterator
		    IntersectionIterator;
		  
		  IntersectionIterator endit = IntersectionIteratorGetter<G,TypeTag>::end(e);
		  for (IntersectionIterator it = IntersectionIteratorGetter<G,TypeTag>::begin(e); 
		       it!=endit; ++it)
			{
			  // if we have a neighbor then we assume there is no boundary (forget interior boundaries)
			  // in level assemble treat non-level neighbors as boundary
			  if (it.neighbor())
				{
				  if (levelBoundaryAsDirichlet && it.outside()->level()==e.level()) 
					continue;
				  if (!levelBoundaryAsDirichlet)
					continue;
				}

			  // determine boundary condition type for this face, initialize with processor boundary
			  FieldVector<typename BoundaryConditions::Flags, m> bctypeface(BoundaryConditions::process);

			  // handle face on exterior boundary, this assumes there are no interior boundaries
			  if (it.boundary())
				{
				  Dune::GeometryType gtface = it.intersectionSelfLocal().type();
				  for (size_t g=0; g<Dune::QuadratureRules<DT,n-1>::rule(gtface,p).size(); ++g)
					{
					  const Dune::FieldVector<DT,n-1>& facelocal = Dune::QuadratureRules<DT,n-1>::rule(gtface,p)[g].position();
					  FieldVector<DT,n> local = it.intersectionSelfLocal().global(facelocal);
					  FieldVector<DT,n> global = it.intersectionGlobal().global(facelocal);
					  bctypeface = problem.bctype(global,e,it,local); // eval bctype


					  if (bctypeface[0]!=BoundaryConditions::neumann) break;

					  VBlockType J = problem.J(global,e,it,local);
					  if (J.two_norm() < 1e-10) 
						  continue;
					  double weightface = Dune::QuadratureRules<DT,n-1>::rule(gtface,p)[g].weight();
					  DT detjacface = it.intersectionGlobal().integrationElement(facelocal);
					  J *= 1.0/(pow(2.0, n-1))*weightface*detjacface;
					  for (int i=0; i<sfs.size(); i++) // loop over test function number
						if (this->bctype[i][0]==BoundaryConditions::neumann)
						  {
							//////////////////////////////////////////////////////////////////////////
							// HACK: piecewise constants with respect to dual grid not implemented yet 
							// works only if exactly one quadrature point is located within each dual 
							// cell boundary (which should be the case for p = 2)
							//////////////////////////////////////////////////////////////////////////
//							if (sfs[i].evaluateFunction(0,local) > 0.5) {
//								J *= weightface*detjacface;
								this->b[i] -= J;
//							}
						  }
					}
				  if (bctypeface[0]==BoundaryConditions::neumann) continue; // was a neumann face, go to next face
				}

			  // If we are here, then it is 
			  // (i)   an exterior boundary face with Dirichlet condition, or
			  // (ii)  a processor boundary (i.e. neither boundary() nor neighbor() was true), or
			  // (iii) a level boundary in case of level-wise assemble
			  // How processor boundaries are handled depends on the processor boundary mode
			  if (bctypeface[0]==BoundaryConditions::process && procBoundaryAsDirichlet==false 
				  && levelBoundaryAsDirichlet==false) 
				continue; // then it acts like homogeneous Neumann

			  // now handle exterior or interior Dirichlet boundary
			  for (int i=0; i<sfs.size(); i++) // loop over test function number
				{
				  if (sfs[i].codim()==0) continue; // skip interior dof
				  if (sfs[i].codim()==1) // handle face dofs
					{
					  if (sfs[i].entity()==it.numberInSelf())
						{
						  if (this->bctype[i][0]<bctypeface[0])
							{
							  this->bctype[i].assign(bctypeface[0]);
							  if (bctypeface[0]==BoundaryConditions::process)
								this->b[i] = 0;
							  if (bctypeface[0]==BoundaryConditions::dirichlet)
								{
								  this->b[i] = 0;
								}
							}
						}
					  continue;
					}
				  // handle subentities of this face
				  for (int j=0; j<ReferenceElements<DT,n>::general(gt).size(it.numberInSelf(),1,sfs[i].codim()); j++)
					if (sfs[i].entity()==ReferenceElements<DT,n>::general(gt).subEntity(it.numberInSelf(),1,j,sfs[i].codim()))
					  {
						if (this->bctype[i][0]<bctypeface[0])
						  {
						    this->bctype[i].assign(bctypeface[0]);
							if (bctypeface[0]==BoundaryConditions::process)
							  this->b[i] = 0;
							if (bctypeface[0]==BoundaryConditions::dirichlet)
							  {
								this->b[i] = 0;
							  }
						  }
					  }
				}
			}
		}

    // parameters given in constructor
    TwoPhaseProblem<G,RT>& problem;
    bool levelBoundaryAsDirichlet;
    bool procBoundaryAsDirichlet;
    const BoxFunction& currentSolution;
    BoxFunction oldSolution;
    VertexMapper vertexMapper;
    std::vector<StaticNodeData> statNData;
    std::vector<VariableNodeData> varNData;
    
  public:
    double dt;
    VBlockType uold[SIZE];
  };

  /** @} */
}
#endif
