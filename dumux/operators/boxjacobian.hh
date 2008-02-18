#ifndef DUNE_BOXJACOBIAN_HH
#define DUNE_BOXJACOBIAN_HH

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
  template<class Imp, class G, class RT, int m, class BoxFunction = LeafP1Function<G, RT, m> >
  class BoxJacobian 
    : public LocalJacobian<Imp,G,RT,m>
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
    typedef typename LocalJacobian<Imp,G,RT,m>::VBlockType VBlockType;
    typedef typename LocalJacobian<Imp,G,RT,m>::MBlockType MBlockType;
 	typedef FVElementGeometry<G> FVElementGeometry;
	typedef MultipleCodimMultipleGeomTypeMapper<G, typename G::Traits::LeafIndexSet, P1Layout> VertexMapper;	
	
  public:
    // define the number of components of your system, this is used outside
    // to allocate the correct size of (dense) blocks with a FieldMatrix
    enum {n=G::dimension};
    enum {SIZE=LagrangeShapeFunctionSetContainer<DT,RT,n>::maxsize};
    
    //! Constructor
    BoxJacobian (bool levelBoundaryAsDirichlet_, const G& grid, 
			      BoxFunction& sol, 
			      bool procBoundaryAsDirichlet_=true)
    : levelBoundaryAsDirichlet(levelBoundaryAsDirichlet_), 
      procBoundaryAsDirichlet(procBoundaryAsDirichlet_), 
      currentSolution(sol), oldSolution(grid), 
      vertexMapper(grid, grid.leafIndexSet()), dt(1)
    {    }
    

    //**********************************************************
    //*																			*
    //*	Computation of the local defect								*
    //*																			*
    //**********************************************************
    
    template<class TypeTag>
    void localDefect (const Entity& e, const FVElementGeometry& fvGeom, const VBlockType* sol)
    {
      computeElementData(e, fvGeom);
      updateVariableData(e, fvGeom, sol);

	  for (int i=0; i < fvGeom.nodes; i++) // begin loop over vertices / sub control volumes
	  {
		  // implicit Euler
		  VBlockType massContrib = computeM(e, fvGeom, sol, i);
		  massContrib -= computeM(e, fvGeom, uold, i);
		  massContrib *= fvGeom.subContVol[i].volume/dt;
		  this->def[i] += massContrib;
		  
		  // get source term 
		  VBlockType q = computeQ(e, fvGeom, sol, i);
		  q *= fvGeom.subContVol[i].volume;
		  this->def[i] -= q;
	  } // end loop over vertices / sub control volumes
	  
	  for (int face = 0; face < e.template count<n-1>(); face++) // begin loop over edges / sub control volume faces
	  {
		  int i = fvGeom.subContVolFace[face].i;
		  int j = fvGeom.subContVolFace[face].j;
		  
		  VBlockType flux = computeA(e, fvGeom, sol, face);
		  
		  // obtain integrated Flux 
		  flux *= fvGeom.subContVolFace[face].area; 
			  
		  // add to defect 
		  this->def[i] -= flux;
		  this->def[j] += flux;
	  } // end loop over edges / sub control volume faces

	  // assemble boundary conditions 
	  assembleBC<TypeTag> (e); 
	  
	  // add to defect 
	  for (int i=0; i < fvGeom.nodes; i++) {
		  this->def[i] -= this->b[i];
	  }

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
  
    VBlockType computeM (const Entity& e, const FVElementGeometry& fvGeom, const VBlockType* sol, int node)
    {
   	 return this->getImp().computeM(e, fvGeom, sol, node);
    };
    
    VBlockType computeQ (const Entity& e, const FVElementGeometry& fvGeom, const VBlockType* sol, int node)
    {
   	 return this->getImp().computeQ(e, fvGeom, sol, node);
    };
    
    VBlockType computeA (const Entity& e, const FVElementGeometry& fvGeom, const VBlockType* sol, int face)
    {
   	 return this->getImp().computeA(e, fvGeom, sol, face);
    };

    void computeElementData (const Entity& e, const FVElementGeometry& fvGeom)
    {
   	 return this->getImp().computeElementData(e, fvGeom);
    };

    // analog to EvalStaticData in MUFTE
    virtual void updateStaticData (const Entity& e, const FVElementGeometry& fvGeom, const VBlockType* sol)
    {
   	 return this->getImp().updateStaticData(e, fvGeom, sol);
    }

    // analog to EvalPrimaryData in MUFTE, uses members of varNData
    virtual void updateVariableData (const Entity& e, const FVElementGeometry& fvGeom, const VBlockType* sol)
    {
   	 return this->getImp().updateVariableData(e, fvGeom, sol);
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
					  bctypeface = this->getImp().problem.bctype(global,e,it,local); // eval bctype


					  if (bctypeface[0]!=BoundaryConditions::neumann) break;

					  VBlockType J = this->getImp().problem.J(global,e,it,local);
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
    bool levelBoundaryAsDirichlet;
    bool procBoundaryAsDirichlet;
    const BoxFunction& currentSolution;
    BoxFunction oldSolution;
    VertexMapper vertexMapper;
    
  public:
    double dt;
    VBlockType uold[SIZE];
  };

  /** @} */
}
#endif
