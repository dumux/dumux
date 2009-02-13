// $Id$

#ifndef DUNE_BOX1P2CJACOBIAN_HH
#define DUNE_BOX1P2CJACOBIAN_HH

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

#include"dumux/operators/boxjacobian.hh"
#include"dumux/1p2c/1p2cproblem.hh"
#include"dumux/io/vtkmultiwriter.hh"


/**
 * @file
 * @brief  compute local jacobian matrix for box scheme for one-phase two-component flow equation
 * @author Karin Erbertseder
 */



namespace Dune
{
  /** @addtogroup DISC_Disc
   *
   * @{
   */
  /**
   * @brief compute local jacobian matrix for the boxfile for one-phase two-component flow equation
   *
   */


  //! Derived class for computing local jacobian matrices
  /*! A class for computing local jacobian matrix for the one-phase two-component flow equation

	    div j = q; j = -K grad u; in Omega

		u = g on Gamma1; j*n = J on Gamma2.

	Uses box scheme with the Lagrange shape functions.
	It should work for all dimensions and element types.
	All the numbering is with respect to the reference element and the
	Lagrange shape functions

	Template parameters are:

	- Grid  a DUNE grid type
	- Scalar    type used for return values

	A class for computing the transport equation with the primary variables x and p.
	The transport equation is coupled with the continuity equation.
  */

	//the influence of the gravity is neglected

  template<class Grid, class Scalar, class BoxFunction = LeafP1Function<Grid, Scalar, 2> >
  class Box1P2CJacobian
    : public BoxJacobian<Box1P2CJacobian<Grid,Scalar,BoxFunction>,Grid,Scalar,2,BoxFunction>
  {
    typedef typename Grid::ctype CoordScalar;
    typedef typename Grid::Traits::template Codim<0>::Entity Element;
    typedef typename Element::Geometry Geometry;
    typedef Box1P2CJacobian<Grid,Scalar,BoxFunction> ThisType;
    typedef typename LocalJacobian<ThisType,Grid,Scalar,2>::VBlockType VBlockType;
    typedef typename LocalJacobian<ThisType,Grid,Scalar,2>::MBlockType MBlockType;
    typedef FVElementGeometry<Grid> FVElementGeometry;

 	enum {konti = 0, transport = 1};	// Solution vector index


  public:
    enum {dim=Grid::dimension};
    enum {numComp=2};					//number of components
    enum {SIZE=LagrangeShapeFunctionSetContainer<CoordScalar,Scalar,dim>::maxsize};
    struct VariableNodeData;

    typedef FieldMatrix<Scalar,dim,dim> FMatrix;
    typedef FieldVector<Scalar,dim> FVector;

    //! Constructor
    Box1P2CJacobian (OnePTwoCProblem<Grid,Scalar>& params, bool levelBoundaryAsDirichlet_, const Grid& grid,
			      BoxFunction& sol, bool procBoundaryAsDirichlet_=true)
    : BoxJacobian<ThisType,Grid,Scalar,2,BoxFunction>(levelBoundaryAsDirichlet_, grid, sol, procBoundaryAsDirichlet_),
      problem(params), vNDat(SIZE), oldVNDat(SIZE)
    {
      this->analytic = false;
    }

	/** @brief compute time dependent term (storage), loop over nodes / subcontrol volumes
	 *  @param e element
	 *  @param sol solution vector
	 *  @param node local node id
	 *  @return storage term
	 */
    virtual VBlockType computeM (const Element& element, const VBlockType* sol,
    		int node, std::vector<VariableNodeData>& varData)
    {

   	 VBlockType result;

   	 result[0] = 0; //continuity equation has no storage term.

   	 result[1] = vNDat[node].porosity * vNDat[node].molefraction; //storage term of the transport equation

   	 //std::cout << result << " " << node << std::endl;
   	 return result;
    };

    virtual VBlockType computeM (const Element& element, const VBlockType* sol, int node, bool old = false)
    {
    	if (old)
    		return computeM(element, sol, node, oldVNDat);
    	else
    		return computeM(element, sol, node, vNDat);
    }

    /** @brief compute diffusive/advective fluxes, loop over subcontrol volume faces
	 *  @param e element
	 *  @param sol solution vector
	 *  @param face face id
	 *  @return flux term
     */
    virtual VBlockType computeA (const Element& element, const VBlockType* sol, int face)
    {
   	 int i = this->fvGeom.subContVolFace[face].i;
 	 int j = this->fvGeom.subContVolFace[face].j;

 	 // normal vector, value of the area of the scvf
	 const FieldVector<Scalar,dim> normal(this->fvGeom.subContVolFace[face].normal);

  	 // get global coordinates of nodes i,j
	 const FieldVector<CoordScalar,dim> global_i = this->fvGeom.subContVol[i].global;
	 const FieldVector<CoordScalar,dim> global_j = this->fvGeom.subContVol[j].global;
	 const FieldVector<CoordScalar,dim> local_i = this->fvGeom.subContVol[i].local;
	 const FieldVector<CoordScalar,dim> local_j = this->fvGeom.subContVol[j].local;

	 VBlockType flux;

	 //calculation of the flux of the continuity equation
	 FieldVector<Scalar,dim> gradP(0);
	 for (int k = 0; k < this->fvGeom.numVertices; k++) {
	     FieldVector<Scalar,dim> grad(this->fvGeom.subContVolFace[face].grad[k]);
	     grad *= sol[k][0];
	     gradP += grad;
	 }

	 FieldVector<Scalar,dim> KGradP(0);
	 FMatrix Ki(0), Kj(0);

	 // calculate harmonic mean of permeabilities of nodes i and j
	 Ki = this->problem.soil().K(global_i,element,local_i);
	 Kj = this->problem.soil().K(global_j,element,local_j);
	 const FMatrix K = harmonicMeanK(Ki, Kj); //harmonic mean of the permeability

	 K.umv(gradP, KGradP);					 //KGradP = KGradP + gradP

	 flux[0] = ( KGradP* this->fvGeom.subContVolFace[face].normal) / vNDat[i].viscosity;

	 //calculation of the flux of the transport equation
	 FieldVector<Scalar,dim> gradX(0);
	 for (int k=0; k<this->fvGeom.numVertices; k++) {
	     FieldVector<Scalar,dim> grad(this->fvGeom.subContVolFace[face].grad[k]);
	     grad *= sol[k][1];
	     gradX += grad;
	 }

	 const Scalar T = arithmeticMeanT(i, j, element);	//arithmetic mean of the tortuosity T

	 const Scalar P = arithmeticMeanP(i, j, element);	//arithmetic mean of the porosity P

	 Scalar outward = KGradP * this->fvGeom.subContVolFace[face].normal;													//calcualte the flux using upwind
	 if (outward < 0)
	    flux[1] = (outward * vNDat[i].molefraction / vNDat[i].viscosity) + (T * P * vNDat[i].diffCoeff * (gradX * this->fvGeom.subContVolFace[face].normal));
	 else
	    flux[1] = (outward * vNDat[j].molefraction / vNDat[j].viscosity) + (T * P * vNDat[i].diffCoeff * (gradX * this->fvGeom.subContVolFace[face].normal));

	 return flux;
	 }


      /** @brief integrate sources / sinks
       *  @param e element
       *  @param sol solution vector
       *  @param node local node id
       *  @return source/sink term
       */
      virtual VBlockType computeQ (const Element& element, const VBlockType* sol, const int& node)
          {
              // ASSUME problem.q already contains \rho.q
              return problem.q(this->fvGeom.subContVol[node].global, element, this->fvGeom.subContVol[node].local);
          }


    // harmonic mean of the permeability computed directly
      virtual FMatrix harmonicMeanK (FMatrix& Ki, const FMatrix& Kj) const
      {
      	double eps = 1e-20;

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

      //arithmetic mean of the tortuosity is computed directly
      virtual Scalar arithmeticMeanT (const int i, const int j, const Element& element) const
      {

      	Scalar Ti, Tj;

      	Ti = this->problem.soil().tortuosity(this->fvGeom.subContVol[i].global, element, this->fvGeom.subContVol[i].local);
      	Tj = this->problem.soil().tortuosity(this->fvGeom.subContVol[j].global, element, this->fvGeom.subContVol[j].local);

      	if (Ti != Tj)
      	    {
      	     Ti = (Ti + Tj)/2;
      	     return Ti;
      	    }
      	else{
      		return Ti;
      	}
      }

      //arithmetic mean of the porosity is computed directly
      virtual Scalar arithmeticMeanP (const int i, const int j, const Element& element) const
      {

       Scalar Pi, Pj;

       Pi = this->problem.soil().porosity(this->fvGeom.subContVol[i].global, element, this->fvGeom.subContVol[i].local);
       Pj = this->problem.soil().porosity(this->fvGeom.subContVol[j].global, element, this->fvGeom.subContVol[j].local);

       if (Pi != Pj)
          {
          Pi = (Pi + Pj)/2;
          return Pi;
          }
          else{
          return Pi;
          }
       }


    virtual void clearVisited ()
    {
    	return;
   	}

    // updates old phase state after each time step
    virtual void updatePhaseState ()
    {
       return;
    }

    virtual void resetPhaseState ()
    {
       return;
    }

	  //*********************************************************
	  //*														*
	  //*	Calculation of Data at Elements (elData) 			*
	  //*						 								*
	  //*														*
	  //*********************************************************

    virtual void computeElementData (const Element& element)
    {
//		 // ASSUMING element-wise constant permeability, evaluate K at the cell center
// 		 elData.K = problem.K(this->fvGeom.cellGlobal, element, this->fvGeom.cellLocal);
//
//		 // ASSUMING element-wise constant porosity
// 		 elData.porosity = problem.porosity(this->fvGeom.cellGlobal, element, this->fvGeom.cellLocal);
   	 return;
    }


	  //*********************************************************
	  //*														*
	  //*	Calculation of Data at Nodes that has to be			*
	  //*	determined only once	(sNDat)						*
	  //*														*
	  //*********************************************************

    // analog to EvalStaticData in MUFTE
    virtual void updateStaticDataVS (const Element& element, VBlockType* sol)
    {
   	 return;
    }

    // for initialization of the Static Data (sets porosity)
    virtual void updateStaticData (const Element& element, VBlockType* sol)
    {
   	 return;
    }


	  //*********************************************************
	  //*														*
	  //*	Calculation of variable Data at Nodes				*
	  //*	(vNDat)												*
	  //*														*
	  //*********************************************************


    struct VariableNodeData
    {
    	Scalar porosity;
    	Scalar molefraction;
    	Scalar viscosity;
    	Scalar tortuosity;
    	Scalar diffCoeff;
    	FMatrix K;
    };

    // analog to EvalPrimaryData in MUFTE, uses members of vNDat
	virtual void updateVariableData(const Element& element, const VBlockType* sol,
			int i, std::vector<VariableNodeData>& varData)
    {	double T=273; 		//Temperatur in Kelvin
    	double p=1e5; 		//Druck in Pa

		varData[i].porosity = problem.soil().porosity(this->fvGeom.subContVol[i].global, element, this->fvGeom.subContVol[i].local);
   		varData[i].viscosity = problem.phase().viscosity(T, p);
   		varData[i].tortuosity = problem.soil().tortuosity(this->fvGeom.subContVol[i].global, element, this->fvGeom.subContVol[i].local);
	   	varData[i].diffCoeff = problem.phase().diffCoeff();
	   	varData[i].molefraction = sol[i][1];
    }


 	virtual void updateVariableData(const Element& element, const VBlockType* sol, int i, bool old = false)
 		{

 			if (old)
 				updateVariableData(element, sol, i, oldVNDat);
 			else
 				updateVariableData(element, sol, i, vNDat);
 		}

	void updateVariableData (const Element& element, const VBlockType* sol, bool old = false)
	{
		int size = this->fvGeom.numVertices;
		vNDat.resize(size);
		oldVNDat.resize(size);

		for (int i = 0; i < size; i++)
			updateVariableData(element, sol, i, old);
	}


	struct StaticNodeData
    {
   	 bool visited;
    };




    // parameters given in constructor
    OnePTwoCProblem<Grid,Scalar>& problem;
    std::vector<StaticNodeData> sNDat;
    std::vector<VariableNodeData> vNDat;
    std::vector<VariableNodeData> oldVNDat;


  };

}
#endif
