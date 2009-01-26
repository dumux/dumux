// $Id$

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
#include"dumux/operators/boxjacobian.hh"
#include"dumux/twophase/twophaseproblem.hh"

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
   * @brief compute local jacobian matrix for the box method for two-phase equation
   *
   */


  //! A class for computing local jacobian matrices
  /*! A class for computing local jacobian matrix for the
    fully coupled two-phase model with Pw and Sn as primary variables

    Uses the box scheme.
    It should work for all dimensions and element types.
    All the numbering is with respect to the reference element and the
    Lagrange shape functions

    Template parameters are:

    - Grid    a DUNE grid type
    - Scalar  type used for return values
  */
  template<class Grid, class Scalar, class BoxFunction = LeafP1FunctionExtended<Grid, Scalar, 2> >
  class BoxPwSnJacobian
    : public BoxJacobian<BoxPwSnJacobian<Grid,Scalar,BoxFunction>,Grid,Scalar,2,BoxFunction>
  {
    typedef typename Grid::Traits::template Codim<0>::Entity Entity;
    typedef typename Entity::Geometry Geometry;
    typedef BoxPwSnJacobian<Grid,Scalar,BoxFunction> ThisType;
    typedef typename LocalJacobian<ThisType,Grid,Scalar,2>::VBlockType VBlockType;
    typedef typename LocalJacobian<ThisType,Grid,Scalar,2>::MBlockType MBlockType;
    enum {pWIdx = 0, satNIdx = 1};
    enum {wPhase = 0, nPhase = 1};


  public:
        // define the number of components of your system, this is used outside
        // to allocate the correct size of (dense) blocks with a FieldMatrix
        enum {dim=Grid::dimension};
        enum {numEq=2};
        enum {SIZE=LagrangeShapeFunctionSetContainer<Scalar,Scalar,dim>::maxsize};
        struct VariableNodeData;
        typedef FieldMatrix<Scalar,dim,dim> FMatrix;
        typedef FieldVector<Scalar,dim> FVector;

         //! Constructor
    BoxPwSnJacobian (TwoPhaseProblem<Grid,Scalar>& params,
                  bool levelBoundaryAsDirichlet_, const Grid& grid,
                  BoxFunction& sol,
                  bool procBoundaryAsDirichlet_=true)
    : BoxJacobian<ThisType,Grid,Scalar,2,BoxFunction>(levelBoundaryAsDirichlet_, grid, sol, procBoundaryAsDirichlet_),
      problem(params),
      statNData(this->vertexMapper.size()), varNData(SIZE), oldVarNData(SIZE)
    {
      this->analytic = false;
    }

    void clearVisited ()
    {
        return;
    }

    // compute storage term
    VBlockType computeM (const Entity& element, const VBlockType* sol,
            int node, const std::vector<VariableNodeData>& varData)
    {
        VBlockType result;
        //std::cout << "rhoW = " << varData[node].density[pWIdx] << ", rhoN = " << varData[node].density[satNIdx] << std::endl;
        result[wPhase] = varData[node].density[wPhase]*elData.porosity*varData[node].satW;//sol[node][satNIdx];
        result[nPhase] = varData[node].density[nPhase]*elData.porosity*varData[node].satN;//*sol[node][satNIdx];

        return result;
    };

    VBlockType computeM (const Entity& element, const VBlockType* sol, int node, bool old = false)
    {
        if (old)
            return computeM(element, sol, node, oldVarNData);
        else
            return computeM(element, sol, node, varNData);
    }

    // compute flux term
    VBlockType computeA (const Entity& element, const VBlockType* sol, int face)
    {
         int i = this->fvGeom.subContVolFace[face].i;
          int j = this->fvGeom.subContVolFace[face].j;

          // permeability in edge direction
          FieldVector<Scalar,dim> Kij(0);
          elData.K.umv(this->fvGeom.subContVolFace[face].normal, Kij);

          VBlockType flux;
          for (int phase = 0; phase < numEq; phase++)
          {
              // calculate FE gradient of the pressure
              FieldVector<Scalar, dim> pGrad(0);
              Scalar densityIJ = 0;
              for (int vert = 0; vert < this->fvGeom.numVertices; vert++)
              {
                  FieldVector<Scalar,dim> grad(this->fvGeom.subContVolFace[face].grad[vert]);
                  if (phase == wPhase)
                      grad *= sol[vert][pWIdx];
                  else if (phase == nPhase)
                      grad *= varNData[vert].pN;
                  pGrad += grad;

                  densityIJ += varNData[vert].density[phase]*this->fvGeom.subContVolFace[face].shapeValue[vert];
              }

              // adjust by gravity
              FieldVector<Scalar, dim> gravity = problem.gravity();
              gravity *= densityIJ;
              pGrad -= gravity;

              // calculate the flux using fully upwind
              Scalar outward = pGrad*Kij;
              if (outward < 0)
                  flux[phase] = varNData[i].density[phase]*varNData[i].mobility[phase]*outward;
              else
                  flux[phase] = varNData[j].density[phase]*varNData[j].mobility[phase]*outward;
          }

          return flux;
   };


    // compute sources/sinks
    VBlockType computeQ (const Entity& element, const VBlockType* sol, const int& node)
    {
        // ASSUME problem.q already contains \rho.q
        return problem.q(this->fvGeom.subContVol[node].global, element, this->fvGeom.subContVol[node].local);
    };


      // *********************************************************
      // *                                                         *
      // *    Calculation of Data at Elements                          *
      // *                                                          *
      // *                                                         *
      // *********************************************************

    void computeElementData (const Entity& element)
    {
         // ASSUMING element-wise constant permeability and porosity, evaluate at the element center
          elData.K = this->problem.soil().K(this->fvGeom.elementGlobal, element, this->fvGeom.elementLocal);
          elData.porosity = this->problem.soil().porosity(this->fvGeom.elementGlobal, element, this->fvGeom.elementLocal);
    };


      // *********************************************************
      // *                                                         *
      // *    Calculation of Data at Nodes that has to be             *
      // *    determined only once    (statNData)                     *
      // *                                                         *
      // *********************************************************

    // analog to EvalStaticData in MUFTE
    void updateStaticData (const Entity& element, const VBlockType* sol)
    {
      return;
    }


      //*********************************************************
      //*                                                        *
      //*    Calculation of variable Data at Nodes                 *
      //*    (varNData)                                             *
      //*                                                         *
      //*********************************************************


    // the members of the struct are defined here
    struct VariableNodeData
    {
        Scalar satN;
        Scalar satW;
        Scalar pC;
        Scalar pW;
        Scalar pN;
        VBlockType mobility;  //Vector with the number of phases
        VBlockType density;
        FieldMatrix<Scalar,dim,dim> K;
    };

    // analog to EvalPrimaryData in MUFTE, uses members of varNData
    void updateVariableData(const Entity& element, const VBlockType* sol,
            int i, std::vector<VariableNodeData>& varData)
    {
           const int globalIdx = this->vertexMapper.template map<dim>(element, i);
        FVector& global = this->fvGeom.subContVol[i].global;
             FVector& local = this->fvGeom.subContVol[i].local;

        varData[i].satN = sol[i][satNIdx];
        varData[i].satW = 1.0 - sol[i][satNIdx];
           varData[i].pC = problem.materialLaw().pC(varData[i].satW, global, element, local);
        varData[i].pW = sol[i][pWIdx];
           varData[i].pN = sol[i][pWIdx] + varData[i].pC;
           varData[i].mobility[wPhase] = problem.materialLaw().mobW(varData[i].satW, global, element, local);
           varData[i].mobility[nPhase] = problem.materialLaw().mobN(varData[i].satN, global, element, local);
           varData[i].density[wPhase] = problem.wettingPhase().density(283.15, varData[i].pW);
           varData[i].density[nPhase] = problem.nonwettingPhase().density(283.15, varData[i].pN);

           // for output
         (*outPressureN)[globalIdx] = varData[i].pN;
         (*outCapillaryP)[globalIdx] = varData[i].pC;
         (*outSaturationW)[globalIdx] = varData[i].satW;
         (*outSaturationN)[globalIdx] = varData[i].satN;
         (*outDensityW)[globalIdx] = varData[i].density[wPhase];
         (*outDensityN)[globalIdx] = varData[i].density[nPhase];
         (*outMobilityW)[globalIdx] = varData[i].mobility[wPhase];
         (*outMobilityN)[globalIdx] = varData[i].mobility[nPhase];

              return;
    }

    void updateVariableData(const Entity& element, const VBlockType* sol, int i, bool old = false)
    {
        if (old) {
            updateVariableData(element, sol, i, oldVarNData);
        }
        else
            updateVariableData(element, sol, i, varNData);
    }

    void updateVariableData(const Entity& element, const VBlockType* sol, bool old = false)
    {
        int size = this->fvGeom.numVertices;

        for (int i = 0; i < size; i++)
                updateVariableData(element, sol, i, old);
    }

    struct StaticNodeData
    {
        bool visited;
    };

    struct ElementData {
        Scalar cellVolume;
     Scalar porosity;
        FieldMatrix<Scalar,dim,dim> K;
        } elData;

     // parameters given in constructor
     TwoPhaseProblem<Grid,Scalar>& problem;
     std::vector<StaticNodeData> statNData;
     std::vector<VariableNodeData> varNData;
     std::vector<VariableNodeData> oldVarNData;

     // for output files
     BlockVector<FieldVector<Scalar, 1> > *outPressureN;
     BlockVector<FieldVector<Scalar, 1> > *outCapillaryP;
     BlockVector<FieldVector<Scalar, 1> > *outSaturationN;
     BlockVector<FieldVector<Scalar, 1> > *outSaturationW;
     BlockVector<FieldVector<Scalar, 1> > *outDensityW;
     BlockVector<FieldVector<Scalar, 1> > *outDensityN;
     BlockVector<FieldVector<Scalar, 1> > *outMobilityW;
     BlockVector<FieldVector<Scalar, 1> > *outMobilityN;


  };

  /** @} */
}
#endif
