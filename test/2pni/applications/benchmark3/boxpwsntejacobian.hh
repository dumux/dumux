// $Id$

#ifndef DUNE_BOXPWSNTEJACOBIAN_HH
#define DUNE_BOXPWSNTEJACOBIAN_HH

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
#include"dumux/2pni/2pniproblem.hh"
#include"dumux/io/vtkmultiwriter.hh"

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
   * @brief compute local jacobian matrix for the box method for nonisothermal two-phase equation
   *
   */


  //! A class for computing local jacobian matrices
  /*! A class for computing local jacobian matrix for the
    fully coupled two-phase model with Pw, Sn and Te as primary variables

    Uses the box scheme.
    It should work for all dimensions and element types.
    All the numbering is with respect to the reference element and the
    Lagrange shape functions

    Template parameters are:

    - Grid  a DUNE grid type
    - RT    type used for return values
  */
  template<class G, class RT, class BoxFunction = LeafP1FunctionExtended<G, RT, 3> >
  class BoxPwSnTeJacobian
    : public BoxJacobian<BoxPwSnTeJacobian<G,RT,BoxFunction>,G,RT,3,BoxFunction>
  {
    typedef typename G::ctype DT;
    typedef typename G::Traits::template Codim<0>::Entity Entity;
    typedef typename Entity::Geometry Geometry;
    typedef BoxPwSnTeJacobian<G,RT,BoxFunction> ThisType;
    typedef typename LocalJacobian<ThisType,G,RT,3>::VBlockType VBlockType;
    typedef typename LocalJacobian<ThisType,G,RT,3>::MBlockType MBlockType;
    enum {pWIdx = 0, satNIdx = 1, teIdx = 2};
    enum {wPhase = 0, nPhase = 1, heat = 2};


  public:
        // define the number of components of your system, this is used outside
        // to allocate the correct size of (dense) blocks with a FieldMatrix
        enum {dim=G::dimension};
        enum {m=3};
        enum {SIZE=LagrangeShapeFunctionSetContainer<DT,RT,dim>::maxsize};
        struct VariableNodeData;
        typedef FieldMatrix<RT,dim,dim> FMatrix;
        typedef FieldVector<RT,dim> FVector;

         //! Constructor
    BoxPwSnTeJacobian (TwoPhaseHeatProblem<G,RT>& params,
                  bool levelBoundaryAsDirichlet_, const G& grid,
                  BoxFunction& sol,
                  bool procBoundaryAsDirichlet_=true)
    : BoxJacobian<ThisType,G,RT,m,BoxFunction>(levelBoundaryAsDirichlet_, grid, sol, procBoundaryAsDirichlet_),
      problem(params),
      sNDat(this->vertexMapper.size()), vNDat(SIZE), oldVNDat(SIZE)
    {
      this->analytic = false;
    }

    /** @brief compute time dependent term (storage), loop over nodes / subcontrol volumes
     *  @param e entity
     *  @param sol solution vector
     *  @param node local node id
     *  @return storage term
     */
    virtual VBlockType computeM (const Entity& e, const VBlockType* sol,
            int node, std::vector<VariableNodeData>& varData)
    {
         GeometryType gt = e.geometry().type();
         const typename LagrangeShapeFunctionSetContainer<DT,RT,dim>::value_type&
                     sfs=LagrangeShapeFunctions<DT,RT,dim>::general(gt,1);

        int globalIdx = this->vertexMapper.template map<dim>(e, sfs[node].entity());

        VBlockType result;
        RT satN = varData[node].satN;
        RT satW = varData[node].satW;

        // storage of component wPhase
        result[wPhase] =
             sNDat[globalIdx].porosity*varData[node].density[wPhase]*satW;
        // storage of component co2
        result[nPhase] =
             sNDat[globalIdx].porosity*varData[node].density[nPhase]*satN;

        // storage term of energy equation
        result[heat] = sNDat[globalIdx].porosity * (varData[node].density[wPhase] * varData[node].intenergy[wPhase] * satW
                    + varData[node].density[nPhase] * varData[node].intenergy[nPhase] * satN)
                    + sNDat[globalIdx].heatCap * varData[node].temperature;
         // soil properties defined at the elements!!!

        return result;
    };

    virtual VBlockType computeM (const Entity& e, const VBlockType* sol, int node, bool old = false)
     {
         if (old)
             return computeM(e, sol, node, oldVNDat);
         else
             return computeM(e, sol, node, vNDat);
     }

    /** @brief compute fluxes and heat conduction, loop over subcontrol volume faces
     *  @param e entity
     *  @param sol solution vector
     *  @param face face id
     *  @return flux term
     */
    VBlockType computeA (const Entity& e, const VBlockType* sol, int face)
    {
         int i = this->fvGeom.subContVolFace[face].i;
          int j = this->fvGeom.subContVolFace[face].j;

      // normal vector, value of the area of the scvf
     const FieldVector<RT,dim> normal(this->fvGeom.subContVolFace[face].normal);
     GeometryType gt = e.geometry().type();
     const typename LagrangeShapeFunctionSetContainer<DT,RT,dim>::value_type&
      sfs=LagrangeShapeFunctions<DT,RT,dim>::general(gt,1);

      // global index of the subcontrolvolume face neighbour nodes in element e
     int globalIdx_i = this->vertexMapper.template map<dim>(e, sfs[i].entity());
       int globalIdx_j = this->vertexMapper.template map<dim>(e, sfs[j].entity());

       // get global coordinates of nodes i,j
     const FieldVector<DT,dim> global_i = this->fvGeom.subContVol[i].global;
     const FieldVector<DT,dim> global_j = this->fvGeom.subContVol[j].global;
     const FieldVector<DT,dim> local_i = this->fvGeom.subContVol[i].local;
     const FieldVector<DT,dim> local_j = this->fvGeom.subContVol[j].local;

     /////////////////////////////////////////////////////////////////////////////
     // AVERAGING to get parameter values at integration points
         // Harmonic Mean:

         // harmonic mean of the permeability
            FMatrix Ki(0.), Kj(0.);
         Ki = this->problem.soil().K(global_i,e,local_i,i);
         Kj = this->problem.soil().K(global_j,e,local_j,j);
         const FMatrix K = harmonicMeanK(Ki, Kj);

         // harmonic mean of the heat conductivity lambda
         RT lambda;
         lambda = 2./((1./vNDat[i].lambda) + (1./vNDat[j].lambda));

           // Calculate arithmetic mean of the densities
             VBlockType avgDensity;
            avgDensity[wPhase] = 0.5*(vNDat[i].density[wPhase] + vNDat[j].density[wPhase]);
            avgDensity[nPhase] = 0.5*(vNDat[i].density[nPhase] + vNDat[j].density[nPhase]);
    //////////////////////////////////////////////////////////////////////////////////////////
    // GRADIENTS

     FieldMatrix<RT,2,dim> pGrad(0.);
     FieldVector<RT,dim> teGrad(0.);
     FieldVector<RT,dim> temp(0.);
     VBlockType flux(0.);

     // calculate FE gradient at subcontrolvolumeface
         for (int k = 0; k < this->fvGeom.numVertices; k++) // loop over adjacent nodes
         {
             // FEGradient at subcontrolvolumeface face
             const FieldVector<DT,dim> feGrad(this->fvGeom.subContVolFace[face].grad[k]);
             FieldVector<RT,2> pressure(0.0);

             pressure[wPhase] = vNDat[k].pW;
             pressure[nPhase] = vNDat[k].pN;

             // compute pressure gradients for each phase at integration point of subcontrolvolumeface face
             for (int phase = 0; phase < 2; phase++)
               {
                   temp = feGrad;
                   temp *= pressure[phase];
                   pGrad[phase] += temp;
               }

             // compute temperature gradient
               temp = feGrad;
              temp *= vNDat[k].temperature;
              teGrad += temp;
         }

      // deduce gravity*density of each phase
         FieldMatrix<RT,2,dim> contribComp(0);
         for (int phase=0; phase<2; phase++)
         {
             contribComp[phase] = problem.gravity();
             contribComp[phase] *= avgDensity[phase];
             pGrad[phase] -= contribComp[phase]; // grad p - rho*g
         }

     // Darcy velocity in normal direction for each phase K*n(grad p -rho*g)
         VBlockType outward(0);
         FieldVector<RT,dim> v_tilde_w(0);
         FieldVector<RT,dim> v_tilde_n(0);

          K.umv(pGrad[wPhase], v_tilde_w);  // v_tilde=K*gradP
         outward[wPhase] = v_tilde_w*normal;
         K.umv(pGrad[nPhase], v_tilde_n);  // v_tilde=K*gradP
         outward[nPhase] = v_tilde_n*normal;

     // Heat conduction
         outward[heat] = teGrad * normal;
           outward[heat] *= lambda;

     // evaluate upwind nodes
         int up_w, dn_w, up_n, dn_n;
         if (outward[wPhase] <= 0) {up_w = i; dn_w = j;}
         else {up_w = j; dn_w = i;};
         if (outward[nPhase] <= 0) {up_n = i; dn_n = j;}
         else {up_n = j; dn_n = i;};

     RT alpha = 1.0;  // Upwind parameter

     ////////////////////////////////////////////////////////////////////////////////////////////////
     // ADVECTIVE TRANSPORT
         // Water conservation
          flux[wPhase] =   (alpha* vNDat[up_w].density[wPhase]*vNDat[up_w].mobility[wPhase]
                     + (1-alpha)* vNDat[dn_w].density[wPhase]*vNDat[dn_w].mobility[wPhase])
                                * outward[wPhase];
          // CO2 conservation
          flux[nPhase]   =   (alpha* vNDat[up_n].density[nPhase]*vNDat[up_n].mobility[nPhase]
                     + (1-alpha)* vNDat[dn_n].density[nPhase]*vNDat[dn_n].mobility[nPhase])
                                * outward[nPhase];
          // Heat conservation
         flux[heat]  =  (alpha* vNDat[up_n].density[nPhase]*vNDat[up_n].mobility[nPhase] * vNDat[up_n].enthalpy[nPhase]
                     + (1-alpha)* vNDat[dn_n].density[nPhase]*vNDat[dn_n].mobility[nPhase] * vNDat[dn_n].enthalpy[nPhase])
                     * outward[nPhase];
          flux[heat]  +=  (alpha* vNDat[up_w].density[wPhase]*vNDat[up_w].mobility[wPhase] * vNDat[up_w].enthalpy[wPhase]
                     + (1-alpha)* vNDat[dn_w].density[wPhase]*vNDat[dn_w].mobility[wPhase] * vNDat[dn_w].enthalpy[wPhase])
                     * outward[wPhase];

      //////////////////////////////////////////////////////////////////////////////////////////////
      // HEAT CONDUCTION

          flux[heat]    +=    outward[heat];

         return flux;
  };

      /** @brief integrate sources / sinks
       *  @param e entity
     *  @param sol solution vector
     *  @param node local node id
     *  @return source/sink term
     */
       virtual VBlockType computeQ (const Entity& e, const VBlockType* sol, const int& node)
       {
           // ASSUME problem.q already contains \rho.q
               VBlockType result(0);
                 result[0] = problem.q(this->fvGeom.subContVol[node].global, e, this->fvGeom.subContVol[node].local)[0];
                 result[1] = problem.q(this->fvGeom.subContVol[node].global, e, this->fvGeom.subContVol[node].local)[1];
                 result[2] = result[0] * vNDat[node].enthalpy[wPhase] + result[1] * vNDat[node].enthalpy[nPhase];
                // ASSUME problem.q already contains \rho.q
                return result;
}

       // harmonic mean computed directly
    virtual FMatrix harmonicMeanK (FMatrix& Ki, const FMatrix& Kj) const
    {
        double eps = 1e-20;
        FMatrix K(0.);
        for (int kx=0; kx<dim; kx++){
            for (int ky=0; ky<dim; ky++){
                if (Ki[kx][ky] != Kj[kx][ky])
                {
                    K[kx][ky] = 2 / (1/(Ki[kx][ky]+eps) + (1/(Kj[kx][ky]+eps)));
                }
                else K = Ki;
            }
        }
        return K;
    }

       virtual void clearVisited ()
    {
        for (int i = 0; i < this->vertexMapper.size(); i++){
           sNDat[i].visited = false;
        }
        return;
       }

          // *********************************************************
      // *                                                         *
      // *    Calculation of Data at Elements                          *
      // *                                                          *
      // *                                                         *
      // *********************************************************
    void computeElementData (const Entity& e)
    {
    };

      // *********************************************************
      // *                                                         *
      // *    Calculation of Data at Nodes that has to be             *
      // *    determined only once    (statNData)                     *
      // *                                                         *
      // *********************************************************

    // analog to EvalStaticData in MUFTE
    void updateStaticData (const Entity& e, const VBlockType* sol)
    {
        // size of the sNDat vector is determined in the constructor

        // local to global id mapping (do not ask vertex mapper repeatedly
        //int localToGlobal[LagrangeShapeFunctionSetContainer<DT,RT,n>::maxsize];

        // get access to shape functions for P1 elements
        GeometryType gt = e.geometry().type();
        const typename LagrangeShapeFunctionSetContainer<DT,RT,dim>::value_type&
        sfs=LagrangeShapeFunctions<DT,RT,dim>::general(gt,1);

        // get local to global id map
        for (int k = 0; k < sfs.size(); k++) {
           const int globalIdx = this->vertexMapper.template map<dim>(e, sfs[k].entity());

           // if nodes are not already visited
          if (!sNDat[globalIdx].visited)
           {
                // ASSUME porosity defined at nodes
                sNDat[globalIdx].porosity = problem.soil().porosity(this->fvGeom.subContVol[k].global, e, this->fvGeom.subContVol[k].local,k);

                sNDat[globalIdx].heatCap = problem.soil().heatCap(this->fvGeom.subContVol[k].global,e, this->fvGeom.subContVol[k].local,k);

                // mark elements that were already visited
                sNDat[globalIdx].visited = true;
           }
        }

      return;
    }

    virtual void initiateStaticData (const Entity& e)
    {
        // get access to shape functions for P1 elements
        GeometryType gt = e.geometry().type();
        const typename LagrangeShapeFunctionSetContainer<DT,RT,dim>::value_type&
        sfs=LagrangeShapeFunctions<DT,RT,dim>::general(gt,1);

        // get local to global id map
        for (int k = 0; k < sfs.size(); k++) {
           const int globalIdx = this->vertexMapper.template map<dim>(e, sfs[k].entity());

           // if nodes are not already visited
           if (!sNDat[globalIdx].visited)
            {
                // ASSUME porosity defined at nodes
                sNDat[globalIdx].porosity = problem.soil().porosity(this->fvGeom.elementGlobal, e, this->fvGeom.elementLocal,k);

                // mark elements that were already visited
                sNDat[globalIdx].visited = true;
            }
        }

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
        RT satN;
     RT satW;
     RT pW;
     RT pC;
     RT pN;
     RT temperature;
     RT lambda;
     RT viscosityCO2;
     RT krCO2;
     RT permeabilityXDir;
     RT porosity;
     FieldVector<RT,2> mobility;  //Vector with the number of phases
     FieldVector<RT,2> density;
     FieldVector<RT,2> enthalpy;
     FieldVector<RT,2> intenergy;
    };

    // analog to EvalPrimaryData in MUFTE, uses members of vNDat
    virtual void updateVariableData(const Entity& e, const VBlockType* sol,
            int i, std::vector<VariableNodeData>& varData)
    {
	   	 const int global = this->vertexMapper.template map<dim>(e, i);
 	   	 GeometryType gt = e.geometry().type();
 	   	 const typename LagrangeShapeFunctionSetContainer<DT,RT,dim>::value_type&
 	   	 sfs=LagrangeShapeFunctions<DT,RT,dim>::general(gt,1);

 	  	 // get global and local coordinates of node i
 		 const FieldVector<DT,dim> globalPos = this->fvGeom.subContVol[i].global;
 		 const FieldVector<DT,dim> localPos = this->fvGeom.subContVol[i].local;


            varData[i].pW = sol[i][pWIdx];
            varData[i].satN = sol[i][satNIdx];
            varData[i].satW = 1.0 - varData[i].satN;
           varData[i].pC = 0.0;// problem.materialLaw().pC(varData[i].satW, this->fvGeom.subContVol[i].global, e, this->fvGeom.subContVol[i].local);
            varData[i].pN = varData[i].pW + varData[i].pC;
            varData[i].temperature = sol[i][teIdx]; // in [K]

            varData[i].density[wPhase] = problem.wettingPhase().density(varData[i].temperature, varData[i].pW);
            varData[i].density[nPhase] = problem.nonwettingPhase().density(varData[i].temperature, varData[i].pN);
            varData[i].mobility[wPhase] = problem.materialLaw().mobW(varData[i].satW, this->fvGeom.subContVol[i].global, e, this->fvGeom.subContVol[i].local, varData[i].temperature, varData[i].pW);
            varData[i].mobility[nPhase] = problem.materialLaw().mobN(varData[i].satN, this->fvGeom.subContVol[i].global, e, this->fvGeom.subContVol[i].local, varData[i].temperature, varData[i].pN);
            varData[i].lambda = problem.soil().heatCond(this->fvGeom.subContVol[i].global, e, this->fvGeom.subContVol[i].local, varData[i].satW, i);
         varData[i].enthalpy[wPhase] = problem.wettingPhase().enthalpy(varData[i].temperature,varData[i].pW);
         varData[i].enthalpy[nPhase] = problem.nonwettingPhase().enthalpy(varData[i].temperature,varData[i].pN);
         varData[i].intenergy[wPhase] = problem.wettingPhase().intEnergy(varData[i].temperature,varData[i].pW);
         varData[i].intenergy[nPhase] = problem.nonwettingPhase().intEnergy(varData[i].temperature,varData[i].pN);

         varData[i].permeabilityXDir = this->problem.soil().K(globalPos,e,localPos, i)[0][0];
         varData[i].porosity = this->problem.soil().porosity(globalPos,e,localPos,i);

         // CONSTANT solubility (for comparison with twophase)
//         varData[i].massfrac[nPhase][wPhase] = 0.0; varData[i].massfrac[wPhase][wPhase] = 1.0;
//         varData[i].massfrac[wPhase][nPhase] = 0.0; varData[i].massfrac[nPhase][nPhase] = 1.0;

         //std::cout << "wPhase in gasphase: " << varData[i].massfrac[wPhase][nPhase] << std::endl;
         //std::cout << "nPhase in waterphase: " << varData[i].massfrac[co2][wPhase] << std::endl;

            // for output
            (*outPressureN)[global] = varData[i].pN;
            (*outCapillaryP)[global] = varData[i].pC;
              (*outSaturationW)[global] = varData[i].satW;
               (*outSaturationN)[global] = varData[i].satN;
               (*outTemperature)[global] = varData[i].temperature;
               (*outDensityW)[global] = varData[i].density[wPhase];
               (*outDensityN)[global] = varData[i].density[nPhase];
               (*outMobilityW)[global] = varData[i].mobility[wPhase];
               (*outMobilityN)[global] = varData[i].mobility[nPhase];
               (*outPermeabilityXDir)[global] = varData[i].permeabilityXDir;
               (*outPorosity)[global] = varData[i].porosity;

               return;
    }

    void updateVariableData(const Entity& e, const VBlockType* sol, int i, bool old = false)
    {
        if (old) {
            updateVariableData(e, sol, i, oldVNDat);
        }
        else
            updateVariableData(e, sol, i, vNDat);
    }

    void updateVariableData(const Entity& e, const VBlockType* sol, bool old = false)
    {
        int size = this->fvGeom.numVertices;

        for (int i = 0; i < size; i++)
                updateVariableData(e, sol, i, old);
    }

    struct StaticNodeData
    {
        bool visited;
        RT cellVolume;
        RT porosity;
        RT heatCap;
    };

    struct ElementData {
//     RT heatCap;
        } elData;


    // parameters given in constructor
       TwoPhaseHeatProblem<G,RT>& problem;
    std::vector<StaticNodeData> sNDat;
    std::vector<VariableNodeData> vNDat;
    std::vector<VariableNodeData> oldVNDat;

    // for output files
    BlockVector<FieldVector<RT, 1> > *outPressureN;
    BlockVector<FieldVector<RT, 1> > *outCapillaryP;
    BlockVector<FieldVector<RT, 1> > *outSaturationN;
    BlockVector<FieldVector<RT, 1> > *outSaturationW;
    BlockVector<FieldVector<RT, 1> > *outTemperature;
    BlockVector<FieldVector<RT, 1> > *outDensityW;
    BlockVector<FieldVector<RT, 1> > *outDensityN;
    BlockVector<FieldVector<RT, 1> > *outMobilityW;
    BlockVector<FieldVector<RT, 1> > *outMobilityN;
    BlockVector<FieldVector<RT, 1> > *outPhaseState;
    BlockVector<FieldVector<RT, 1> > *outPermeabilityXDir;
    BlockVector<FieldVector<RT, 1> > *outPorosity;

  };

  /** @} */
}
#endif
