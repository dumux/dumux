// $Id: boxco2jacobian.hh 972 2009-01-12 10:15:57Z lauser $

#ifndef DUNE_BOXCO2JACOBIAN_HH
#define DUNE_BOXCO2JACOBIAN_HH

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
#include"dumux/2p2cni/2p2cniproblem.hh"
//#include "varswitch.hh"
#include"dumux/io/vtkmultiwriter.hh"

/**
 * @file
 * @brief  compute local jacobian matrix for box scheme for two-phase two-component flow equation
 * @author Bernd Flemisch, Klaus Mosthaf, Melanie Darcis
 */



namespace Dune
{
  /** @addtogroup DISC_Disc
   *
   * @{
   */
  /**
   * @brief compute local jacobian matrix for the boxfile for two-phase two-component flow equation
   *
   */


  //! Derived class for computing local jacobian matrices
  /*! A class for computing local jacobian matrix for the two-phase two-component flow equation

        div j = q; j = -K grad u; in Omega

        u = g on Gamma1; j*n = J on Gamma2.

    Uses box scheme with the Lagrange shape functions.
    It should work for all dimensions and element types.
    All the numbering is with respect to the reference element and the
    Lagrange shape functions

    Template parameters are:

    - Grid  a DUNE grid type
    - RT    type used for return values
  */
  template<class G, class RT, class BoxFunction = LeafP1FunctionExtended<G, RT, 3> >
  class BoxCO2Jacobian
    : public BoxJacobian<BoxCO2Jacobian<G,RT,BoxFunction>,G,RT,3,BoxFunction>
  {
    typedef typename G::ctype DT;
    typedef typename G::Traits::template Codim<0>::Entity Entity;
    typedef typename Entity::Geometry Geometry;
    typedef BoxCO2Jacobian<G,RT,BoxFunction> ThisType;
    typedef typename LocalJacobian<ThisType,G,RT,3>::VBlockType VBlockType;
    typedef typename LocalJacobian<ThisType,G,RT,3>::MBlockType MBlockType;

     enum {pWIdx = 0, switchIdx = 1, teIdx=2};    // Solution vector index
    enum {wPhase = 0, nPhase = 1};                    // Phase index
    enum {gasPhase = 0, waterPhase = 1, bothPhases = 2};    // Phase state
    enum {water = 0, co2 = 1, heat = 2};                // Component index

  public:
    // define the number of phases (m) and components (c) of your system, this is used outside
    // to allocate the correct size of (dense) blocks with a FieldMatrix
    enum {dim=G::dimension};
    enum {m=3, c=2};
    enum {SIZE=LagrangeShapeFunctionSetContainer<DT,RT,dim>::maxsize};
    struct VariableNodeData;

    typedef FieldMatrix<RT,dim,dim> FMatrix;
    typedef FieldVector<RT,dim> FVector;

    //! Constructor
    BoxCO2Jacobian (TwoPTwoCNIProblem<G,RT>& params,
                  bool levelBoundaryAsDirichlet_, const G& grid,
                  BoxFunction& sol,
                  bool procBoundaryAsDirichlet_=true)
    : BoxJacobian<ThisType,G,RT,m,BoxFunction>(levelBoundaryAsDirichlet_, grid, sol, procBoundaryAsDirichlet_),
      problem(params),
      sNDat(this->vertexMapper.size()), vNDat(SIZE), oldVNDat(SIZE), switched(false), switchBreak(false), switchedGlobal(false), primVarSet(false)
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

        // storage of component water
        result[water] =
             sNDat[globalIdx].porosity*(varData[node].density[wPhase]*satW*varData[node].massfrac[water][wPhase]
                            +varData[node].density[nPhase]*satN*varData[node].massfrac[water][nPhase]);
        // storage of component co2
        result[co2] =
             sNDat[globalIdx].porosity*(varData[node].density[nPhase]*satN*varData[node].massfrac[co2][nPhase]
                           +varData[node].density[wPhase]*satW*varData[node].massfrac[co2][wPhase]);

        // storage term of energy equation
        result[heat] = sNDat[globalIdx].porosity * (varData[node].density[wPhase] * varData[node].intenergy[wPhase] * satW
                    + varData[node].density[nPhase] * varData[node].intenergy[nPhase] * satN)
                    + sNDat[globalIdx].heatCap * varData[node].temperature;
         // soil properties defined at the elements!!!

        //DEBUG

        for(int j=0; j<3; j++)
        {
            if(isinf(result[j]))
            {
                std::cout<<"INF in ComputeM \n" << "Coordinates:X="<< this->fvGeom.subContVol[node].global[0] <<" Y="<< this->fvGeom.subContVol[node].global[1] <<
                " Z="<<  this->fvGeom.subContVol[node].global[2]<<"\n"
                <<"water pressure: " << varData[node].pW << "water saturation: "<< satW << "temperature: "<< varData[node].temperature << std::endl;
            }
        }

        return result;
    };

    virtual VBlockType computeM (const Entity& e, const VBlockType* sol, int node, bool old = false)
     {
         if (old)
             return computeM(e, sol, node, oldVNDat);
         else
             return computeM(e, sol, node, vNDat);
     }

    /** @brief compute diffusive/advective fluxes, loop over subcontrol volume faces
     *  @param e entity
     *  @param sol solution vector
     *  @param face face id
     *  @return flux term
     */
     virtual VBlockType computeA (const Entity& e, const VBlockType* sol, int face)
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

         // Arithmetic Mean:

         // calculate tortuosity at the nodes i and j needed for porous media diffusion coefficient
             RT tauW_i, tauW_j, tauN_i, tauN_j; // tortuosity of wetting and nonwetting phase
             tauW_i = pow(sNDat[globalIdx_i].porosity * vNDat[i].satW,(7/3))/
                 (sNDat[globalIdx_i].porosity*sNDat[globalIdx_i].porosity);
             tauW_j = pow(sNDat[globalIdx_j].porosity * vNDat[j].satW,(7/3))/
                  (sNDat[globalIdx_j].porosity*sNDat[globalIdx_j].porosity);
             tauN_i = pow(sNDat[globalIdx_i].porosity * vNDat[i].satN,(7/3))/
                 (sNDat[globalIdx_i].porosity*sNDat[globalIdx_i].porosity);
             tauN_j = pow(sNDat[globalIdx_j].porosity * vNDat[j].satN,(7/3))/
                 (sNDat[globalIdx_j].porosity*sNDat[globalIdx_j].porosity);

             // arithmetic mean of porous media diffusion coefficient
              RT Dwg, Daw;
              Dwg = 0.5*(sNDat[globalIdx_i].porosity * vNDat[i].satN * tauN_i * vNDat[i].diff[nPhase] +
                  sNDat[globalIdx_j].porosity * vNDat[j].satN * tauN_j * vNDat[j].diff[nPhase]);
              Daw = 0.5*(sNDat[globalIdx_i].porosity * vNDat[i].satW * tauW_i * vNDat[i].diff[wPhase] +
                 sNDat[globalIdx_j].porosity * vNDat[j].satW * tauW_j * vNDat[j].diff[wPhase]);

             // arithmetic mean of phase enthalpies (dissolved components neglected)
              RT enthW;
              RT enthCO2;
              enthW = 0.5*(vNDat[i].enthalpy[pWIdx] + vNDat[j].enthalpy[pWIdx]);
              enthCO2 = 0.5*(vNDat[i].enthalpy[switchIdx] + vNDat[j].enthalpy[switchIdx]);

            // Calculate arithmetic mean of the densities
             VBlockType avgDensity;
            avgDensity[wPhase] = 0.5*(vNDat[i].density[wPhase] + vNDat[j].density[wPhase]);
            avgDensity[nPhase] = 0.5*(vNDat[i].density[nPhase] + vNDat[j].density[nPhase]);

    //////////////////////////////////////////////////////////////////////////////////////////
    // GRADIENTS

     FieldMatrix<RT,2,dim> pGrad(0.), xGrad(0.);
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

              // compute concentration gradients
              // for diffusion of co2 in wetting phase
               temp = feGrad;
              temp *= vNDat[k].massfrac[co2][wPhase];
              xGrad[wPhase] += temp;

               // for diffusion of water in nonwetting phase
              temp = feGrad;
              temp *= vNDat[k].massfrac[water][nPhase];
              xGrad[nPhase] += temp;
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
         FieldVector<RT,dim> v_tilde_w(0.);
         FieldVector<RT,dim> v_tilde_n(0.);

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
          flux[water] =   (alpha* vNDat[up_w].density[wPhase]*vNDat[up_w].mobility[wPhase]
                                * vNDat[up_w].massfrac[water][wPhase]
                     + (1-alpha)* vNDat[dn_w].density[wPhase]*vNDat[dn_w].mobility[wPhase]
                                * vNDat[dn_w].massfrac[water][wPhase])
                                * outward[wPhase];
          flux[water] +=  (alpha* vNDat[up_n].density[nPhase]*vNDat[up_n].mobility[nPhase]
                                * vNDat[up_n].massfrac[water][nPhase]
                     + (1-alpha)* vNDat[dn_n].density[nPhase]*vNDat[dn_n].mobility[nPhase]
                                * vNDat[dn_n].massfrac[water][nPhase])
                                * outward[nPhase];
          // CO2 conservation
          flux[co2]   =   (alpha* vNDat[up_n].density[nPhase]*vNDat[up_n].mobility[nPhase]
                                * vNDat[up_n].massfrac[co2][nPhase]
                     + (1-alpha)* vNDat[dn_n].density[nPhase]*vNDat[dn_n].mobility[nPhase]
                                * vNDat[dn_n].massfrac[co2][nPhase])
                                * outward[nPhase];
          flux[co2]  +=   (alpha* vNDat[up_w].density[wPhase]*vNDat[up_w].mobility[wPhase]
                                * vNDat[up_w].massfrac[co2][wPhase]
                     + (1-alpha)* vNDat[dn_w].density[wPhase]*vNDat[dn_w].mobility[wPhase]
                                * vNDat[dn_w].massfrac[co2][wPhase])
                                * outward[wPhase];
          // Heat conservation
         flux[heat]  =  (alpha* vNDat[up_n].density[nPhase]*vNDat[up_n].mobility[nPhase] * vNDat[up_n].enthalpy[nPhase]
                     + (1-alpha)* vNDat[dn_n].density[nPhase]*vNDat[dn_n].mobility[nPhase] * vNDat[dn_n].enthalpy[nPhase])
                     * outward[nPhase];
          flux[heat]  +=  (alpha* vNDat[up_w].density[wPhase]*vNDat[up_w].mobility[wPhase] * vNDat[up_w].enthalpy[wPhase]
                     + (1-alpha)* vNDat[dn_w].density[wPhase]*vNDat[dn_w].mobility[wPhase] * vNDat[dn_w].enthalpy[wPhase])
                     * outward[wPhase];

      /////////////////////////////////////////////////////////////////////////////////////////////////////////
      // DIFFUSIVE TRANSPORT

          VBlockType normDiffGrad;

           // get local to global id map
         int state_i = sNDat[globalIdx_i].phaseState;
         int state_j = sNDat[globalIdx_j].phaseState;

         RT diffusionAW(0.0), diffusionWW(0.0), diffusionWN(0.0), diffusionAN(0.0);

         normDiffGrad[wPhase] = xGrad[wPhase]*normal;
         normDiffGrad[nPhase] = xGrad[nPhase]*normal;

         if (state_i == bothPhases && state_j == bothPhases)
         {
             diffusionAW = Daw * avgDensity[wPhase] * normDiffGrad[wPhase];
             diffusionWW = - diffusionAW;
             diffusionWN = Dwg * avgDensity[nPhase] * normDiffGrad[nPhase];
             diffusionAN = - diffusionWN;
         }
         else if (state_i == waterPhase || state_j == waterPhase)
         {
             diffusionAW = Daw * avgDensity[wPhase] * normDiffGrad[wPhase];
             diffusionWW = - diffusionAW;
         }
         else if (state_i == gasPhase || state_j == gasPhase)
         {
             diffusionWN = Dwg * avgDensity[nPhase] * normDiffGrad[nPhase];
             diffusionAN = - diffusionWN;
         }

         // Water conservation
         flux[water] += (diffusionWW + diffusionWN);

         // CO2 conservation
         flux[co2] += (diffusionAN + diffusionAW);

//         // Heat conservation
//          flux[heat] += diffusionWW*enthW + diffusionAW*enthCO2; diffusive transport of enthalpy neglected here

      //////////////////////////////////////////////////////////////////////////////////////////////
      // HEAT CONDUCTION

          flux[heat]    +=    outward[heat];

//    BEGIN DEBUG
//          for(int j=0; j<3; j++)
//             {
//                 if(isinf(flux[j]))
//                 {
//                     std::cout<<"INF in ComputeA \n" << "Coordinates upw_node:X="<< this->fvGeom.subContVol[up_w].global[0] <<" Y="<< this->fvGeom.subContVol[up_w].global[1] <<
//                     " Z="<<  this->fvGeom.subContVol[up_w].global[2]<<"\n"
//                     <<"water pressure: " << vNDat[up_w].pW << "water saturation: "<< vNDat[up_w].satW
//                     << " \n temperature: "<< vNDat[up_w].temperature << "Xaw" << vNDat[up_w].massfrac[co2][wPhase] << std::endl;
//                 }
//             }
//     END DEBUG
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

      /** @brief perform variable switch
       *  @param global global node id
     *  @param sol solution vector
     *  @param local local node id
     */
       virtual void primaryVarSwitch (const Entity& e, int global, VBlockType* sol, int local)
    {

        if(!switchBreak)
         {
           switched = false;
           int state = sNDat[global].phaseState;

        RT pW = sol[local][pWIdx];
        RT satW = 0.0;
        if (state == bothPhases) satW = 1.0-sol[local][switchIdx];
          if (state == waterPhase) satW = 1.0;
          if (state == gasPhase) satW = 0.0;

        RT pC = problem.materialLaw().pC(satW, this->fvGeom.subContVol[local].global, e, this->fvGeom.subContVol[local].local);
          RT pN = pW + pC;

        FVector Coordinates = this->fvGeom.subContVol[local].global;

        switch(state)
        {
        case gasPhase :
            RT xWNmass;
            xWNmass = sol[local][switchIdx];

            if (xWNmass > 0.001 && switched == false)
            {
                // appearance of water phase
                std::cout << "Water appears at node " << global << "  Coordinates: " << Coordinates << std::endl;
                sNDat[global].phaseState = bothPhases;
                sol[local][switchIdx] = 1.0 - 1.e-3; // initialize solution vector
                switched = true;
            }
            break;

        case waterPhase :
            RT xAWmax, xAWmass;
             xAWmass = sol[local][switchIdx];
             xAWmax = problem.multicomp().xAW(pN, sol[local][teIdx]);

            if (xAWmass > xAWmax && switched == false)
            {
                // appearance of gas phase
                std::cout << "Gas appears at node " << global << "  Coordinates: " << Coordinates << std::endl;
                sNDat[global].phaseState = bothPhases;
                sol[local][switchIdx] = 1.e-3; // initialize solution vector
                switched = true;
            }
            break;

        case bothPhases :
            RT satN = sol[local][switchIdx];
                if (satN < 0.0  && switched == false)
                {
                    // disappearance of gas phase
                    std::cout << "Gas disappears at node " << global << "  Coordinates: " << Coordinates << std::endl;
                    sNDat[global].phaseState = waterPhase;
                    sol[local][switchIdx] = 1e-3; // initialize solution vector
                    switched = true;
            }
                else if (satW < -1.e-5  && switched == false)
                {
                    // disappearance of water phase
                    std::cout << "Water disappears at node " << global << "  Coordinates: " << Coordinates << std::endl;
                    sNDat[global].phaseState = gasPhase;
                    sol[local][switchIdx] = 1e-8; // initialize solution vector
                switched = true;
                }
                break;
        }

        if(switched)
         {
                 // fill global solution vector with initialised local values
                 BoxJacobian<ThisType,G,RT,m,BoxFunction>::localToGlobal(e, sol);
            switchedGlobal = true;
         }
         }
       return;
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

    // updates old phase state after each time step
    virtual void updatePhaseState ()
    {
          for (int i = 0; i < this->vertexMapper.size(); i++){
               sNDat[i].oldPhaseState = sNDat[i].phaseState;
           }
       return;
    }

    virtual void resetPhaseState ()
    {
          for (int i = 0; i < this->vertexMapper.size(); i++){
               sNDat[i].phaseState = sNDat[i].oldPhaseState;
           }
       return;
    }
      //*********************************************************
      //*                                                        *
      //*    Calculation of Data at Elements (elData)             *
      //*                                                         *
      //*                                                        *
      //*********************************************************

    virtual void computeElementData (const Entity& e)
    {


//       // ASSUME element-wise constant parameters for the material law
//          elData.parameters = problem.materialLawParameters
//          (this->fvGeom.cellGlobal, e, this->fvGeom.cellLocal);
//
//         // ASSUMING element-wise constant permeability, evaluate K at the cell center
//          elData.K = problem.K(this->fvGeom.cellGlobal, e, this->fvGeom.cellLocal);
//
//         // ASSUMING element-wise constant porosity
//          elData.porosity = problem.porosity(this->fvGeom.cellGlobal, e, this->fvGeom.cellLocal);
        return;
    }


      //*********************************************************
      //*                                                        *
      //*    Calculation of Data at Nodes that has to be            *
      //*    determined only once    (sNDat)                        *
      //*                                                        *
      //*********************************************************

    // analog to EvalStaticData in MUFTE
    virtual void updateStaticData (const Entity& e, VBlockType* sol)
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

                sNDat[globalIdx].heatCap = problem.soil().heatCap(this->fvGeom.elementGlobal,e, this->fvGeom.elementLocal,k);

               // global coordinates
               FieldVector<DT,dim> global_i = this->fvGeom.subContVol[k].global;

                 // evaluate primary variable switch
               primaryVarSwitch(e, globalIdx, sol, k);

                // mark elements that were already visited
                sNDat[globalIdx].visited = true;
           }
        }

      return;
    }


   // for initialization of the Static Data (sets porosity)
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
      //*    Calculation of variable Data at Nodes                *
      //*    (vNDat)                                                *
      //*                                                        *
      //*********************************************************

    // for output files
//    BlockVector<FieldVector<RT, 1> > *hackyMassFracCO2;
//    BlockVector<FieldVector<RT, 1> > *hackyMassFracWater;
//    BlockVector<FieldVector<RT, 1> > *hackySaturationN;

//    void printVariableData()
//    {
//        for (int i = 0; i < 4; i++)
//        {
//            std::cout << "new: i = " << i << ": satN = " << vNDat[i].satN << ", satW = " << vNDat[i].satW
//                << ", pW = " << vNDat[i].pW << ", pC = " << vNDat[i].pC << ", pN = " << vNDat[i].pN
//                << ", T = " << vNDat[i].temperature << ", lambda = " << vNDat[i].lambda << std::endl;
//            std::cout << "old: i = " << i << ": satN = " << oldVNDat[i].satN << ", satW = " << oldVNDat[i].satW
//                << ", pW = " << oldVNDat[i].pW << ", pC = " << oldVNDat[i].pC << ", pN = " << oldVNDat[i].pN
//                << ", T = " << oldVNDat[i].temperature << ", lambda = " << oldVNDat[i].lambda << std::endl;
//        }
//    }



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
     FieldMatrix<RT,c,2> massfrac;
     FieldVector<RT,2> enthalpy;
     FieldVector<RT,2> intenergy;
     FieldVector<RT,2> diff;
     int phasestate;
    };

    // analog to EvalPrimaryData in MUFTE, uses members of vNDat
    virtual void updateVariableData(const Entity& e, const VBlockType* sol,
            int i, std::vector<VariableNodeData>& varData, int state)
    {
         const int global = this->vertexMapper.template map<dim>(e, i);
         GeometryType gt = e.geometry().type();
         const typename LagrangeShapeFunctionSetContainer<DT,RT,dim>::value_type&
         sfs=LagrangeShapeFunctions<DT,RT,dim>::general(gt,1);

         // get global and local coordinates of node i
         const FieldVector<DT,dim> globalPos = this->fvGeom.subContVol[i].global;
         const FieldVector<DT,dim> localPos = this->fvGeom.subContVol[i].local;

            varData[i].pW = sol[i][pWIdx];
            if (state == bothPhases) varData[i].satN = sol[i][switchIdx];
            if (state == waterPhase) varData[i].satN = 0.0;
            if (state == gasPhase) varData[i].satN = 1.0;

            varData[i].satW = 1.0 - varData[i].satN;


           varData[i].pC = problem.materialLaw().pC(varData[i].satW, this->fvGeom.subContVol[i].global, e, this->fvGeom.subContVol[i].local);
            varData[i].pN = varData[i].pW + varData[i].pC;
            varData[i].temperature = sol[i][teIdx]; // in [K]

            // Solubilities of components in phases
            if (state == bothPhases){
                       varData[i].massfrac[co2][wPhase] = problem.multicomp().xAW(varData[i].pN, varData[i].temperature);
                          varData[i].massfrac[water][nPhase] = problem.multicomp().xWN(varData[i].pN, varData[i].temperature);
            }
            if (state == waterPhase){
                   varData[i].massfrac[water][nPhase] = 0.0;
                   varData[i].massfrac[co2][wPhase] =  sol[i][switchIdx];
            }
            if (state == gasPhase){
                   varData[i].massfrac[water][nPhase] = sol[i][switchIdx];
                   varData[i].massfrac[co2][wPhase] = 0.0;
            }
               varData[i].massfrac[water][wPhase] = 1.0 - varData[i].massfrac[co2][wPhase];
               varData[i].massfrac[co2][nPhase] = 1.0 - varData[i].massfrac[water][nPhase];
               // for output
//               (*hackySaturationN)[global] = varData[i].satN;
//               (*hackyMassFracCO2)[global] = varData[i].massfrac[co2][wPhase];
//               (*hackyMassFracWater)[global] = varData[i].massfrac[water][nPhase];
               // Diffusion coefficients
            // Mobilities & densities

           varData[i].density[wPhase] = problem.wettingPhase().density(varData[i].temperature, varData[i].pW, varData[i].massfrac[co2][wPhase]);
        varData[i].density[nPhase] = problem.nonwettingPhase().density(varData[i].temperature, varData[i].pN);
        varData[i].mobility[wPhase] = problem.materialLaw().mobW(varData[i].satW, this->fvGeom.subContVol[i].global, e, this->fvGeom.subContVol[i].local, varData[i].temperature, varData[i].pW);
        varData[i].mobility[nPhase] = problem.materialLaw().mobN(varData[i].satN, this->fvGeom.subContVol[i].global, e, this->fvGeom.subContVol[i].local, varData[i].temperature, varData[i].pN);
        varData[i].lambda = problem.soil().heatCond(this->fvGeom.subContVol[i].global, e, this->fvGeom.subContVol[i].local, varData[i].satW, i);
        varData[i].enthalpy[wPhase] = problem.wettingPhase().enthalpy(varData[i].temperature,varData[i].pW);
        varData[i].enthalpy[nPhase] = problem.nonwettingPhase().enthalpy(varData[i].temperature,varData[i].pN);
        varData[i].intenergy[wPhase] = problem.wettingPhase().intEnergy(varData[i].temperature,varData[i].pW);
        varData[i].intenergy[nPhase] = problem.nonwettingPhase().intEnergy(varData[i].temperature,varData[i].pN);
        varData[i].diff[wPhase] = problem.wettingPhase().diffCoeff();
        varData[i].diff[nPhase] = problem.nonwettingPhase().diffCoeff();

         varData[i].permeabilityXDir = this->problem.soil().K(globalPos,e,localPos, i)[0][0];
         varData[i].porosity = this->problem.soil().porosity(globalPos,e,localPos,i);

         // CONSTANT solubility (for comparison with twophase)
//         varData[i].massfrac[co2][wPhase] = 0.0; varData[i].massfrac[water][wPhase] = 1.0;
//         varData[i].massfrac[water][nPhase] = 0.0; varData[i].massfrac[co2][nPhase] = 1.0;

         //std::cout << "water in gasphase: " << varData[i].massfrac[water][nPhase] << std::endl;
         //std::cout << "co2 in waterphase: " << varData[i].massfrac[co2][wPhase] << std::endl;

            // for output
            (*outPressureN)[global] = varData[i].pN;
            (*outCapillaryP)[global] = varData[i].pC;
              (*outSaturationW)[global] = varData[i].satW;
               (*outSaturationN)[global] = varData[i].satN;
               (*outTemperature)[global] = varData[i].temperature;
               (*outMassFracAir)[global] = varData[i].massfrac[co2][wPhase];
               (*outMassFracWater)[global] = varData[i].massfrac[water][nPhase];
               (*outDensityW)[global] = varData[i].density[wPhase];
               (*outDensityN)[global] = varData[i].density[nPhase];
               (*outMobilityW)[global] = varData[i].mobility[wPhase];
               (*outMobilityN)[global] = varData[i].mobility[nPhase];
               (*outPermeabilityXDir)[global] = varData[i].permeabilityXDir;
               (*outPorosity)[global] = varData[i].porosity;
               (*outPhaseState)[global] = state;

               return;
    }

    virtual void updateVariableData(const Entity& e, const VBlockType* sol, int i, bool old = false)
    {
        int state;
        const int global = this->vertexMapper.template map<dim>(e, i);
        if (old)
        {
                  state = sNDat[global].oldPhaseState;
            updateVariableData(e, sol, i, oldVNDat, state);
        }
        else
        {
            state = sNDat[global].phaseState;
            updateVariableData(e, sol, i, vNDat, state);
        }
    }

    void updateVariableData(const Entity& e, const VBlockType* sol, bool old = false)
    {
        int size = this->fvGeom.numVertices;

        for (int i = 0; i < size; i++)
                updateVariableData(e, sol, i, old);
    }


    struct StaticNodeData
    {
        bool visited, primVarSet;
        int phaseState;
        int oldPhaseState;
        RT cellVolume;
        RT porosity;
        RT heatCap;
        FieldVector<RT, 4> parameters;
        FMatrix K;
    };

     struct StaticIPData
     {
         bool visited;
         FMatrix K;
     };


    struct ElementData {
//     RT heatCap;

//        RT cellVolume;
//          RT porosity;
//        RT gravity;
//        FieldVector<RT, 4> parameters;
//        FieldMatrix<RT,dim,dim> K;
        } elData;


    // parameters given in constructor
       TwoPTwoCNIProblem<G,RT>& problem;
    CBrineCO2 multicomp;
    std::vector<StaticNodeData> sNDat;
    std::vector<StaticIPData> sIPDat;
    std::vector<VariableNodeData> vNDat;
    std::vector<VariableNodeData> oldVNDat;
    bool switched, switchBreak, switchedGlobal, primVarSet;

    // for output files
    BlockVector<FieldVector<RT, 1> > *outPressureN;
    BlockVector<FieldVector<RT, 1> > *outCapillaryP;
    BlockVector<FieldVector<RT, 1> > *outSaturationN;
    BlockVector<FieldVector<RT, 1> > *outSaturationW;
    BlockVector<FieldVector<RT, 1> > *outTemperature;
    BlockVector<FieldVector<RT, 1> > *outMassFracAir;
    BlockVector<FieldVector<RT, 1> > *outMassFracWater;
    BlockVector<FieldVector<RT, 1> > *outDensityW;
    BlockVector<FieldVector<RT, 1> > *outDensityN;
    BlockVector<FieldVector<RT, 1> > *outMobilityW;
    BlockVector<FieldVector<RT, 1> > *outMobilityN;
    BlockVector<FieldVector<RT, 1> > *outPhaseState;
    BlockVector<FieldVector<RT, 1> > *outPermeabilityXDir;
    BlockVector<FieldVector<RT, 1> > *outPorosity;

  };

}
#endif

//    // average permeability from the staticNode Vector
//    virtual FMatrix harmonicMeanK (const Entity& e, int k) const
//    {
//     FMatrix Ki, Kj;
//     const RT eps = 1e-20;
//
//     int i = this->fvGeom.subContVolFace[k].i;
//     int j = this->fvGeom.subContVolFace[k].j;
//
//     int global_i = this->vertexMapper.template map<dim>(e, i);
//     int global_j = this->vertexMapper.template map<dim>(e, j);
//
//
//         Ki = sNDat[global_i].K;
//         Kj = sNDat[global_j].K;
//
//         for (int kx=0; kx<dim; kx++){
//             for (int ky=0; ky<dim; ky++){
//                 if (Ki[kx][ky] != Kj[kx][ky])
//                 {
//                     Ki[kx][ky] = 2 / (1/(Ki[kx][ky]+eps) + (1/(Kj[kx][ky]+eps)));
//                 }
//             }
//         }
//         return Ki;
//    }

