// $Id$

#ifndef DUNE_BOX2P2CNIJACOBIAN_HH
#define DUNE_BOX2P2CNIJACOBIAN_HH

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
 * @author Bernd Flemisch, Melanie Darcis, Klaus Mosthaf
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
    - Scalar    type used for return values
  */
  template<class Grid, class Scalar, class BoxFunction = LeafP1Function<Grid, Scalar, 3> >
  class Box2P2CNIJacobian
    : public BoxJacobian<Box2P2CNIJacobian<Grid,Scalar,BoxFunction>,Grid,Scalar,3,BoxFunction>
  {
    typedef typename Grid::Traits::template Codim<0>::Entity Element;
    typedef typename Element::Geometry Geometry;
    typedef Box2P2CNIJacobian<Grid,Scalar,BoxFunction> ThisType;
    typedef typename LocalJacobian<ThisType,Grid,Scalar,3>::VBlockType VBlockType;
    typedef typename LocalJacobian<ThisType,Grid,Scalar,3>::MBlockType MBlockType;

     enum {pWIdx = 0, switchIdx = 1, teIdx=2, numberOfComponents = 2};    // Solution vector index
    enum {wPhase = 0, nPhase = 1};                    // Phase index
    enum {gasPhase = 0, waterPhase = 1, bothPhases = 2};    // Phase state
    enum {water = 0, air = 1, heat =2};                // Component index

  public:
    // define the number of phases (numPhases) and components (numComp) of your system, this is used outside
    // to allocate the correct size of (dense) blocks with a FieldMatrix
    enum {dim=Grid::dimension};
    enum {numEq=3, numComp=2};
    enum {SIZE=LagrangeShapeFunctionSetContainer<Scalar,Scalar,dim>::maxsize};
    struct VariableNodeData;

    typedef FieldMatrix<Scalar,dim,dim> FMatrix;
    typedef FieldVector<Scalar,dim> FVector;

    //! Constructor
    Box2P2CNIJacobian (TwoPTwoCNIProblem<Grid,Scalar>& params,
                  bool levelBoundaryAsDirichlet_, const Grid& grid,
                  BoxFunction& sol,
                  bool procBoundaryAsDirichlet_=true)
    : BoxJacobian<ThisType,Grid,Scalar,numEq,BoxFunction>(levelBoundaryAsDirichlet_, grid, sol, procBoundaryAsDirichlet_),
      problem(params),
      sNDat(this->vertexMapper.size()),
      vNDat(SIZE),
      oldVNDat(SIZE)
      {
      this->analytic = false;
    }

    /** @brief compute time dependent term (storage), loop over nodes / subcontrol volumes
     *  @param element entity
     *  @param sol solution vector
     *  @param node local node id
     *  @return storage term
     */
    virtual VBlockType computeM (const Element& element, const VBlockType* sol,
            int idx, std::vector<VariableNodeData>& varData)
    {
         GeometryType gt = element.geometry().type();
         const typename LagrangeShapeFunctionSetContainer<Scalar,Scalar,dim>::value_type&
          sfs=LagrangeShapeFunctions<Scalar,Scalar,dim>::general(gt,1);

        int globalIdx = this->vertexMapper.template map<dim>(element, sfs[idx].entity());

        VBlockType result;
        Scalar satN = varData[idx].satN;
        Scalar satW = varData[idx].satW;

        // storage of component water
        result[water] =
             sNDat[globalIdx].porosity*(varData[idx].density[wPhase]*satW*varData[idx].massfrac[water][wPhase]
                            +varData[idx].density[nPhase]*satN*varData[idx].massfrac[water][nPhase]);
        // storage of component air
        result[air] =
             sNDat[globalIdx].porosity*(varData[idx].density[nPhase]*satN*varData[idx].massfrac[air][nPhase]
                           +varData[idx].density[wPhase]*satW*varData[idx].massfrac[air][wPhase]);

        // storage term of energy equation
        result[heat] = sNDat[globalIdx].porosity * (varData[idx].density[wPhase] * varData[idx].intenergy[wPhase] * satW
                    + varData[idx].density[nPhase] * varData[idx].intenergy[nPhase] * satN)
                    + elData.heatCap * varData[idx].temperature;
         // soil properties defined at the elements!!!

        //std::cout << result << " " << idx << std::endl;
        return result;
    };

    virtual VBlockType computeM (const Element& element, const VBlockType* sol, int idx, bool old = false)
     {
         if (old)
             return computeM(element, sol, idx, oldVNDat);
         else
             return computeM(element, sol, idx, vNDat);
     }

    /** @brief compute diffusive/advective fluxes, loop over subcontrol volume faces
     *  @param element entity
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
     const FieldVector<Scalar,dim> global_i = this->fvGeom.subContVol[i].global;
     const FieldVector<Scalar,dim> global_j = this->fvGeom.subContVol[j].global;
     const FieldVector<Scalar,dim> local_i = this->fvGeom.subContVol[i].local;
     const FieldVector<Scalar,dim> local_j = this->fvGeom.subContVol[j].local;

     FieldMatrix<Scalar,2,dim> pGrad(0.), xGrad(0.);
     FMatrix Ki, Kj;
     FieldVector<Scalar,dim> teGrad(0.);
     FieldVector<Scalar,dim> temp(0.);
//         FieldVector<Scalar,dim> Kij(0);
          VBlockType flux(0.);

         // effective permeability in edge direction
          // Scalar Kij = sIPDat[global_j].K_eff[face];

          // calculate harmonic mean of permeabilities of nodes i and j
          Ki = this->problem.soil().K(global_i,element,local_i);
          Kj = this->problem.soil().K(global_j,element,local_j);
         const FMatrix K = harmonicMeanK(Ki, Kj);
         //const FMatrix K = harmonicMeanK(element, face);
          //K.umv(normal, Kij);  // Kij=K*n

     // Harmonic mean:
        // Heat Conductivity
         Scalar lambda;
         lambda = 2./((1./vNDat[i].lambda) + (1./vNDat[j].lambda));

     // Arithmetic mean:
         // Enthalpy
         Scalar avgEnthW;
         Scalar avgEnthN;
         avgEnthW = (vNDat[i].enthalpy[wPhase] + vNDat[j].enthalpy[wPhase]) / 2.;
         avgEnthN = (vNDat[i].enthalpy[nPhase] + vNDat[j].enthalpy[nPhase]) / 2.;



         // calculate FE gradient (grad p for each phase)
         for (int k = 0; k < this->fvGeom.numVertices; k++) // loop over adjacent nodes
         {
             // FEGradient at idx k
             const FieldVector<Scalar,dim> feGrad(this->fvGeom.subContVolFace[face].grad[k]);
             FieldVector<Scalar,2> pressure(0.0), massfrac(0.0);

             pressure[wPhase] = vNDat[k].pW;
             pressure[nPhase] = vNDat[k].pN;

               // compute sum of pressure gradients for each phase
               for (int phase = 0; phase < 2; phase++)
               {
                   temp = feGrad;
                   temp *= pressure[phase];
                   pGrad[phase] += temp;
               }
               // Temperature gradient
               temp = feGrad;
              temp *= vNDat[k].temperature;
              teGrad += temp;

               // for diffusion of air in wetting phase
               temp = feGrad;
              temp *= vNDat[k].massfrac[air][wPhase];
              xGrad[wPhase] += temp;

               // for diffusion of water in nonwetting phase
              temp = feGrad;
              temp *= vNDat[k].massfrac[water][nPhase];
              xGrad[nPhase] += temp;
         }

       // deduce gravity*density of each phase
         FieldMatrix<Scalar,2,dim> contribComp(0);
         for (int phase=0; phase<2; phase++)
         {
             contribComp[phase] = problem.gravity();
             contribComp[phase] *= vNDat[i].density[phase];
             pGrad[phase] -= contribComp[phase]; // grad p - rho*g
         }

         VBlockType outward(0);  // Darcy velocity of each phase
         FieldVector<Scalar,dim> v_tilde(0);

         // calculate the advective flux using upwind: K*n(grad p -rho*g)
         for (int phase=0; phase<2; phase++)
             {
              K.umv(pGrad[phase], v_tilde);  // v_tilde=K*gradP
              outward[phase] = v_tilde*normal;
             }
             outward[heat] = teGrad * normal;
              outward[heat] *= lambda;  // schau mer mal

         // evaluate upwind nodes
         int up_w, dn_w, up_n, dn_n;
         if (outward[wPhase] <= 0) {up_w = i; dn_w = j;}
         else {up_w = j; dn_w = i;};
         if (outward[nPhase] <= 0) {up_n = i; dn_n = j;}
         else {up_n = j; dn_n = i;};

         Scalar alpha = 1.0;  // Upwind parameter

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
          // Air conservation
          flux[air]   =   (alpha* vNDat[up_n].density[nPhase]*vNDat[up_n].mobility[nPhase]
                                * vNDat[up_n].massfrac[air][nPhase]
                     + (1-alpha)* vNDat[dn_n].density[nPhase]*vNDat[dn_n].mobility[nPhase]
                                * vNDat[dn_n].massfrac[air][nPhase])
                                * outward[nPhase];
          flux[air]  +=   (alpha* vNDat[up_w].density[wPhase]*vNDat[up_w].mobility[wPhase]
                                * vNDat[up_w].massfrac[air][wPhase]
                     + (1-alpha)* vNDat[dn_w].density[wPhase]*vNDat[dn_w].mobility[wPhase]
                                * vNDat[dn_w].massfrac[air][wPhase])
                                * outward[wPhase];

          // Heat conservation
                // flux term of energy balance equation
          flux[heat]   =  (alpha* vNDat[up_n].density[nPhase]*vNDat[up_n].mobility[nPhase]
                     + (1-alpha)* vNDat[dn_n].density[nPhase]*vNDat[dn_n].mobility[nPhase])
                                * avgEnthN * outward[nPhase];
          flux[heat]  +=  (alpha* vNDat[up_w].density[wPhase]*vNDat[up_w].mobility[wPhase]
                     + (1-alpha)* vNDat[dn_w].density[wPhase]*vNDat[dn_w].mobility[wPhase])
                                * avgEnthW * outward[wPhase];
          flux[heat]    +=    outward[heat];


         // DIFFUSION
         VBlockType normDiffGrad;

            // get local to global id map
         int state_i = vNDat[i].phasestate;
         int state_j = vNDat[j].phasestate;

           Scalar diffusionWW(0.0), diffusionWN(0.0); // diffusion of water
         Scalar diffusionAW(0.0), diffusionAN(0.0); // diffusion of air
         VBlockType avgDensity, avgDpm;
         avgDpm[wPhase]=1e-9; // needs to be changed !!!
         avgDpm[nPhase]=1e-5; // water in the gasphase

         normDiffGrad[wPhase] = xGrad[wPhase]*normal;
         normDiffGrad[nPhase] = xGrad[nPhase]*normal;

         // calculate the arithmetic mean of densities
         avgDensity[wPhase] = 0.5*(vNDat[i].density[wPhase] + vNDat[j].density[wPhase]);
         avgDensity[nPhase] = 0.5*(vNDat[i].density[nPhase] + vNDat[j].density[nPhase]);


         if (state_i==2 && state_j==2)
         {
             diffusionAW = avgDpm[wPhase] * avgDensity[wPhase] * normDiffGrad[wPhase];
             diffusionWW = - diffusionAW;
             diffusionWN = avgDpm[nPhase] * avgDensity[nPhase] * normDiffGrad[nPhase];
             diffusionAN = - diffusionWN;
         }
         else if ((state_i == 1 || state_j == 1) || (state_i == 1 && state_j == 1))
         {
             diffusionAW = avgDpm[wPhase] * avgDensity[wPhase] * normDiffGrad[wPhase];
             diffusionWW = - diffusionAW;
         }
         else if ((state_i == 0 || state_j == 0) || (state_i == 0 && state_j == 0))
         {
             diffusionWN = avgDpm[nPhase] * avgDensity[nPhase] * normDiffGrad[nPhase];
             diffusionAN = - diffusionWN;
         }

          /////////////////////////////////////////////////////////////////////////////////////////////////////////
          // DIFFUSIVE TRANSPORT


         //TODO: diffusion must be checked again!
         // add water diffusion to flux
         flux[water] += (diffusionWW + diffusionWN);
//         std::cout << "Water Flux: " << flux[water] << std::endl;
         // air diffusion not implemeted
         flux[air] += (diffusionAN + diffusionAW);
//         std::cout << "Air Flux: " << flux[air] << std::endl;

         //TODO: average enthalpy correct here?
         // Heat conservation
         flux[heat] += diffusionWW*avgEnthW + diffusionAW*avgEnthN;

          //////////////////////////////////////////////////////////////////////////////////////////////
          // HEAT CONVECTION

          flux[heat]    += outward[heat];

         return flux;
  };

      /** @brief integrate sources / sinks
       *  @param element entity
     *  @param sol solution vector
     *  @param idx local node id
     *  @return source/sink term
     */
       virtual VBlockType computeQ (const Element& element, const VBlockType* sol, const int& idx)
       {
           // ASSUME problem.q already contains \rho.q
           return problem.q(this->fvGeom.subContVol[idx].global, element, this->fvGeom.subContVol[idx].local);
       }

      /** @brief perform variable switch
       *  @param global global vertex id
     *  @param sol solution vector
     *  @param local local vertex id
     */
       virtual void primaryVarSwitch (const Element& element, int global, VBlockType* sol, int local)
    {

        bool switched = false;
           int state = sNDat[global].phaseState;

        Scalar pW = sol[local][pWIdx];
        Scalar satW = 0.0;
        if (state == bothPhases) satW = 1.0-sol[local][switchIdx];
          if (state == waterPhase) satW = 1.0;
          if (state == gasPhase) satW = 0.0;

        Scalar pC = problem.materialLaw().pC(satW, this->fvGeom.subContVol[local].global, element,
                    this->fvGeom.subContVol[local].local);
          Scalar pN = pW + pC;
          Scalar temperature = sol[local][teIdx];

        FVector& Coordinates = this->fvGeom.subContVol[local].global;

        switch(state)
        {
        case gasPhase :
            Scalar xWNmass, xWNmolar, pwn, pWSat;
            xWNmass = sol[local][switchIdx];
            xWNmolar = problem.multicomp().convertMassToMoleFraction(xWNmass, gasPhase);
               pwn = xWNmolar * pN;
            pWSat = problem.multicomp().vaporPressure(temperature);

            if (pwn > 1.01*pWSat && switched == false)
            {
                // appearance of water phase
                std::cout << "Water appears at node " << global << "  Coordinates: " << Coordinates << std::endl;
                sNDat[global].phaseState = bothPhases;
                sol[local][switchIdx] = 1.0 - 1.e-6; // initialize solution vector
                switched = true;
            }
            break;

        case waterPhase :
            Scalar pbub, xAWmass, xAWmolar, henryInv;
               xAWmass = sol[local][switchIdx];
             xAWmolar = problem.multicomp().convertMassToMoleFraction(xAWmass, waterPhase);
            henryInv = problem.multicomp().henry(sol[local][teIdx]);
            pWSat = problem.multicomp().vaporPressure(sol[local][teIdx]);
            pbub = pWSat + xAWmolar/henryInv;

            if (pbub > pN && switched == false)
            {
                // appearance of gas phase
                std::cout << "Gas appears at node " << global << "  Coordinates: " << Coordinates << std::endl;
                sNDat[global].phaseState = bothPhases;
                sol[local][switchIdx] = 2e-5; // initialize solution vector
                switched = true;
            }
            break;

        case bothPhases :
            Scalar satN = sol[local][switchIdx];

            if (satN < -1e-4  && switched == false)
                {
                    // disappearance of gas phase
                    std::cout << "Gas disappears at node " << global << "  Coordinates: " << Coordinates << std::endl;
                    sNDat[global].phaseState = waterPhase;
                    sol[local][switchIdx] = 2e-5; // initialize solution vector
                    switched = true;
            }
            else if (satW < -1e-4  && switched == false)
                {
                    // disappearance of water phase
                    std::cout << "Water disappears at node " << global << "  Coordinates: " << Coordinates << std::endl;
                    sNDat[global].phaseState = gasPhase;
                    sol[local][switchIdx] = 2e-5; // initialize solution vector
                switched = true;
                }
                break;

        }
        if (switched == true) updateVariableData(element, sol, local, vNDat, sNDat[global].phaseState);

       return;
    }

    // harmonic mean computed directly
    virtual FMatrix harmonicMeanK (FMatrix& Ki, const FMatrix& Kj) const
    {
        double eps = 1e-20;
//        FMatrix K;

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


    virtual void clearVisited ()
    {
        for (int i = 0; i < this->vertexMapper.size(); i++){
           sNDat[i].visited = false;
//            sNDat[i].switched = false;
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

      //*********************************************************
      //*                                                        *
      //*    Calculation of Data at Elements (elData)             *
      //*                                                         *
      //*                                                        *
      //*********************************************************

    virtual void computeElementData (const Element& element)
    {


          elData.heatCap = problem.soil().heatCap(this->fvGeom.elementGlobal, element, this->fvGeom.elementLocal);


//       // ASSUME element-wise constant parameters for the material law
//          elData.parameters = problem.materialLawParameters
//          (this->fvGeom.elementGlobal, element, this->fvGeom.elementLocal);
//
//         // ASSUMING element-wise constant permeability, evaluate K at the cell center
//          elData.K = problem.K(this->fvGeom.elementGlobal, element, this->fvGeom.elementLocal);
//
//         // ASSUMING element-wise constant porosity
//          elData.porosity = problem.porosity(this->fvGeom.elementGlobal, element, this->fvGeom.elementLocal);
        return;
    }


      //*********************************************************
      //*                                                        *
      //*    Calculation of Data at Nodes that has to be            *
      //*    determined only once    (sNDat)                        *
      //*                                                        *
      //*********************************************************

    // analog to EvalStaticData in MUFTE
    virtual void updateStaticData (const Element& element, VBlockType* sol)
    {
        // size of the sNDat vector is determined in the constructor

        // local to global id mapping (do not ask vertex mapper repeatedly
        //int localToGlobal[LagrangeShapeFunctionSetContainer<Scalar,Scalar,n>::maxsize];

        // get access to shape functions for P1 elements
        GeometryType gt = element.geometry().type();
        const typename LagrangeShapeFunctionSetContainer<Scalar,Scalar,dim>::value_type&
        sfs=LagrangeShapeFunctions<Scalar,Scalar,dim>::general(gt,1);

        // get local to global id map
        for (int k = 0; k < sfs.size(); k++) {
           const int globalIdx = this->vertexMapper.template map<dim>(element, sfs[k].entity());

           // if nodes are not already visited
           if (!sNDat[globalIdx].visited)
            {
                // ASSUME porosity defined at nodes
                sNDat[globalIdx].porosity = problem.soil().porosity(this->fvGeom.elementGlobal, element, this->fvGeom.elementLocal);

               // global coordinates
               FieldVector<Scalar,dim> global_i = this->fvGeom.subContVol[k].global;

                 // evaluate primary variable switch
               primaryVarSwitch(element, globalIdx, sol, k);

                // mark elements that were already visited
                sNDat[globalIdx].visited = true;
            }
        }

      return;
    }

    // for initialization of the Static Data (sets porosity)
    virtual void initiateStaticData (const Element& element)
    {
        // get access to shape functions for P1 elements
        GeometryType gt = element.geometry().type();
        const typename LagrangeShapeFunctionSetContainer<Scalar,Scalar,dim>::value_type&
        sfs=LagrangeShapeFunctions<Scalar,Scalar,dim>::general(gt,1);

        // get local to global id map
        for (int k = 0; k < sfs.size(); k++) {
           const int globalIdx = this->vertexMapper.template map<dim>(element, sfs[k].entity());

           // if nodes are not already visited
           if (!sNDat[globalIdx].visited)
            {
                // ASSUME porosity defined at nodes
                sNDat[globalIdx].porosity = problem.soil().porosity(this->fvGeom.elementGlobal, element, this->fvGeom.elementLocal);

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
        Scalar satN;
     Scalar satW;
     Scalar pW;
     Scalar pC;
     Scalar pN;
     Scalar temperature;
     Scalar lambda;
     FieldVector<Scalar,2> mobility;  //Vector with the number of phases
     FieldVector<Scalar,2> density;
     FieldMatrix<Scalar,numComp,2> massfrac;
     FieldVector<Scalar,2> enthalpy;
     FieldVector<Scalar,2> intenergy;
     int phasestate;
    };

    // analog to EvalPrimaryData in MUFTE, uses members of vNDat
    virtual void updateVariableData(const Element& element, const VBlockType* sol,
            int i, std::vector<VariableNodeData>& varData, int state)
    {
               const int globalIdx = this->vertexMapper.template map<dim>(element, i);
               FVector& global = this->fvGeom.subContVol[i].global;
               FVector& local = this->fvGeom.subContVol[i].local;

            varData[i].pW = sol[i][pWIdx];
            if (state == bothPhases) varData[i].satN = sol[i][switchIdx];
            if (state == waterPhase) varData[i].satN = 0.0;
            if (state == gasPhase) varData[i].satN = 1.0;

            varData[i].satW = 1.0 - varData[i].satN;

            varData[i].pC = problem.materialLaw().pC(varData[i].satW, global, element, local);
            varData[i].pN = varData[i].pW + varData[i].pC;
            varData[i].temperature = sol[i][teIdx]; // in [K]

            // Solubilities of components in phases
            if (state == bothPhases){
                   varData[i].massfrac[air][wPhase] = problem.multicomp().xAW(varData[i].pN, varData[i].temperature);
                   varData[i].massfrac[water][nPhase] = problem.multicomp().xWN(varData[i].pN, varData[i].temperature);
            }
            if (state == waterPhase){
                   varData[i].massfrac[water][nPhase] = 0.0;
                   varData[i].massfrac[air][wPhase] =  sol[i][switchIdx];
            }
            if (state == gasPhase){
                   varData[i].massfrac[water][nPhase] = sol[i][switchIdx];
                   varData[i].massfrac[air][wPhase] = 0.0;
            }
               varData[i].massfrac[water][wPhase] = 1.0 - varData[i].massfrac[air][wPhase];
               varData[i].massfrac[air][nPhase] = 1.0 - varData[i].massfrac[water][nPhase];
               varData[i].phasestate = state;

            // Mobilities & densities
            varData[i].mobility[wPhase] = problem.materialLaw().mobW(varData[i].satW, global, element, local, varData[i].temperature, varData[i].pW);
            varData[i].mobility[nPhase] = problem.materialLaw().mobN(varData[i].satN, global, element, local, varData[i].temperature, varData[i].pN);
            varData[i].density[wPhase] = problem.wettingPhase().density(varData[i].temperature, varData[i].pN);
            varData[i].density[nPhase] = problem.nonwettingPhase().density(varData[i].temperature, varData[i].pN,
                    varData[i].massfrac[air][nPhase]);
         varData[i].lambda = problem.soil().heatCond(global, element, local, varData[i].satW);
         varData[i].enthalpy[pWIdx] = problem.wettingPhase().enthalpy(varData[i].temperature,varData[i].pW);
         varData[i].enthalpy[switchIdx] = problem.nonwettingPhase().enthalpy(varData[i].temperature,varData[i].pN);
         varData[i].intenergy[pWIdx] = problem.wettingPhase().intEnergy(varData[i].temperature,varData[i].pW);
         varData[i].intenergy[switchIdx] = problem.nonwettingPhase().intEnergy(varData[i].temperature,varData[i].pN);

         // CONSTANT solubility (for comparison with twophase)
//         varData[i].massfrac[air][wPhase] = 0.0; varData[i].massfrac[water][wPhase] = 1.0;
//         varData[i].massfrac[water][nPhase] = 0.0; varData[i].massfrac[air][nPhase] = 1.0;

         //std::cout << "water in gasphase: " << varData[i].massfrac[water][nPhase] << std::endl;
         //std::cout << "air in waterphase: " << varData[i].massfrac[air][wPhase] << std::endl;

            // for output
            (*outPressureN)[globalIdx] = varData[i].pN;
            (*outCapillaryP)[globalIdx] = varData[i].pC;
              (*outSaturationW)[globalIdx] = varData[i].satW;
               (*outSaturationN)[globalIdx] = varData[i].satN;
               (*outTemperature)[globalIdx] = varData[i].temperature;
               (*outMassFracAir)[globalIdx] = varData[i].massfrac[air][wPhase];
               (*outMassFracWater)[globalIdx] = varData[i].massfrac[water][nPhase];
               (*outDensityW)[globalIdx] = varData[i].density[wPhase];
               (*outDensityN)[globalIdx] = varData[i].density[nPhase];
               (*outMobilityW)[globalIdx] = varData[i].mobility[wPhase];
               (*outMobilityN)[globalIdx] = varData[i].mobility[nPhase];
               (*outPhaseState)[globalIdx] = varData[i].phasestate;

               return;
    }

    virtual void updateVariableData(const Element& element, const VBlockType* sol, int i, bool old = false)
    {
        int state;
        const int global = this->vertexMapper.template map<dim>(element, i);
        if (old)
        {
                  state = sNDat[global].oldPhaseState;
            updateVariableData(element, sol, i, oldVNDat, state);
        }
        else
        {
            state = sNDat[global].phaseState;
            updateVariableData(element, sol, i, vNDat, state);
        }
    }

    void updateVariableData(const Element& element, const VBlockType* sol, bool old = false)
    {
        int size = this->fvGeom.numVertices;

        for (int i = 0; i < size; i++)
                updateVariableData(element, sol, i, old);
    }


    struct StaticNodeData
    {
        bool visited;
//        bool switched;
        int phaseState;
        int oldPhaseState;
        Scalar cellVolume;
        Scalar porosity;
        FieldVector<Scalar, 4> parameters;
        FMatrix K;
    };

     struct StaticIPData
     {
         bool visited;
         FMatrix K;
     };


    struct ElementData {
     Scalar heatCap;

//        Scalar cellVolume;
//          Scalar porosity;
//        Scalar gravity;
//        FieldVector<Scalar, 4> parameters;
//        FieldMatrix<Scalar,dim,dim> K;
        } elData;


    // parameters given in constructor
       TwoPTwoCNIProblem<Grid,Scalar>& problem;
    CWaterAir multicomp;
    std::vector<StaticNodeData> sNDat;
    std::vector<StaticIPData> sIPDat;
    std::vector<VariableNodeData> vNDat;
    std::vector<VariableNodeData> oldVNDat;

    // for output files
    BlockVector<FieldVector<Scalar, 1> > *outPressureN;
    BlockVector<FieldVector<Scalar, 1> > *outCapillaryP;
    BlockVector<FieldVector<Scalar, 1> > *outSaturationN;
    BlockVector<FieldVector<Scalar, 1> > *outSaturationW;
    BlockVector<FieldVector<Scalar, 1> > *outTemperature;
    BlockVector<FieldVector<Scalar, 1> > *outMassFracAir;
    BlockVector<FieldVector<Scalar, 1> > *outMassFracWater;
    BlockVector<FieldVector<Scalar, 1> > *outDensityW;
    BlockVector<FieldVector<Scalar, 1> > *outDensityN;
    BlockVector<FieldVector<Scalar, 1> > *outMobilityW;
    BlockVector<FieldVector<Scalar, 1> > *outMobilityN;
    BlockVector<FieldVector<Scalar, 1> > *outPhaseState;

  };

}
#endif
