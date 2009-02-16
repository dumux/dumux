// $Id$

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
#include<dune/grid/utility/intersectiongetter.hh>
#include<dune/disc/shapefunctions/lagrangeshapefunctions.hh>
#include<dune/disc/operators/boundaryconditions.hh>
#include"dumux/operators/boxjacobian.hh"
#include"dumux/2p2c/2p2cproblem.hh"
#include"dumux/io/vtkmultiwriter.hh"

/**
 * @file
 * @brief  compute local jacobian matrix for box scheme for two-phase two-component flow equation
 * @author Bernd Flemisch, Klaus Mosthaf
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
  template<class Grid, class Scalar, class BoxFunction = LeafP1Function<Grid, Scalar, 2> >
  class Box2P2CJacobian
    : public BoxJacobian<Box2P2CJacobian<Grid,Scalar,BoxFunction>,Grid,Scalar,2,BoxFunction>
  {
    typedef typename Grid::Traits::template Codim<0>::Entity Element;
    typedef typename Element::Geometry Geometry;
    typedef Box2P2CJacobian<Grid,Scalar,BoxFunction> ThisType;
    typedef typename LocalJacobian<ThisType,Grid,Scalar,2>::VBlockType SolutionVector;
    typedef typename LocalJacobian<ThisType,Grid,Scalar,2>::MBlockType MBlockType;

     enum {pWIdx = 0, switchIdx = 1, numberOfComponents = 2};    // Solution vector index
    enum {wPhase = 0, nPhase = 1};                                    // Phase index
    enum {gasPhase = 0, waterPhase = 1, bothPhases = 2};        // Phase state
    enum {water = 0, air = 1};                                        // Component index

  public:
    // define the number of phases (numEq) and components (numComp) of your system, this is used outside
    // to allocate the correct size of (dense) blocks with a FieldMatrix
    enum {dim=Grid::dimension};
    enum {numEq=2, numComp=2};
    enum {SIZE=LagrangeShapeFunctionSetContainer<Scalar,Scalar,dim>::maxsize};
    struct VariableNodeData;

    typedef FieldMatrix<Scalar,dim,dim> FMatrix;
    typedef FieldVector<Scalar,dim> FVector;

    //! Constructor
    Box2P2CJacobian (TwoPTwoCProblem<Grid,Scalar>& params, bool levelBoundaryAsDirichlet_, const Grid& grid,
                  BoxFunction& sol, bool procBoundaryAsDirichlet_=true)
    : BoxJacobian<ThisType,Grid,Scalar,2,BoxFunction>(levelBoundaryAsDirichlet_, grid, sol, procBoundaryAsDirichlet_),
      problem(params), sNDat(this->vertexMapper.size()), vNDat(SIZE), oldVNDat(SIZE), switchFlag(false)
    {
      this->analytic = false;
      switchFlag = false;
      switchFlagLocal = false;
      temperature = 283.15;
    }

    /** @brief compute time dependent term (storage), loop over nodes / subcontrol volumes
     *  @param e entity
     *  @param sol solution vector
     *  @param node local node id
     *  @return storage term
     */
    virtual SolutionVector computeM (const Element& e, const SolutionVector* sol,
            int node, std::vector<VariableNodeData>& varData)
    {
         GeometryType gt = e.geometry().type();
         const typename LagrangeShapeFunctionSetContainer<Scalar,Scalar,dim>::value_type&
         sfs=LagrangeShapeFunctions<Scalar,Scalar,dim>::general(gt,1);

        int globalIdx = this->vertexMapper.template map<dim>(e, sfs[node].entity());

        SolutionVector result;
        Scalar satN = varData[node].satN;
        Scalar satW = varData[node].satW;

        // storage of component water
        result[water] =
            sNDat[globalIdx].porosity*(varData[node].density[wPhase]*satW*varData[node].massfrac[water][wPhase]
                           +varData[node].density[nPhase]*satN*varData[node].massfrac[water][nPhase]);
        // storage of component air
        result[air] =
            sNDat[globalIdx].porosity*(varData[node].density[nPhase]*satN*varData[node].massfrac[air][nPhase]
                           +varData[node].density[wPhase]*satW*varData[node].massfrac[air][wPhase]);

        //std::cout << result << " " << node << std::endl;
        return result;
    };

    virtual SolutionVector computeM (const Element& e, const SolutionVector* sol, int node, bool old = false)
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
    virtual SolutionVector computeA (const Element& e, const SolutionVector* sol, int face)
    {
      int i = this->fvGeom.subContVolFace[face].i;
      int j = this->fvGeom.subContVolFace[face].j;

      // normal vector, value of the area of the scvf
      const FieldVector<Scalar,dim> normal(this->fvGeom.subContVolFace[face].normal);
      GeometryType gt = e.geometry().type();
      const typename LagrangeShapeFunctionSetContainer<Scalar,Scalar,dim>::value_type&
      sfs=LagrangeShapeFunctions<Scalar,Scalar,dim>::general(gt,1);

      // global index of the subcontrolvolume face neighbor nodes in element e
      int globalIdx_i = this->vertexMapper.template map<dim>(e, sfs[i].entity());
      int globalIdx_j = this->vertexMapper.template map<dim>(e, sfs[j].entity());

       // get global coordinates of nodes i,j
     const FieldVector<Scalar,dim> global_i = this->fvGeom.subContVol[i].global;
     const FieldVector<Scalar,dim> global_j = this->fvGeom.subContVol[j].global;
     const FieldVector<Scalar,dim> local_i = this->fvGeom.subContVol[i].local;
     const FieldVector<Scalar,dim> local_j = this->fvGeom.subContVol[j].local;

     FieldMatrix<Scalar,numEq,dim> pGrad(0.), xGrad(0.);
     FieldVector<Scalar,dim> temp(0.);
     SolutionVector flux(0.);
     FieldVector<Scalar,numEq> densityIJ(0.);
     FMatrix Ki(0), Kj(0);

      // calculate harmonic mean of permeabilities of nodes i and j
      Ki = this->problem.soil().K(global_i,e,local_i);
      Kj = this->problem.soil().K(global_j,e,local_j);
      const FMatrix K = harmonicMeanK(Ki, Kj);

     // calculate FE gradient (grad p for each phase)
     for (int k = 0; k < this->fvGeom.numVertices; k++) // loop over adjacent nodes
     {
         // FEGradient at node k
         const FieldVector<Scalar,dim> feGrad(this->fvGeom.subContVolFace[face].grad[k]);
         FieldVector<Scalar,numEq> pressure(0.0), massfrac(0.0);

         pressure[wPhase] = vNDat[k].pW;
         pressure[nPhase] = vNDat[k].pN;

         // compute sum of pressure gradients for each phase
         for (int phase = 0; phase < numEq; phase++)
         {
               temp = feGrad;
               temp *= pressure[phase];
               pGrad[phase] += temp;

               densityIJ[phase] += vNDat[k].density[phase]*this->fvGeom.subContVolFace[face].shapeValue[k];
          }
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
     FieldMatrix<Scalar,numEq,dim> contribComp(0);
     for (int phase=0; phase<numEq; phase++)
     {
         contribComp[phase] = problem.gravity();
         contribComp[phase] *= densityIJ[phase];
         pGrad[phase] -= contribComp[phase]; // grad p - rho*g
     }

     SolutionVector outward(0);  // Darcy velocity of each phase

     // calculate the advective flux using upwind: K*n(grad p -rho*g)
     for (int phase=0; phase<numEq; phase++)
     {
         FieldVector<Scalar,dim> v_tilde(0);
         K.mv(pGrad[phase], v_tilde);  // v_tilde=K*gradP
         outward[phase] = v_tilde*normal;
     }

     // evaluate upwind nodes
     int up_w, dn_w, up_n, dn_n;

     if (outward[wPhase] <= 0) {up_w = i; dn_w = j;}
     else {up_w = j; dn_w = i;};
     if (outward[nPhase] <= 0) {up_n = i; dn_n = j;}
     else {up_n = j; dn_n = i;};

     Scalar alpha = 1.0;  // Upwind parameter

     // water conservation
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
     // air conservation
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

         return flux;

     // DIFFUSION
     SolutionVector normDiffGrad;

     // get local to global id map
     int state_i = vNDat[i].phasestate;
     int state_j = vNDat[j].phasestate;

     Scalar diffusionWW(0.0), diffusionWN(0.0); // diffusion of water
     Scalar diffusionAW(0.0), diffusionAN(0.0); // diffusion of air
     SolutionVector avgDensity, avgDpm;

         // calculate tortuosity at the nodes i and j needed for porous media diffusion coefficient
     Scalar tauW_i, tauW_j, tauN_i, tauN_j; // tortuosity of wetting and nonwetting phase
     tauW_i = pow(sNDat[globalIdx_i].porosity * vNDat[i].satW,(7/3))/
             (sNDat[globalIdx_i].porosity*sNDat[globalIdx_i].porosity);
     tauW_j = pow(sNDat[globalIdx_j].porosity * vNDat[j].satW,(7/3))/
             (sNDat[globalIdx_j].porosity*sNDat[globalIdx_j].porosity);
     tauN_i = pow(sNDat[globalIdx_i].porosity * vNDat[i].satN,(7/3))/
             (sNDat[globalIdx_i].porosity*sNDat[globalIdx_i].porosity);
     tauN_j = pow(sNDat[globalIdx_j].porosity * vNDat[j].satN,(7/3))/
             (sNDat[globalIdx_j].porosity*sNDat[globalIdx_j].porosity);

     // arithmetic mean of porous media diffusion coefficient
     Scalar Dwn, Daw;
     Dwn = (sNDat[globalIdx_i].porosity * vNDat[i].satN * tauN_i * vNDat[i].diff[nPhase] +
                sNDat[globalIdx_j].porosity * vNDat[j].satN * tauN_j * vNDat[j].diff[nPhase])/2;
     Daw = (sNDat[globalIdx_i].porosity * vNDat[i].satW * tauW_i * vNDat[i].diff[wPhase] +
                sNDat[globalIdx_j].porosity * vNDat[j].satW * tauW_j * vNDat[j].diff[wPhase])/2;

     // adapt the diffusion coefficent according to the phase state.
     // TODO: make this continuously dependent on the phase saturations
     if (state_i == gasPhase || state_j == gasPhase) {
             // one element is only gas -> no diffusion in water phase
             avgDpm[wPhase] = 0;
     }
     if (state_i == waterPhase || state_j == waterPhase) {
             // one element is only water -> no diffusion in gas phase
             avgDpm[nPhase] = 0;
     }

     normDiffGrad[wPhase] = xGrad[wPhase]*normal;
     normDiffGrad[nPhase] = xGrad[nPhase]*normal;

     // calculate the arithmetic mean of densities
     avgDensity[wPhase] = 0.5*(vNDat[i].density[wPhase] + vNDat[j].density[wPhase]);
     avgDensity[nPhase] = 0.5*(vNDat[i].density[nPhase] + vNDat[j].density[nPhase]);

     // diffusion in the wetting phase
     diffusionAW = Daw * avgDensity[wPhase] * normDiffGrad[wPhase];
     diffusionWW = - diffusionAW; // air must be replaced by water
     diffusionWN = Dwn * avgDensity[nPhase] * normDiffGrad[nPhase];
     diffusionAN = - diffusionWN;

     // add diffusion of water to flux
     flux[water] += (diffusionWW + diffusionWN);
     //    std::cout << "Water Flux: " << flux[water] << std::endl;

     // add diffusion of air to flux
     flux[air] += (diffusionAN + diffusionAW);
     // std::cout << "Air Flux: " << flux[air] << std::endl;


     return flux;
    };

      /** @brief integrate sources / sinks
       *  @param e entity
       *  @param sol solution vector
       *  @param node local node id
       *  @return source/sink term
       */
      virtual SolutionVector computeQ (const Element& e, const SolutionVector* sol, const int node)
          {
              // ASSUME problem.q already contains \rho.q
              return problem.q(this->fvGeom.subContVol[node].global, e, this->fvGeom.subContVol[node].local);
          }

      /** @brief perform variable switch
       *  @param global global node id
       *  @param sol solution vector
       *  @param local local node id
       */
      virtual void primaryVarSwitch (const Element& e, int globalIdx, SolutionVector* sol, int localIdx)
          {
        bool switched = false;
        const FVector global = this->fvGeom.subContVol[localIdx].global;
        const FVector local = this->fvGeom.subContVol[localIdx].local;
        int state = sNDat[globalIdx].phaseState;
//        int switch_counter = sNDat[globalIdx].switched;

        // Evaluate saturation and pressures first
        Scalar pW = sol[localIdx][pWIdx];
        Scalar satW = 0.0;
        if (state == bothPhases)
            satW = 1.0-sol[localIdx][switchIdx];
        if (state == waterPhase)
            satW = 1.0;
        if (state == gasPhase)
            satW = 0.0;
        Scalar pC = problem.materialLaw().pC(satW, global, e, local);
        Scalar pN = pW + pC;

        switch(state)
        {
            case gasPhase :
            Scalar xWNmass, xWNmolar, pwn, pWSat; // auxiliary variables

            xWNmass = sol[localIdx][switchIdx];
            xWNmolar = problem.multicomp().convertMassToMoleFraction(xWNmass, gasPhase);
            pwn = xWNmolar * pN;
            pWSat = problem.multicomp().vaporPressure(temperature);

            if (pwn > pWSat && !switched)// && switch_counter < 3)
            {
                // appearance of water phase
                std::cout << "Water appears at node " << globalIdx << "  Coordinates: " << global << std::endl;
                sNDat[globalIdx].phaseState = bothPhases;
                sol[localIdx][switchIdx] = 1.0 - 2e-5; // initialize solution vector
                sNDat[globalIdx].switched += 1;
                switched = true;
            }
            break;

            case waterPhase :
            Scalar xAWmass, xAWmolar, henryInv, pbub; // auxiliary variables

            xAWmass = sol[localIdx][switchIdx];
            xAWmolar = problem.multicomp().convertMassToMoleFraction(xAWmass, waterPhase);

            henryInv = problem.multicomp().henry(temperature);
            pWSat = problem.multicomp().vaporPressure(temperature);
            pbub = pWSat + xAWmolar/henryInv; // pWSat + pAW

            if (pN < pbub && !switched)// && switch_counter < 3)
            {
                // appearance of gas phase
                std::cout << "Gas appears at node " << globalIdx << ",  Coordinates: " << global << std::endl;
                sNDat[globalIdx].phaseState = bothPhases;
                sol[localIdx][switchIdx] = 2e-5; // initialize solution vector
                sNDat[globalIdx].switched += 1;
                switched = true;
            }
            break;

            case bothPhases:
            Scalar satN = sol[localIdx][switchIdx];

            if (satN < -1e-5  && !switched)// && switch_counter < 3)
                {
                // disappearance of gas phase
                std::cout << "Gas disappears at node " << globalIdx << "  Coordinates: " << global << std::endl;
                sNDat[globalIdx].phaseState = waterPhase;
                sol[localIdx][switchIdx] = problem.multicomp().xAW(pN); // initialize solution vector
                sNDat[globalIdx].switched += 1;
                switched = true;
            }
            else if (satW < -1e-5 && !switched)// && switch_counter < 3)
                {
                // disappearance of water phase
                std::cout << "Water disappears at node " << globalIdx << "  Coordinates: " << global << std::endl;
                sNDat[globalIdx].phaseState = gasPhase;
                sol[localIdx][switchIdx] = problem.multicomp().xWN(pN); // initialize solution vector
                sNDat[globalIdx].switched += 1;
                switched = true;
                }
                break;

        }
        if (switched){
            updateVariableData(e, sol, localIdx, vNDat, sNDat[globalIdx].phaseState);
            BoxJacobian<ThisType,Grid,Scalar,2,BoxFunction>::localToGlobal(e,sol);
            setSwitchedLocal(); // if switch is triggered at any node, switchFlagLocal is set
        }

       return;
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

    virtual void computeElementData (const Element& e)
    {
//         // ASSUMING element-wise constant permeability, evaluate K at the element center
//          elData.K = problem.K(this->fvGeom.elementGlobal, e, this->fvGeom.elementLocal);
//
//         // ASSUMING element-wise constant porosity
//          elData.porosity = problem.porosity(this->fvGeom.elementGlobal, e, this->fvGeom.elementLocal);
        return;
    }


      //*********************************************************
      //*                                                        *
      //*    Calculation of Data at Nodes that has to be            *
      //*    determined only once    (sNDat)                        *
      //*                                                        *
      //*********************************************************

    // analog to EvalStaticData in MUFTE
    virtual void updateStaticDataVS (const Element& e, SolutionVector* sol)
    {
        // size of the sNDat vector is determined in the constructor

        // get access to shape functions for P1 elements
        GeometryType gt = e.geometry().type();
        const typename LagrangeShapeFunctionSetContainer<Scalar,Scalar,dim>::value_type&
        sfs=LagrangeShapeFunctions<Scalar,Scalar,dim>::general(gt,1);

        // get local to global id map
        for (int k = 0; k < sfs.size(); k++)
        {
           const int globalIdx = this->vertexMapper.template map<dim>(e, sfs[k].entity());

           // if nodes are not already visited
           if (!sNDat[globalIdx].visited)
            {
                 // evaluate primary variable switch
               primaryVarSwitch(e, globalIdx, sol, k);

                // mark elements that were already visited
                sNDat[globalIdx].visited = true;
            }
        }

      return;
    }

    // for initialization of the Static Data (sets porosity)
    virtual void updateStaticData (const Element& e, SolutionVector* sol)
    {
        // get access to shape functions for P1 elements
        GeometryType gt = e.geometry().type();
        const typename LagrangeShapeFunctionSetContainer<Scalar,Scalar,dim>::value_type&
        sfs=LagrangeShapeFunctions<Scalar,Scalar,dim>::general(gt,1);

        // get local to global id map
        for (int k = 0; k < sfs.size(); k++)
        {
            const int globalIdx = this->vertexMapper.template map<dim>(e, sfs[k].entity());

           // if nodes are not already visited
           if (!sNDat[globalIdx].visited)
           {
               // ASSUME porosity defined at nodes
               sNDat[globalIdx].porosity = problem.soil().porosity(this->fvGeom.elementGlobal, e, this->fvGeom.elementLocal);

               // set counter for variable switch to zero
               sNDat[globalIdx].switched = 0;

//               if (!checkSwitched())
//                  {
             primaryVarSwitch(e, globalIdx, sol, k);
//                  }

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


    struct VariableNodeData
    {
        Scalar satN;
     Scalar satW;
     Scalar pW;
     Scalar pC;
     Scalar pN;
     Scalar temperature;
     SolutionVector mobility;  //Vector with the number of phases
     SolutionVector density;
     FieldMatrix<Scalar,numComp,numEq> massfrac;
     int phasestate;
     SolutionVector diff;
    };

    // analog to EvalPrimaryData in MUFTE, uses members of vNDat
    virtual void updateVariableData(const Element& e, const SolutionVector* sol,
            int i, std::vector<VariableNodeData>& varData, int state)
    {
               const int globalIdx = this->vertexMapper.template map<dim>(e, i);
               FVector& global = this->fvGeom.subContVol[i].global;
               FVector& local = this->fvGeom.subContVol[i].local;

            varData[i].pW = sol[i][pWIdx];
            if (state == bothPhases) varData[i].satN = sol[i][switchIdx];
            if (state == waterPhase) varData[i].satN = 0.0;
            if (state == gasPhase) varData[i].satN = 1.0;

            varData[i].satW = 1.0 - varData[i].satN;

            varData[i].pC = problem.materialLaw().pC(varData[i].satW, global, e, local);
// For computation of a second pc-Sw law
//           const std::vector<double>& param = problem.soil().paramRelPerm2(global, e, local);
//           varData[i].pC = problem.materialLaw().pC(varData[i].satW, global, e, local, param);

            varData[i].pN = varData[i].pW + varData[i].pC;
            varData[i].temperature = temperature; // in [K], constant

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
            varData[i].mobility[wPhase] = problem.materialLaw().mobW(varData[i].satW, global, e, local, varData[i].temperature, varData[i].pW);
            varData[i].mobility[nPhase] = problem.materialLaw().mobN(varData[i].satN, global, e, local, varData[i].temperature, varData[i].pN);
            // Density of Water is set constant here!
            varData[i].density[wPhase] = 1000;//problem.wettingPhase().density(varData[i].temperature, varData[i].pN);
            varData[i].density[nPhase] = problem.nonwettingPhase().density(varData[i].temperature, varData[i].pN,
                    varData[i].massfrac[water][nPhase]);

         varData[i].diff[wPhase] = problem.wettingPhase().diffCoeff();
         varData[i].diff[nPhase] = problem.nonwettingPhase().diffCoeff();

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
               (*outMassFracAir)[globalIdx] = varData[i].massfrac[air][wPhase];
               (*outMassFracWater)[globalIdx] = varData[i].massfrac[water][nPhase];
               (*outDensityW)[globalIdx] = varData[i].density[wPhase];
               (*outDensityN)[globalIdx] = varData[i].density[nPhase];
               (*outMobilityW)[globalIdx] = varData[i].mobility[wPhase];
               (*outMobilityN)[globalIdx] = varData[i].mobility[nPhase];
               (*outPhaseState)[globalIdx] = varData[i].phasestate;

               return;
    }

    virtual void updateVariableData(const Element& e, const SolutionVector* sol, int i, bool old = false)
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

    void updateVariableData(const Element& e, const SolutionVector* sol, bool old = false)
    {
        int size = this->fvGeom.numVertices;

        for (int i = 0; i < size; i++)
                updateVariableData(e, sol, i, old);
    }

    bool checkSwitched()
    {
        return switchFlag;
    }

    bool checkSwitchedLocal()
    {
        return switchFlagLocal;
    }

    void setSwitchedLocalToGlobal()
    {
        switchFlag = switchFlagLocal;
        return;
    }

    void setSwitchedLocal()
    {
        switchFlagLocal = true;
        return;
    }

    void resetSwitched()
    {
        switchFlag = false;
        return;
    }

    void resetSwitchedLocal()
    {
        switchFlagLocal = false;
        return;
    }

    struct StaticNodeData
    {
        bool visited;
        int switched;
        int phaseState;
        int oldPhaseState;
        Scalar elementVolume;
        Scalar porosity;
        FMatrix K;
    };

    struct ElementData {
//        Scalar elementVolume;
//          Scalar porosity;
//        Scalar gravity;
        } elData;


    // parameters given in constructor
    TwoPTwoCProblem<Grid,Scalar>& problem;
    CWaterAir multicomp;
    std::vector<StaticNodeData> sNDat;
    std::vector<VariableNodeData> vNDat;
    std::vector<VariableNodeData> oldVNDat;

    // for output files
    BlockVector<FieldVector<Scalar, 1> > *outPressureN;
    BlockVector<FieldVector<Scalar, 1> > *outCapillaryP;
    BlockVector<FieldVector<Scalar, 1> > *outSaturationN;
    BlockVector<FieldVector<Scalar, 1> > *outSaturationW;
    BlockVector<FieldVector<Scalar, 1> > *outMassFracAir;
    BlockVector<FieldVector<Scalar, 1> > *outMassFracWater;
    BlockVector<FieldVector<Scalar, 1> > *outDensityW;
    BlockVector<FieldVector<Scalar, 1> > *outDensityN;
    BlockVector<FieldVector<Scalar, 1> > *outMobilityW;
    BlockVector<FieldVector<Scalar, 1> > *outMobilityN;
    BlockVector<FieldVector<Scalar, 1> > *outPhaseState;
//    BlockVector<FieldVector<Scalar, 1> > *outPermeability;

  protected:
        bool switchFlag;
        bool switchFlagLocal;
        double temperature;
  };

}
#endif
