// $Id: box2p2cjacobian.hh 1514 2009-03-31 10:08:24Z klaus $

#ifndef DUNE_BOXDARCYTRANSPORTJACOBIAN_HH
#define DUNE_BOXDARCYTRANSPORTJACOBIAN_HH

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
class BoxDarcyTransportJacobian
    : public BoxJacobian<BoxDarcyTransportJacobian<Grid,Scalar,BoxFunction>,Grid,Scalar,2,BoxFunction>
{
    enum {dim=Grid::dimension};
    enum {numEq=2, numComp=2};
    typedef typename Grid::Traits::template Codim<0>::Entity Element;
    typedef typename Element::Geometry Geometry;
    typedef BoxDarcyTransportJacobian<Grid,Scalar,BoxFunction> ThisType;
    typedef BoxJacobian<ThisType,Grid,Scalar,numEq,BoxFunction> BoxJacobianType;
    typedef typename LocalJacobian<ThisType,Grid,Scalar,numEq>::VBlockType SolutionVector;
    typedef typename LocalJacobian<ThisType,Grid,Scalar,numEq>::MBlockType MBlockType;
    typedef typename Grid::template Codim<0>::LeafIntersectionIterator IntersectionIterator;

    enum {pWIdx = 0, switchIdx = 1, numberOfComponents = 2};    // Solution vector index
    enum {wPhase = 0, nPhase = 1};                              // Phase index
    enum {gasPhase = 0, waterPhase = 1, bothPhases = 2};        // Phase state
    enum {water = 0, air = 1};                                  // Component index

public:
    // define the number of phases (numEq) and components (numComp) of your system, this is used outside
    // to allocate the correct size of (dense) blocks with a FieldMatrix
    enum {SIZE=LagrangeShapeFunctionSetContainer<Scalar,Scalar,dim>::maxsize};
    struct VariableNodeData;

    typedef FieldMatrix<Scalar,dim,dim> FMatrix;
    typedef FieldVector<Scalar,dim> FVector;

    //! Constructor
    BoxDarcyTransportJacobian (TwoPTwoCProblem<Grid,Scalar>& params, bool levelBoundaryAsDirichlet_, const Grid& grid,
                     BoxFunction& sol, bool procBoundaryAsDirichlet_=true)
        : BoxJacobian<ThisType,Grid,Scalar,2,BoxFunction>(levelBoundaryAsDirichlet_, grid, sol, procBoundaryAsDirichlet_),
          problem(params), sNDat(this->vertexMapper.size()), vNDat(SIZE), oldVNDat(SIZE), switchFlag(false)
    {
        this->analytic = false;
        switchFlag = false;
        switchFlagLocal = false;
        temperature = 283.15;
    }

    void localDefect(const Element& element, const SolutionVector* sol, bool withBC = true)
    {
        BoxJacobianType::localDefect(element, sol, withBC);

        this->assembleBoundaryCondition(element);

        Dune::GeometryType gt = element.geometry().type();
        const typename ReferenceElementContainer<Scalar,dim>::value_type& referenceElement = ReferenceElements<Scalar, dim>::general(gt);

        for (int vert=0; vert < this->fvGeom.numVertices; vert++) // begin loop over vertices / sub control volumes
            if (!this->fvGeom.subContVol[vert].inner)
            {
                FieldVector<Scalar,dim> averagedNormal(0);
                int faces = 0;
                IntersectionIterator endit = element.ileafend();
                for (IntersectionIterator it = element.ileafbegin(); it!=endit; ++it)
                {
                    if (it->boundary()) {
                        int faceIdx = it->indexInInside();
                        int numVerticesOfFace = referenceElement.size(faceIdx, 1, dim);
                        for (int nodeInFace = 0; nodeInFace < numVerticesOfFace; nodeInFace++) {
                            int nodeInElement = referenceElement.subEntity(faceIdx, 1, nodeInFace, dim);
                            if (nodeInElement != vert)
                                continue;

                            faces++;
                        }
                    }
                }

                if (faces == 2 && this->fvGeom.numVertices == 4)
                {
                    this->def[vert][pWIdx] = sol[0][pWIdx] + sol[3][pWIdx] - sol[1][pWIdx] - sol[2][pWIdx];
                }
            }

        return;
    }

    /** @brief compute time dependent term (storage), loop over nodes / subcontrol volumes
     *  @param element entity
     *  @param sol solution vector
     *  @param node local node id
     *  @return storage term
     */
    SolutionVector computeM (const Element& element, const SolutionVector* sol,
                                     int node, std::vector<VariableNodeData>& varData)
    {
        GeometryType gt = element.geometry().type();
        const typename LagrangeShapeFunctionSetContainer<Scalar,Scalar,dim>::value_type&
            sfs=LagrangeShapeFunctions<Scalar,Scalar,dim>::general(gt,1);

        int globalIdx = this->vertexMapper.template map<dim>(element, sfs[node].entity());

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

    SolutionVector computeM (const Element& element, const SolutionVector* sol, int node, bool old = false)
    {
        if (old)
            return computeM(element, sol, node, oldVNDat);
        else
            return computeM(element, sol, node, vNDat);
    }

    /** @brief compute diffusive/advective fluxes, loop over subcontrol volume faces
     *  @param element entity
     *  @param sol solution vector
     *  @param face face id
     *  @return flux term
     */
    SolutionVector computeA (const Element& element, const SolutionVector* sol, int face)
    {
        int i = this->fvGeom.subContVolFace[face].i;
        int j = this->fvGeom.subContVolFace[face].j;

        // normal vector, value of the area of the scvf
        const FieldVector<Scalar,dim> normal(this->fvGeom.subContVolFace[face].normal);
        GeometryType gt = element.geometry().type();
        const typename LagrangeShapeFunctionSetContainer<Scalar,Scalar,dim>::value_type&
            sfs=LagrangeShapeFunctions<Scalar,Scalar,dim>::general(gt,1);

        // global index of the subcontrolvolume face neighbor nodes in element element
        int globalIdx_i = this->vertexMapper.template map<dim>(element, sfs[i].entity());
        int globalIdx_j = this->vertexMapper.template map<dim>(element, sfs[j].entity());

        // get global coordinates of nodes i,j
        const FieldVector<Scalar,dim> global_i = this->fvGeom.subContVol[i].global;
        const FieldVector<Scalar,dim> global_j = this->fvGeom.subContVol[j].global;
        const FieldVector<Scalar,dim> local_i = this->fvGeom.subContVol[i].local;
        const FieldVector<Scalar,dim> local_j = this->fvGeom.subContVol[j].local;

        FieldMatrix<Scalar,numEq,dim> pGrad(0.), xGrad(0.);
        FieldVector<Scalar,dim> temp(0.);
        SolutionVector flux(0.);
        FieldVector<Scalar,numEq> densityIneumann(0.);
        FMatrix Ki(0), Kj(0);

        // calculate harmonic mean of permeabilities of nodes i and j
        Ki = this->problem.soil().K(global_i,element,local_i);
        Kj = this->problem.soil().K(global_j,element,local_j);
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
     *  @param element entity
     *  @param sol solution vector
     *  @param node local node id
     *  @return source/sink term
     */
    SolutionVector computeQ (const Element& element, const SolutionVector* sol, const int node)
    {
        // ASSUME problem.q already contains \rho.q
        return problem.q(this->fvGeom.subContVol[node].global, element, this->fvGeom.subContVol[node].local);
    }

    /** @brief perform variable switch
     *  @param global global node id
     *  @param sol solution vector
     *  @param local local node id
     */
    virtual void primaryVarSwitch (const Element& element, int globalIdx, SolutionVector* sol, int localIdx)
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
        Scalar pC = problem.materialLaw().pC(satW, global, element, local);
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
            updateVariableData(element, sol, localIdx, vNDat, sNDat[globalIdx].phaseState);
            BoxJacobian<ThisType,Grid,Scalar,2,BoxFunction>::localToGlobal(element,sol);
            setSwitchedLocal(); // if switch is triggered at any node, switchFlagLocal is set
        }

        return;
    }

    // harmonic mean of the permeability computed directly
    virtual FMatrix harmonicMeanK (FMatrix& Ki, const FMatrix& Kj) const
    {
        for (int kx=0; kx<dim; kx++){
            for (int ky=0; ky<dim; ky++){
                if (Ki[kx][ky] != Kj[kx][ky])
                {
		    Ki[kx][ky] = 2*Ki[kx][ky]*Kj[kx][ky] / (Ki[kx][ky]+Kj[kx][ky]);
                }
            }
        }
        return Ki;
    }


    virtual void clearVisited ()
    {
        for (int vertex = 0; vertex < this->vertexMapper.size(); vertex++){
            sNDat[vertex].visited = false;
            //            sNDat[vertex].switched = false;
        }
        return;
    }

    // updates old phase state after each time step
    virtual void updatePhaseState ()
    {
        for (int vertex = 0; vertex < this->vertexMapper.size(); vertex++){
            sNDat[vertex].oldPhaseState = sNDat[vertex].phaseState;
        }
        return;
    }

    virtual void resetPhaseState ()
    {
        for (int vertex = 0; vertex < this->vertexMapper.size(); vertex++){
            sNDat[vertex].phaseState = sNDat[vertex].oldPhaseState;
        }
        return;
    }

    //*********************************************************
    //*                                                        *
    //*    Calculation of Data at Elements (elData)             *
    //*                                                         *
    //*                                                        *
    //*********************************************************

    void computeElementData (const Element& element)
    {
        //         // ASSUMING element-wise constant permeability, evaluate K at the element center
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
    void updateStaticDataVS (const Element& element, SolutionVector* sol)
    {
        // size of the sNDat vector is determined in the constructor

        // get access to shape functions for P1 elements
        GeometryType gt = element.geometry().type();
        const typename LagrangeShapeFunctionSetContainer<Scalar,Scalar,dim>::value_type&
            sfs=LagrangeShapeFunctions<Scalar,Scalar,dim>::general(gt,1);

        // get local to global id map
        for (int k = 0; k < sfs.size(); k++)
        {
            const int globalIdx = this->vertexMapper.template map<dim>(element, sfs[k].entity());

            // if nodes are not already visited
            if (!sNDat[globalIdx].visited)
            {
                // evaluate primary variable switch
//                primaryVarSwitch(element, globalIdx, sol, k);

                // mark elements that were already visited
                sNDat[globalIdx].visited = true;
            }
        }

        return;
    }

    // for initialization of the Static Data (sets porosity)
    void updateStaticData (const Element& element, SolutionVector* sol)
    {
        // get access to shape functions for P1 elements
        GeometryType gt = element.geometry().type();
        const typename LagrangeShapeFunctionSetContainer<Scalar,Scalar,dim>::value_type&
            sfs=LagrangeShapeFunctions<Scalar,Scalar,dim>::general(gt,1);

        // get local to global id map
        for (int k = 0; k < sfs.size(); k++)
        {
            const int globalIdx = this->vertexMapper.template map<dim>(element, sfs[k].entity());

            // if nodes are not already visited
            if (!sNDat[globalIdx].visited)
            {
                // ASSUME porosity defined at nodes
                sNDat[globalIdx].porosity = problem.soil().porosity(this->fvGeom.elementGlobal, element, this->fvGeom.elementLocal);

                // set counter for variable switch to zero
                sNDat[globalIdx].switched = 0;

                //               if (!checkSwitched())
                //                  {
//                primaryVarSwitch(element, globalIdx, sol, k);
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
    void updateVariableData(const Element& element, const SolutionVector* sol,
                                    int vertex, std::vector<VariableNodeData>& varData, int state)
    {
//        const int globalIdx = this->vertexMapper.template map<dim>(element, vertex);
//        FVector& global = this->fvGeom.subContVol[vertex].global;
//        FVector& local = this->fvGeom.subContVol[vertex].local;

        varData[vertex].pW = sol[vertex][pWIdx];
        if (state == bothPhases) varData[vertex].satN = sol[vertex][switchIdx];
        if (state == waterPhase) varData[vertex].satN = 0.0;
        if (state == gasPhase) varData[vertex].satN = 1.0;

        varData[vertex].satW = 1.0 - varData[vertex].satN;

        varData[vertex].pC = 0;//problem.materialLaw().pC(varData[vertex].satW, global, element, local);
        // For computation of a second pc-Sw law
        //           const std::vector<Scalar>& param = problem.soil().paramRelPerm2(global, element, local);
        //           varData[vertex].pC = problem.materialLaw().pC(varData[vertex].satW, global, element, local, param);

        varData[vertex].pN = varData[vertex].pW + varData[vertex].pC;
        varData[vertex].temperature = temperature; // in [K], constant

        // Solubilities of components in phases
        if (state == bothPhases){
            varData[vertex].massfrac[air][wPhase] = problem.multicomp().xAW(varData[vertex].pN, varData[vertex].temperature);
            varData[vertex].massfrac[water][nPhase] = problem.multicomp().xWN(varData[vertex].pN, varData[vertex].temperature);
        }
        if (state == waterPhase){
            varData[vertex].massfrac[water][nPhase] = 0.0;
            varData[vertex].massfrac[air][wPhase] =  sol[vertex][switchIdx];
        }
        if (state == gasPhase){
            varData[vertex].massfrac[water][nPhase] = sol[vertex][switchIdx];
            varData[vertex].massfrac[air][wPhase] = 0.0;
        }
        varData[vertex].massfrac[water][wPhase] = 1.0 - varData[vertex].massfrac[air][wPhase];
        varData[vertex].massfrac[air][nPhase] = 1.0 - varData[vertex].massfrac[water][nPhase];
        varData[vertex].phasestate = state;

        // Mobilities & densities
        varData[vertex].mobility[wPhase] = 0;//problem.materialLaw().mobW(varData[vertex].satW, global, element, local, varData[vertex].temperature, varData[vertex].pW);
        varData[vertex].mobility[nPhase] = 100;//problem.materialLaw().mobN(varData[vertex].satN, global, element, local, varData[vertex].temperature, varData[vertex].pN);
        // Density of Water is set constant here!
        varData[vertex].density[wPhase] = 1000;//problem.wettingPhase().density(varData[vertex].temperature, varData[vertex].pN);
        varData[vertex].density[nPhase] = 1.23;//problem.nonwettingPhase().density(varData[vertex].temperature, varData[vertex].pN,
                                                                            //varData[vertex].massfrac[water][nPhase]);

        varData[vertex].diff[wPhase] = problem.wettingPhase().diffCoeff();
        varData[vertex].diff[nPhase] = problem.nonwettingPhase().diffCoeff();

        // CONSTANT solubility (for comparison with twophase)
        //         varData[vertex].massfrac[air][wPhase] = 0.0; varData[vertex].massfrac[water][wPhase] = 1.0;
        //         varData[vertex].massfrac[water][nPhase] = 0.0; varData[vertex].massfrac[air][nPhase] = 1.0;

        //std::cout << "water in gasphase: " << varData[vertex].massfrac[water][nPhase] << std::endl;
        //std::cout << "air in waterphase: " << varData[vertex].massfrac[air][wPhase] << std::endl;

//        // for output
//        (*outPressureN)[globalIdx] = varData[vertex].pN;
//        (*outCapillaryP)[globalIdx] = varData[vertex].pC;
//        (*outSaturationW)[globalIdx] = varData[vertex].satW;
//        (*outSaturationN)[globalIdx] = varData[vertex].satN;
//        (*outMassFracAir)[globalIdx] = varData[vertex].massfrac[air][wPhase];
//        (*outMassFracWater)[globalIdx] = varData[vertex].massfrac[water][nPhase];
//        (*outDensityW)[globalIdx] = varData[vertex].density[wPhase];
//        (*outDensityN)[globalIdx] = varData[vertex].density[nPhase];
//        (*outMobilityW)[globalIdx] = varData[vertex].mobility[wPhase];
//        (*outMobilityN)[globalIdx] = varData[vertex].mobility[nPhase];
//        (*outPhaseState)[globalIdx] = varData[vertex].phasestate;

        return;
    }

    void updateVariableData(const Element& element, const SolutionVector* sol, int vertex, bool old = false)
    {
        int state;
        const int global = this->vertexMapper.template map<dim>(element, vertex);
        if (old)
        {
            state = sNDat[global].oldPhaseState;
            updateVariableData(element, sol, vertex, oldVNDat, state);
        }
        else
        {
            state = sNDat[global].phaseState;
            updateVariableData(element, sol, vertex, vNDat, state);
        }
    }

    void updateVariableData(const Element& element, const SolutionVector* sol, bool old = false)
    {
        int size = this->fvGeom.numVertices;

        for (int vertex = 0; vertex < size; vertex++)
            updateVariableData(element, sol, vertex, old);
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
//    BlockVector<FieldVector<Scalar, 1> > *outPressureN;
//    BlockVector<FieldVector<Scalar, 1> > *outCapillaryP;
//    BlockVector<FieldVector<Scalar, 1> > *outSaturationN;
//    BlockVector<FieldVector<Scalar, 1> > *outSaturationW;
//    BlockVector<FieldVector<Scalar, 1> > *outMassFracAir;
//    BlockVector<FieldVector<Scalar, 1> > *outMassFracWater;
//    BlockVector<FieldVector<Scalar, 1> > *outDensityW;
//    BlockVector<FieldVector<Scalar, 1> > *outDensityN;
//    BlockVector<FieldVector<Scalar, 1> > *outMobilityW;
//    BlockVector<FieldVector<Scalar, 1> > *outMobilityN;
//    BlockVector<FieldVector<Scalar, 1> > *outPhaseState;
//    //    BlockVector<FieldVector<Scalar, 1> > *outPermeability;

protected:
    bool switchFlag;
    bool switchFlagLocal;
    Scalar temperature;
};

}
#endif
