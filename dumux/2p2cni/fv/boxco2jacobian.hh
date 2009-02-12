// $Id$

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
 - Scalar    type used for return values
 */
template<class Grid, class Scalar, class BoxFunction = LeafP1Function<
        Grid, Scalar, 3> >
class BoxCO2Jacobian: public BoxJacobian<BoxCO2Jacobian<Grid, Scalar,
        BoxFunction> , Grid, Scalar, 3, BoxFunction>
{
typedef    typename Grid::ctype CoordScalar;
    typedef typename Grid::Traits::template Codim<0>::Entity Element;
    typedef typename Element::Geometry Geometry;
    typedef BoxCO2Jacobian<Grid,Scalar,BoxFunction> ThisType;
    typedef typename LocalJacobian<ThisType,Grid,Scalar,3>::VBlockType VBlockType;
    typedef typename LocalJacobian<ThisType,Grid,Scalar,3>::MBlockType MBlockType;

    enum
    {   pWIdx = 0, switchIdx = 1, teIdx=2}; // Solution vector index
    enum
    {   wPhase = 0, nPhase = 1, temp = 2}; // Phase index
    enum
    {   gasPhase = 0, waterPhase = 1, bothPhases = 2}; // Phase state
    enum
    {   wComp = 0, nComp = 1, heat = 2}; // Component index

public:
    // define the number of phases (numEq) and components (numComp) of your system, this is used outside
    // to allocate the correct size of (dense) blocks with a FieldMatrix
    enum
    {   dim=Grid::dimension};
    enum
    {   numEq=3, numComp=2};
    enum
    {   SIZE=LagrangeShapeFunctionSetContainer<CoordScalar,Scalar,dim>::maxsize};
    struct VariableNodeData;

    typedef FieldMatrix<Scalar,dim,dim> FMatrix;
    typedef FieldVector<Scalar,dim> FVector;

    //! Constructor
    BoxCO2Jacobian (TwoPTwoCNIProblem<Grid,Scalar>& params,
            bool levelBoundaryAsDirichlet_, const Grid& grid,
            BoxFunction& sol,
            bool procBoundaryAsDirichlet_=true)
    : BoxJacobian<ThisType,Grid,Scalar,numEq,BoxFunction>(levelBoundaryAsDirichlet_, grid, sol, procBoundaryAsDirichlet_),
    problem(params),
    sNDat(this->vertexMapper.size()), vNDat(SIZE), oldVNDat(SIZE), switched(false), switchBreak(false), switchedGlobal(false), primVarSet(false)
    {
        this->analytic = false;
    }

    /** @brief compute time dependent term (storage), loop over nodes / subcontrol volumes
     *  @param element entity
     *  @param sol solution vector
     *  @param idx local node id
     *  @return storage term
     */
    virtual VBlockType computeM (const Element& element, const VBlockType* sol,
            int idx, std::vector<VariableNodeData>& varData)
    {
        GeometryType gt = element.geometry().type();
        const typename LagrangeShapeFunctionSetContainer<CoordScalar,Scalar,dim>::value_type&
        sfs=LagrangeShapeFunctions<CoordScalar,Scalar,dim>::general(gt,1);

        int globalIdx = this->vertexMapper.template map<dim>(element, sfs[idx].entity());

        VBlockType result;
        Scalar satN = varData[idx].satN;
        Scalar satW = varData[idx].satW;

        // storage of main component wetting phase
        result[wComp] =
        sNDat[globalIdx].porosity*(varData[idx].density[wPhase]*satW*varData[idx].massfrac[wComp][wPhase]
                +varData[idx].density[nPhase]*satN*varData[idx].massfrac[wComp][nPhase]);
        // storage of main component in nonwetting phase
        result[nComp] =
        sNDat[globalIdx].porosity*(varData[idx].density[nPhase]*satN*varData[idx].massfrac[nComp][nPhase]
                +varData[idx].density[wPhase]*satW*varData[idx].massfrac[nComp][wPhase]);

        // storage term of energy equation
        result[heat] = sNDat[globalIdx].porosity * (varData[idx].density[wPhase] * varData[idx].intenergy[wPhase] * satW
                + varData[idx].density[nPhase] * varData[idx].intenergy[nPhase] * satN)
        + elData.heatCap * varData[idx].temperature;
        // soil properties defined at the elements!!!

        //DEBUG
//        for(int j=0; j<3; j++)
//        {
//            if(isinf(result[j]))
//            {
//                std::cout<<"INF in ComputeM \n" << "Coordinates:X="<< this->fvGeom.subContVol[idx].global[0] <<" Y="<< this->fvGeom.subContVol[idx].global[1] <<
//                " Z="<< this->fvGeom.subContVol[idx].global[2]<<"\n"
//                <<"wetting phase pressure: " << varData[idx].pW << "wetting phase saturation: "<< satW << "temperature: "<< varData[idx].temperature << std::endl;
//            }
//        }

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
        int idx_i = this->fvGeom.subContVolFace[face].i;
        int idx_j = this->fvGeom.subContVolFace[face].j;

        // normal vector, value of the area of the scvf
        const FieldVector<Scalar,dim> normal(this->fvGeom.subContVolFace[face].normal);
        GeometryType gt = element.geometry().type();
        const typename LagrangeShapeFunctionSetContainer<CoordScalar,Scalar,dim>::value_type&
        sfs=LagrangeShapeFunctions<CoordScalar,Scalar,dim>::general(gt,1);

        // global index of the subcontrolvolume face neighbour nodes in element element
        int globalIdx_i = this->vertexMapper.template map<dim>(element, sfs[idx_i].entity());
        int globalIdx_j = this->vertexMapper.template map<dim>(element, sfs[idx_j].entity());

        // get global coordinates of nodes idx_i,idx_j
        const FieldVector<CoordScalar,dim> globalPos_i = this->fvGeom.subContVol[idx_i].global;
        const FieldVector<CoordScalar,dim> globalPos_j = this->fvGeom.subContVol[idx_j].global;
        const FieldVector<CoordScalar,dim> localPos_i = this->fvGeom.subContVol[idx_i].local;
        const FieldVector<CoordScalar,dim> localPos_j = this->fvGeom.subContVol[idx_j].local;

        /////////////////////////////////////////////////////////////////////////////
        // AVERAGING to get parameter values at integration points
        // Harmonic Mean:

        // harmonic mean of the permeability
        FMatrix Ki(0.), Kj(0.);
        Ki = this->problem.soil().K(globalPos_i,element,localPos_i);
        Kj = this->problem.soil().K(globalPos_j,element,localPos_j);
        const FMatrix K = harmonicMeanK(Ki, Kj);

        // harmonic mean of the heat conductivity lambda
        Scalar lambda;
        lambda = 2./((1./vNDat[idx_i].lambda) + (1./vNDat[idx_j].lambda));

        // Arithmetic Mean:

        // calculate tortuosity at the nodes idx_i and idx_j needed for porous media diffusion coefficient
        Scalar tauW_i, tauW_j, tauN_i, tauN_j; // tortuosity of wetting and nonwetting phase
        tauW_i = pow(sNDat[globalIdx_i].porosity * vNDat[idx_i].satW,(7/3))/
        (sNDat[globalIdx_i].porosity*sNDat[globalIdx_i].porosity);
        tauW_j = pow(sNDat[globalIdx_j].porosity * vNDat[idx_j].satW,(7/3))/
        (sNDat[globalIdx_j].porosity*sNDat[globalIdx_j].porosity);
        tauN_i = pow(sNDat[globalIdx_i].porosity * vNDat[idx_i].satN,(7/3))/
        (sNDat[globalIdx_i].porosity*sNDat[globalIdx_i].porosity);
        tauN_j = pow(sNDat[globalIdx_j].porosity * vNDat[idx_j].satN,(7/3))/
        (sNDat[globalIdx_j].porosity*sNDat[globalIdx_j].porosity);

        // arithmetic mean of porous media diffusion coefficient
        Scalar Dwg, Daw;
        Dwg = 0.5*(sNDat[globalIdx_i].porosity * vNDat[idx_i].satN * tauN_i * vNDat[idx_i].diff[nPhase] +
                sNDat[globalIdx_j].porosity * vNDat[idx_j].satN * tauN_j * vNDat[idx_j].diff[nPhase]);
        Daw = 0.5*(sNDat[globalIdx_i].porosity * vNDat[idx_i].satW * tauW_i * vNDat[idx_i].diff[wPhase] +
                sNDat[globalIdx_j].porosity * vNDat[idx_j].satW * tauW_j * vNDat[idx_j].diff[wPhase]);

        // arithmetic mean of phase enthalpies (dissolved components neglected)
        Scalar enthW;
        Scalar enthCO2;
        enthW = 0.5*(vNDat[idx_i].enthalpy[wPhase] + vNDat[idx_j].enthalpy[wPhase]);
        enthCO2 = 0.5*(vNDat[idx_i].enthalpy[nPhase] + vNDat[idx_j].enthalpy[nPhase]);

        // Calculate arithmetic mean of the densities
        VBlockType avgDensity;
        avgDensity[wPhase] = 0.5*(vNDat[idx_i].density[wPhase] + vNDat[idx_j].density[wPhase]);
        avgDensity[nPhase] = 0.5*(vNDat[idx_i].density[nPhase] + vNDat[idx_j].density[nPhase]);

        //////////////////////////////////////////////////////////////////////////////////////////
        // GRADIENTS

        FieldMatrix<Scalar,2,dim> pGrad(0.), xGrad(0.);
        FieldVector<Scalar,dim> teGrad(0.);
        FieldVector<Scalar,dim> tmp(0.);
        VBlockType flux(0.);

        // calculate FE gradient at subcontrolvolumeface
        for (int k = 0; k < this->fvGeom.numVertices; k++) // loop over adjacent nodes

        {
            // FEGradient at subcontrolvolumeface face
            const FieldVector<CoordScalar,dim> feGrad(this->fvGeom.subContVolFace[face].grad[k]);
            FieldVector<Scalar,2> pressure(0.0);

            pressure[wPhase] = vNDat[k].pW;
            pressure[nPhase] = vNDat[k].pN;

            // compute pressure gradients for each phase at integration point of subcontrolvolumeface face
            for (int phase = 0; phase < 2; phase++)
            {
                tmp = feGrad;
                tmp *= pressure[phase];
                pGrad[phase] += tmp;
            }

            // compute temperature gradient
            tmp = feGrad;
            tmp *= vNDat[k].temperature;
            teGrad += tmp;

            // compute concentration gradients
            // for diffusion of nComp in wetting phase
            tmp = feGrad;
            tmp *= vNDat[k].massfrac[nComp][wPhase];
            xGrad[wPhase] += tmp;

            // for diffusion of wetting phase in nonwetting phase
            tmp = feGrad;
            tmp *= vNDat[k].massfrac[wComp][nPhase];
            xGrad[nPhase] += tmp;
        }

        // deduce gravity*density of each phase
        FieldMatrix<Scalar,2,dim> contribComp(0);
        for (int phase=0; phase<2; phase++)
        {
            contribComp[phase] = problem.gravity();
            contribComp[phase] *= avgDensity[phase];
            pGrad[phase] -= contribComp[phase]; // grad p - rho*g
        }

        // Darcy velocity in normal direction for each phase K*n(grad p -rho*g)
        VBlockType outward(0);
        FieldVector<Scalar,dim> v_tilde_w(0.);
        FieldVector<Scalar,dim> v_tilde_n(0.);

        K.umv(pGrad[wPhase], v_tilde_w); // v_tilde=K*gradP
        outward[wPhase] = v_tilde_w*normal;
        K.umv(pGrad[nPhase], v_tilde_n); // v_tilde=K*gradP
        outward[nPhase] = v_tilde_n*normal;

        // Heat conduction
        outward[heat] = teGrad * normal;
        outward[heat] *= lambda;

        // evaluate upwind nodes
        int up_w, dn_w, up_n, dn_n;
        if (outward[wPhase] <= 0)
        {   up_w = idx_i; dn_w = idx_j;}
        else
        {   up_w = idx_j; dn_w = idx_i;};
        if (outward[nPhase] <= 0)
        {   up_n = idx_i; dn_n = idx_j;}
        else
        {   up_n = idx_j; dn_n = idx_i;};

        Scalar alpha = 1.0; // Upwind parameter

        ////////////////////////////////////////////////////////////////////////////////////////////////
        // ADVECTIVE TRANSPORT
        // Water conservation
        flux[wComp] = (alpha* vNDat[up_w].density[wPhase]*vNDat[up_w].mobility[wPhase]
                * vNDat[up_w].massfrac[wComp][wPhase]
                + (1-alpha)* vNDat[dn_w].density[wPhase]*vNDat[dn_w].mobility[wPhase]
                * vNDat[dn_w].massfrac[wComp][wPhase])
        * outward[wPhase];
        flux[wComp] += (alpha* vNDat[up_n].density[nPhase]*vNDat[up_n].mobility[nPhase]
                * vNDat[up_n].massfrac[wComp][nPhase]
                + (1-alpha)* vNDat[dn_n].density[nPhase]*vNDat[dn_n].mobility[nPhase]
                * vNDat[dn_n].massfrac[wComp][nPhase])
        * outward[nPhase];
        // CO2 conservation
        flux[nComp] = (alpha* vNDat[up_n].density[nPhase]*vNDat[up_n].mobility[nPhase]
                * vNDat[up_n].massfrac[nComp][nPhase]
                + (1-alpha)* vNDat[dn_n].density[nPhase]*vNDat[dn_n].mobility[nPhase]
                * vNDat[dn_n].massfrac[nComp][nPhase])
        * outward[nPhase];
        flux[nComp] += (alpha* vNDat[up_w].density[wPhase]*vNDat[up_w].mobility[wPhase]
                * vNDat[up_w].massfrac[nComp][wPhase]
                + (1-alpha)* vNDat[dn_w].density[wPhase]*vNDat[dn_w].mobility[wPhase]
                * vNDat[dn_w].massfrac[nComp][wPhase])
        * outward[wPhase];
        // Heat conservation
        flux[heat] = (alpha* vNDat[up_n].density[nPhase]*vNDat[up_n].mobility[nPhase] * vNDat[up_n].enthalpy[nPhase]
                + (1-alpha)* vNDat[dn_n].density[nPhase]*vNDat[dn_n].mobility[nPhase] * vNDat[dn_n].enthalpy[nPhase])
        * outward[nPhase];
        flux[heat] += (alpha* vNDat[up_w].density[wPhase]*vNDat[up_w].mobility[wPhase] * vNDat[up_w].enthalpy[wPhase]
                + (1-alpha)* vNDat[dn_w].density[wPhase]*vNDat[dn_w].mobility[wPhase] * vNDat[dn_w].enthalpy[wPhase])
        * outward[wPhase];

        /////////////////////////////////////////////////////////////////////////////////////////////////////////
        // DIFFUSIVE TRANSPORT

        VBlockType normDiffGrad;

        // get local to global id map
        int state_i = sNDat[globalIdx_i].phaseState;
        int state_j = sNDat[globalIdx_j].phaseState;

        Scalar diffusionAW(0.0), diffusionWW(0.0), diffusionWN(0.0), diffusionAN(0.0);

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
        flux[wComp] += (diffusionWW + diffusionWN);

        // CO2 conservation
        flux[nComp] += (diffusionAN + diffusionAW);

        //         // Heat conservation
        //          flux[heat] += diffusionWW*enthW + diffusionAW*enthCO2; diffusive transport of enthalpy neglected here

        //////////////////////////////////////////////////////////////////////////////////////////////
        // HEAT CONDUCTION

        flux[heat] += outward[heat];

        //    BEGIN DEBUG
        //          for(int j=0; j<3; j++)
        //             {
        //                 if(isinf(flux[j]))
        //                 {
        //                     std::cout<<"INF in ComputeA \n" << "Coordinates upw_node:X="<< this->fvGeom.subContVol[up_w].global[0] <<" Y="<< this->fvGeom.subContVol[up_w].global[1] <<
        //                     " Z="<<  this->fvGeom.subContVol[up_w].global[2]<<"\n"
        //                     <<"wetting phase pressure: " << vNDat[up_w].pW << "wetting phase saturation: "<< vNDat[up_w].satW
        //                     << " \n temperature: "<< vNDat[up_w].temperature << "Xaw" << vNDat[up_w].massfrac[nComp][wPhase] << std::endl;
        //                 }
        //             }
        //     END DEBUG
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
     *  @param global global node id
     *  @param sol solution vector
     *  @param local local node id
     */
    virtual void primaryVarSwitch (const Element& element, int globalIdx, VBlockType* sol, int idx)
    {

        if(!switchBreak)
        {
            switched = false;
            int state = sNDat[globalIdx].phaseState;

            Scalar pW = sol[idx][wPhase];
            Scalar satW = 0.0;
            if (state == bothPhases) satW = 1.0-sol[idx][nPhase];
            if (state == waterPhase) satW = 1.0;
            if (state == gasPhase) satW = 0.0;

            Scalar pC = problem.materialLaw().pC(satW, this->fvGeom.subContVol[idx].global, element, this->fvGeom.subContVol[idx].local);
            Scalar pN = pW + pC;

            FVector Coordinates = this->fvGeom.subContVol[idx].global;

            switch(state)
            {
                case gasPhase :
                Scalar xWNmass;
                xWNmass = sol[idx][switchIdx];

                if (xWNmass> 0.001 && switched == false)
                {
                    // appearance of wetting phase phase
                    std::cout << "Water appears at node " << globalIdx << "  Coordinates: " << Coordinates << std::endl;
                    sNDat[globalIdx].phaseState = bothPhases;
                    sol[idx][nPhase] = 1.0 - 1.e-3; // initialize solution vector
                    switched = true;
                }
                break;

                case waterPhase :
                Scalar xAWmax, xAWmass;
                xAWmass = sol[idx][switchIdx];
                xAWmax = problem.multicomp().xAW(pN, sol[idx][temp]);

                if (xAWmass> xAWmax && switched == false)
                {
                    // appearance of gas phase
                    std::cout << "Gas appears at node " << globalIdx << "  Coordinates: " << Coordinates << std::endl;
                    sNDat[globalIdx].phaseState = bothPhases;
                    sol[idx][nPhase] = 1.e-3; // initialize solution vector
                    switched = true;
                }
                break;

                case bothPhases :
                Scalar satN = sol[idx][nPhase];
                if (satN < 0.0 && switched == false)
                {
                    // disappearance of gas phase
                    std::cout << "Gas disappears at node " << globalIdx << "  Coordinates: " << Coordinates << std::endl;
                    sNDat[globalIdx].phaseState = waterPhase;
                    sol[idx][switchIdx] = 1e-3; // initialize solution vector
                    switched = true;
                }
                else if (satW < -1.e-5 && switched == false)
                {
                    // disappearance of wetting phase phase
                    std::cout << "Water disappears at node " << globalIdx << "  Coordinates: " << Coordinates << std::endl;
                    sNDat[globalIdx].phaseState = gasPhase;
                    sol[idx][switchIdx] = 1e-8; // initialize solution vector
                    switched = true;
                }
                break;
            }

            if(switched)
            {
                // fill global solution vector with initialised local values
                BoxJacobian<ThisType,Grid,Scalar,numEq,BoxFunction>::localToGlobal(element, sol);
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
        for (int kx=0; kx<dim; kx++)
        {
            for (int ky=0; ky<dim; ky++)
            {
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
        for (int globalIdx = 0; globalIdx < this->vertexMapper.size(); globalIdx++)
        {
            sNDat[globalIdx].visited = false;
        }
        return;
    }

    // updates old phase state after each time step
    virtual void updatePhaseState ()
    {
        for (int globalIdx = 0; globalIdx < this->vertexMapper.size(); globalIdx++)
        {
            sNDat[globalIdx].oldPhaseState = sNDat[globalIdx].phaseState;
        }
        return;
    }

    virtual void resetPhaseState ()
    {
        for (int globalIdx = 0; globalIdx < this->vertexMapper.size(); globalIdx++)
        {
            sNDat[globalIdx].phaseState = sNDat[globalIdx].oldPhaseState;
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
        //          (this->fvGeom.cellGlobal, element, this->fvGeom.cellLocal);
        //
        //         // ASSUMING element-wise constant permeability, evaluate K at the cell center
        //          elData.K = problem.K(this->fvGeom.cellGlobal, element, this->fvGeom.cellLocal);
        //
        //         // ASSUMING element-wise constant porosity
        //          elData.porosity = problem.porosity(this->fvGeom.cellGlobal, element, this->fvGeom.cellLocal);
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
        //int localToGlobal[LagrangeShapeFunctionSetContainer<CoordScalar,Scalar,n>::maxsize];

        // get access to shape functions for P1 elements
        GeometryType gt = element.geometry().type();
        const typename LagrangeShapeFunctionSetContainer<CoordScalar,Scalar,dim>::value_type&
        sfs=LagrangeShapeFunctions<CoordScalar,Scalar,dim>::general(gt,1);

        // get local to global id map
        for (int k = 0; k < sfs.size(); k++)
        {
            const int globalIdx = this->vertexMapper.template map<dim>(element, sfs[k].entity());

            // if nodes are not already visited
            if (!sNDat[globalIdx].visited)
            {
                // ASSUME porosity defined at nodes
                sNDat[globalIdx].porosity = problem.soil().porosity(this->fvGeom.subContVol[k].global, element, this->fvGeom.subContVol[k].local);

                // global coordinates
                FieldVector<CoordScalar,dim> globalPos_i = this->fvGeom.subContVol[k].global;

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
        const typename LagrangeShapeFunctionSetContainer<CoordScalar,Scalar,dim>::value_type&
        sfs=LagrangeShapeFunctions<CoordScalar,Scalar,dim>::general(gt,1);

        // get local to global id map
        for (int k = 0; k < sfs.size(); k++)
        {
            const int globalIdx = this->vertexMapper.template map<dim>(element, sfs[k].entity());

            // if nodes are not already visited
            if (!sNDat[globalIdx].visited)
            {
                // ASSUME porosity defined at nodes
                sNDat[globalIdx].porosity = problem.soil().porosity(this->fvGeom.subContVol[k].global, element, this->fvGeom.subContVol[k].local);

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
    //    BlockVector<FieldVector<Scalar, 1> > *hackyMassFracCO2;
    //    BlockVector<FieldVector<Scalar, 1> > *hackyMassFracWater;
    //    BlockVector<FieldVector<Scalar, 1> > *hackySaturationN;

    //    void printVariableData()
    //    {
    //        for (int idx = 0; idx < 4; idx++)
    //        {
    //            std::cout << "new: idx = " << idx << ": satN = " << vNDat[idx].satN << ", satW = " << vNDat[idx].satW
    //                << ", pW = " << vNDat[idx].pW << ", pC = " << vNDat[idx].pC << ", pN = " << vNDat[idx].pN
    //                << ", T = " << vNDat[idx].temperature << ", lambda = " << vNDat[idx].lambda << std::endl;
    //            std::cout << "old: idx = " << idx << ": satN = " << oldVNDat[idx].satN << ", satW = " << oldVNDat[idx].satW
    //                << ", pW = " << oldVNDat[idx].pW << ", pC = " << oldVNDat[idx].pC << ", pN = " << oldVNDat[idx].pN
    //                << ", T = " << oldVNDat[idx].temperature << ", lambda = " << oldVNDat[idx].lambda << std::endl;
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
        Scalar viscosityCO2;
        Scalar krCO2;
        FieldVector<Scalar,2> mobility; //Vector with the number of phases
        FieldVector<Scalar,2> density;
        FieldMatrix<Scalar,numComp,2> massfrac;
        FieldVector<Scalar,2> enthalpy;
        FieldVector<Scalar,2> intenergy;
        FieldVector<Scalar,2> diff;
        int phasestate;
    };

    // analog to EvalPrimaryData in MUFTE, uses members of vNDat
    virtual void updateVariableData(const Element& element, const VBlockType* sol,
            int idx, std::vector<VariableNodeData>& varData, int state)
    {
        const int globalIdx = this->vertexMapper.template map<dim>(element, idx);

        varData[idx].pW = sol[idx][wPhase];
        if (state == bothPhases) varData[idx].satN = sol[idx][nPhase];
        if (state == waterPhase) varData[idx].satN = 0.0;
        if (state == gasPhase) varData[idx].satN = 1.0;

        varData[idx].satW = 1.0 - varData[idx].satN;

        varData[idx].pC = problem.materialLaw().pC(varData[idx].satW, this->fvGeom.subContVol[idx].global, element, this->fvGeom.subContVol[idx].local);
        varData[idx].pN = varData[idx].pW + varData[idx].pC;
        varData[idx].temperature = sol[idx][temp]; // in [K]

        // Solubilities of components in phases
        if (state == bothPhases)
        {
            varData[idx].massfrac[nComp][wPhase] = problem.multicomp().xAW(varData[idx].pN, varData[idx].temperature);
            varData[idx].massfrac[wComp][nPhase] = problem.multicomp().xWN(varData[idx].pN, varData[idx].temperature);
        }
        if (state == waterPhase)
        {
            varData[idx].massfrac[wComp][nPhase] = 0.0;
            varData[idx].massfrac[nComp][wPhase] = sol[idx][switchIdx];
        }
        if (state == gasPhase)
        {
            varData[idx].massfrac[wComp][nPhase] = sol[idx][switchIdx];
            varData[idx].massfrac[nComp][wPhase] = 0.0;
        }
        varData[idx].massfrac[wComp][wPhase] = 1.0 - varData[idx].massfrac[nComp][wPhase];
        varData[idx].massfrac[nComp][nPhase] = 1.0 - varData[idx].massfrac[wComp][nPhase];
        // for output
        //               (*hackySaturationN)[global] = varData[idx].satN;
        //               (*hackyMassFracCO2)[global] = varData[idx].massfrac[nComp][wPhase];
        //               (*hackyMassFracWater)[global] = varData[idx].massfrac[wComp][nPhase];
        // Diffusion coefficients
        // Mobilities & densities

        varData[idx].density[wPhase] = problem.wettingPhase().density(varData[idx].temperature, varData[idx].pW, varData[idx].massfrac[nComp][wPhase]);
        varData[idx].density[nPhase] = problem.nonwettingPhase().density(varData[idx].temperature, varData[idx].pN);
        varData[idx].mobility[wPhase] = problem.materialLaw().mobW(varData[idx].satW, this->fvGeom.subContVol[idx].global, element, this->fvGeom.subContVol[idx].local, varData[idx].temperature, varData[idx].pW);
        varData[idx].mobility[nPhase] = problem.materialLaw().mobN(varData[idx].satN, this->fvGeom.subContVol[idx].global, element, this->fvGeom.subContVol[idx].local, varData[idx].temperature, varData[idx].pN);
        varData[idx].lambda = problem.soil().heatCond(this->fvGeom.subContVol[idx].global, element, this->fvGeom.subContVol[idx].local, varData[idx].satW);
        varData[idx].enthalpy[wPhase] = problem.wettingPhase().enthalpy(varData[idx].temperature,varData[idx].pW);
        varData[idx].enthalpy[nPhase] = problem.nonwettingPhase().enthalpy(varData[idx].temperature,varData[idx].pN);
        varData[idx].intenergy[wPhase] = problem.wettingPhase().intEnergy(varData[idx].temperature,varData[idx].pW);
        varData[idx].intenergy[nPhase] = problem.nonwettingPhase().intEnergy(varData[idx].temperature,varData[idx].pN);
        varData[idx].diff[wPhase] = problem.wettingPhase().diffCoeff();
        varData[idx].diff[nPhase] = problem.nonwettingPhase().diffCoeff();

        // CONSTANT solubility (for comparison with twophase)
        //         varData[idx].massfrac[nComp][wPhase] = 0.0; varData[idx].massfrac[wComp][wPhase] = 1.0;
        //         varData[idx].massfrac[wComp][nPhase] = 0.0; varData[idx].massfrac[nComp][nPhase] = 1.0;

        //std::cout << "wComp in gasphase: " << varData[idx].massfrac[wComp][nPhase] << std::endl;
        //std::cout << "nComp in waterphase: " << varData[idx].massfrac[nComp][wPhase] << std::endl;

        // for output
        (*outPressureN)[globalIdx] = varData[idx].pN;
        (*outCapillaryP)[globalIdx] = varData[idx].pC;
        (*outSaturationW)[globalIdx] = varData[idx].satW;
        (*outSaturationN)[globalIdx] = varData[idx].satN;
        (*outTemperature)[globalIdx] = varData[idx].temperature;
        (*outMassFracAir)[globalIdx] = varData[idx].massfrac[nComp][wPhase];
        (*outMassFracWater)[globalIdx] = varData[idx].massfrac[wComp][nPhase];
        (*outDensityW)[globalIdx] = varData[idx].density[wPhase];
        (*outDensityN)[globalIdx] = varData[idx].density[nPhase];
        (*outMobilityW)[globalIdx] = varData[idx].mobility[wPhase];
        (*outMobilityN)[globalIdx] = varData[idx].mobility[nPhase];
        (*outPhaseState)[globalIdx] = state;

        return;
    }

    virtual void updateVariableData(const Element& element, const VBlockType* sol, int idx, bool old = false)
    {
        int state;
        const int globalIdx = this->vertexMapper.template map<dim>(element, idx);
        if (old)
        {
            state = sNDat[globalIdx].oldPhaseState;
            updateVariableData(element, sol, idx, oldVNDat, state);
        }
        else
        {
            state = sNDat[globalIdx].phaseState;
            updateVariableData(element, sol, idx, vNDat, state);
        }
    }

    void updateVariableData(const Element& element, const VBlockType* sol, bool old = false)
    {
        int size = this->fvGeom.numVertices;

        for (int idx = 0; idx < size; idx++)
        updateVariableData(element, sol, idx, old);
    }

    struct StaticNodeData
    {
        bool visited, primVarSet;
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

    struct ElementData
    {
        Scalar heatCap;

        //        Scalar cellVolume;
        //          Scalar porosity;
        //        Scalar gravity;
        //        FieldVector<Scalar, 4> parameters;
        //        FieldMatrix<Scalar,dim,dim> K;
    }elData;

    // parameters given in constructor
    TwoPTwoCNIProblem<Grid,Scalar>& problem;
    CBrineCO2 multicomp;
    std::vector<StaticNodeData> sNDat;
    std::vector<StaticIPData> sIPDat;
    std::vector<VariableNodeData> vNDat;
    std::vector<VariableNodeData> oldVNDat;
    bool switched, switchBreak, switchedGlobal, primVarSet;

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

//    // average permeability from the staticNode Vector
//    virtual FMatrix harmonicMeanK (const Element& element, int k) const
//    {
//     FMatrix Ki, Kj;
//     const Scalar eps = 1e-20;
//
//     int idx = this->fvGeom.subContVolFace[k].idx;
//     int j = this->fvGeom.subContVolFace[k].j;
//
//     int global_i = this->vertexMapper.template map<dim>(element, idx);
//     int global_j = this->vertexMapper.template map<dim>(element, j);
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

