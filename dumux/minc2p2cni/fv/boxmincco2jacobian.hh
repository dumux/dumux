// $Id: boxco2jacobian.hh 767 2008-11-05 13:11:14Z melanie $

#ifndef DUNE_BOXMINCCO2JACOBIAN_HH
#define DUNE_BOXMINCCO2JACOBIAN_HH

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
#include"dumux/minc2p2cni/minc2p2cniproblem.hh"
//#include "varswitch.hh"
#include"dumux/io/vtkmultiwriter.hh"

/**
 * @file
 * @brief  compute local jacobian matrix for box scheme for two-phase two-component flow equation
 * @author Bernd Flemisch, Klaus Mosthaf, Melanie Darcis, Alex Tatomir
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
template<class G, class RT, int m, class BoxFunction = LeafP1Function<G, RT, m> >
class BoxMincCO2Jacobian
    : public BoxJacobian<BoxMincCO2Jacobian<G,RT,m,BoxFunction>,G,RT,m,BoxFunction>
{
    typedef typename G::ctype DT;
    typedef typename G::Traits::template Codim<0>::Entity Entity;
    typedef typename Entity::Geometry Geometry;
    typedef BoxMincCO2Jacobian<G,RT,m,BoxFunction> ThisType;
    typedef typename LocalJacobian<ThisType,G,RT,m>::VBlockType VBlockType;
    typedef typename LocalJacobian<ThisType,G,RT,m>::MBlockType MBlockType;

    enum {pWIdx = 0, switchIdx = 1, teIdx=2, numberOfComponents = 2};    // Solution vector index
    enum {wPhase = 0, nPhase = 1};                                        // Phase index
    enum {gasPhase = 0, waterPhase = 1, bothPhases = 2};                // Phase state
    enum {waterF = 0, co2F = 1, heatF = 2};                                // Component index
    enum {water = 0, co2 = 1, heat = 2};                                // Component index
    enum {nPhases = 2};                                                    // Number of phases
    enum {nPorosities = 2};                                                // Number of porosities
    enum {nPermeabilities =2};                                            // Number of permeabilities
    //    enum {waterF = 0, co2F = 1 , heatF = 2 , waterM = 3, co2M = 4, heatM = 5};  // Minc Component index
    //    enum {waterF, co2F, heatF, waterM, co2M, heatM};  // Minc Component index

public:
    // define the number of phases (m) and components (c) of your system, this is used outside
    // to allocate the correct size of (dense) blocks with a FieldMatrix
    enum {dim=G::dimension};
    enum {c=2};
    //number of interacting continua
    enum {nCont = m/3};
    enum {SIZE=LagrangeShapeFunctionSetContainer<DT,RT,dim>::maxsize};
    struct VariableNodeData;

    typedef FieldMatrix<RT,dim,dim> FMatrix;
    typedef FieldVector<RT,dim> FVector;

    //! Constructor
    BoxMincCO2Jacobian (TwoPTwoCNIProblem<G,RT,m>& params,
                        bool levelBoundaryAsDirichlet_, const G& grid,
                        BoxFunction& sol,
                        bool procBoundaryAsDirichlet_=true)
        : BoxJacobian<ThisType,G,RT,m,BoxFunction>(levelBoundaryAsDirichlet_, grid, sol, procBoundaryAsDirichlet_),
          problem(params),
          sNDat(this->vertexMapper.size()), vNDat(SIZE), oldVNDat(SIZE), switched(false), switchBreak(false)
    {
        this->analytic = false;
    }


    virtual FieldVector<RT, nCont> AreaP (const Entity& e, int face)
    {
        FieldVector<RT,nCont> x(0.0);
        FieldVector<RT,nCont> y(0.0);
        FieldVector<RT,nCont> L(0.0);
        double length_x =0;      //the length of the subcontrol face on x direction
        double length_y =0;      //the length of the subcontrol face on y direction
        double ExtLength =0;        //Exterior length of the subcontrol face
        for (int face =0; face < this->fvGeom.numVertices; face++){
            int it = this->fvGeom.subContVolFace[face].i;
            int jt = this->fvGeom.subContVolFace[face].j;
            if (it == 0)
            {
                FieldVector<DT,dim> ht;
                ht = this->fvGeom.subContVolFace[face].normal;
                if (ht[0]<0) {ht[0]*=(-1);}
                if (ht[1]<0) {ht[1]*=(-1);}
                if (ht[0]==0){ht[0]=ht[1];}
                if (ht[1]==0){ht[1]=ht[0];}
                length_y = ht[0];
            }
            else if (jt==0)   {
                FieldVector<DT,dim> ht;
                ht = this->fvGeom.subContVolFace[face].normal;
                if (ht[0]<0) {ht[0]*=(-1);}
                if (ht[1]<0) {ht[1]*=(-1);}
                if (ht[0]==0){ht[0]=ht[1];}
                if (ht[1]==0){ht[1]=ht[0];}
                length_x = ht[0];
            }
            ExtLength = length_x+length_y;

        }
        //              double f = 1.0 / nCont; // volume fraction
        //              // Calculation of the distance between the nested volume elements delta
        //              // delta in this case divides the surfaces equidistantly (delta is const)
        elData.delta = 0.0;
        if (length_x > length_y){
            elData.delta = 2*length_y / (2.0 * (nCont-1) + 1.0);
        }
        else {
            elData.delta = 2*length_x / (2.0 * (nCont-1) + 1.0);
        }

        double TotalLength = 4*ExtLength;

        L[0]  = TotalLength;
        x[0]  = 2*length_x;    //length of the unitary element (CV) on x direction
        y[0]  = 2*length_y;    //length of the unitary element (CV) on y direction
        for (int nC=1; nC<nCont; nC++){
            x[nC] = x[nC-1] - 2 * elData.delta;
            y[nC] = y[nC-1] - 2 * elData.delta;
            L[nC] = 2* (x[nC]+y[nC]);
        }
        return L;

        //              return TotalLength;
    }



    /***********************************************************************************************************/
    /*the harmonicMeanKMinc function computes the harmonic mean of the perameabilities of fractures and rock at
     * the same global coordinate " i". It is used afterwards in the computation of the interporosity flux */
    virtual FMatrix harmonicMeanKMinc (FMatrix& KF, const FMatrix& KM) const
    {
        double eps = 1e-20;
        FieldMatrix<RT,dim, dim> K(0.0);
        for (int kx=0; kx<nPermeabilities; kx++){
            for (int ky=0; ky<nPermeabilities; ky++){
                if (KF[kx][ky] != KM[kx][ky])
                {
                    K[kx][ky] = 2 / (1/(KF[kx][ky]+eps) + (1/(KM[kx][ky]+eps)));
                }
                else
                {
                    K = KF;
                }
            }
        }
        return K;
    }
    //    /*the harmonicMeanK function computes the harmonic mean of the perameabilities between the two nodes of different
    //     *  permeabilities in the fracture domain */


    // harmonic mean computed directly
    virtual FMatrix harmonicMeanK (FMatrix& Ki, const FMatrix& Kj) const
    {
        double eps = 1e-20;
        FMatrix K(0.);
        for (int kx=0; kx<nPermeabilities; kx++){
            for (int ky=0; ky<nPermeabilities; ky++){
                if (Ki[kx][ky] != Kj[kx][ky])
                {
                    K[kx][ky] = 2 / (1/(Ki[kx][ky]+eps) + (1/(Kj[kx][ky]+eps)));
                }
                else K = Ki;
            }
        }
        return K;
    }

    ///*************************************************************************************************************/



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


        VBlockType AccumulationTerm(0);
        MBlockType result(0);
        RT satNFracture = varData[node].satNFracture;
        RT satWFracture = varData[node].satWFracture;
        //                 RT satNMatrix = varData[node].satNMatrix;
        //                 RT satWMatrix = varData[node].satWMatrix;
        RT satNFractureT = varData[node].saturation[nPhase][0];
        RT satWFractureT = varData[node].saturation[wPhase][0];


        for (int nC = 0; nC < nCont; nC++){
            //// In the fracture domain
            if (nC==0){

                result[waterF][nC] =
                    sNDat[globalIdx].porosityFracture*(varData[node].densityF[wPhase]*varData[node].saturation[wPhase][nC]*varData[node].massfrac[waterF][wPhase]
                                                       +varData[node].densityF[nPhase]*varData[node].saturation[nPhase][nC]*varData[node].massfrac[waterF][nPhase]);
                result[co2F][nC] =
                    sNDat[globalIdx].porosityFracture*(varData[node].densityF[nPhase]*varData[node].saturation[nPhase][nC]*varData[node].massfrac[co2F][nPhase]
                                                       +varData[node].densityF[wPhase]*varData[node].saturation[wPhase][nC]*varData[node].massfrac[co2F][wPhase]);
                result[heatF][nC] = sNDat[globalIdx].porosityFracture * (varData[node].densityF[wPhase] * varData[node].intenergy[wPhase] * varData[node].saturation[wPhase][nC]
                                                                         + varData[node].densityF[nPhase] * varData[node].intenergy[nPhase] * varData[node].saturation[nPhase][nC])
                    + elData.heatCap * varData[node].temperatureF;
            }
            /// In the rock domain
            else {
                result[waterF][nC] =
                    sNDat[globalIdx].porosityMatrix*(varData[node].densityF[wPhase]*varData[node].saturation[water][nC]*varData[node].massfrac[waterF][wPhase]
                                                     +varData[node].densityF[nPhase]*varData[node].saturation[co2][nC]*varData[node].massfrac[waterF][nPhase]);
                result[co2F][nC] =
                    sNDat[globalIdx].porosityMatrix*(varData[node].densityF[nPhase]*varData[node].saturation[co2][nC]*varData[node].massfrac[co2F][nPhase]
                                                     +varData[node].densityF[wPhase]*varData[node].saturation[water][nC]*varData[node].massfrac[co2F][wPhase]);
                result[heatF][nC] = sNDat[globalIdx].porosityMatrix * (varData[node].densityF[wPhase] * varData[node].intenergy[wPhase] * varData[node].saturation[water][nC]
                                                                       + varData[node].densityF[nPhase] * varData[node].intenergy[nPhase] * varData[node].saturation[co2][nC])
                    + elData.heatCap * varData[node].temperatureF;
            }

        }

        for (int nC = 0; nC < nCont ;nC++ ){
            //3 is the number of equations (water, CO2, heat) and nC - the number of Continua
            int count = 3*nC;
            AccumulationTerm[count]=result[water][nC];
            AccumulationTerm[count+1]=result[co2][nC];
            AccumulationTerm[count+2]=result[heat][nC];
        }

        return AccumulationTerm;

    };

    virtual VBlockType computeM (const Entity& e, const VBlockType* sol, int node, bool old = false)
    {
        if (old)
            return computeM(e, sol, node, oldVNDat);
        else
            return computeM(e, sol, node, vNDat);
    }
    //************************************************************************************************************************************//
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
        FMatrix Ki(0.0), Kj(0.0);
        Ki = this->problem.soil().KFracture(global_i,e,local_i);
        Kj = this->problem.soil().KFracture(global_j,e,local_j);
        const FMatrix K = harmonicMeanK(Ki, Kj);

        // harmonic mean of the heat conductivity lambda
        RT lambda;
        lambda = 2./((1./vNDat[i].lambda) + (1./vNDat[j].lambda));

        // Arithmetic Mean:

        // calculate tortuosity at the nodes i and j needed for porous media diffusion coefficient
        RT tauW_i, tauW_j, tauN_i, tauN_j; // tortuosity of wetting and nonwetting phase
        tauW_i = pow(sNDat[globalIdx_i].porosityFracture * vNDat[i].satWFracture,(7/3))/
            (sNDat[globalIdx_i].porosityFracture*sNDat[globalIdx_i].porosityFracture);
        tauW_j = pow(sNDat[globalIdx_j].porosityFracture * vNDat[j].satWFracture,(7/3))/
            (sNDat[globalIdx_j].porosityFracture*sNDat[globalIdx_j].porosityFracture);
        tauN_i = pow(sNDat[globalIdx_i].porosityFracture * vNDat[i].satNFracture,(7/3))/
            (sNDat[globalIdx_i].porosityFracture*sNDat[globalIdx_i].porosityFracture);
        tauN_j = pow(sNDat[globalIdx_j].porosityFracture * vNDat[j].satNFracture,(7/3))/
            (sNDat[globalIdx_j].porosityFracture*sNDat[globalIdx_j].porosityFracture);

        // arithmetic mean of porous media diffusion coefficient
        RT Dwg, Daw;
        Dwg = (sNDat[globalIdx_i].porosityFracture * vNDat[i].satNFracture * tauN_i * vNDat[i].diff[nPhase] +
               sNDat[globalIdx_j].porosityFracture * vNDat[j].satNFracture * tauN_j * vNDat[j].diff[nPhase])/2;
        Daw = (sNDat[globalIdx_i].porosityFracture * vNDat[i].satWFracture * tauW_i * vNDat[i].diff[wPhase] +
               sNDat[globalIdx_j].porosityFracture * vNDat[j].satWFracture * tauW_j * vNDat[j].diff[wPhase])/2;

        // arithmetic mean of phase enthalpies (dissolved components neglected)
        RT enthW;
        RT enthCO2;
        enthW = (vNDat[i].enthalpy[pWIdx] + vNDat[j].enthalpy[pWIdx]) / 2.;
        enthCO2 = (vNDat[i].enthalpy[switchIdx] + vNDat[j].enthalpy[switchIdx]) / 2.;

        // Calculate arithmetic mean of the densities
        VBlockType avgDensity;
        avgDensity[wPhase] = 0.5*(vNDat[i].densityF[wPhase] + vNDat[j].densityF[wPhase]);
        avgDensity[nPhase] = 0.5*(vNDat[i].densityF[nPhase] + vNDat[j].densityF[nPhase]);

        //////////////////////////////////////////////////////////////////////////////////////////
        // GRADIENTS ----- in the fracture domain F

        FieldMatrix<RT,2,dim> pGradF(0.0), xGradF(0.0);
        FieldVector<RT,dim> teGradF(0.0);
        FieldVector<RT,dim> tempF(0.0);
        VBlockType flux(0.0);

        // calculate FE gradient at subcontrolvolumeface
        for (int k = 0; k < this->fvGeom.numVertices; k++) // loop over adjacent nodes
        {
            // FEGradient at subcontrolvolumeface face
            const FieldVector<DT,dim> feGrad(this->fvGeom.subContVolFace[face].grad[k]);
            FieldVector<RT,2> pressureF(0.0);
            FieldVector<RT,2> pressureM(0.0);

            pressureF[wPhase] = vNDat[k].pWFracture;
            pressureF[nPhase] = vNDat[k].pNFracture;
            pressureM[wPhase] = vNDat[k].pWMatrix;
            pressureM[nPhase] = vNDat[k].pNMatrix;



            // compute pressure gradients for each phase at integration point of subcontrolvolumeface face
            for (int phase = 0; phase < 2; phase++)
            {
                tempF = feGrad;
                tempF *= pressureF[phase];
                pGradF[phase] += tempF;
            }

            // compute temperature gradient
            tempF = feGrad;
            tempF *= vNDat[k].temperatureF;
            teGradF += tempF;

            // compute concentration gradients
            // for diffusion of co2 in wetting phase
            tempF = feGrad;
            tempF *= vNDat[k].massfrac[co2F][wPhase];
            xGradF[wPhase] += tempF;

            // for diffusion of water in nonwetting phase
            tempF = feGrad;
            tempF *= vNDat[k].massfrac[waterF][nPhase];
            xGradF[nPhase] += tempF;
        }

        // deduce gravity*density of each phase
        FieldMatrix<RT,2,dim> contribComp(0);
        for (int phase=0; phase<2; phase++)
        {
            contribComp[phase] = problem.gravity();
            contribComp[phase] *= vNDat[i].densityF[phase];
            pGradF[phase] -= contribComp[phase]; // grad p - rho*g
        }

        // Darcy velocity in normal direction for each phase K*n(grad p -rho*g)
        VBlockType outward(0);
        FieldVector<RT,dim> v_tilde_w(0);
        FieldVector<RT,dim> v_tilde_n(0);

        K.umv(pGradF[wPhase], v_tilde_w);  // v_tilde=K*gradP
        outward[wPhase] = v_tilde_w*normal;
        K.umv(pGradF[nPhase], v_tilde_n);  // v_tilde=K*gradP
        outward[nPhase] = v_tilde_n*normal;

        // Heat conduction
        outward[heatF] = teGradF * normal;
        outward[heatF] *= lambda;

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
        flux[waterF] =   (alpha* vNDat[up_w].densityF[wPhase]*vNDat[up_w].mobilityF[wPhase]
                          * vNDat[up_w].massfrac[waterF][wPhase]
                          + (1-alpha)* vNDat[dn_w].densityF[wPhase]*vNDat[dn_w].mobilityF[wPhase]
                          * vNDat[dn_w].massfrac[waterF][wPhase])
            * outward[wPhase];
        flux[waterF] +=  (alpha* vNDat[up_n].densityF[nPhase]*vNDat[up_n].mobilityF[nPhase]
                          * vNDat[up_n].massfrac[waterF][nPhase]
                          + (1-alpha)* vNDat[dn_n].densityF[nPhase]*vNDat[dn_n].mobilityF[nPhase]
                          * vNDat[dn_n].massfrac[waterF][nPhase])
            * outward[nPhase];
        // CO2 conservation
        flux[co2F]   =   (alpha* vNDat[up_n].densityF[nPhase]*vNDat[up_n].mobilityF[nPhase]
                          * vNDat[up_n].massfrac[co2F][nPhase]
                          + (1-alpha)* vNDat[dn_n].densityF[nPhase]*vNDat[dn_n].mobilityF[nPhase]
                          * vNDat[dn_n].massfrac[co2F][nPhase])
            * outward[nPhase];
        flux[co2F]  +=   (alpha* vNDat[up_w].densityF[wPhase]*vNDat[up_w].mobilityF[wPhase]
                          * vNDat[up_w].massfrac[co2F][wPhase]
                          + (1-alpha)* vNDat[dn_w].densityF[wPhase]*vNDat[dn_w].mobilityF[wPhase]
                          * vNDat[dn_w].massfrac[co2F][wPhase])
            * outward[wPhase];
        // Heat conservation
        flux[heatF]  =  (alpha* vNDat[up_n].densityF[nPhase]*vNDat[up_n].mobilityF[nPhase] * vNDat[up_n].enthalpy[nPhase]
                         + (1-alpha)* vNDat[dn_n].densityF[nPhase]*vNDat[dn_n].mobilityF[nPhase] * vNDat[dn_n].enthalpy[nPhase])
            * outward[nPhase];
        flux[heatF]  +=  (alpha* vNDat[up_w].densityF[wPhase]*vNDat[up_w].mobilityF[wPhase] * vNDat[up_w].enthalpy[wPhase]
                          + (1-alpha)* vNDat[dn_w].densityF[wPhase]*vNDat[dn_w].mobilityF[wPhase] * vNDat[dn_w].enthalpy[wPhase])
            * outward[wPhase];

        /////////////////////////////////////////////////////////////////////////////////////////////////////////
        // DIFFUSIVE TRANSPORT

        VBlockType normDiffGrad;

        // get local to global id map
        int state_i = sNDat[globalIdx_i].phaseState;
        int state_j = sNDat[globalIdx_j].phaseState;

        RT diffusionAW(0.0), diffusionWW(0.0), diffusionWN(0.0), diffusionAN(0.0);

        normDiffGrad[wPhase] = xGradF[wPhase]*normal;
        normDiffGrad[nPhase] = xGradF[nPhase]*normal;

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
        flux[waterF] += (diffusionWW + diffusionWN);

        // CO2 conservation
        flux[co2F] += (diffusionAN + diffusionAW);

        // Heat conservation
        flux[heatF] += diffusionWW*enthW + diffusionAW*enthCO2;

        //////////////////////////////////////////////////////////////////////////////////////////////
        // HEAT CONDUCTION

        flux[heatF]    +=    outward[heatF];

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

    ///**************************************************************************************************************************///
    /** @brief integrate sources / sinks
     *  @param e entity
     *  @param sol solution vector
     *  @param node local node id
     *  @return source/sink term
     */
    /////////////////////////////////////////////////////////////////////////////////////////
    //    INTERPOROSITY FLUX

    virtual VBlockType computeQ (const Entity& e, const VBlockType* sol, const int& node)
    {
        // ASSUME problem.q already contains \rho.q
        VBlockType q(0.0);
        q = problem.q(this->fvGeom.subContVol[node].global, e, this->fvGeom.subContVol[node].local);
        const FieldVector<DT,dim> global_i=this->fvGeom.subContVol[node].global;
        const FieldVector<DT,dim> local_i = this->fvGeom.subContVol[node].local;

        FieldVector <RT, nCont> A(0.0);
        A = AreaP(e, node);
        FieldMatrix<RT,dim, dim> KF(0.);
        FieldMatrix<RT,dim, dim> KM(0.);
        KF =  this->problem.soil().KFracture(global_i,e,local_i);
        KM =  this->problem.soil().KMatrix(global_i,e,local_i);
        const FMatrix Kharmonic = harmonicMeanKMinc(KF,KM);
        int i = node;
        FieldVector<RT, nCont> KA_nestVol(0.0);
        ////  Product of permeability with the nestedvolume area
        for (int nC=0; nC<nCont; nC++)
        {
            KA_nestVol[nC] = Kharmonic[0][0] * A[nC];
        }

        /***********************************************************************************/
        //// MINC Flux Matrix
        FieldMatrix<RT, nPhases, nCont> IPFlux(0.0);
        FieldMatrix<RT, nPhases, nCont> PGrad(0.0);
        FieldVector<RT, nCont> inward(0.0);
        //// Here should be added the normal and what sign it gives

        for (int nC = 0; nC < nCont-1; nC++)
        {
            PGrad[wPhase][nC] += vNDat[i].p[wPhase][nC];
            PGrad[wPhase][nC] -= vNDat[i].p[wPhase][nC+1];
            PGrad[nPhase][nC] += vNDat[i].p[nPhase][nC];
            PGrad[nPhase][nC] -= vNDat[i].p[nPhase][nC+1];

            ////       //if the normal of the face gives negative numbers it multiplies with -1
            inward[wPhase] = PGrad[wPhase][nC] * KA_nestVol[nC] / elData.delta;
            if (inward[wPhase] < 0){
                inward[wPhase]*=(-1);
            }
            inward[nPhase] = PGrad[wPhase][nC] * KA_nestVol[nC] / elData.delta;
            if (inward[nPhase] < 0){
                inward[nPhase]*=(-1);
            }
            //*******************************************************************//
            ////      the position of the node
            //      for (int node = 0; node < numVertices; node++) {
            //                          fvGeom.subContVol[node].local  = referenceElement.position(node, dim);
            //                          subContVol[node].global = geometry.global(subContVol[node].local);
            //                      }
            //*******************************************************************//

            //// If the pressure in the fracture is bigger than the one in the matrix the flow occurs from fracture to matrix
            if (vNDat[i].p[wPhase][nC]>vNDat[i].p[wPhase][nC+1])
            {
                // wetting interporosity flux calculated with the wetting mobility of the Fracture;
                IPFlux[wPhase][nC]=vNDat[i].density[wPhase][nC]*vNDat[i].mobility[wPhase][nC]*inward[wPhase];
                q[3*nC]-=IPFlux[wPhase][nC]; //pWFIdx
                q[3*(nC+1)]+=IPFlux[wPhase][nC];//pWMIdx
            }
            else
            {
                // wetting interporosity flux calculated with the wetting mobility of the Matrix;
                IPFlux[wPhase][nC]=vNDat[i].density[wPhase][nC+1]*vNDat[i].mobility[wPhase][nC+1]*inward[wPhase];
                q[3*nC]+=IPFlux[wPhase][nC];//pWFIdx
                q[3*(nC+1)]-=IPFlux[wPhase][nC];//pWMIdx
            }
            if (vNDat[i].p[nPhase][nC]>vNDat[i].p[nPhase][nC+1])
            {
                IPFlux[nPhase][nC]=vNDat[i].density[nPhase][nC]*vNDat[i].mobility[nPhase][nC]*inward[nPhase];
                q[3*nC+1]-=IPFlux[nPhase][0];//satNFIdx
                q[3*(nC+1)+1]+=IPFlux[nPhase][0];//satNMIdx
            }
            else
            {
                IPFlux[nPhase][nC]=vNDat[i].density[nPhase][nC+1]*vNDat[i].mobility[nPhase][nC+1]*inward[nPhase];
                q[3*nC+1]+=IPFlux[nPhase][nC];//satNFIdx
                q[3*(nC+1)+1]-=IPFlux[nPhase][nC];//satNMIdx
            }
        }

        return q;
    }





    /** @brief perform variable switch
     *  @param global global node id
     *  @param sol solution vector
     *  @param local local node id
     */
    virtual void primaryVarSwitch (const Entity& e, int global, VBlockType* sol, int local)
    {

        //       if(!switchBreak)
        //        {
        switched = false;
        int state = sNDat[global].phaseState;

        //        RT pW = sol[local][pWIdx];
        RT pWFracture = sol[local][pWIdx];
        RT pWMatrix = sol[local][pWIdx];
        RT satWFracture = 0.0;
        RT satWMatrix = 0.0;
        if (state == bothPhases) satWFracture = 1.0-sol[local][switchIdx];
        if (state == waterPhase) satWFracture = 1.0;
        if (state == gasPhase)      satWFracture = 0.0;

        if (state == bothPhases) satWMatrix = 1.0-sol[local][switchIdx];
        if (state == waterPhase) satWMatrix = 1.0;
        if (state == gasPhase)      satWMatrix = 0.0;


        // capillary pressure neglected
        RT pCFracture = 0.0; //problem.materialLaw().pC(satW, this->fvGeom.subContVol[local].global, e, this->fvGeom.subContVol[local].local);
        RT pCMatrix   = 0.0; //problem.materialLaw().pC(satW, this->fvGeom.subContVol[local].global, e, this->fvGeom.subContVol[local].local);
        RT pNFracture = pWFracture + pCFracture;
        RT pNMatrix   = pWMatrix + pCMatrix;

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
                sol[local][switchIdx] = 1.0 - 1.e-6; // initialize solution vector
                switched = true;
            }
            break;

        case waterPhase :
            RT xAWmax, xAWmass;
            xAWmass = sol[local][switchIdx];
            xAWmax = problem.multicomp().xAW(pNFracture, sol[local][teIdx]);

            if (xAWmass > xAWmax && switched == false)
            {
                // appearance of gas phase
                std::cout << "Gas appears at node " << global << "  Coordinates: " << Coordinates << std::endl;
                sNDat[global].phaseState = bothPhases;
                sol[local][switchIdx] = 1.e-6; // initialize solution vector
                switched = true;
            }
            break;

        case bothPhases :
            RT satNFracture = sol[local][switchIdx];
            if (satNFracture < 0.0  && switched == false)
            {
                // disappearance of gas phase
                std::cout << "Gas disappears at node " << global << "  Coordinates: " << Coordinates << std::endl;
                sNDat[global].phaseState = waterPhase;
                sol[local][switchIdx] = 1e-6; // initialize solution vector
                switched = true;
            }
            else if (satWFracture < 0.0  && switched == false)
            {
                // disappearance of water phase
                std::cout << "Water disappears at node " << global << "  Coordinates: " << Coordinates << std::endl;
                sNDat[global].phaseState = gasPhase;
                sol[local][switchIdx] = 1e-6; // initialize solution vector
                switched = true;
            }
            break;
        }
        //        }
        return;
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
        elData.heatCap = problem.soil().heatCap(this->fvGeom.elementGlobal, e, this->fvGeom.elementLocal);

        //       // ASSUME element-wise constant parameters for the material law
        //          elData.parameters = problem.materialLawParameters
        //          (this->fvGeom.elementGlobal, e, this->fvGeom.elementLocal);
        //
        //         // ASSUMING element-wise constant permeability, evaluate K at the cell center
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
                sNDat[globalIdx].porosityFracture = problem.soil().porosity(this->fvGeom.subContVol[k].global, e, this->fvGeom.subContVol[k].local);
                sNDat[globalIdx].porosityMatrix = problem.soil().porosity(this->fvGeom.subContVol[k].global, e, this->fvGeom.subContVol[k].local);

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
                sNDat[globalIdx].porosityFracture = problem.soil().porosity(this->fvGeom.elementGlobal, e, this->fvGeom.elementLocal);

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
        RT satNFracture;
        RT satNMatrix;
        RT satWFracture;
        RT satWMatrix;
        RT pWFracture;
        RT pWMatrix;
        RT pCFracture;
        RT pCMatrix;
        RT pNFracture;
        RT pNMatrix;
        RT temperatureF;
        RT lambda;
        RT viscosityCO2;
        RT krCO2;
        FieldVector<RT,2> mobilityF;  //Vector with the number of phases
        FieldVector<RT,2> densityF;
        FieldMatrix<RT,c,2> massfrac;
        FieldVector<RT,2> enthalpy;
        FieldVector<RT,2> intenergy;
        FieldVector<RT,2> diff;

        FieldVector<RT, nCont> pC;                    //capillary pressure
        FieldVector<RT, nCont> temperature;         //temperature (minc)

        FieldMatrix<RT, nPhases, nCont> p;
        FieldMatrix<RT, nPhases, nCont> saturation;
        FieldMatrix<RT, nPhases, nCont> density;     //density in the (M)atrix form
        FieldMatrix<RT, nPhases, nCont> mobility;     //density in the (M)atrix form


        int phasestate;
    };

    // analog to EvalPrimaryData in MUFTE, uses members of vNDat
    virtual void updateVariableData(const Entity& e, const VBlockType* sol,
                                    int i, std::vector<VariableNodeData>& varData, int state)
    {
        ///////////////////////////////////////////////////////////////////////////////

        const int global = this->vertexMapper.template map<dim>(e, i);

        for (int nC=0; nC<nCont; nC++){

            varData[i].p[wPhase][nC] = sol[i][pWIdx];

            if (state == bothPhases) varData[i].saturation[nPhase][nC] = sol[i][switchIdx];    //S nPhase F
            if (state == waterPhase) varData[i].saturation[nPhase][nC] = 0.0;
            if (state == gasPhase)   varData[i].saturation[nPhase][nC] = 1.0;

            varData[i].saturation[wPhase][nC]  = 1.0 - varData[i].saturation[nPhase][nC];        //S wPhase F
            //Capillary pressure
            varData[i].pC[nC] = 0;
            varData[i].p[nPhase][nC] = varData[i].p[wPhase][nC]+varData[i].pC[nC];
            varData[i].temperature[nC] = sol[i][teIdx];

        }




        ///////////////////////////////////////////////////////////////////////////////

        varData[i].pWFracture = sol[i][pWIdx];
        if (state == bothPhases) varData[i].satNFracture = sol[i][switchIdx];
        if (state == waterPhase) varData[i].satNFracture = 0.0;
        if (state == gasPhase)   varData[i].satNFracture = 1.0;

        varData[i].satWFracture = 1.0 - varData[i].satNFracture;

        // Capillary pressure neglected (may produce negative pw values)
        varData[i].pCFracture = 0.0; //problem.materialLaw().pC(varData[i].satW, this->fvGeom.subContVol[i].global, e, this->fvGeom.subContVol[i].local);
        varData[i].pNFracture = varData[i].pWFracture + varData[i].pCFracture;
        varData[i].temperatureF = sol[i][teIdx]; // in [K]

        // Solubilities of components in phases
        if (state == bothPhases){
            varData[i].massfrac[co2F][wPhase] = problem.multicomp().xAW(varData[i].pNFracture, varData[i].temperatureF);
            varData[i].massfrac[waterF][nPhase] = problem.multicomp().xWN(varData[i].pNFracture, varData[i].temperatureF);
        }
        if (state == waterPhase){
            varData[i].massfrac[waterF][nPhase] = 0.0;
            varData[i].massfrac[co2F][wPhase] =  sol[i][switchIdx];
        }
        if (state == gasPhase){
            varData[i].massfrac[waterF][nPhase] = sol[i][switchIdx];
            varData[i].massfrac[co2F][wPhase] = 0.0;
        }
        varData[i].massfrac[waterF][wPhase] = 1.0 - varData[i].massfrac[co2F][wPhase];
        varData[i].massfrac[co2F][nPhase] = 1.0 - varData[i].massfrac[waterF][nPhase];
        // for output
        //               (*hackySaturationN)[global] = varData[i].satN;
        //               (*hackyMassFracCO2)[global] = varData[i].massfrac[co2][wPhase];
        //               (*hackyMassFracWater)[global] = varData[i].massfrac[water][nPhase];
        // Diffusion coefficients
        // Mobilities & densities

        varData[i].densityF[wPhase] = problem.wettingPhase().density(varData[i].temperatureF, varData[i].pWFracture, varData[i].massfrac[co2F][wPhase]);
        varData[i].densityF[nPhase] = problem.nonwettingPhase().density(varData[i].temperatureF, varData[i].pNFracture);
        varData[i].mobilityF[wPhase] = problem.materialLaw().mobW(varData[i].satWFracture, this->fvGeom.subContVol[i].global, e, this->fvGeom.subContVol[i].local, varData[i].temperatureF, varData[i].pWFracture);
        varData[i].lambda = problem.soil().heatCond(this->fvGeom.subContVol[i].global, e, this->fvGeom.subContVol[i].local, varData[i].satWFracture);
        varData[i].enthalpy[wPhase] = problem.wettingPhase().enthalpy(varData[i].temperatureF,varData[i].pWFracture);
        varData[i].enthalpy[nPhase] = problem.nonwettingPhase().enthalpy(varData[i].temperatureF,varData[i].pNFracture);
        varData[i].intenergy[wPhase] = problem.wettingPhase().intEnergy(varData[i].temperatureF,varData[i].pWFracture);
        varData[i].intenergy[nPhase] = problem.nonwettingPhase().intEnergy(varData[i].temperatureF,varData[i].pNFracture);
        varData[i].diff[wPhase] = problem.wettingPhase().diffCoeff();
        varData[i].diff[nPhase] = problem.nonwettingPhase().diffCoeff();
        // calculate mobility of CO2
        varData[i].viscosityCO2 = problem.nonwettingPhase().viscosityCO2(varData[i].temperatureF, varData[i].pNFracture, varData[i].densityF[nPhase]);
        varData[i].krCO2 = problem.materialLaw().krn(varData[i].satNFracture, this->fvGeom.subContVol[i].global, e, this->fvGeom.subContVol[i].local);
        varData[i].mobilityF[nPhase] = varData[i].krCO2/varData[i].viscosityCO2;

        for (int nC = 0; nC<nCont; nC++){
            varData[i].density[wPhase][nC] = problem.wettingPhase().density(varData[i].temperatureF, varData[i].pWFracture, varData[i].massfrac[co2F][wPhase]);
            varData[i].density[nPhase][nC] = problem.nonwettingPhase().density(varData[i].temperatureF, varData[i].pNFracture);
            varData[i].mobility[wPhase][nC] = problem.materialLaw().mobW(varData[i].satWFracture, this->fvGeom.subContVol[i].global, e, this->fvGeom.subContVol[i].local, varData[i].temperatureF, varData[i].pWFracture);
        }


        // CONSTANT solubility (for comparison with twophase)
        //         varData[i].massfrac[co2][wPhase] = 0.0; varData[i].massfrac[water][wPhase] = 1.0;
        //         varData[i].massfrac[water][nPhase] = 0.0; varData[i].massfrac[co2][nPhase] = 1.0;

        //std::cout << "water in gasphase: " << varData[i].massfrac[water][nPhase] << std::endl;
        //std::cout << "co2 in waterphase: " << varData[i].massfrac[co2][wPhase] << std::endl;

        // for output
        (*outPressureN)[global] = varData[i].pNFracture;

        (*outCapillaryP)[global] = varData[i].pCFracture;
        (*outSaturationW)[global] = varData[i].satWFracture;
        (*outSaturationN)[global] = varData[i].satNFracture;
        (*outTemperature)[global] = varData[i].temperatureF;
        (*outMassFracAir)[global] = varData[i].massfrac[co2F][wPhase];
        (*outMassFracWater)[global] = varData[i].massfrac[waterF][nPhase];
        (*outDensityW)[global] = varData[i].densityF[wPhase];
        (*outDensityN)[global] = varData[i].densityF[nPhase];
        (*outMobilityW)[global] = varData[i].mobilityF[wPhase];
        (*outMobilityN)[global] = varData[i].mobilityF[nPhase];
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
        bool visited;
        int phaseState;
        int oldPhaseState;
        RT cellVolume;
        RT porosity;
        RT porosityFracture;
        RT porosityMatrix;


        FieldVector<RT, 4> parameters;
        FMatrix K;
    };

    struct StaticIPData
    {
        bool visited;
        FMatrix K;
        FMatrix KFracture, KMatrix;
    };


    struct ElementData {
        RT heatCap;
        double delta; //distance between nested volume faces (kept constant - so all faces are equally distanced)
        //        RT cellVolume;
        //          RT porosity;
        //        RT gravity;
        //        FieldVector<RT, 4> parameters;
        //        FieldMatrix<RT,dim,dim> K;
    } elData;


    // parameters given in constructor
    TwoPTwoCNIProblem<G,RT,m>& problem;
    CBrineCO2 multicomp;
    std::vector<StaticNodeData> sNDat;
    std::vector<StaticIPData> sIPDat;
    std::vector<VariableNodeData> vNDat;
    std::vector<VariableNodeData> oldVNDat;
    bool switched, switchBreak;

    // for output files
    BlockVector<FieldVector<RT, 1> > *outPressureN;
    BlockVector<FieldVector<RT, 1> > *outPressureNFracture;
    BlockVector<FieldVector<RT, 1> > *outPressureNMatrix;
    BlockVector<FieldVector<RT, 1> > *outCapillaryP;
    BlockVector<FieldVector<RT, 1> > *outCapillaryPFracture;
    BlockVector<FieldVector<RT, 1> > *outCapillaryPMatrix;
    BlockVector<FieldVector<RT, 1> > *outSaturationN;
    BlockVector<FieldVector<RT, 1> > *outSaturationNFracture;
    BlockVector<FieldVector<RT, 1> > *outSaturationNMatrix;
    BlockVector<FieldVector<RT, 1> > *outSaturationW;
    BlockVector<FieldVector<RT, 1> > *outSaturationWFracture;
    BlockVector<FieldVector<RT, 1> > *outSaturationWMatrix;
    BlockVector<FieldVector<RT, 1> > *outTemperature;
    BlockVector<FieldVector<RT, 1> > *outTemperatureFracture;
    BlockVector<FieldVector<RT, 1> > *outTemperatureMatrix;
    BlockVector<FieldVector<RT, 1> > *outMassFracAir;
    BlockVector<FieldVector<RT, 1> > *outMassFracAirFracture;
    BlockVector<FieldVector<RT, 1> > *outMassFracAirMatrix;
    BlockVector<FieldVector<RT, 1> > *outMassFracWater;
    BlockVector<FieldVector<RT, 1> > *outDensityW;
    BlockVector<FieldVector<RT, 1> > *outDensityN;
    BlockVector<FieldVector<RT, 1> > *outMobilityW;
    BlockVector<FieldVector<RT, 1> > *outMobilityN;
    BlockVector<FieldVector<RT, 1> > *outPhaseState;

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

