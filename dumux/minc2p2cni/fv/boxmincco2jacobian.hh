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

namespace Dune {
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
template<class Grid, class Scalar, int numEq,
         class BoxFunction =LeafP1Function<Grid, Scalar, numEq> > class BoxMincCO2Jacobian :
        public BoxJacobian<BoxMincCO2Jacobian<Grid,Scalar,numEq,BoxFunction>,Grid,Scalar,numEq,BoxFunction> {
    typedef typename Grid::ctype DT;
    typedef typename Grid::Traits::template Codim<0>::Entity Entity;
    typedef typename Entity::Geometry Geometry;
    typedef BoxMincCO2Jacobian<Grid,Scalar,numEq,BoxFunction> ThisType;
    typedef typename LocalJacobian<ThisType,Grid,Scalar,numEq>::VBlockType VBlockType;
    typedef typename LocalJacobian<ThisType,Grid,Scalar,numEq>::MBlockType MBlockType;

    enum {pWIdx = 0, switchIdx = 1, teIdx=2, numberOfComponents = 2}; // Solution vector index
    enum {fracture = 0};
    enum {wPhase = 0, nPhase = 1}; // Phase index
    enum {gasPhase = 0, waterPhase = 1, bothPhases = 2}; // Phase state
    enum {waterF = 0, co2F = 1, heatF = 2}; // Component index
    enum {water = 0, co2 = 1, heat = 2}; // Component index
    enum {numPhases = 2}; // Number of phases
    enum {nPorosities = 2}; // Number of porosities
    enum {nPermeabilities =2}; // Number of permeabilities

public:
    // define the number of phases (numEq) and components (c) of your system, this is used outside
    // to allocate the correct size of (dense) blocks with a FieldMatrix
    enum {dim=Grid::dimension};
    enum {c=2};
    //number of interacting continua
    enum {numCont = numEq/3};
    enum {SIZE=LagrangeShapeFunctionSetContainer<DT,Scalar,dim>::maxsize};
    struct VariableNodeData;

    typedef FieldMatrix<Scalar,dim,dim> FMatrix;
    typedef FieldVector<Scalar,dim> FVector;

    //! Constructor
    BoxMincCO2Jacobian(TwoPTwoCNIProblem<Grid,Scalar,numEq>& params,
                       bool levelBoundaryAsDirichlet_, const Grid& grid, BoxFunction& sol,
                       bool procBoundaryAsDirichlet_=true) :
        BoxJacobian<ThisType, Grid, Scalar, numEq, BoxFunction>(levelBoundaryAsDirichlet_,
                                                                grid, sol, procBoundaryAsDirichlet_), problem(params),
        sNDat(this->vertexMapper.size()), vNDat(SIZE), oldVNDat(SIZE),
        switched(false), switchBreak(false) {
        this->analytic = false;
    }

    virtual FieldVector<Scalar, numCont> AreaP(const Entity& element, int face) {
        FieldVector<Scalar,numCont> x(0.0);
        FieldVector<Scalar,numCont> y(0.0);
        FieldVector<Scalar,numCont> L(0.0);
        double length_x =0; //the length of the subcontrol face on x direction
        double length_y =0; //the length of the subcontrol face on y direction
        double ExtLength =0; //Exterior length of the subcontrol face
        for (int face =0; face < this->fvGeom.numVertices; face++) {
            int it = this->fvGeom.subContVolFace[face].i;
            int jt = this->fvGeom.subContVolFace[face].j;
            if (it == 0) {
                FieldVector<DT,dim> ht;
                ht = this->fvGeom.subContVolFace[face].normal;
                if (ht[0]<0) {
                    ht[0]*=(-1);
                }
                if (ht[1]<0) {
                    ht[1]*=(-1);
                }
                if (ht[0]==0) {
                    ht[0]=ht[1];
                }
                if (ht[1]==0) {
                    ht[1]=ht[0];
                }
                length_y = ht[0];
            } else if (jt==0) {
                FieldVector<DT,dim> ht;
                ht = this->fvGeom.subContVolFace[face].normal;
                if (ht[0]<0) {
                    ht[0]*=(-1);
                }
                if (ht[1]<0) {
                    ht[1]*=(-1);
                }
                if (ht[0]==0) {
                    ht[0]=ht[1];
                }
                if (ht[1]==0) {
                    ht[1]=ht[0];
                }
                length_x = ht[0];
            }
            ExtLength = length_x+length_y;

        }
        //              double f = 1.0 / numCont; // volume fraction
        //              // Calculation of the distance between the nested volume elements delta
        //              // delta in this case divides the surfaces equidistantly (delta is const)
        elData.delta = 0.0;
        if (length_x > length_y) {
            elData.delta = 2*length_y / (2.0 * (numCont-1) + 1.0);
        } else {
            elData.delta = 2*length_x / (2.0 * (numCont-1) + 1.0);
        }

        double TotalLength = 4*ExtLength;

        L[0] = TotalLength;
        x[0] = 2*length_x; //length of the unitary element (CV) on x direction
        y[0] = 2*length_y; //length of the unitary element (CV) on y direction
        for (int nC=1; nC < numCont; nC++) {
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
    virtual FMatrix harmonicMeanKMinc(FMatrix& KF, const FMatrix& KM) const {
        double eps = 1e-20;
        FieldMatrix<Scalar,dim, dim> K(0.0);
        for (int kx=0; kx<nPermeabilities; kx++) {
            for (int ky=0; ky<nPermeabilities; ky++) {
                if (KF[kx][ky] != KM[kx][ky]) {
                    K[kx][ky] = 2 / (1/(KF[kx][ky]+eps) + (1/(KM[kx][ky]+eps)));
                } else {
                    K = KF;
                }
            }
        }
        return K;
    }
    //    /*the harmonicMeanK function computes the harmonic mean of the perameabilities between the two nodes of different
    //     *  permeabilities in the fracture domain */


    // harmonic mean computed directly
    virtual FMatrix harmonicMeanK(FMatrix& Ki, const FMatrix& Kj) const {
        double eps = 1e-20;
        FMatrix K(0.);
        for (int kx=0; kx<nPermeabilities; kx++) {
            for (int ky=0; ky<nPermeabilities; ky++) {
                if (Ki[kx][ky] != Kj[kx][ky]) {
                    K[kx][ky] = 2 / (1/(Ki[kx][ky]+eps) + (1/(Kj[kx][ky]+eps)));
                } else
                    K = Ki;
            }
        }
        return K;
    }

    ///*************************************************************************************************************/


    /** @brief compute time dependent term (storage), loop over nodes / subcontrol volumes
     *  @param element entity
     *  @param sol solution vector
     *  @param node local node id
     *  @return storage term
     */
    virtual VBlockType computeM(const Entity& element, const VBlockType* sol,
                                int idx, std::vector<VariableNodeData>& varData) {
        GeometryType gt = element.geometry().type();
        const typename LagrangeShapeFunctionSetContainer<DT,Scalar,dim>::value_type
            & sfs=LagrangeShapeFunctions<DT, Scalar, dim>::general(gt, 1);

        int globalIdx = this->vertexMapper.template map<dim>(element,
                                                             sfs[idx].entity());

        VBlockType AccumulationTerm(0.);
        MBlockType result(0.);

        for (int nC = 0; nC < numCont; nC++) {
            result[water][nC]
                = sNDat[globalIdx].porosity[nC]
                *(varData[idx].density[wPhase][nC]
                  *varData[idx].saturation[water][nC]
                  *varData[idx].massfracMincW[water][nC]
                  +varData[idx].density[nPhase][nC]
                  *varData[idx].saturation[co2][nC]
                  *varData[idx].massfracMincN[water][nC]);
            result[co2][nC] = sNDat[globalIdx].porosity[nC]
                *(varData[idx].density[nPhase][nC]
                  *varData[idx].saturation[co2][nC]
                  *varData[idx].massfracMincN[co2F][nC]
                  +varData[idx].density[wPhase][nC]
                  *varData[idx].saturation[water][nC]
                  *varData[idx].massfracMincW[co2][nC]);
            result[heat][nC] = sNDat[globalIdx].porosity[nC]
                * (varData[idx].density[wPhase][nC]
                   * varData[idx].intenergy[wPhase][nC]
                   * varData[idx].saturation[water][nC]
                   + varData[idx].density[nPhase][nC]
                   * varData[idx].intenergy[nPhase][nC]
                   * varData[idx].saturation[co2][nC])
                + elData.heatCap * varData[idx].temperature[nC];
        }

        for (int nC = 0; nC < numCont; nC++) {
            //3 is the number of equations (water, CO2, heat) and nC - the number of Continua
            int count = 3*nC;
            AccumulationTerm[count]=result[water][nC];
            AccumulationTerm[count+1]=result[co2][nC];
            AccumulationTerm[count+2]=result[heat][nC];
        }

        return AccumulationTerm;

    };

    virtual VBlockType computeM(const Entity& element, const VBlockType* sol,
                                int idx, bool old = false) {
        if (old)
            return computeM(element, sol, idx, oldVNDat);
        else
            return computeM(element, sol, idx, vNDat);
    }
    //************************************************************************************************************************************//
    /** @brief compute diffusive/advective fluxes, loop over subcontrol volume faces
     *  @param element entity
     *  @param sol solution vector
     *  @param face face id
     *  @return flux term
     */
    virtual VBlockType computeA(const Entity& element, const VBlockType* sol, int face) {
        int i = this->fvGeom.subContVolFace[face].i;
        int j = this->fvGeom.subContVolFace[face].j;

        // normal vector, value of the area of the scvf
        const FieldVector<Scalar,dim>
            normal(this->fvGeom.subContVolFace[face].normal);
        GeometryType gt = element.geometry().type();
        const typename LagrangeShapeFunctionSetContainer<DT,Scalar,dim>::value_type
            & sfs=LagrangeShapeFunctions<DT, Scalar, dim>::general(gt, 1);

        // global index of the subcontrolvolume face neighbour nodes in element element
        int globalIdx_i = this->vertexMapper.template map<dim>(element,
                                                               sfs[i].entity());
        int globalIdx_j = this->vertexMapper.template map<dim>(element,
                                                               sfs[j].entity());

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
        Ki = this->problem.soil().KFracture(global_i, element, local_i);
        Kj = this->problem.soil().KFracture(global_j, element, local_j);
        const FMatrix K = harmonicMeanK(Ki, Kj);

        // harmonic mean of the heat conductivity lambda
        Scalar lambdaF;
        lambdaF = 2./((1./vNDat[i].lambda[fracture]) + (1./vNDat[j].lambda[fracture]));

        // Arithmetic Mean:

        // calculate tortuosity at the nodes i and j needed for porous media diffusion coefficient
        Scalar tauW_i, tauW_j, tauN_i, tauN_j; // tortuosity of wetting and nonwetting phase
        tauW_i = pow(sNDat[globalIdx_i].porosity[fracture]
                     * vNDat[i].saturation[wPhase][fracture], (7/3))
            / (sNDat[globalIdx_i].porosity[fracture]
               *sNDat[globalIdx_i].porosity[fracture]);
        tauW_j = pow(sNDat[globalIdx_j].porosity[fracture]
                     * vNDat[j].saturation[wPhase][fracture], (7/3))
            / (sNDat[globalIdx_j].porosity[fracture]
               *sNDat[globalIdx_j].porosity[fracture]);
        tauN_i = pow(sNDat[globalIdx_i].porosity[fracture]
                     * vNDat[i].saturation[nPhase][fracture], (7/3))
            / (sNDat[globalIdx_i].porosity[fracture]
               *sNDat[globalIdx_i].porosity[fracture]);
        tauN_j = pow(sNDat[globalIdx_j].porosity[fracture]
                     * vNDat[j].saturation[nPhase][fracture], (7/3))
            / (sNDat[globalIdx_j].porosity[fracture]
               *sNDat[globalIdx_j].porosity[fracture]);

        // arithmetic mean of porous media diffusion coefficient
        Scalar Dwg, Daw;
        Dwg = (sNDat[globalIdx_i].porosity[fracture] * vNDat[i].saturation[nPhase][fracture]
               * tauN_i * vNDat[i].diff[nPhase]
               + sNDat[globalIdx_j].porosity[fracture] * vNDat[j].saturation[nPhase][fracture]
               * tauN_j * vNDat[j].diff[nPhase])/2;
        Daw = (sNDat[globalIdx_i].porosity[fracture] * vNDat[i].saturation[wPhase][fracture]
               * tauW_i * vNDat[i].diff[wPhase]
               + sNDat[globalIdx_j].porosity[fracture] * vNDat[j].saturation[wPhase][fracture]
               * tauW_j * vNDat[j].diff[wPhase])/2;

        // arithmetic mean of phase enthalpies (dissolved components neglected)
        Scalar enthW;
        Scalar enthCO2;
        enthW = (vNDat[i].enthalpy[wPhase][fracture] + vNDat[j].enthalpy[wPhase][fracture]) / 2.;
        enthCO2 = (vNDat[i].enthalpy[nPhase][fracture] + vNDat[j].enthalpy[nPhase][fracture]) / 2.;
        // Calculate arithmetic mean of the densities
        VBlockType avgDensity;
        avgDensity[wPhase] = 0.5*(vNDat[i].density[wPhase][fracture] + vNDat[j].density[wPhase][fracture]);
        avgDensity[nPhase] = 0.5*(vNDat[i].density[nPhase][fracture] + vNDat[j].density[nPhase][fracture]);


        //////////////////////////////////////////////////////////////////////////////////////////
        // GRADIENTS ----- in the fracture domain F

        FieldMatrix<Scalar,2,numPhases> pGradF(0.0), xGradF(0.0);
        FieldVector<Scalar,dim> teGradF(0.0);
        FieldVector<Scalar,dim> tempF(0.0);
        VBlockType flux(0.0);

        // calculate FE gradient at subcontrolvolumeface
        for (int k = 0; k < this->fvGeom.numVertices; k++) // loop over adjacent nodes
        {
            // FEGradient at subcontrolvolumeface face
            const FieldVector<DT,dim> feGrad(this->fvGeom.subContVolFace[face].grad[k]);
            FieldVector<Scalar,numPhases> pressureF(0.0);

            pressureF[wPhase] = vNDat[k].p[wPhase][fracture];
            pressureF[nPhase] = vNDat[k].p[nPhase][fracture];

            // compute pressure gradients for each phase at integration point of subcontrolvolumeface face
            for (int phase = 0; phase < 2; phase++) {
                tempF = feGrad;
                tempF *= pressureF[phase];
                pGradF[phase] += tempF;
            }

            // compute temperature gradient
            tempF = feGrad;
            tempF *= vNDat[k].temperature[fracture];
            teGradF += tempF;

            // compute concentration gradients
            // for diffusion of co2 in wetting phase
            tempF = feGrad;
            tempF *= vNDat[k].massfracMincW[co2][fracture];
            xGradF[wPhase] += tempF;

            // for diffusion of water in nonwetting phase
            tempF = feGrad;
            tempF *= vNDat[k].massfracMincN[water][fracture];
            xGradF[nPhase] += tempF;
        }

        // deduce gravity*density of each phase
        FieldMatrix<Scalar,2,dim> contribComp(0);
        for (int phase=0; phase<2; phase++) {
            contribComp[phase] = problem.gravity();
            contribComp[phase] *= vNDat[i].density[phase][fracture];
            pGradF[phase] -= contribComp[phase]; // grad p - rho*g
        }

        // Darcy velocity in normal direction for each phase K*n(grad p -rho*g)
        VBlockType outward(0);
        FieldVector<Scalar,dim> v_tilde_w(0);
        FieldVector<Scalar,dim> v_tilde_n(0);

        K.umv(pGradF[wPhase], v_tilde_w); // v_tilde+=K*gradP
        outward[wPhase] = v_tilde_w*normal;
        K.umv(pGradF[nPhase], v_tilde_n); // v_tilde+=K*gradP
        outward[nPhase] = v_tilde_n*normal;

        // Heat conduction
        outward[heatF] = teGradF * normal;
        outward[heatF] *= lambdaF;

        // evaluate upwind nodes
        int up_w, dn_w, up_n, dn_n;
        if (outward[wPhase] <= 0) {
            up_w = i;
            dn_w = j;
        } else {
            up_w = j;
            dn_w = i;
        };
        if (outward[nPhase] <= 0) {
            up_n = i;
            dn_n = j;
        } else {
            up_n = j;
            dn_n = i;
        };

        Scalar alpha = 1.0; // Upwind parameter

        ////////////////////////////////////////////////////////////////////////////////////////////////
        // ADVECTIVE TRANSPORT
        // Water conservation
        flux[waterF] = (alpha* vNDat[up_w].density[wPhase][fracture]
                        * vNDat[up_w].mobility[wPhase][fracture]
                        * vNDat[up_w].massfracMincW[water][fracture] + (1-alpha)
                        * vNDat[dn_w].density[wPhase][fracture]*vNDat[dn_w].mobility[wPhase][fracture]
                        * vNDat[dn_w].massfracMincW[water][fracture]) * outward[wPhase];
        flux[waterF] += (alpha* vNDat[up_n].density[nPhase][fracture]
                         * vNDat[up_n].mobility[nPhase][fracture]
                         * vNDat[up_n].massfracMincN[water][fracture] + (1-alpha)
                         * vNDat[dn_n].density[nPhase][fracture]*vNDat[dn_n].mobility[nPhase][fracture]
                         * vNDat[dn_n].massfracMincN[water][fracture]) * outward[nPhase];
        // CO2 conservation
        flux[co2F] = (alpha* vNDat[up_n].density[nPhase][fracture]
                      * vNDat[up_n].mobility[nPhase][fracture]
                      * vNDat[up_n].massfracMincN[co2][fracture] + (1-alpha)
                      * vNDat[dn_n].density[nPhase][fracture]*vNDat[dn_n].mobility[nPhase][fracture]
                      * vNDat[dn_n].massfracMincN[co2][fracture]) * outward[nPhase];
        flux[co2F] += (alpha* vNDat[up_w].density[wPhase][fracture]
                       * vNDat[up_w].mobility[wPhase][fracture]
                       * vNDat[up_w].massfracMincW[co2][fracture] + (1-alpha)
                       * vNDat[dn_w].density[wPhase][fracture]*vNDat[dn_w].mobility[wPhase][fracture]
                       * vNDat[dn_w].massfracMincW[co2][fracture]) * outward[wPhase];

        // Heat conservation
        flux[heatF] = (alpha* vNDat[up_n].density[nPhase][fracture]
                       * vNDat[up_n].mobility[nPhase][fracture] * vNDat[up_n].enthalpy[nPhase][fracture]
                       + (1-alpha)* vNDat[dn_n].density[nPhase][fracture]
                       * vNDat[dn_n].mobility[nPhase][fracture]
                       * vNDat[dn_n].enthalpy[nPhase][fracture]) * outward[nPhase];
        flux[heatF] += (alpha* vNDat[up_w].density[wPhase][fracture]
                        * vNDat[up_w].mobility[wPhase][fracture] * vNDat[up_w].enthalpy[wPhase][fracture]
                        + (1-alpha)* vNDat[dn_w].density[wPhase][fracture]
                        *vNDat[dn_w].mobility[wPhase][fracture]
                        * vNDat[dn_w].enthalpy[wPhase][fracture]) * outward[wPhase];

        /////////////////////////////////////////////////////////////////////////////////////////////////////////
        // DIFFUSIVE TRANSPORT

        VBlockType normDiffGrad;

        // get local to global id map
        int state_i = sNDat[globalIdx_i].phaseStateF;
        int state_j = sNDat[globalIdx_j].phaseStateF;

        Scalar diffusionAW(0.0), diffusionWW(0.0), diffusionWN(0.0), diffusionAN(0.0);

        normDiffGrad[wPhase] = xGradF[wPhase]*normal;
        normDiffGrad[nPhase] = xGradF[nPhase]*normal;

        if (state_i == bothPhases && state_j == bothPhases) {
            diffusionAW = Daw * avgDensity[wPhase] * normDiffGrad[wPhase];
            diffusionWW = -diffusionAW;
            diffusionWN = Dwg * avgDensity[nPhase] * normDiffGrad[nPhase];
            diffusionAN = -diffusionWN;
        } else if (state_i == waterPhase || state_j == waterPhase) {
            diffusionAW = Daw * avgDensity[wPhase] * normDiffGrad[wPhase];
            diffusionWW = -diffusionAW;
        } else if (state_i == gasPhase || state_j == gasPhase) {
            diffusionWN = Dwg * avgDensity[nPhase] * normDiffGrad[nPhase];
            diffusionAN = -diffusionWN;
        }

        // Water conservation
        flux[waterF] += (diffusionWW + diffusionWN);

        // CO2 conservation
        flux[co2F] += (diffusionAN + diffusionAW);

        // Heat conservation
        flux[heatF] += diffusionWW*enthW + diffusionAW*enthCO2;

        //////////////////////////////////////////////////////////////////////////////////////////////
        // HEAT CONDUCTION

        flux[heatF] += outward[heatF];

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
    }
    ;

    ///**************************************************************************************************************************///
    /** @brief integrate sources / sinks
     *  @param element entity
     *  @param sol solution vector
     *  @param idx local node id
     *  @return source/sink term
     */
    /////////////////////////////////////////////////////////////////////////////////////////
    //    INTERPOROSITY FLUX

    virtual VBlockType computeQ(const Entity& element, const VBlockType* sol,
                                const int& idx) {
        // ASSUME problem.q already contains \rho.q
        VBlockType q(0.0);
        q = problem.q(this->fvGeom.subContVol[idx].global, element,
                      this->fvGeom.subContVol[idx].local);
        const FieldVector<DT,dim> globalIdx_i=this->fvGeom.subContVol[idx].global;
        const FieldVector<DT,dim> local_i = this->fvGeom.subContVol[idx].local;

        FieldVector <Scalar, numCont> A(0.0);
        A = AreaP(element, idx);
        FieldMatrix<Scalar,dim, dim> KF(0.);
        FieldMatrix<Scalar,dim, dim> KM(0.);
        KF = this->problem.soil().KFracture(globalIdx_i, element, local_i);
        KM = this->problem.soil().KMatrix(globalIdx_i, element, local_i);
        const FMatrix Kharmonic = harmonicMeanKMinc(KF, KM);
        int i = idx;
        FieldVector<Scalar, numCont> KA_nestVol(0.0);
        ////  Product of permeability with the nestedvolume area
        for (int nC=0; nC<numCont; nC++) {
            KA_nestVol[nC] = Kharmonic[0][0] * A[nC];
        }

        /***********************************************************************************/
        //// MINC Mass + Heat Flux
        /***********************************************************************************/
        /***********************************************************************************/
        FieldMatrix<Scalar, numPhases, numCont> InterPorosityFlux(0.0);
        FieldMatrix<Scalar, numPhases, numCont> PGrad(0.0);
        FieldVector<Scalar, numCont> inward(0.0);
        FieldVector<Scalar, numCont> tauWMinc(0.0);
        FieldVector<Scalar, numCont> tauNMinc(0.0); // tortuosity of wetting and nonwetting phase extrapolated to the MINC concept
        FieldVector<Scalar, numCont> DwgMinc(0.0);
        FieldVector<Scalar, numCont> DawMinc(0.0);

        for (int nC = 0; nC < numCont-1; nC++) {
            PGrad[wPhase][nC] += vNDat[i].p[wPhase][nC];
            PGrad[wPhase][nC] -= vNDat[i].p[wPhase][nC+1];
            PGrad[nPhase][nC] += vNDat[i].p[nPhase][nC];
            PGrad[nPhase][nC] -= vNDat[i].p[nPhase][nC+1];


            ///************************************************************************// Pasted from compute A and should be modified
            //// calculate tortuosity at the nodes i and j needed for porous media diffusion coefficient
            //
            //            tauWMinc[nC] = pow(sNDat[i].porosity[nC]
            //                    * vNDat[i].saturation[wPhase][nC], (7/3))
            //                    / (sNDat[i].porosity[nC]
            //                            *sNDat[i].porosity[nC]);
            //            tauWMinc[nC+1] = pow(sNDat[i].porosity[nC+1]
            //                    * vNDat[i].saturation[wPhase][nC+1], (7/3))
            //                    / (sNDat[i].porosity[nC+1]
            //                    *sNDat[i].porosity[nC+1]);
            //
            //            tauNMinc[nC] = pow(sNDat[i].porosity[nC]
            //                    * vNDat[i].saturation[nPhase][nC], (7/3))
            //                    / (sNDat[i].porosity[nC]
            //                    *sNDat[i].porosity[nC]);
            //            tauNMinc[nC+1] = pow(sNDat[i].porosity[nC+1]
            //                    * vNDat[i].saturation[nPhase][nC+1], (7/3))
            //                    / (sNDat[i].porosity[nC+1]
            //                    *sNDat[i].porosity[nC+1]);
            //
            ////            // arithmetic mean of porous media diffusion coefficient
            //            DwgMinc[nC] = (sNDat[i].porosity[nC] * vNDat[i].saturation[nPhase][nC]
            //                    * tauNMinc[nC] * vNDat[i].diff[nPhase]
            //                    + sNDat[i].porosity[nC+1] * vNDat[i].saturation[nPhase][nC+1]
            //                            * tauNMinc[nC+1] * vNDat[i].diff[nPhase])/2;
            //            DawMinc[nC] = (sNDat[i].porosity[nC] * vNDat[i].saturation[wPhase][nC]
            //                    * tauWMinc[nC] * vNDat[i].diff[wPhase]
            //                    + sNDat[i].porosity[nC+1] * vNDat[i].saturation[wPhase][nC+1]
            //                            * tauWMinc[nC+1] * vNDat[i].diff[wPhase])/2;
            //            // arithmetic mean of phase enthalpies (dissolved components neglected)
            FieldVector <Scalar, numCont> enthWMinc(0);
            FieldVector <Scalar, numCont> enthCO2Minc(0);
            enthWMinc[nC] = (vNDat[i].enthalpy[wPhase][nC] + vNDat[i].enthalpy[wPhase][nC+1]) / 2.;
            enthCO2Minc[nC] = (vNDat[i].enthalpy[nPhase][nC] + vNDat[i].enthalpy[nPhase][nC+1]) / 2.;

            //            // Calculate arithmetic mean of the densities
            FieldMatrix <Scalar, numPhases, numCont-1> avgDensityMinc(0);
            avgDensityMinc[wPhase][nC] = 0.5*(vNDat[i].density[wPhase][nC]
                                              + vNDat[i].density[wPhase][nC+1]);
            avgDensityMinc[nPhase][nC] = 0.5*(vNDat[i].density[nPhase][nC]
                                              + vNDat[i].density[nPhase][nC+1]);
            //Calculate arithmetic mean of enthalpies
            FieldMatrix<Scalar, numPhases, numCont-1> avgEnthMinc(0.0);
            avgEnthMinc[wPhase][nC] = 0.5 * (vNDat[i].enthalpy[wPhase][nC]+vNDat[i].enthalpy[wPhase][nC+1]);
            avgEnthMinc[nPhase][nC] = 0.5 * (vNDat[i].enthalpy[nPhase][nC]+vNDat[i].enthalpy[nPhase][nC+1]);//////////////????? what's the difference enthWMinc

            inward[wPhase] = PGrad[wPhase][nC] * KA_nestVol[nC] / elData.delta;

            int up_w, dn_w, up_n, dn_n;
            if (inward[wPhase] > 0) {up_w = nC; dn_w = nC+1;}
            else {up_w = nC+1; dn_w = nC;}
            inward[nPhase] = PGrad[wPhase][nC] * KA_nestVol[nC] / elData.delta;
            if (inward[nPhase] > 0) {up_n = nC; dn_n = nC+1;}
            else {up_n = nC+1; dn_w = nC;}

            InterPorosityFlux[wPhase][nC]=vNDat[i].density[wPhase][up_w]*vNDat[i].mobility[wPhase][up_w]
                * vNDat[i].massfracMincW[water][up_w]
                * inward[wPhase];

            q[3*nC]-= InterPorosityFlux[wPhase][nC];        //wetting flux nC    continua
            q[3*(nC+1)]+= InterPorosityFlux[wPhase][nC];     //wetting flux nC+1 continua
            q[3*nC + 1]-= InterPorosityFlux[nPhase][nC];    //non-wetting flux nC continua
            q[3*(nC+1) + 1]+= InterPorosityFlux[nPhase][nC];//non-wetting flux nC+1 continua
            // Heat InterPorosity Flux for the wetting flux
            q[3*nC+2]    -=  vNDat[i].density[wPhase][up_w] * vNDat[i].mobility[wPhase][up_w]
                * inward[wPhase]
                * avgEnthMinc[wPhase][nC];
            q[3*(nC+1)+2]+=  vNDat[i].density[wPhase][up_w] * vNDat[i].mobility[wPhase][up_w]
                * inward[wPhase]
                * avgEnthMinc[wPhase][nC];
            q[3*nC+2]    -=  vNDat[i].density[nPhase][up_n] * vNDat[i].mobility[nPhase][up_n]
                * inward[nPhase]
                * avgEnthMinc[nPhase][nC];
            q[3*(nC+1)+2]+=  vNDat[i].density[nPhase][up_n] * vNDat[i].mobility[nPhase][up_n]
                * inward[nPhase]
                * avgEnthMinc[nPhase][nC];


        }
        /***********************************************************************************/
        //// MINC Temperature Flux
        /***********************************************************************************/
        /***********************************************************************************/
        FieldMatrix<Scalar, numPhases, numCont> InterPorosityFluxTemperature(0.0);
        FieldVector<Scalar, numCont-1> TempGrad(0.0);
        Scalar inwardTemp(0.0);

        for (int nC = 0; nC < numCont-1; nC++) {
            TempGrad[nC] += vNDat[i].temperature[nC];
            TempGrad[nC] -= vNDat[i].temperature[nC+1];

            // harmonic mean of the heat conductivity lambda
            Scalar lambda[nC];
            lambda[nC] = 2./((1./vNDat[i].lambda[nC]) + (1./vNDat[i].lambda[nC+1]));
            ////       //if the normal of the face gives negative numbers it multiplies with -1
            inwardTemp = TempGrad[nC] * A[nC] / elData.delta;
            if (inwardTemp < 0) {
                inward *=(-1);
            }
            inwardTemp *= lambda[nC];

            //// If the temperature in the fracture is bigger than the one in the matrix the flow occurs from fracture to matrix
            if (vNDat[i].temperature[nC]>vNDat[i].temperature[nC+1]) {
                // temperature interporosity flux calculated;
                InterPorosityFluxTemperature[wPhase][nC]= inwardTemp;
                q[3*nC+2]-=InterPorosityFluxTemperature[wPhase][nC]; // taken from the model of pressure
                q[3*(nC+1)+2]+=InterPorosityFluxTemperature[wPhase][nC];
            }
            else {
                // wetting interporosity flux calculated with the wetting mobility of the Matrix;
                InterPorosityFluxTemperature[wPhase][nC]= inwardTemp;
                q[3*nC+2]+=InterPorosityFluxTemperature[wPhase][nC];
                q[3*(nC+1)+2]-=InterPorosityFluxTemperature[wPhase][nC];
            }

            //                if (vNDat[i].p[nPhase][nC]>vNDat[i].p[nPhase][nC+1]) {
            //
            //                InterPorosityFluxTemperature[nPhase][nC]=inwardTemp;
            //                InterPorosityFluxTemperature[nPhase][nC]=vNDat[i].density[nPhase][nC]
            //                        *vNDat[i].mobility[nPhase][nC]*inward[nPhase];
            //                Scalar test_pressure_n = vNDat[i].p[nPhase][nC];
            //                Scalar test_pressure_n1 = vNDat[i].p[nPhase][nC+1];
            //                Scalar test_mobility = vNDat[i].mobility[nPhase][nC];
            //                Scalar test_Flux = InterPorosityFluxTemperature[nPhase][nC];
            //                q[3*nC+1]-=InterPorosityFluxTemperature[nPhase][nC];//satNFIdx
            //                q[3*(nC+1)+1]+=InterPorosityFluxTemperature[nPhase][nC];//satNMIdx
            //            }
            //            else {
            //                InterPorosityFluxTemperature[nPhase][nC]=inwardTemp;
            ////                IPFluxTemperature[nPhase][nC]=vNDat[i].density[nPhase][nC+1]
            ////                        *vNDat[i].mobility[nPhase][nC+1]*inward[nPhase];
            //                q[3*nC+1]+=InterPorosityFluxTemperature[nPhase][nC];//satNFIdx
            //                q[3*(nC+1)+1]-=InterPorosityFluxTemperature[nPhase][nC];//satNMIdx
            //            }

        }
        return q;
    }




    /** @brief perform variable switch
     *  @param global global node id
     *  @param sol solution vector
     *  @param local local node id
     */
    virtual void primaryVarSwitch(const Entity& element, int global, VBlockType* sol,
                                  int idx) {

        //               if(!switchBreak)
        //                {
        switched = false;
        FieldVector <Scalar, numCont> state;
        FVector Coordinates = this->fvGeom.subContVol[idx].global;

        for (int nC=0; nC<numCont; nC++){
            int state = sNDat[global].phaseState[nC];
            FieldVector<Scalar, numCont> pW(0.), pC(0.), pN(0.), satW(0.), satN(0.);
            pW[nC] = sol[idx][3*nC];
            if (state == bothPhases){satW[nC]=1.0 - sol[idx][3*nC+1];}
            if (state == waterPhase){satW[nC]=1.0;}
            if (state == gasPhase){satW[nC]=0.0;}

            //TODO introduce the pN
            pC[nC] =  problem.materialLaw().pC(satW[nC], this->fvGeom.subContVol[global].global, element,
                                               this->fvGeom.subContVol[global].local);
            pN[nC] = pW[nC]+pC[nC];
            switch (state) {
            case gasPhase:
                Scalar xWNmass;
                xWNmass = sol[idx][3*nC+1];

                if (xWNmass > 0.001 && switched == false) {
                    // appearance of water phase
                    std::cout << "Water appears at node " << global
                              << "  Coordinates: " << Coordinates << std::endl;
                    sNDat[global].phaseState[nC] = bothPhases;
                    sol[idx][3*nC+1] = 1.0 - 1.e-6; // initialize solution vector
                    switched = true;
                }
                break;

            case waterPhase:
                Scalar xAWmax, xAWmass;
                xAWmass = sol[idx][3*nC+1];
                xAWmax = problem.multicomp().xAW(pN[nC], sol[idx][3*nC+2]);

                if (xAWmass > xAWmax && switched == false) {
                    // appearance of gas phase
                    std::cout << "Gas appears at node " << global
                              << "  Coordinates: " << Coordinates << std::endl;
                    sNDat[global].phaseState[nC] = bothPhases;
                    sol[idx][3*nC+1] = 1.e-6; // initialize solution vector
                    switched = true;
                }
                break;

            case bothPhases:
                Scalar satNFracture = sol[idx][3*nC+1];
                if (satNFracture < 0.0 && switched == false) {
                    // disappearance of gas phase
                    std::cout << "Gas disappears at node " << global
                              << "  Coordinates: " << Coordinates << std::endl;
                    sNDat[global].phaseState[nC] = waterPhase;
                    sol[idx][3*nC+1] = 1e-6; // initialize solution vector
                    switched = true;
                }
                else if (satW[nC] < 0.0 && switched == false) {
                    // disappearance of water phase
                    std::cout << "Water disappears at node " << global
                              << "  Coordinates: " << Coordinates << std::endl;
                    sNDat[global].phaseState[nC] = gasPhase;
                    sol[idx][3*nC+1] = 1e-6; // initialize solution vector
                    switched = true;
                }
                break;
            }
        }
        //        ////************************************************************/////
        //          } //End switchBreak
        return;
    }

    virtual void clearVisited() {
        for (int i = 0; i < this->vertexMapper.size(); i++) {
            sNDat[i].visited = false;
        }
        return;
    }

    // updates old phase state after each time step
    virtual void updatePhaseState() {
        for (int i = 0; i < this->vertexMapper.size(); i++) {
            ////MINC
            for (int nC=0; nC<numCont; nC++){
                sNDat[i].oldPhaseState[nC] = sNDat[i].phaseState[nC];
            }
        }

        return;
    }

    virtual void resetPhaseState() {
        for (int i = 0; i < this->vertexMapper.size(); i++) {
            ////MINC
            for (int nC=0; nC<numCont; nC++){
                sNDat[i].phaseState[nC] = sNDat[i].oldPhaseState[nC];
            }
        }
        return;
    }
    //*********************************************************
    //*                                                        *
    //*    Calculation of Data at Elements (elData)             *
    //*                                                         *
    //*                                                        *
    //*********************************************************

    virtual void computeElementData(const Entity& element) {
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
    virtual void updateStaticData(const Entity& element, VBlockType* sol) {
        // size of the sNDat vector is determined in the constructor

        // local to global id mapping (do not ask vertex mapper repeatedly
        //int localToGlobal[LagrangeShapeFunctionSetContainer<DT,Scalar,n>::maxsize];

        // get access to shape functions for P1 elements
        GeometryType gt = element.geometry().type();
        const typename LagrangeShapeFunctionSetContainer<DT,Scalar,dim>::value_type
            & sfs=LagrangeShapeFunctions<DT, Scalar, dim>::general(gt, 1);



        // get local to global id map
        for (int k = 0; k < sfs.size(); k++) {
            const int globalIdx = this->vertexMapper.template map<dim>(element,
                                                                       sfs[k].entity());

            for (int nC=0; nC < numCont; nC++){
                if (nC==0){
                    sNDat[globalIdx].porosity[nC] = problem.soil().porosityFracture(this->fvGeom.subContVol[k].global, element,
                                                                                    this->fvGeom.subContVol[k].local);
                }
                else {
                    sNDat[globalIdx].porosity[nC] = problem.soil().porosityMatrix(this->fvGeom.subContVol[k].global, element,
                                                                                  this->fvGeom.subContVol[k].local);
                }

                // if nodes are not already visited
                if (!sNDat[globalIdx].visited) {
                    // global coordinates
                    FieldVector<DT,dim> globalIdx_i = this->fvGeom.subContVol[k].global;

                    // evaluate primary variable switch
                    primaryVarSwitch(element, globalIdx, sol, k);

                    // mark elements that were already visited
                    sNDat[globalIdx].visited = true;
                }
            }
        }

        return;
    }

    // for initialization of the Static Data (sets porosity)
    virtual void initiateStaticData(const Entity& element) {
        // get access to shape functions for P1 elements
        GeometryType gt = element.geometry().type();
        const typename LagrangeShapeFunctionSetContainer<DT,Scalar,dim>::value_type
            & sfs=LagrangeShapeFunctions<DT, Scalar, dim>::general(gt, 1);

        // get local to global id map
        for (int k = 0; k < sfs.size(); k++) {
            const int globalIdx = this->vertexMapper.template map<dim>(element,
                                                                       sfs[k].entity());
            // if nodes are not already visited
            if (!sNDat[globalIdx].visited) {
                // ASSUME porosity defined at nodes
                sNDat[globalIdx].porosityFracture = problem.soil().porosityFracture(this->fvGeom.subContVol[k].global, element,
                                                                                    this->fvGeom.subContVol[k].local);
                sNDat[globalIdx].porosityMatrix = problem.soil().porosityMatrix(this->fvGeom.subContVol[k].global, element,
                                                                                this->fvGeom.subContVol[k].local);
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


    struct VariableNodeData {
        FieldVector<Scalar,2> mobilityF; //Vector with the number of phases
        FieldVector<Scalar,2> densityF;
        FieldMatrix<Scalar,c,2> massfrac;
        FieldMatrix<FieldVector<Scalar,numPhases>, c, numCont > massMinc;
        FieldMatrix<FieldVector<Scalar,numPhases>, c, numCont > massfracMinc;

        FieldVector<Scalar,2> enthalpyF;
        FieldVector<Scalar,2> intenergyF;
        FieldVector<Scalar,2> diff;

        //MINC
        FieldVector<Scalar, numCont> viscosityCO2;
        FieldVector<Scalar, numCont> krCO2;

        FieldVector<Scalar, numCont> lambda;
        FieldVector<Scalar, numCont> pC; //capillary pressure
        FieldVector<Scalar, numCont> temperature; //temperature (minc)


        FieldMatrix<Scalar, numPhases, numCont> p; //pressure
        FieldMatrix<Scalar, numPhases, numCont> saturation;
        FieldMatrix<Scalar, numPhases, numCont> density; //density in the (M)atrix form
        FieldMatrix<Scalar, numPhases, numCont> mobility; //density in the (M)atrix form
        FieldMatrix<Scalar, numPhases, numCont> enthalpy;
        FieldMatrix<Scalar, numPhases, numCont> intenergy;
        FieldMatrix<Scalar, c, numCont > massfracMincW;
        FieldMatrix<Scalar, c, numCont > massfracMincN;
        FieldVector<Scalar, numCont> phasestate;

        int phasestateF;
    };

    // analog to EvalPrimaryData in MUFTE, uses members of vNDat
    virtual void updateVariableData(const Entity& element, const VBlockType* sol,
                                    int i, std::vector<VariableNodeData>& varData, FieldVector<Scalar, numCont> state) {
        ///////////////////////////////////////////////////////////////////////////////

        const int global = this->vertexMapper.template map<dim>(element, i);
        // Initializing all variables to 0
        varData[i].temperature = 0;
        varData[i].p = 0;
        varData[i].massfracMincW = 0;
        varData[i].massfracMincN = 0;
        for (int nEq=0; nEq<numEq; nEq++) {
            int a = nEq%6;
            int pos = nEq/3; //the number of continuum the variable is located
            // Pressure
            if (a==0) {
                varData[i].p[wPhase][pos] = sol[i][nEq];
            }
            if (a==3) {
                varData[i].p[wPhase][pos] = sol[i][nEq];
            }
            //    Saturation / Mass Fraction
            if (a==1) {
                if (state[pos] == bothPhases) {
                    varData[i].saturation[nPhase][pos] = sol[i][nEq];
                }
                if (state[pos] == waterPhase) {
                    varData[i].saturation[nPhase][pos] = 0.0;
                    varData[i].massfracMincN[water][pos]= 0.0;
                    varData[i].massfracMincW[co2][pos]= sol[i][nEq];
                }
                if (state[pos] == gasPhase)   {
                    varData[i].saturation[nPhase][pos] = 1.0;
                    varData[i].massfracMincN[water]=sol[i][nEq];
                    varData[i].massfracMincW[co2][pos]=0.0;
                }
                varData[i].saturation[wPhase][pos] = 1.0 - varData[i].saturation[nPhase][pos]; //S nPhase

            }
            if (a==4) {
                if (state[pos] == bothPhases) {
                    varData[i].saturation[nPhase][pos] = sol[i][nEq];
                }
                if (state[pos] == waterPhase) {
                    varData[i].saturation[nPhase][pos] = 0.0;
                    varData[i].massfracMincN[water][pos]= 0.0;
                    varData[i].massfracMincW[co2][pos]= sol[i][nEq];
                }
                if (state[pos] == gasPhase) {
                    varData[i].saturation[nPhase][pos] = 1.0;
                    varData[i].massfracMincN[water]=sol[i][nEq];
                    varData[i].massfracMincW[co2][pos]=0.0;
                }
                varData[i].saturation[wPhase][pos] = 1.0 - varData[i].saturation[nPhase][pos]; //S nPhase
            }

            // Temperature
            if (a==2){ varData[i].temperature[pos]=sol[i][nEq];
            }
            if (a==5){ varData[i].temperature[pos]=sol[i][nEq];
            }
        }

        ////// MINC implementation of the Variable Data
        ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
        for (int nC = 0; nC < numCont; nC++) {
            //Capilary pressure and Non-wetting pressure (MINC)
            varData[i].pC[nC] = problem.materialLaw().pC(varData[i].saturation[wPhase][nC],
                                                         this->fvGeom.subContVol[i].global, element, this->fvGeom.subContVol[i].local); //Capillary pressure
            varData[i].p[nPhase][nC] = varData[i].p[wPhase][nC] + varData[i].pC[nC];
            // Mobilities and densities (MINC)
            varData[i].density[wPhase][nC] = problem.wettingPhase().density(varData[i].temperature[nC], varData[i].p[wPhase][nC],
                                                                            varData[i].massfrac[co2F][wPhase]);
            varData[i].density[nPhase][nC] = problem.nonwettingPhase().density(varData[i].temperature[nC], varData[i].p[nPhase][nC]);
            //            std::cout << "temperature: " << varData[i].temperature[nC]<<"pressure :"<< varData[i].p[nPhase][nC]<<"node "<<i<<"\n";
            varData[i].mobility[wPhase][nC] = problem.materialLaw().mobW(varData[i].saturation[wPhase][nC],
                                                                         this->fvGeom.subContVol[i].global, element,
                                                                         this->fvGeom.subContVol[i].local, varData[i].temperature[nC],
                                                                         varData[i].p[wPhase][nC]);

            varData[i].lambda[nC] = problem.soil().heatCond(this->fvGeom.subContVol[i].global, element,
                                                            this->fvGeom.subContVol[i].local, varData[i].saturation[wPhase][nC]);

            varData[i].enthalpy[wPhase][nC]  = problem.wettingPhase().enthalpy(varData[i].temperature[nC], varData[i].p[wPhase][nC]);
            varData[i].enthalpy[nPhase][nC]  = problem.nonwettingPhase().enthalpy(varData[i].temperature[nC], varData[i].p[nPhase][nC]);
            varData[i].intenergy[wPhase][nC] = problem.wettingPhase().intEnergy(varData[i].temperature[nC], varData[i].p[wPhase][nC]);
            varData[i].intenergy[nPhase][nC] = problem.nonwettingPhase().intEnergy(varData[i].temperature[nC], varData[i].p[nPhase][nC]);
            //            varData[idx].mobility[wPhase] = problem.materialLaw().mobW(varData[idx].satW,
            //                    this->fvGeom.subContVol[idx].global, element,
            //                    this->fvGeom.subContVol[idx].local, varData[idx].temperature,
            //                    varData[idx].pW);
            varData[i].mobility[nPhase][nC] = problem.materialLaw().mobN(varData[i].saturation[nPhase][nC],
                                                                         this->fvGeom.subContVol[i].global, element,
                                                                         this->fvGeom.subContVol[i].local, varData[i].temperature[nC],
                                                                         varData[i].p[nPhase][nC]);

            // calculate mobility of CO2 (MINC)

            //            varData[i].viscosityCO2[nC] = problem.nonwettingPhase().viscosityCO2(varData[i].temperature[nC], varData[i].p[nPhase][nC],
            //                    varData[i].density[nPhase][nC]);
            //            if (varData[i].viscosityCO2[nC]) { int x = 5; }
            //
            //            varData[i].krCO2[nC] = problem.materialLaw().krn(varData[i].saturation[nPhase][nC],
            //                    this->fvGeom.subContVol[i].global, element,
            //                    this->fvGeom.subContVol[i].local);
            //            if (varData[i].krCO2[nC]) { int x = 5; }
            //            varData[i].mobility[nPhase][nC] = varData[i].krCO2[nC]/varData[i].viscosityCO2[nC];
            //            if (varData[i].mobility[nPhase][nC]) { int x = 5; }

            ////////******************************************************************////////////////
            // Solubilities of components in phases (MINC)
            if (state[nC] == bothPhases) {
                varData[i].massfracMincW[co2][nC] = problem.multicomp().xAW(varData[i].p[nPhase][nC], varData[i].temperature[nC]);
                varData[i].massfracMincN[water][nC] = problem.multicomp().xWN(varData[i].p[nPhase][nC], varData[i].temperature[nC]);
            }
            varData[i].massfracMincW[water][nC] = 1.0 - varData[i].massfracMincW[co2][nC];
            varData[i].massfracMincN[co2][nC] = 1.0 - varData[i].massfracMincN[water][nC];
        }
        ////////******************************************************************////////////////


        ///////////////////////////////////////////////////////////////////////////////

        // For the fracture domain  used by ComputeA only//
        ///////////////////////////////////////////////////////////////////////////////
        varData[i].diff[wPhase] = problem.wettingPhase().diffCoeff();
        varData[i].diff[nPhase] = problem.nonwettingPhase().diffCoeff();
        // calculate mobility of CO2
        //

        (*outPressureN)[global] = varData[i].p[nPhase][fracture];
        (*outCapillaryP)[global]  = varData[i].pC[fracture];
        (*outSaturationW)[global] = varData[i].saturation[wPhase][fracture];
        (*outSaturationN)[global] = varData[i].saturation[nPhase][fracture];
        (*outTemperature)[global] = varData[i].temperature[fracture];
        (*outMassFracAir)[global] = varData[i].massfrac[co2F][wPhase];
        (*outMassFracWater)[global] = varData[i].massfracMincN[water][fracture];
        (*outDensityW)[global] = varData[i].density[wPhase][fracture];
        (*outDensityN)[global] = varData[i].density[nPhase][fracture];
        (*outMobilityW)[global] = varData[i].mobility[wPhase][fracture];
        (*outMobilityN)[global] = varData[i].mobility[nPhase][fracture];


        return;
    }

    virtual void updateVariableData(const Entity& element, const VBlockType* sol,
                                    int i, bool old = false) {
        FieldVector<Scalar, numCont>state(0.);
        const int global = this->vertexMapper.template map<dim>(element, i);

        if (old){
            state=sNDat[global].oldPhaseState;
            updateVariableData(element, sol, i, oldVNDat, state);
        }
        else {
            state = sNDat[global].phaseState;
            updateVariableData(element, sol, i, vNDat, state);
        }
    }

    void updateVariableData(const Entity& element, const VBlockType* sol,
                            bool old = false) {
        int size = this->fvGeom.numVertices;

        for (int i = 0; i < size; i++)
            updateVariableData(element, sol, i, old);
    }

    struct StaticNodeData {
        bool visited;
        int phaseStateF;
        int oldPhaseStateF;
        Scalar cellVolume;
        Scalar porosityFracture;
        Scalar porosityMatrix;

        FieldVector<Scalar, numCont>phaseState;
        FieldVector<Scalar, numCont>oldPhaseState;



        FieldVector<Scalar, numCont> porosity;
        FieldVector<Scalar, 4> parameters;
        FMatrix K;
    };

    struct StaticIPData {
        bool visited;
        FMatrix K;
        FMatrix KFracture, KMatrix;
    };

    struct ElementData {
        Scalar heatCap;
        double delta; //distance between nested volume faces (kept constant - so all faces are equally distanced)
        //        Scalar cellVolume;
        //          Scalar porosity;
        //        Scalar gravity;
        //        FieldVector<Scalar, 4> parameters;
        //        FieldMatrix<Scalar,dim,dim> K;
    } elData;

    // parameters given in constructor
    TwoPTwoCNIProblem<Grid,Scalar,numEq>& problem;
    CBrineCO2 multicomp;
    std::vector<StaticNodeData> sNDat;
    std::vector<StaticIPData> sIPDat;
    std::vector<VariableNodeData> vNDat;
    std::vector<VariableNodeData> oldVNDat;
    bool switched, switchBreak;

    // for output files
    BlockVector<FieldVector<Scalar, 1> > *outPressureN;
    BlockVector<FieldVector<Scalar, 1> > *outPressureNFracture;
    BlockVector<FieldVector<Scalar, 1> > *outPressureNMatrix;
    BlockVector<FieldVector<Scalar, 1> > *outCapillaryP;
    BlockVector<FieldVector<Scalar, 1> > *outCapillaryPFracture;
    BlockVector<FieldVector<Scalar, 1> > *outCapillaryPMatrix;
    BlockVector<FieldVector<Scalar, 1> > *outSaturationN;
    BlockVector<FieldVector<Scalar, 1> > *outSaturationNFracture;
    BlockVector<FieldVector<Scalar, 1> > *outSaturationNMatrix;
    BlockVector<FieldVector<Scalar, 1> > *outSaturationW;
    BlockVector<FieldVector<Scalar, 1> > *outSaturationWFracture;
    BlockVector<FieldVector<Scalar, 1> > *outSaturationWMatrix;
    BlockVector<FieldVector<Scalar, 1> > *outTemperature;
    BlockVector<FieldVector<Scalar, 1> > *outTemperatureFracture;
    BlockVector<FieldVector<Scalar, 1> > *outTemperatureMatrix;
    BlockVector<FieldVector<Scalar, 1> > *outMassFracAir;
    BlockVector<FieldVector<Scalar, 1> > *outMassFracAirFracture;
    BlockVector<FieldVector<Scalar, 1> > *outMassFracAirMatrix;
    BlockVector<FieldVector<Scalar, 1> > *outMassFracWater;
    BlockVector<FieldVector<Scalar, 1> > *outDensityW;
    BlockVector<FieldVector<Scalar, 1> > *outDensityN;
    BlockVector<FieldVector<Scalar, 1> > *outMobilityW;
    BlockVector<FieldVector<Scalar, 1> > *outMobilityN;
    BlockVector<FieldVector<Scalar, 1> > *outPhaseState;

};

}
#endif


