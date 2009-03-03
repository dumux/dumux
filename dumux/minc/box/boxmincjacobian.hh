

#ifndef DUNE_BOXMINCJACOBIAN_HH
#define DUNE_BOXMINCJACOBIAN_HH

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
#include"dumux/minc/mincproblem.hh"

/**
 * @file
 * @brief  compute local jacobian matrix for conforming finite elements for diffusion equation
 * @author Peter Bastian
 * @modified Alex Tatomir
 */

namespace Dune
{
/** @addtogroup DISC_Disc
 *
 * @{
 */
/**
 * @brief compute local jacobian matrix for conforming finite elements for diffusion equation
 *
 */


//! A class for computing local jacobian matrices
/*! A class for computing local jacobian matrix for the
  diffusion equation

  div j = q; j = -K grad u; in Omega

  u = g on Gamma1; j*n = J on Gamma2.

  Uses conforming finite elements with the Lagrange shape functions.
  It should work for all dimensions and element types.
  All the numbering is with respect to the reference element and the
  Lagrange shape functions

  Template parameters are:

  - Grid  a DUNE grid type
  - Scalar type used for return values
*/
template<class Grid, class Scalar, int numEq, class BoxFunction = LeafP1Function<Grid, Scalar, numEq> > class BoxMincJacobian
    : public BoxJacobian<BoxMincJacobian<Grid,Scalar, numEq,BoxFunction>,Grid,Scalar, numEq,BoxFunction>
{
    typedef typename Grid::ctype DT;
    typedef typename Grid::Traits::template Codim<0>::Entity Entity;
    typedef typename Entity::Geometry Geometry;
    typedef BoxMincJacobian<Grid,Scalar,numEq,BoxFunction> ThisType;
    typedef typename LocalJacobian<ThisType,Grid,Scalar,numEq>::VBlockType VBlockType;
    typedef typename LocalJacobian<ThisType,Grid,Scalar,numEq>::MBlockType MBlockType;
    //     typedef Dune::FVElementGeometry<G> FVElementGeometry;
    // number of phases

    enum {wPhase = 0, nPhase = 1};
    enum {F = 0, M = 1};
    enum {pWFIdx = 0, satNFIdx = 1, pWMIdx = 2, satNMIdx = 3};
    enum {WF = 0, NF = 1, WM = 2, NM = 3};
    enum {nPermeabilities =2}; // Number of permeabilities

public:
    // define the number of components of your system, this is used outside
    // to allocate the correct size of (dense) blocks with a FieldMatrix
    // m represents the number of equations
    enum {dim = Grid::dimension};
    enum {nPhases = 2}; //number of phases
    // number of interacting continua
    enum {numCont = numEq/2};
    enum {numPhases = 2};
    enum {SIZE=LagrangeShapeFunctionSetContainer<DT,Scalar,dim>::maxsize};
    struct VariableNodeData;

    typedef FieldMatrix<Scalar,dim,dim> FMatrix;
    typedef FieldVector<Scalar,dim> FVector;
    //! Constructor
    BoxMincJacobian (MincProblem<Grid,Scalar, numEq>& params,
                     bool levelBoundaryAsDirichlet_, const Grid& grid,
                     BoxFunction& sol,
                     bool procBoundaryAsDirichlet_=true)
        : BoxJacobian<ThisType,Grid,Scalar,numEq,BoxFunction>(levelBoundaryAsDirichlet_, grid, sol, procBoundaryAsDirichlet_),
          problem(params),
          statNData(this->vertexMapper.size()),
          vNData(SIZE), oldVNData(SIZE)
    {
        this->analytic = false;
    }
    /***********************************************************************************************************/
    ////  Calculating the perimeter of the control volume that will be used as Aj,j+1
    virtual FieldVector<Scalar, numCont> AreaP (const Entity& e, int face)
    //    virtual double AreaP (const Entity& e, int face)
    {
        FieldVector<Scalar,numCont> x(0.0);
        FieldVector<Scalar,numCont> y(0.0);
        FieldVector<Scalar,numCont> L(0.0);
        double length_x =0;      //the length of the subcontrol face on x direction
        double length_y =0;      //the length of the subcontrol face on y direction
        double ExtLength =0;        //Exterior length of the subcontrol face
        for (int face =0; face < this->fvGeom.numVertices; face++){
            int it =this->fvGeom.subContVolFace[face].i;
            int jt =this->fvGeom.subContVolFace[face].j;
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
        //              double f = 1.0 / numCont; // volume fraction
        //              // Calculation of the distance between the nested volume elements delta
        //              // delta in this case divides the surfaces equidistantly (delta is const)
        elData.delta = 0.0;
        if (length_x > length_y){
            elData.delta = 2*length_y / (2.0 * (numCont-1) + 1.0);
        }
        else {
            elData.delta = 2*length_x / (2.0 * (numCont-1) + 1.0);
        }

        double TotalLength = 4*ExtLength;

        L[0]  = TotalLength;
        x[0]  = 2*length_x;    //length of the unitary element (CV) on x direction
        y[0]  = 2*length_y;    //length of the unitary element (CV) on y direction
        for (int nC=1; nC<numCont; nC++){
            x[nC] = x[nC-1] - 2 * elData.delta;
            y[nC] = y[nC-1] - 2 * elData.delta;
            L[nC] = 2* (x[nC]+y[nC]);
        }
        return L;

        //              return TotalLength;
    }

    ///***********************************************************************************************************/
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



    /*the harmonicMeanK function computes the harmonic mean of the perameabilities between the two nodes of different
     *  permeabilities in the fracture domain */

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

    /*************************************************************************************************************/
    virtual void clearVisited ()
    {
        return;
    }

    virtual VBlockType computeM (const Entity& e, const VBlockType* sol,
                                 int node, std::vector<VariableNodeData>& vNData)
    {
        VBlockType AccumulationTerm;
        MBlockType result;

        for (int nC = 0; nC < numCont; nC++){

            result[wPhase][nC] = - vNData[node].density[wPhase][nC]* elData.porosity[nC] * vNData[node].saturation[nPhase][nC];
            result[nPhase][nC] = vNData[node].density[nPhase][nC]* elData.porosity[nC] * vNData[node].saturation[nPhase][nC];
        }

        for (int nC = 0; nC < numCont ;nC++ )
        {
            //2 is the number of phases (wetting and non-wetting) and nC - the number of Continua
            int count = 2*nC;
            AccumulationTerm[count]=result[wPhase][nC];
            AccumulationTerm[count+1]=result[nPhase][nC];
        }

        return AccumulationTerm;
    };

    virtual VBlockType computeM (const Entity& e, const VBlockType* sol, int node, bool old = false)
    {
        if (old)
            return computeM(e, sol, node, oldVNData);
        else
            return computeM(e, sol, node, vNData);
    }

    //***********************************************************************************************//
    virtual VBlockType computeA (const Entity& e, const VBlockType* sol, int face)
    {
        int i = this->fvGeom.subContVolFace[face].i;
        int j = this->fvGeom.subContVolFace[face].j;

        // permeability in edge direction
        FieldVector<DT,dim> KFractureij(0.);
        elData.KFracture.umv(this->fvGeom.subContVolFace[face].normal, KFractureij);
        FieldVector<DT,dim> h;
        h = this->fvGeom.subContVolFace[face].normal;
        // permeability of the Interacting Continua
        FieldVector<DT,dim> KMatrixij(0.);
        elData.KMatrix.umv(this->fvGeom.subContVolFace[face].normal, KMatrixij);
        const FieldVector<DT,dim> global_i = this-> fvGeom.subContVol[i].global;
        const FieldVector<DT,dim> global_j = this-> fvGeom.subContVol[j].global;
        const FieldVector<DT,dim> local_i = this->fvGeom.subContVol[i].local;
        const FieldVector<DT,dim> local_j = this->fvGeom.subContVol[j].local;
        FMatrix Ki(0.0), Kj(0.0);
        Ki = this->problem.soil().KFracture(global_i, e, local_i);
        Kj = this->problem.soil().KFracture(global_j, e, local_j);

        const FMatrix K = harmonicMeanK(Ki, Kj);
        const FieldVector<DT,dim> normal(this->fvGeom.subContVolFace[face].normal);
        VBlockType flux(0.0);
        FieldVector<Scalar, dim> pGrad(0.0);
        FieldVector<Scalar, dim> pWFGrad(0.0);
        FieldVector<Scalar, dim> pNFGrad(0.0);
        FieldVector<Scalar, dim> pWMGrad(0.0);
        FieldVector<Scalar, dim> pNMGrad(0.0);
        for (int comp = 0; comp < 4; comp++)
        {
            FieldVector<Scalar, dim> gravity = problem.gravity();
            switch (comp)    {
            case WF:
                for (int k = 0; k < this->fvGeom.numVertices; k++) {
                    FieldVector<DT,dim> grad(this->fvGeom.subContVolFace[face].grad[k]);
                    //                  grad *= sol[k][pWFIdx];
                    grad *= vNData[k].p[wPhase][F];
                    pWFGrad += grad;
                }
                // adjust by gravity
                gravity *= vNData[i].density[wPhase][0];
                pWFGrad -= gravity;
                break;
            case NF:
                for (int k = 0; k < this->fvGeom.numVertices; k++) {
                    FieldVector<DT,dim> grad(this->fvGeom.subContVolFace[face].grad[k]);
                    grad *= vNData[k].p[nPhase][F];
                    pNFGrad += grad;
                }
                // adjust by gravity
                gravity *= vNData[i].density[nPhase][0];
                pNFGrad -= gravity;
                break;
            case WM:
                break;
            case NM:
                break;

            }

            FieldVector<Scalar,dim>Kij(0.0);
            K.umv(normal, Kij);

            Scalar out_WF = pWFGrad*Kij;
            Scalar out_NF = pNFGrad*Kij;

            int up_WF, down_WF, up_NF, down_NF;
            if (out_WF<0){
                up_WF = i; down_WF=j;}
            else {
                up_WF = j; down_WF=i;}
            if (out_NF<0){
                up_NF = i; down_NF=j;}
            else {
                up_NF = j; down_NF=i;}
            // assigns the fully upwind mobility
            if (out_WF<0)
                flux[pWFIdx]=vNData[i].density[wPhase][F]*vNData[i].mobility[wPhase][F]*out_WF;
            else
                flux[pWFIdx]=vNData[j].density[wPhase][F]*vNData[j].mobility[wPhase][F]*out_WF;
            if (out_NF<0)
                flux[satNFIdx]=vNData[i].density[nPhase][F]*vNData[i].mobility[nPhase][F]*out_NF;
            else
                flux[satNFIdx]=vNData[j].density[nPhase][F]*vNData[j].mobility[nPhase][F]*out_NF;
        }

        return flux;
    };

    //******************************************************************************************//
    //******************************************************************************************//

    virtual VBlockType computeQ (const Entity& e, const VBlockType* sol, const int& node)
    {
        // ASSUME problem.q already contains \rho.q
        VBlockType q(0);

        ////      std::ofstream file;
        ////      file.open("file.txt", std::ios::app);


        q = problem.q(this->fvGeom.subContVol[node].global, e, this->fvGeom.subContVol[node].local);
        const FieldVector<DT,dim> globalIdx_i=this->fvGeom.subContVol[node].global;
        const FieldVector<DT,dim> localIdx_i = this->fvGeom.subContVol[node].local;
        FieldVector <Scalar, numCont> A(0.0);
        A = AreaP(e, node);
        FieldMatrix<Scalar,dim, dim> KF(0.);
        FieldMatrix<Scalar,dim, dim> KM(0.);
        KF = this->problem.soil().KFracture(globalIdx_i, e, localIdx_i);
        KM = this->problem.soil().KMatrix(globalIdx_i, e, localIdx_i);
        const FMatrix Kharmonic = harmonicMeanKMinc(KF, KM);
        //      const FMatrix Kharmonic = harmonicMeanKMinc (global_i);
        int i = node;
        //      double K_used = Kharmonic[0][0];
        ////  Product of permeability with the nestedvolume area
        FieldVector<Scalar, numCont> KA_nestVol(0.0);
        for (int nC=0; nC<numCont; nC++)
        {
            KA_nestVol[nC] = Kharmonic[0][0] * A[nC];
        }
        /***********************************************************************************/
        //// MINC Flux Matrix
        FieldMatrix <Scalar, dim, numCont> IPFlux(0.);
        FieldMatrix <Scalar, dim, numCont> PGrad(0.);
        FieldVector <Scalar, numCont> inward(0.);
        for (int nC = 0; nC < numCont-1; nC++)
        {
            PGrad[wPhase][nC] += vNData[i].p[wPhase][nC];
            PGrad[wPhase][nC] -= vNData[i].p[wPhase][nC+1];
            PGrad[nPhase][nC] += vNData[i].p[nPhase][nC];
            PGrad[nPhase][nC] -= vNData[i].p[nPhase][nC+1];

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
            if (vNData[i].p[wPhase][nC]>vNData[i].p[wPhase][nC+1])
            {
                // wetting interporosity flux calculated with the wetting mobility of the Fracture;
                IPFlux[wPhase][nC]=vNData[i].density[wPhase][nC]*vNData[i].mobility[wPhase][nC]*inward[wPhase];
                q[2*nC]-=IPFlux[wPhase][nC]; //pWFIdx
                q[2*(nC+1)]+=IPFlux[wPhase][nC];//pWMIdx
            }
            else
            {
                // wetting interporosity flux calculated with the wetting mobility of the Matrix;
                IPFlux[wPhase][nC]=vNData[i].density[wPhase][nC+1]*vNData[i].mobility[wPhase][nC+1]*inward[wPhase];
                q[2*nC]+=IPFlux[wPhase][nC];//pWFIdx
                q[2*(nC+1)]-=IPFlux[wPhase][nC];//pWMIdx
            }
            if (vNData[i].p[nPhase][nC]>vNData[i].p[nPhase][nC+1])
            {
                IPFlux[nPhase][nC]=vNData[i].density[nPhase][nC]*vNData[i].mobility[nPhase][nC]*inward[nPhase];
                q[2*nC+1]-=IPFlux[nPhase][0];//satNFIdx
                q[2*(nC+1)+1]+=IPFlux[nPhase][0];//satNMIdx
            }
            else
            {
                IPFlux[nPhase][nC]=vNData[i].density[nPhase][nC+1]*vNData[i].mobility[nPhase][nC+1]*inward[nPhase];
                q[2*nC+1]+=IPFlux[nPhase][nC];//satNFIdx
                q[2*(nC+1)+1]-=IPFlux[nPhase][nC];//satNMIdx
            }
        }
        /************************************************************************************/
        //      VBlockType q_previous(0);
        //      q_previous = q;
        //      for (int c=0; c<4; c++){
        //          file << q[c]<<"\t\t";
        //              }
        //      file << std::endl;
        //      file.close();

        return q;
    };



    //*********************************************************
    //*                                                                *
    //*    Calculation of Data at Elements                                 *
    //*                                                                 *
    //*                                                                 *
    //*********************************************************

    virtual void computeElementData (const Entity& e)
    {
        // ASSUME element-wise constant parameters for the material law
        //          elData.parametersFracture = problem.materialLawParametersFracture(this->fvGeom.elementGlobal, e, this->fvGeom.elementLocal);
        //          elData.parametersMatrix = problem.materialLawParametersMatrix(this->fvGeom.elementGlobal, e, this->fvGeom.elementLocal);
        // ASSUMING element-wise constant permeability, evaluate K at the cell center
        elData.KFracture = problem.soil().KFracture(this->fvGeom.elementGlobal, e, this->fvGeom.elementLocal);
        elData.KMatrix = problem.soil().KMatrix(this->fvGeom.elementGlobal, e, this->fvGeom.elementLocal);

        for (int nC=0; nC<numCont; nC++){
            if (nC==0){
                elData.porosity[nC] = problem.soil().porosityFracture(this->fvGeom.elementGlobal, e, this->fvGeom.elementLocal);
            }
            else{
                elData.porosity[nC] = problem.soil().porosityMatrix(this->fvGeom.elementGlobal, e, this->fvGeom.elementLocal);
            }
        }
    };


    //*********************************************************
    //*                                                            *
    //*    Calculation of Data at Nodes that has to be                 *
    //*    determined only once    (statNData)                        *
    //*                                                            *
    //*********************************************************

    // analog to EvalStaticData in MUFTE
    virtual void updateStaticData (const Entity& e, const VBlockType* sol)
    {
        return;
    }


    //*********************************************************
    //*                                                        *
    //*    Calculation of variable Data at Nodes                 *
    //*    (vNData)                                             *
    //*                                                         *
    //*********************************************************


    // analog to EvalPrimaryData in MUFTE, uses members of varNData
    virtual void updateVariableData(const Entity& e, const VBlockType* sol,
                                    int i, std::vector<VariableNodeData>& vNData)
    {
        vNData.resize(this->fvGeom.numVertices);
        int size = vNData.size();

        for (int i = 0; i < size; i++) {

            for (int nEq=0; nEq<numEq; nEq++)
            {
                int a = nEq%4;
                int pos = nEq/2;
                if (a==1){
                    vNData[i].saturation[wPhase][pos]  = 1.0 - sol[i][nEq];        //S wPhase F
                    vNData[i].saturation[nPhase][pos]  = sol[i][nEq];            //S nPhase F
                }
                if (a==3){
                    vNData[i].saturation[wPhase][pos]= 1.0 - sol[i][nEq];         //S wPhase M
                    vNData[i].saturation[nPhase][pos]= sol[i][nEq];            //S nPhase M
                }
            }

            vNData[i].pCFracture = problem.materialLaw().pC(vNData[i].saturation[wPhase][F],
                                                            this->fvGeom.subContVol[i].global, e, this->fvGeom.subContVol[i].local);
            vNData[i].pCMatrix   = problem.materialLaw().pC(vNData[i].saturation[wPhase][M],
                                                            this->fvGeom.subContVol[i].global, e, this->fvGeom.subContVol[i].local);
            for (int nC=0; nC<numCont; nC++){
                if (nC==0){
                    vNData[i].pC[nC]=problem.materialLaw().pC(vNData[i].saturation[wPhase][F],
                                                              this->fvGeom.subContVol[i].global, e, this->fvGeom.subContVol[i].local);
                }
                else {
                    vNData[i].pC[nC]=problem.materialLaw().pC(vNData[i].saturation[wPhase][nC],
                                                              this->fvGeom.subContVol[i].global, e, this->fvGeom.subContVol[i].local);
                }
            }
            // it assigns the values of the solution vector which has 4 components to the vNData Matrix (2x2) for 2 continua.
            for (int nEq = 0; nEq < numEq; nEq++)
            {
                int a = nEq%4;
                int pos = nEq/2;
                if (a==0) {
                    vNData[i].p[wPhase][pos]   = sol[i][nEq];                             //P wPhase F
                    vNData[i].p[nPhase][pos]   = sol[i][nEq] + vNData[i].pCFracture;     //P nPhase F
                }
                if (a==2) {
                    vNData[i].p[wPhase][pos] = sol[i][nEq];                             //P wPhase M
                    vNData[i].p[nPhase][pos] = sol[i][nEq] + vNData[i].pCMatrix;         //P nPhase M
                }
            }

            //Mobilities
            for (int nC = 0; nC < numCont ;nC++ )
            {
                vNData[i].mobility[wPhase][nC] = problem.materialLaw().mobW(vNData[i].saturation[wPhase][nC],
                                                                            this->fvGeom.subContVol[i].global, e, this->fvGeom.subContVol[i].local);
                vNData[i].mobility[nPhase][nC] = problem.materialLaw().mobN(vNData[i].saturation[nPhase][nC],
                                                                            this->fvGeom.subContVol[i].global, e, this->fvGeom.subContVol[i].local);
                vNData[i].density[wPhase][nC]  = problem.materialLaw().wettingPhase.density();
                vNData[i].density[nPhase][nC]  = problem.materialLaw().nonwettingPhase.density();
            }
        }

    }
    //*******************************************************************
    virtual void updateVariableData(const Entity& e, const VBlockType* sol, int i, bool old = false)
    {
        if (old)
            updateVariableData(e, sol, i, oldVNData);
        else
            updateVariableData(e, sol, i, vNData);
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
        FieldMatrix<DT,dim,dim> K;
        FieldMatrix<DT,dim,dim> K1Fracture;
        FieldMatrix<DT,dim,dim> K1Matrix;
    };

    // the members of the struct are defined here
    struct VariableNodeData
    {

        Scalar pCFracture;
        Scalar pCMatrix;
        FieldVector<DT, numCont> pC;
        FieldMatrix<DT, numPhases, numCont> p;
        FieldMatrix<DT, numPhases, numCont> saturation;
        FieldMatrix<DT, numPhases, numCont> mobility;
        FieldMatrix<DT, numPhases, numCont> density;
    };

    struct ElementData {
        Scalar cellVolume;
        FieldVector<Scalar, numCont> porosity;
        Scalar gravity;

        FieldMatrix<DT,dim,dim> KFracture;
        FieldMatrix<DT,dim,dim> KMatrix;

        double delta; //distance between nested volume faces (kept constant - so all faces are equally distanced)
    } elData;


    // parameters given in constructor
    MincProblem<Grid,Scalar, numEq>& problem;
    std::vector<StaticNodeData> statNData;
    std::vector<VariableNodeData> vNData;
    std::vector<VariableNodeData> oldVNData;
};

/** @} */
}
#endif
