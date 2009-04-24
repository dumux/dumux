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
  - Scalar    type used for return values
*/
template<class Grid, class Scalar, class BoxFunktion = LeafP1Function<
                                       Grid, Scalar, 3> >
class BoxPwSnTeJacobian: public BoxJacobian<BoxPwSnTeJacobian<Grid, Scalar,
                                                              BoxFunktion> , Grid, Scalar, 3, BoxFunktion>
{
    typedef    typename Grid::ctype CoordScalar;
    typedef typename Grid::Traits::template Codim<0>::Entity Element;
    typedef typename Element::Geometry Geometry;
    typedef BoxPwSnTeJacobian<Grid,Scalar,BoxFunktion> ThisType;
    typedef typename LocalJacobian<ThisType,Grid,Scalar,3>::VBlockType VBlockType;
    typedef typename LocalJacobian<ThisType,Grid,Scalar,3>::VBlockType SolutionVector;
    typedef typename LocalJacobian<ThisType,Grid,Scalar,3>::MBlockType MBlockType;
    enum
        {   wPhase = 0, nPhase = 1, temp = 2};

public:
    // define the number of components of your system, this is used outside
    // to allocate the correct size of (dense) blocks with a FieldMatrix
    enum
        {   dim=Grid::dimension};
    enum
        {   numEq=3};
    enum
        {   SIZE=LagrangeShapeFunctionSetContainer<CoordScalar,Scalar,dim>::maxsize};
    struct VariableNodeData;
    typedef FieldMatrix<Scalar,dim,dim> FMatrix;
    typedef FieldVector<Scalar,dim> FVector;

    //! Constructor
    BoxPwSnTeJacobian (TwoPhaseHeatProblem<Grid,Scalar>& params,
                       bool levelBoundaryAsDirichlet_, const Grid& grid,
                       BoxFunktion& sol,
                       bool procBoundaryAsDirichlet_=true)
        : BoxJacobian<ThisType,Grid,Scalar,numEq,BoxFunktion>(levelBoundaryAsDirichlet_, grid, sol, procBoundaryAsDirichlet_),
          problem(params),
          sNDat(this->vertexMapper.size()), vNDat(SIZE), oldVNDat(SIZE)
    {
        this->analytic = false;
    }

    /** @brief compute time dependent term (storage), loop over nodes / subcontrol volumes
     *  @param element entity
     *  @param sol solution vector
     *  @param idx local node id
     *  @return storage term
     */
    virtual VBlockType computeStorage (const Element& element, const SolutionVector* sol,
                                 int idx, std::vector<VariableNodeData>& varData)
    {
        GeometryType gt = element.geometry().type();
        const typename LagrangeShapeFunctionSetContainer<CoordScalar,Scalar,dim>::value_type&
            sfs=LagrangeShapeFunctions<CoordScalar,Scalar,dim>::general(gt,1);

        int globalIdx = this->vertexMapper.template map<dim>(element, sfs[idx].entity());

        VBlockType result;
        Scalar satN = varData[idx].satN;
        Scalar satW = varData[idx].satW;

        // storage of component wPhase
        result[wPhase] =
            sNDat[globalIdx].porosity*varData[idx].density[wPhase]*satW;
        // storage of component co2
        result[nPhase] =
            sNDat[globalIdx].porosity*varData[idx].density[nPhase]*satN;

        // storage term of energy equation
        result[temp] = sNDat[globalIdx].porosity * (varData[idx].density[wPhase] * varData[idx].intenergy[wPhase] * satW
                                                    + varData[idx].density[nPhase] * varData[idx].intenergy[nPhase] * satN)
            + elData.heatCap * varData[idx].temperature;
        // soil properties defined at the elements!!!

        return result;
    };

    virtual VBlockType computeStorage (const Element& element, const SolutionVector* sol, int idx, bool old = false)
    {
        if (old)
            return computeStorage(element, sol, idx, oldVNDat);
        else
            return computeStorage(element, sol, idx, vNDat);
    }

    /** @brief compute fluxes and heat conduction, loop over subcontrol volume faces
     *  @param element entity
     *  @param sol solution vector
     *  @param face face id
     *  @return flux term
     */
    VBlockType computeFlux (const Element& element, const SolutionVector* sol, int face)
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

        // Calculate arithmetic mean of the densities
        VBlockType avgDensity;
        avgDensity[wPhase] = 0.5*(vNDat[idx_i].density[wPhase] + vNDat[idx_j].density[wPhase]);
        avgDensity[nPhase] = 0.5*(vNDat[idx_i].density[nPhase] + vNDat[idx_j].density[nPhase]);

        //////////////////////////////////////////////////////////////////////////////////////////
        // GRADIENTS

        FieldMatrix<Scalar,2,dim> pGrad(0.);
        FieldVector<Scalar,dim> teGrad(0.);
        FieldVector<Scalar,dim> tmp(0.);
        VBlockType flux(0.);
        FieldVector<Scalar,2> densityIP(0.0);

        // calculate FE gradient at subcontrolvolumeface
        for (int idx = 0; idx < this->fvGeom.numVertices; idx++) // loop over adjacent nodes

        {
            // FEGradient at subcontrolvolumeface face
            const FieldVector<CoordScalar,dim> feGrad(this->fvGeom.subContVolFace[face].grad[idx]);
            FieldVector<Scalar,2> pressure(0.0);

            pressure[wPhase] = vNDat[idx].pW;
            pressure[nPhase] = vNDat[idx].pN;

            // compute pressure gradients for each phase at integration point of subcontrolvolumeface face
            for (int phase = 0; phase < 2; phase++)
            {
                tmp = feGrad;
                tmp *= pressure[phase];
                pGrad[phase] += tmp;
            }

            // compute temperature gradient
            tmp = feGrad;
            tmp *= vNDat[idx].temperature;
            teGrad += tmp;

            densityIP[wPhase] += vNDat[idx].density[wPhase] * this->fvGeom.subContVolFace[face].shapeValue[idx];
            densityIP[nPhase] += vNDat[idx].density[nPhase] * this->fvGeom.subContVolFace[face].shapeValue[idx];
        }

        // deduce gravity*density of each phase
        FieldMatrix<Scalar,2,dim> contribComp(0);
        for (int phase=0; phase<2; phase++)
        {
            contribComp[phase] = problem.gravity();
            contribComp[phase] *= densityIP[phase];
            pGrad[phase] -= contribComp[phase]; // grad p - rho*g
        }

        // Darcy velocity in normal direction for each phase K*n(grad p -rho*g)
        VBlockType outward(0);
        FieldVector<Scalar,dim> v_tilde_w(0);
        FieldVector<Scalar,dim> v_tilde_n(0);

        K.umv(pGrad[wPhase], v_tilde_w); // v_tilde=K*gradP
        outward[wPhase] = v_tilde_w*normal;
        K.umv(pGrad[nPhase], v_tilde_n); // v_tilde=K*gradP
        outward[nPhase] = v_tilde_n*normal;

        // Heat conduction
        outward[temp] = teGrad * normal;
        outward[temp] *= lambda;

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
        flux[wPhase] = (alpha* vNDat[up_w].density[wPhase]*vNDat[up_w].mobility[wPhase]
                        + (1-alpha)* vNDat[dn_w].density[wPhase]*vNDat[dn_w].mobility[wPhase])
            * outward[wPhase];
        // CO2 conservation
        flux[nPhase] = (alpha* vNDat[up_n].density[nPhase]*vNDat[up_n].mobility[nPhase]
                        + (1-alpha)* vNDat[dn_n].density[nPhase]*vNDat[dn_n].mobility[nPhase])
            * outward[nPhase];
        // Heat conservation
        flux[temp] = (alpha* vNDat[up_n].density[nPhase]*vNDat[up_n].mobility[nPhase] * vNDat[up_n].enthalpy[nPhase]
                      + (1-alpha)* vNDat[dn_n].density[nPhase]*vNDat[dn_n].mobility[nPhase] * vNDat[dn_n].enthalpy[nPhase])
            * outward[nPhase];
        flux[temp] += (alpha* vNDat[up_w].density[wPhase]*vNDat[up_w].mobility[wPhase] * vNDat[up_w].enthalpy[wPhase]
                       + (1-alpha)* vNDat[dn_w].density[wPhase]*vNDat[dn_w].mobility[wPhase] * vNDat[dn_w].enthalpy[wPhase])
            * outward[wPhase];

        //////////////////////////////////////////////////////////////////////////////////////////////
        // HEAT CONDUCTION

        flux[temp] += outward[temp];

        return flux;
    };

    /** @brief integrate sources / sinks
     *  @param element entity
     *  @param sol solution vector
     *  @param idx local node id
     *  @return source/sink term
     */
    virtual VBlockType computeSource (const Element& element, const SolutionVector* sol, const int& idx)
    {
        // ASSUME problem.q already contains \rho.q
        return problem.q(this->fvGeom.subContVol[idx].global, element, this->fvGeom.subContVol[idx].local);
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

    // *********************************************************
    // *                                                         *
    // *    Calculation of Data at Elements                          *
    // *                                                          *
    // *                                                         *
    // *********************************************************
    void computeElementData (const Element& element)
    {
        elData.heatCap = problem.soil().heatCap(this->fvGeom.elementGlobal, element, this->fvGeom.elementLocal);
    };

    // *********************************************************
    // *                                                         *
    // *    Calculation of Data at Nodes that has to be             *
    // *    determined only once    (statNData)                     *
    // *                                                         *
    // *********************************************************

    // analog to EvalStaticData in MUFTE
    void updateStaticData (const Element& element, const SolutionVector* sol)
    {
        // size of the sNDat vector is determined in the constructor

        // local to global id mapping (do not ask vertex mapper repeatedly
        //int localToGlobal[LagrangeShapeFunctionSetContainer<CoordScalar,Scalar,n>::maxsize];

        // get access to shape functions for P1 elements
        GeometryType gt = element.geometry().type();
        const typename LagrangeShapeFunctionSetContainer<CoordScalar,Scalar,dim>::value_type&
            sfs=LagrangeShapeFunctions<CoordScalar,Scalar,dim>::general(gt,1);

        // get local to global id map
        for (int idx = 0; idx < sfs.size(); idx++)
        {
            const int globalIdx = this->vertexMapper.template map<dim>(element, sfs[idx].entity());

            // if nodes are not already visited
            if (!sNDat[globalIdx].visited)
            {
                // ASSUME porosity defined at nodes
                sNDat[globalIdx].porosity = problem.soil().porosity(this->fvGeom.subContVol[idx].global, element, this->fvGeom.subContVol[idx].local);

                // mark elements that were already visited
                sNDat[globalIdx].visited = true;
            }
        }

        return;
    }

    virtual void initiateStaticData (const Element& element)
    {
        // get access to shape functions for P1 elements
        GeometryType gt = element.geometry().type();
        const typename LagrangeShapeFunctionSetContainer<CoordScalar,Scalar,dim>::value_type&
            sfs=LagrangeShapeFunctions<CoordScalar,Scalar,dim>::general(gt,1);

        // get local to global id map
        for (int idx = 0; idx < sfs.size(); idx++)
        {
            const int globalIdx = this->vertexMapper.template map<dim>(element, sfs[idx].entity());

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
    //*    Calculation of variable Data at Nodes                 *
    //*    (varNData)                                             *
    //*                                                         *
    //*********************************************************


    // the members of the struct are defined here
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
        FieldVector<Scalar,2> enthalpy;
        FieldVector<Scalar,2> intenergy;
    };

    // analog to EvalPrimaryData in MUFTE, uses members of vNDat
    virtual void updateVariableData(const Element& element, const SolutionVector* sol,
                                    int idx, std::vector<VariableNodeData>& varData)
    {
        const int global = this->vertexMapper.template map<dim>(element, idx);

        varData[idx].pW = sol[idx][wPhase];
        varData[idx].satN = sol[idx][nPhase];
        varData[idx].satW = 1.0 - varData[idx].satN;
        varData[idx].pC = problem.materialLaw().pC(varData[idx].satW, this->fvGeom.subContVol[idx].global, element, this->fvGeom.subContVol[idx].local);
        varData[idx].pN = varData[idx].pW + varData[idx].pC;
        varData[idx].temperature = sol[idx][temp]; // in [K]

        // for output
        varData[idx].density[wPhase] = problem.wettingPhase().density(varData[idx].temperature, varData[idx].pW);
        varData[idx].density[nPhase] = problem.nonwettingPhase().density(varData[idx].temperature, varData[idx].pN);
        varData[idx].mobility[wPhase] = problem.materialLaw().mobW(varData[idx].satW, this->fvGeom.subContVol[idx].global, element, this->fvGeom.subContVol[idx].local, varData[idx].temperature, varData[idx].pW);
        varData[idx].mobility[nPhase] = problem.materialLaw().mobN(varData[idx].satN, this->fvGeom.subContVol[idx].global, element, this->fvGeom.subContVol[idx].local, varData[idx].temperature, varData[idx].pN);
        varData[idx].lambda = problem.soil().heatCond(this->fvGeom.subContVol[idx].global, element, this->fvGeom.subContVol[idx].local, varData[idx].satW);
        varData[idx].enthalpy[wPhase] = problem.wettingPhase().enthalpy(varData[idx].temperature,varData[idx].pW);
        varData[idx].enthalpy[nPhase] = problem.nonwettingPhase().enthalpy(varData[idx].temperature,varData[idx].pN);
        varData[idx].intenergy[wPhase] = problem.wettingPhase().intEnergy(varData[idx].temperature,varData[idx].pW);
        varData[idx].intenergy[nPhase] = problem.nonwettingPhase().intEnergy(varData[idx].temperature,varData[idx].pN);

        // CONSTANT solubility (for comparison with twophase)
        //         varData[idx].massfrac[nPhase][wPhase] = 0.0; varData[i].massfrac[wPhase][wPhase] = 1.0;
        //         varData[i].massfrac[wPhase][nPhase] = 0.0; varData[i].massfrac[nPhase][nPhase] = 1.0;

        //std::cout << "wPhase in gasphase: " << varData[i].massfrac[wPhase][nPhase] << std::endl;
        //std::cout << "nPhase in waterphase: " << varData[i].massfrac[co2][wPhase] << std::endl;

        // for output
        (*outPressureN)[global] = varData[idx].pN;
        (*outCapillaryP)[global] = varData[idx].pC;
        (*outSaturationW)[global] = varData[idx].satW;
        (*outSaturationN)[global] = varData[idx].satN;
        (*outTemperature)[global] = varData[idx].temperature;
        (*outDensityW)[global] = varData[idx].density[wPhase];
        (*outDensityN)[global] = varData[idx].density[nPhase];
        (*outMobilityW)[global] = varData[idx].mobility[wPhase];
        (*outMobilityN)[global] = varData[idx].mobility[nPhase];

        return;
    }

    void updateVariableData(const Element& element, const SolutionVector* sol, int idx, bool old = false)
    {
        if (old)
        {
            updateVariableData(element, sol, idx, oldVNDat);
        }
        else
            updateVariableData(element, sol, idx, vNDat);
    }

    void updateVariableData(const Element& element, const SolutionVector* sol, bool old = false)
    {
        int size = this->fvGeom.numVertices;

        for (int idx = 0; idx < size; idx++)
            updateVariableData(element, sol, idx, old);
    }

    struct StaticNodeData
    {
        bool visited;
        Scalar cellVolume;
        Scalar porosity;
    };

    struct ElementData
    {
        Scalar heatCap;
    }elData;

    // parameters given in constructor
    TwoPhaseHeatProblem<Grid,Scalar>& problem;
    std::vector<StaticNodeData> sNDat;
    std::vector<VariableNodeData> vNDat;
    std::vector<VariableNodeData> oldVNDat;

    // for output files
    BlockVector<FieldVector<Scalar, 1> > *outPressureN;
    BlockVector<FieldVector<Scalar, 1> > *outCapillaryP;
    BlockVector<FieldVector<Scalar, 1> > *outSaturationN;
    BlockVector<FieldVector<Scalar, 1> > *outSaturationW;
    BlockVector<FieldVector<Scalar, 1> > *outTemperature;
    BlockVector<FieldVector<Scalar, 1> > *outDensityW;
    BlockVector<FieldVector<Scalar, 1> > *outDensityN;
    BlockVector<FieldVector<Scalar, 1> > *outMobilityW;
    BlockVector<FieldVector<Scalar, 1> > *outMobilityN;
    BlockVector<FieldVector<Scalar, 1> > *outPhaseState;

};

/** @} */
}
#endif
