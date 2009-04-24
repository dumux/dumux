#ifndef DUNE_BOXSTOKESTRENJACOBIAN_HH
#define DUNE_BOXSTOKESTRENJACOBIAN_HH

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

#include<dumux/operators/boxjacobian.hh>
#include "dumux/stokes/stokestrenproblem.hh"

namespace Dune
{
template<class Grid, class Scalar, class BoxFunction = LeafP1Function<Grid, Scalar, Grid::dimension+3> >
class BoxStokesTrEnJacobian
    : public BoxJacobian<BoxStokesTrEnJacobian<Grid,Scalar,BoxFunction>,Grid,Scalar,Grid::dimension+3,BoxFunction>
{
    enum {dim=Grid::dimension};
    enum {numEq = dim+3};
    enum {SIZE=LagrangeShapeFunctionSetContainer<Scalar,Scalar,dim>::maxsize};

    typedef typename Grid::Traits::template Codim<0>::Entity Element;
    typedef typename Element::Geometry Geometry;
    typedef BoxStokesTrEnJacobian<Grid,Scalar,BoxFunction> ThisType;
    typedef typename LocalJacobian<ThisType,Grid,Scalar,numEq>::VBlockType SolutionVector;
    typedef BoxJacobian<ThisType,Grid,Scalar,numEq,BoxFunction> BoxJacobianType;
    typedef Dune::FVElementGeometry<Grid> FVElementGeometry;

    enum {nPhase = 0};

public:
    struct VariableNodeData;
    typedef FieldVector<Scalar,dim> FVector;
    typedef FieldMatrix<Scalar,dim,dim> FMatrix;

    //! Constructor
    BoxStokesTrEnJacobian (StokesTrEnProblem<Grid,Scalar>& params,
                           bool levelBoundaryAsDirichlet_, const Grid& grid,
                           BoxFunction& sol,
                           bool procBoundaryAsDirichlet_=true)
        : BoxJacobianType(levelBoundaryAsDirichlet_, grid, sol, procBoundaryAsDirichlet_),
          problem(params), varNData(SIZE), oldVarNData(SIZE)
    {
        alpha = -1.0;
        this->analytic = false;
    }

    void clearVisited ()
    {
        return;
    }

    void localDefect(const Element& element, const SolutionVector* sol, bool withBC = true) {
        BoxJacobianType::localDefect(element, sol, withBC);

        assembleBoundaryCondition(element);

        Dune::GeometryType gt = element.geometry().type();
        const typename ReferenceElementContainer<Scalar,dim>::value_type& referenceElement = ReferenceElements<Scalar, dim>::general(gt);

        for (int vert=0; vert < this->fvGeom.numVertices; vert++) // begin loop over vertices / sub control volumes
            if (!this->fvGeom.subContVol[vert].inner)
            {
                typedef typename Grid::LeafGridView::IntersectionIterator IntersectionIterator;

                FieldVector<Scalar,dim> averagedNormal(0);
                int faces = 0;
                IntersectionIterator endit = element.ileafend();
                for (IntersectionIterator it = element.ileafbegin(); it!=endit; ++it)
                {
                    // handle face on exterior boundary, this assumes there are no interior boundaries
                    if (it->boundary()) {
                        // get geometry type of face
                        GeometryType faceGT = it->geometryInInside().type();

                        // center in face's reference element
                        const FieldVector<Scalar,dim-1>& faceLocal = ReferenceElements<Scalar,dim-1>::general(faceGT).position(0,0);

                        int faceIdx = it->indexInInside();
                        int numVerticesOfFace = referenceElement.size(faceIdx, 1, dim);
                        for (int nodeInFace = 0; nodeInFace < numVerticesOfFace; nodeInFace++) {
                            int nodeInElement = referenceElement.subEntity(faceIdx, 1, nodeInFace, dim);
                            if (nodeInElement != vert)
                                continue;
                            int bfIdx = this->fvGeom.boundaryFaceIndex(faceIdx,    nodeInFace);
                            FieldVector<Scalar,dim> normal = it->unitOuterNormal(faceLocal);
                            normal *= this->fvGeom.boundaryFace[bfIdx].area;
                            averagedNormal += normal;
                            faces++;
                        }
                    }
                }
                Scalar defect = 0;
                if (averagedNormal.two_norm())
                    averagedNormal /= averagedNormal.two_norm();
                for (int k = 0; k < dim; k++)
                    defect += this->def[vert][k]*averagedNormal[k];
                //mass balance equation
                if (faces == 2 && this->fvGeom.numVertices == 4)
                    this->def[vert][dim+2] = sol[0][dim+2] + sol[3][dim+2] - sol[1][dim+2] - sol[2][dim+2];
                else if (this->bctype[vert][0] == BoundaryConditions::dirichlet)
                    this->def[vert][dim+2] = defect;
                else // de-stabilize
                {
                    for (int face = 0; face < this->fvGeom.numEdges; face++)
                    {
                        int i = this->fvGeom.subContVolFace[face].i;
                        int j = this->fvGeom.subContVolFace[face].j;

                        if (i != vert && j != vert)
                            continue;

                        FieldVector<Scalar, dim> pressGradient(0);
                        for (int k = 0; k < this->fvGeom.numVertices; k++) {
                            FieldVector<Scalar,dim> grad(this->fvGeom.subContVolFace[face].grad[k]);
                            grad *= sol[k][dim+2];
                            pressGradient += grad;
                        }

                        Scalar alphaH2 = 0.5*alpha*(this->fvGeom.subContVol[i].volume + this->fvGeom.subContVol[j].volume);
                        pressGradient *= alphaH2;
                        if (i == vert)
                            this->def[vert][dim+2] += pressGradient*this->fvGeom.subContVolFace[face].normal;
                        else
                            this->def[vert][dim+2] -= pressGradient*this->fvGeom.subContVolFace[face].normal;
                    }
                }

                //std::cout << "node " << i << ", coord = " << this->fvGeom.subContVol[i].global << ", N = " << averagedNormal << std::endl;
            }

        return;
    }

    // harmonic mean of the permeability computed directly
    virtual FMatrix harmonicMeanK (FMatrix& Ki, const FMatrix& Kj) const
    {
        double eps = 1e-20;

        for (int kx=0; kx<dim; kx++)
            for (int ky=0; ky<dim; ky++)
                if (Ki[kx][ky] != Kj[kx][ky])
                    Ki[kx][ky] = 2 / (1/(Ki[kx][ky]+eps) + (1/(Kj[kx][ky]+eps)));

        return Ki;
    }


    // compute storage term
    SolutionVector computeM (const Element& element, const SolutionVector* sol, int node, const std::vector<VariableNodeData>& varData)
    {
        SolutionVector result(0);

        //velocity u
        result[0] = -1*varData[node].density*sol[node][0];
        //velocity v
        result[1] = -1*varData[node].density*sol[node][1];
        //partial density
        result[2] = -1*sol[node][2];
        //temperature T
        result[3] = -1*varData[node].density*varData[node].internalenergy;
        //pressure p
        result[4] = -1*varData[node].density;

        //          std::cout << "node " <<  node << " time dep = " << result << std::endl;

        return result;
    };

    SolutionVector computeM (const Element& element, const SolutionVector* sol, int node, bool old = false)
    {
        if (old)
            return computeM(element, sol, node, oldVarNData);
        else
            return computeM(element, sol, node, varNData);
    }

    SolutionVector computeQ (const Element& element, const SolutionVector* sol, const int& node)
    {
        SolutionVector result = problem.q(this->fvGeom.subContVol[node].global, element, this->fvGeom.subContVol[node].local);
        result *= -1.0;

        Scalar gravity = problem.gravity()*varNData[node].density;
        result[dim-1] += gravity;

        Scalar alphaH2 = alpha*this->fvGeom.subContVol[node].volume;
        //        result[dim+2] *= alphaH2;

        Scalar MassST = problem.Qg(this->fvGeom.subContVol[node].global, element, this->fvGeom.subContVol[node].local);
        MassST *= -1;
        MassST *= alphaH2;

        result[dim+2] += MassST;

        SolutionVector flux = boundaryFlux(element, sol, node);

        flux /= this->fvGeom.subContVol[node].volume;
        result -= flux;

        return (result);
    }

    SolutionVector computeA (const Element& element, const SolutionVector* sol, int face)
    {
        SolutionVector flux(0);

        Scalar pressValue = 0;
        FieldVector<Scalar, dim> velocityValue(0);
        FieldVector<Scalar, dim> xV(0);
        FieldVector<Scalar, dim> rhoV(0);
        FieldVector<Scalar, dim> gradX(0);
        FieldVector<Scalar, dim> gradP(0);
        FieldVector<Scalar, dim> gradT(0);
        FieldVector<Scalar, dim> rhV(0);

        // get global coordinates of nodes i,j
        int i = this->fvGeom.subContVolFace[face].i;
        int j = this->fvGeom.subContVolFace[face].j;

        const FieldVector<Scalar,dim> global_i = this->fvGeom.subContVol[i].global;
        const FieldVector<Scalar,dim> global_j = this->fvGeom.subContVol[j].global;

        const FieldVector<Scalar,dim> local_i = this->fvGeom.subContVol[i].local;
        const FieldVector<Scalar,dim> local_j = this->fvGeom.subContVol[j].local;

        FMatrix Di(0), Dj(0);

        // calculate harmonic mean of permeabilities of nodes i and j
        Di = this->problem.D(global_i,element,local_i);
        Dj = this->problem.D(global_j,element,local_j);

        const FMatrix D = harmonicMeanK(Di, Dj);

        for (int k = 0; k < this->fvGeom.numVertices; k++)
        {
            pressValue += sol[k][dim+2]*this->fvGeom.subContVolFace[face].shapeValue[k];

            FieldVector<Scalar,dim> grad(this->fvGeom.subContVolFace[face].grad[k]);
            grad *= sol[k][dim];
            gradX += grad;

            grad = this->fvGeom.subContVolFace[face].grad[k];
            grad *= sol[k][dim+2];
            gradP += grad;

            grad = this->fvGeom.subContVolFace[face].grad[k];
            grad *= sol[k][dim+1];
            grad *= varNData[k].heatconductivity;
            gradT -= grad;

            for (int comp = 0; comp < dim; comp++)
            {
                velocityValue[comp] += sol[k][comp]*this->fvGeom.subContVolFace[face].shapeValue[k];
                rhoV[comp] += varNData[k].density*sol[k][comp]*this->fvGeom.subContVolFace[face].shapeValue[k];
            }
        }

        Scalar xValue;
        Scalar rhValue;
        Scalar outward = velocityValue*this->fvGeom.subContVolFace[face].normal;
        if (outward > 0)
        {
            xValue = sol[i][dim];
            rhValue = varNData[i].density*varNData[i].enthalpy;
        }
        else
        {
            xValue = sol[j][dim];
            rhValue = varNData[j].density*varNData[j].enthalpy;
        }

        for (int comp = 0; comp < dim; comp++)
        {
            xV[comp] = velocityValue[comp]*xValue;
            rhV[comp] = velocityValue[comp]*rhValue;
        }

        FieldVector<Scalar,dim> DgradX(0);
        D.mv(gradX, DgradX);  // DgradX=D*gradX

        // momentum balance:
        for (int comp = 0; comp < dim; comp++)
        {
            FieldVector<Scalar, dim> gradVComp(0);
            for (int k = 0; k < this->fvGeom.numVertices; k++)
            {
                FieldVector<Scalar,dim> grad(this->fvGeom.subContVolFace[face].grad[k]);
                grad *= sol[k][comp];
                grad *= varNData[k].viscosity;
                gradVComp -= grad;
            }

            FieldVector<Scalar, dim> pComp(0);
            pComp[comp] = pressValue;
            gradVComp += pComp;
            flux[comp] = gradVComp*this->fvGeom.subContVolFace[face].normal;
        }

        //transport
        //flux[dim] = (DgradX - xV)*this->fvGeom.subContVolFace[face].normal;
        flux[dim] = (xV - DgradX)*this->fvGeom.subContVolFace[face].normal;

        // mass balance:
        Scalar alphaH2 = 0.5*alpha*(this->fvGeom.subContVol[i].volume + this->fvGeom.subContVol[j].volume);
        gradP *= alphaH2;
        rhoV += gradP;
        flux[dim+2] = rhoV*this->fvGeom.subContVolFace[face].normal;

        flux[dim+1] = (gradT + rhV)*this->fvGeom.subContVolFace[face].normal;

        return flux;
    }

    void computeElementData (const Element& element)
    {
        // ASSUMING element-wise constant viscosity, evaluate mu at the cell center
        // elData.mu = problem.mu(this->fvGeom.elementGlobal, element, this->fvGeom.elementLocal);
    };

    // the members of the struct are defined here
    struct VariableNodeData
    {
        Scalar pN;
        Scalar T;
        Scalar Xg;
        Scalar viscosity;
        Scalar density;
        Scalar massfrac;
        Scalar internalenergy;
        Scalar enthalpy;
        Scalar heatconductivity;
        FieldMatrix<Scalar,dim,dim> D;
    };

    void updateVariableData(const Element& element, const SolutionVector* sol, int i, std::vector<VariableNodeData>& varData)
    {
        //analytic problem
        varData[i].pN = sol[i][4];
        varData[i].density = problem.density(this->fvGeom.elementGlobal, element, this->fvGeom.elementLocal);
        varData[i].viscosity = 1.0;
        varData[i].internalenergy = problem.internalenergy(this->fvGeom.elementGlobal, element, this->fvGeom.elementLocal);
        varData[i].enthalpy = problem.enthalpy(this->fvGeom.elementGlobal, element, this->fvGeom.elementLocal);
        varData[i].heatconductivity = problem.heatConductivity(this->fvGeom.elementGlobal, element, this->fvGeom.elementLocal);

        /*
        //real problem
        varData[i].pN = sol[i][4];
        varData[i].T = sol[i][dim+1];
        varData[i].Xg = sol[i][dim];
        varData[i].density = problem.gasPhase().density(varData[i].T, varData[i].pN, varData[i].Xg);
        varData[i].viscosity = problem.gasPhase().viscosity(varData[i].T, varData[i].pN, varData[i].Xg);
        varData[i].internalenergy = problem.gasPhase().intEnergy(varData[i].T, varData[i].pN, varData[i].Xg);
        varData[i].enthalpy = problem.gasPhase().enthalpy(varData[i].T, varData[i].pN, varData[i].Xg);
        varData[i].heatconductivity = problem.heatConductivity(this->fvGeom.elementGlobal, element, this->fvGeom.elementLocal);
        */

        //         std::cout << "enthalpy, " << i << "   = " << varData[i].enthalpy << std::endl;
        //         std::cout << "internal energy, " << i << "   = " << varData[i].internalenergy << std::endl;
        // for output
        //         (*outDensityW)[globalIdx] = varData[i].density[wPhase];
        //         (*outDensityN)[globalIdx] = varData[i].density[nPhase];

        return;
    }

    void updateVariableData(const Element& element, const SolutionVector* sol, int i, bool old = false)
    {
        if (old)
            updateVariableData(element, sol, i, oldVarNData);
        else
            updateVariableData(element, sol, i, varNData);
    }

    void updateVariableData(const Element& element, const SolutionVector* sol, bool old = false)
    {
        int size = this->fvGeom.numVertices;

        for (int i = 0; i < size; i++)
            updateVariableData(element, sol, i, old);
    }

    virtual void updateStaticData (const Element& element, const SolutionVector* sol)
    {
        return;
    }

    //    SolutionVector boundaryFlux(const Element& element, const SolutionVector* sol, int node, std::vector<VariableNodeData>& varData) {
    SolutionVector boundaryFlux(const Element& element, const SolutionVector* sol, int node) {
        SolutionVector  result(0);

        Dune::GeometryType gt = element.geometry().type();
        const typename Dune::LagrangeShapeFunctionSetContainer<Scalar,Scalar,dim>::value_type
            &sfs=Dune::LagrangeShapeFunctions<Scalar, Scalar, dim>::general(gt, 1);
        setcurrentsize(sfs.size());
        this->fvGeom.update(element);

        const typename ReferenceElementContainer<Scalar,dim>::value_type
            &referenceElement = ReferenceElements<Scalar, dim>::general(gt);

        // evaluate boundary conditions via intersection iterator
        typedef typename Grid::LeafGridView::IntersectionIterator IntersectionIterator;

        IntersectionIterator endit = element.ileafend();
        for (IntersectionIterator it = element.ileafbegin(); it!=endit; ++it)
        {
            // if we have a neighbor then we assume there is no boundary (forget interior boundaries)
            // in level assemble treat non-level neighbors as boundary
            if (it->neighbor()) {
                if (this->levelBoundaryAsDirichlet && it->outside()->level()==element.level())
                    continue;
                if (!this->levelBoundaryAsDirichlet)
                    continue;
            }

            // handle face on exterior boundary, this assumes there are no interior boundaries
            if (it->boundary()) {
                int faceIdx = it->indexInInside();

                int numVerticesOfFace = referenceElement.size(faceIdx, 1, dim);
                for (int nodeInFace = 0; nodeInFace < numVerticesOfFace; nodeInFace++) {
                    int nodeInElement = referenceElement.subEntity(faceIdx, 1, nodeInFace, dim);
                    if (node != nodeInElement)
                        continue;

                    int bfIdx = this->fvGeom.boundaryFaceIndex(faceIdx,    nodeInFace);

                    // get geometry type of face
                    GeometryType faceGT = it->geometryInInside().type();

                    // center in face's reference element
                    const FieldVector<Scalar,dim-1>& faceLocal = ReferenceElements<Scalar,dim-1>::general(faceGT).position(0,0);

                    Scalar pressureValue = 0;
                    FieldVector<Scalar,dim>  velocityValue(0);
                    FieldVector<Scalar,dim>  rhoV(0);
                    FieldMatrix<Scalar,dim,dim> velocityGradient(0);

                    for (int vert = 0; vert < this->fvGeom.numVertices; vert++)
                    {
                        pressureValue += sol[vert][dim+2]*this->fvGeom.boundaryFace[bfIdx].shapeValue[vert];
                        FieldVector<Scalar,dim> grad(this->fvGeom.boundaryFace[bfIdx].grad[vert]);
                        //                        grad *= sol[vert][dim+2];
                        for (int comp = 0; comp < dim; comp++)
                        {
                            velocityValue[comp] += sol[vert][comp]*this->fvGeom.boundaryFace[bfIdx].shapeValue[vert];
                            rhoV[comp] += varNData[vert].density*sol[vert][comp]*this->fvGeom.boundaryFace[bfIdx].shapeValue[vert];
                            FieldVector<Scalar,dim> gradVComp = grad;
                            gradVComp *= sol[vert][comp];
                            gradVComp *= varNData[vert].viscosity;
                            velocityGradient[comp] += gradVComp;
                        }
                    }
                    FieldVector<Scalar, dim> massResidual = rhoV;
                    Scalar alphaH2 = alpha*this->fvGeom.subContVol[node].volume;
                    SolutionVector source = problem.q(this->fvGeom.boundaryFace[bfIdx].ipGlobal, element, this->fvGeom.boundaryFace[bfIdx].ipLocal);
                    FieldVector<Scalar, dim> dimSource;
                    for (int comp = 0; comp < dim; comp++)
                        dimSource[comp] = source[comp];
                    dimSource *= alphaH2;
                    massResidual += dimSource;

                    result[dim+2] -= massResidual*it->unitOuterNormal(faceLocal)*this->fvGeom.boundaryFace[bfIdx].area;

                    FieldVector<Scalar,dim> gradVN(0);
                    velocityGradient.umv(it->unitOuterNormal(faceLocal), gradVN);
                    gradVN *= this->fvGeom.boundaryFace[bfIdx].area;

                    for (int comp = 0; comp < dim; comp++)
                    {
                        FieldVector<Scalar,dim>  pressVector(0);
                        pressVector[comp] = -pressureValue;
                        result[comp] += pressVector*it->unitOuterNormal(faceLocal)*this->fvGeom.boundaryFace[bfIdx].area;
                        if (alpha == 0)
                            result[comp] += gradVN[comp];
                    }

                    //transport
                    Scalar pressValue = 0;
                    FieldVector<Scalar, dim> xV(0);
                    FieldVector<Scalar, dim> gradX(0);

                    const FieldVector<Scalar,dim> local = this->fvGeom.boundaryFace[bfIdx].ipLocal;
                    const FieldVector<Scalar,dim> global = this->fvGeom.boundaryFace[bfIdx].ipGlobal;

                    const FMatrix D = this->problem.D(global,element,local);

                    for (int k = 0; k < this->fvGeom.numVertices; k++)
                    {
                        pressValue += sol[k][dim+2]*this->fvGeom.boundaryFace[bfIdx].shapeValue[k];
                        FieldVector<Scalar,dim> grad(this->fvGeom.boundaryFace[bfIdx].grad[k]);
                        grad *= sol[k][dim];
                        gradX += grad;
                    }

                    // NOT CORRECT YET
                    Scalar xValue;
                    Scalar outward = velocityValue*this->fvGeom.boundaryFace[bfIdx].normal;
                    if (outward > 0)
                        //if (outward < 0)
                        xValue = sol[node][dim];
                    else
                        xValue = sol[node][dim];

                    for (int comp = 0; comp < dim; comp++)
                        xV[comp] = velocityValue[comp]*xValue;

                    FieldVector<Scalar,dim> DgradX(0);
                    D.mv(gradX, DgradX);  // DgradX=D*gradX

                    //result[dim] = (DgradX - xV)*it->unitOuterNormal(faceLocal)*this->fvGeom.boundaryFace[bfIdx].area;
                    result[dim] = (xV - DgradX)*it->unitOuterNormal(faceLocal)*this->fvGeom.boundaryFace[bfIdx].area;

                    //energy
                    FieldVector<Scalar, dim> rhV(0);
                    FieldVector<Scalar, dim> gradT(0);

                    for (int k = 0; k < this->fvGeom.numVertices; k++)
                    {
                        FieldVector<Scalar,dim> grad(this->fvGeom.boundaryFace[bfIdx].grad[k]);
                        grad *= sol[k][dim+1];
                        gradT -= grad;
                        gradT *= varNData[k].heatconductivity;
                    }

                    // NOT CORRECT YET
                    Scalar rhValue = varNData[node].density*varNData[node].enthalpy;;

                    for (int comp = 0; comp < dim; comp++)
                        xV[comp] = velocityValue[comp]*rhValue;

                    result[dim] = (gradT + rhV)*it->unitOuterNormal(faceLocal)*this->fvGeom.boundaryFace[bfIdx].area;



                }
            }
        }
        return result;
    }


    void assembleBoundaryCondition(const Element& element, int k = 1) {
        Dune::GeometryType gt = element.geometry().type();
        const typename Dune::LagrangeShapeFunctionSetContainer<Scalar,Scalar,dim>::value_type
            &sfs=Dune::LagrangeShapeFunctions<Scalar, Scalar, dim>::general(gt, 1);
        setcurrentsize(sfs.size());
        this->fvGeom.update(element);

        const typename ReferenceElementContainer<Scalar,dim>::value_type
            &referenceElement = ReferenceElements<Scalar, dim>::general(gt);

        for (int i = 0; i < sfs.size(); i++) {
            this->bctype[i].assign(BoundaryConditions::neumann);
            this->b[i] = 0;
        }

        // evaluate boundary conditions via intersection iterator
        typedef typename Grid::LeafGridView::IntersectionIterator IntersectionIterator;

        IntersectionIterator endit = element.ileafend();
        for (IntersectionIterator it = element.ileafbegin(); it!=endit; ++it)
        {
            // if we have a neighbor then we assume there is no boundary (forget interior boundaries)
            // in level assemble treat non-level neighbors as boundary
            if (it->neighbor()) {
                if (this->levelBoundaryAsDirichlet && it->outside()->level()==element.level())
                    continue;
                if (!this->levelBoundaryAsDirichlet)
                    continue;
            }

            // determine boundary condition type for this face, initialize with processor boundary
            FieldVector<typename BoundaryConditions::Flags, numEq> bctypeface(BoundaryConditions::process);
            FieldVector<int,numEq> dirichletIdx(0);

            // handle face on exterior boundary, this assumes there are no interior boundaries
            if (it->boundary()) {
                int faceIdx = it->indexInInside();
                int numVerticesOfFace = referenceElement.size(faceIdx, 1, dim);
                for (int nodeInFace = 0; nodeInFace < numVerticesOfFace; nodeInFace++) {
                    int nodeInElement = referenceElement.subEntity(faceIdx, 1, nodeInFace, dim);
                    for (int equationNumber = 0; equationNumber < numEq; equationNumber++) {
                        if (this->bctype[nodeInElement][equationNumber] == BoundaryConditions::neumann) {
                            int bfIdx = this->fvGeom.boundaryFaceIndex(faceIdx,    nodeInFace);
                            FieldVector<Scalar,dim> local = this->fvGeom.boundaryFace[bfIdx].ipLocal;
                            FieldVector<Scalar,dim> global = this->fvGeom.boundaryFace[bfIdx].ipGlobal;
                            bctypeface = this->getImp().problem.bctype(global, element, it, local); // eval bctype
                            this->getImp().problem.dirichletIndex(global, element, it, local, dirichletIdx); // eval bctype
                            if (bctypeface[equationNumber]!=BoundaryConditions::neumann)
                                break;
                            FieldVector<Scalar,dim+3> J = this->getImp().problem.neumann(global, element, it, local);
                            if (equationNumber < dim+2) {
                                J[equationNumber] *= this->fvGeom.boundaryFace[bfIdx].area;
                                this->b[nodeInElement][equationNumber] += J[equationNumber];
                            }
                        }
                    }
                }

                bool nface(true); // check if face is a neumann face
                for(int i=0; i<numEq; i++)
                {
                    if(bctypeface[i] != BoundaryConditions::neumann)
                        nface = false; // was not a neumann face
                }
                if(nface == true)
                    continue; // was a neumann face, go to next face
            }

            // If we are here, then it is
            // (i)   an exterior boundary face with Dirichlet condition, or
            // (ii)  a processor boundary (i.element. neither boundary() nor neighbor() was true), or
            // (iii) a level boundary in case of level-wise assemble
            // How processor boundaries are handled depends on the processor boundary mode

            bool pface(false);  // check if face is a process boundary
            for(int i=0; i<numEq; i++)
            {
                if (bctypeface[i]==BoundaryConditions::process
                    && this->procBoundaryAsDirichlet==false
                    && this->levelBoundaryAsDirichlet==false)
                {
                    pface = true;
                    break;
                }
            }
            if(pface == true)
                continue;   // if face is a process boundary it acts like homogeneous Neumann


            for (int equationNumber=0; equationNumber<dim+2; equationNumber++) {
                for (int i=0; i<sfs.size(); i++) // loop over test function number
                {
                    //this->dirichletIndex[i][equationNumber] = equationNumber;

                    //std::cout<<"i = "<<i<<std::endl;
                    if (sfs[i].codim()==0)
                        continue; // skip interior dof
                    if (sfs[i].codim()==1) // handle face dofs
                    {
                        if (sfs[i].entity()==it->indexInInside()) {
                            if (this->bctype[i][equationNumber] < bctypeface[equationNumber]) {
                                this->bctype[i][equationNumber] = bctypeface[equationNumber];

                                if (bctypeface[equationNumber] == BoundaryConditions::process)
                                    this->b[i][equationNumber] = 0;
                                if (bctypeface[equationNumber] == BoundaryConditions::dirichlet) {
                                    this->b[i][equationNumber] = 0;
                                }
                            }
                        }
                        continue;
                    }
                    // handle subentities of this face
                    for (int j=0; j<ReferenceElements<Scalar,dim>::general(gt).size(it->indexInInside(), 1, sfs[i].codim()); j++)
                        if (sfs[i].entity()==ReferenceElements<Scalar,dim>::general(gt).subEntity(it->indexInInside(), 1, j, sfs[i].codim()))
                        {
                            if (this->bctype[i][equationNumber] < bctypeface[equationNumber]) {
                                this->bctype[i][equationNumber] = bctypeface[equationNumber];
                               if (bctypeface[equationNumber] == BoundaryConditions::process)
                                    this->b[i][equationNumber] = 0;
                                if (bctypeface[equationNumber] == BoundaryConditions::dirichlet) {
                                    this->b[i][equationNumber] = 0;
                                }
                            }
                        }
                }
            }
        }

    }

    struct StaticNodeData
    {
        FMatrix K;
    };


    struct ElementData {
        Scalar mu;
        Scalar rho;
        Scalar gravity;
    };

    ElementData elData;
    StokesTrEnProblem<Grid,Scalar>& problem;
    double alpha;

    std::vector<VariableNodeData> varNData;
    std::vector<VariableNodeData> oldVarNData;

};
}
#endif
