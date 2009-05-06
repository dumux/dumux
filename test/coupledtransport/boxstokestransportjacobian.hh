#ifndef DUNE_BOXSTOKESTRANSPORTJACOBIAN_HH
#define DUNE_BOXSTOKESTRANSPORTJACOBIAN_HH

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
#include "dumux/stokes/stokestransportproblem.hh"

namespace Dune
{
template<class Grid, class Scalar, class BoxFunction = LeafP1Function<Grid, Scalar, Grid::dimension+2> >
class BoxStokesTransportJacobian
    : public BoxJacobian<BoxStokesTransportJacobian<Grid,Scalar,BoxFunction>,Grid,Scalar,Grid::dimension+2,BoxFunction>
{
    enum {dim=Grid::dimension};
    enum {numEq = dim+2};
    enum {SIZE=LagrangeShapeFunctionSetContainer<Scalar,Scalar,dim>::maxsize};
    enum {velocityXIdx=0, velocityYIdx=1, velocityZIdx=2, massFracIdx=dim, pressureIdx=dim+1};

    typedef typename Grid::Traits::template Codim<0>::Entity Element;
    typedef typename Grid::LeafGridView::IntersectionIterator IntersectionIterator;
    typedef typename Element::Geometry Geometry;
    typedef BoxStokesTransportJacobian<Grid,Scalar,BoxFunction> ThisType;
    typedef typename LocalJacobian<ThisType,Grid,Scalar,numEq>::VBlockType SolutionVector;
    typedef BoxJacobian<ThisType,Grid,Scalar,numEq,BoxFunction> BoxJacobianType;
    typedef Dune::FVElementGeometry<Grid> FVElementGeometry;

    enum {nPhase = 0};

public:
    struct VariableNodeData;
    typedef FieldVector<Scalar,dim> FVector;
    typedef FieldMatrix<Scalar,dim,dim> FMatrix;

    //! Constructor
    BoxStokesTransportJacobian (StokesTransportProblem<Grid,Scalar>& params,
                                bool levelBoundaryAsDirichlet_, const Grid& grid,
                                BoxFunction& sol,
                                bool procBoundaryAsDirichlet_=true)
        : BoxJacobianType(levelBoundaryAsDirichlet_, grid, sol, procBoundaryAsDirichlet_),
          problem(params), varNData(SIZE), oldVarNData(SIZE)
    {
	// factor for stabilization
	// set to zero, when no stabilization is neccessary
        alpha = 0; // -1e3;//-1.0;
        this->analytic = false;
    }

    void clearVisited ()
    {
        return;
    }

    void localDefect(const Element& element, const SolutionVector* sol, bool withBC = true)
    {
        BoxJacobianType::localDefect(element, sol, withBC);

        this->getImp().assembleBoundaryCondition(element);

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
                    // handle face on exterior boundary, this assumes there are no interior boundaries
                    if (it->boundary()) {
                        // get geometry type of face
                        GeometryType faceGT = it->geometryInInside().type();

                        // center in face's reference element
                        const FieldVector<Scalar,dim-1>& faceLocal = ReferenceElements<Scalar,dim-1>::general(faceGT).position(0,0);

                        int faceIdx = it->indexInInside();
                        int numVerticesOfFace = referenceElement.size(faceIdx, 1, dim);
                        for (int nodeInFace = 0; nodeInFace < numVerticesOfFace; nodeInFace++)
                        {
                            int nodeInElement = referenceElement.subEntity(faceIdx, 1, nodeInFace, dim);
                            if (nodeInElement != vert)
                                continue;
                            int bfIdx = this->fvGeom.boundaryFaceIndex(faceIdx,    nodeInFace);
                            FieldVector<Scalar,dim> local = this->fvGeom.boundaryFace[bfIdx].ipLocal;
                            FieldVector<Scalar,dim> global = this->fvGeom.boundaryFace[bfIdx].ipGlobal;
//                            FieldVector<BoundaryConditions::Flags,numEq> bctypeface = this->getImp().problem.bctype(global, element, it, local);

                            // Checks if boundary conditions of ALL equations are DIRICHLET
                            bool onlyDirichlet = true;
                            FieldVector<BoundaryConditions::Flags, numEq>  bctypeNode = this->getImp().problem.bctype(global, element, it, local);
                            for (int i=0; i<dim; ++i)
                            {
                            	if (bctypeNode[i] != BoundaryConditions::neumann)
                            		onlyDirichlet = false;
                            }

                            // TODO: is this always valid?? Checks only bctype of first primary variable
                            if (onlyDirichlet)
//                            if (bctypeface[velocityXIdx] == BoundaryConditions::dirichlet)
                            {
                            	FieldVector<Scalar,dim> normal = it->unitOuterNormal(faceLocal);
                            	normal *= this->fvGeom.boundaryFace[bfIdx].area;
                            	averagedNormal += normal;
                           	}
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
                {
//                    this->def[vert][velocityXIdx] = sol[0][velocityXIdx] + sol[3][velocityXIdx] - sol[1][velocityXIdx] - sol[2][velocityXIdx];
//                    this->def[vert][velocityYIdx] = sol[0][velocityYIdx] + sol[3][velocityYIdx] - sol[1][velocityYIdx] - sol[2][velocityYIdx];
//                    this->def[vert][massFracIdx] = sol[0][massFracIdx] + sol[3][massFracIdx] - sol[1][massFracIdx] - sol[2][massFracIdx];
                    this->def[vert][pressureIdx] = sol[0][pressureIdx] + sol[3][pressureIdx] - sol[1][pressureIdx] - sol[2][pressureIdx];
                }
                else if (this->bctype[vert][0] == BoundaryConditions::dirichlet)
                    this->def[vert][pressureIdx] = defect;
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
                            grad *= sol[k][pressureIdx];
                            pressGradient += grad;
                        }

                        Scalar alphaH2 = 0.5*alpha*(this->fvGeom.subContVol[i].volume + this->fvGeom.subContVol[j].volume);
                        pressGradient *= alphaH2;
                        if (i == vert)
                            this->def[vert][pressureIdx] -= pressGradient*this->fvGeom.subContVolFace[face].normal;
                        else
                            this->def[vert][pressureIdx] += pressGradient*this->fvGeom.subContVolFace[face].normal;
                    }
                }

                //std::cout << "node " << i << ", coord = " << this->fvGeom.subContVol[i].global << ", N = " << averagedNormal << std::endl;
            }

        return;
    }

    // harmonic mean of the permeability computed directly
    virtual FMatrix harmonicMeanK (FMatrix& Ki, const FMatrix& Kj) const
    {
        for (int kx=0; kx<dim; kx++)
            for (int ky=0; ky<dim; ky++)
                if (Ki[kx][ky] != Kj[kx][ky])
		    Ki[kx][ky] = 2*Ki[kx][ky]*Kj[kx][ky] / (Ki[kx][ky]+Kj[kx][ky]);
        return Ki;
    }


    // compute storage term
    SolutionVector computeStorage (const Element& element, const SolutionVector* sol, int node, const std::vector<VariableNodeData>& varData)
    {
        SolutionVector result(0);

        //velocity u
        result[velocityXIdx] = -1*varData[node].density*sol[node][velocityXIdx];
        //velocity v
        result[velocityYIdx] = -1*varData[node].density*sol[node][velocityYIdx];
        //partial density
        result[massFracIdx] = -1*sol[node][massFracIdx];
        //pressure p
        result[pressureIdx] = -1*varData[node].density;

        //          std::cout << "node " <<  node << " time dep = " << result << std::endl;

        return result;
    };

    SolutionVector computeStorage (const Element& element, const SolutionVector* sol, int node, bool old = false)
    {
        if (old)
            return computeStorage(element, sol, node, oldVarNData);
        else
            return computeStorage(element, sol, node, varNData);
    }

    SolutionVector computeSource (const Element& element, const SolutionVector* sol, const int& node)
    {
        SolutionVector result = problem.q(this->fvGeom.subContVol[node].global, element, this->fvGeom.subContVol[node].local);
        result *= -1.0;

        FieldVector<Scalar,dim> gravityVector(problem.gravity(this->fvGeom.subContVol[node].global));

        gravityVector[dim-1] *= varNData[node].density;

        result[dim-1] -= gravityVector[dim-1]; //sign???

        Scalar alphaH2 = alpha*this->fvGeom.subContVol[node].volume;
//        result[dim+1] *= alphaH2;

        Scalar MassST = problem.Qg(this->fvGeom.subContVol[node].global, element, this->fvGeom.subContVol[node].local);
        MassST *= -1;
        MassST *= alphaH2;

        result[pressureIdx] += MassST;

        SolutionVector flux = boundaryFlux(element, sol, node);

        flux /= this->fvGeom.subContVol[node].volume;
        result -= flux;

        return (result);
    }

    SolutionVector computeFlux (const Element& element, const SolutionVector* sol, int face)
    {
        SolutionVector flux(0);

        Scalar pressureValue = 0;
        FieldVector<Scalar, dim> velocityValue(0);
        FieldVector<Scalar, dim> xV(0);
        FieldVector<Scalar, dim> rhoV(0);
        FieldVector<Scalar, dim> gradRhoX(0);
        FieldVector<Scalar, dim> gradP(0);

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
            pressureValue += sol[k][pressureIdx]*this->fvGeom.subContVolFace[face].shapeValue[k];
            FieldVector<Scalar,dim> grad(this->fvGeom.subContVolFace[face].grad[k]);
            grad *= sol[k][massFracIdx];
            grad *= varNData[k].density;
            gradRhoX += grad;
            grad = this->fvGeom.subContVolFace[face].grad[k];
            grad *= sol[k][pressureIdx];
            gradP += grad;

            for (int comp = 0; comp < dim; comp++)
            {
                velocityValue[comp] += sol[k][comp]*this->fvGeom.subContVolFace[face].shapeValue[k];
                rhoV[comp] += varNData[k].density*sol[k][comp]*this->fvGeom.subContVolFace[face].shapeValue[k];
            }
        }

        Scalar xValue;
        Scalar outward = velocityValue*this->fvGeom.subContVolFace[face].normal;
        if (outward > 0)
            xValue = sol[i][massFracIdx];
        else
            xValue = sol[j][massFracIdx];

        for (int comp = 0; comp < dim; comp++)
            xV[comp] = rhoV[comp]*xValue;

        FieldVector<Scalar,dim> DgradX(0);
        D.mv(gradRhoX, DgradX);  // DgradX=D*gradX

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
            pComp[comp] = pressureValue;
            gradVComp += pComp;
            flux[comp] = gradVComp*this->fvGeom.subContVolFace[face].normal;
        }

        //transport
        flux[massFracIdx] = (xV - DgradX)*this->fvGeom.subContVolFace[face].normal;

        // mass balance:
        Scalar alphaH2 = 0.5*alpha*(this->fvGeom.subContVol[i].volume + this->fvGeom.subContVol[j].volume);
        gradP *= alphaH2;
        rhoV += gradP;
        flux[pressureIdx] = rhoV*this->fvGeom.subContVolFace[face].normal;

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
        Scalar viscosity;
        Scalar density;
        Scalar massfrac;
        Scalar partialpressure;
        FieldMatrix<Scalar,dim,dim> D;
        Scalar Xg;
    };

    void updateVariableData(const Element& element, const SolutionVector* sol, int i, std::vector<VariableNodeData>& varData)
    {
    	 varData[i].pN = sol[i][3];
    	 varData[i].Xg = sol[i][2];
    	 varData[i].density = problem.density(this->fvGeom.elementGlobal, element, this->fvGeom.elementLocal);
    	 varData[i].viscosity = problem.viscosity(this->fvGeom.elementGlobal, element, this->fvGeom.elementLocal);

        //         std::cout << "node " << i << "   " << varData[i].density << std::endl;
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
    SolutionVector boundaryFlux(const Element& element, const SolutionVector* sol, int node)
    {
        SolutionVector  result(0);

        Dune::GeometryType gt = element.geometry().type();
        const typename Dune::LagrangeShapeFunctionSetContainer<Scalar,Scalar,dim>::value_type
            &sfs=Dune::LagrangeShapeFunctions<Scalar, Scalar, dim>::general(gt, 1);
        setcurrentsize(sfs.size());
        this->fvGeom.update(element);

        const typename ReferenceElementContainer<Scalar,dim>::value_type
            &referenceElement = ReferenceElements<Scalar, dim>::general(gt);

        // evaluate boundary conditions via intersection iterator
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

                    FieldVector<typename BoundaryConditions::Flags, numEq> bctypeface(BoundaryConditions::process);
                    bctypeface = this->getImp().problem.bctype(this->fvGeom.boundaryFace[bfIdx].ipGlobal, element, it, this->fvGeom.boundaryFace[bfIdx].ipLocal);

                    // get geometry type of face
                    GeometryType faceGT = it->geometryInInside().type();

                    // center in face's reference element
                    const FieldVector<Scalar,dim-1>& faceLocal = ReferenceElements<Scalar,dim-1>::general(faceGT).position(0,0);

                    Scalar pressureValue = 0;
                    FieldVector<Scalar,dim>  velocityValue(0);
                    FieldVector<Scalar,dim>  xV(0);
                    FieldVector<Scalar,dim>  rhoV(0);
                    FieldVector<Scalar,dim>  gradRhoX(0);
                    FieldVector<Scalar,dim>  gradP(0);

                    const FieldVector<Scalar,dim> local = this->fvGeom.boundaryFace[bfIdx].ipLocal;
                    const FieldVector<Scalar,dim> global = this->fvGeom.boundaryFace[bfIdx].ipGlobal;
                    const FMatrix D = this->problem.D(global,element,local);

                    FMatrix velocityGradient(0);

                    for (int k = 0; k < this->fvGeom.numVertices; k++)
                    {
                        pressureValue += sol[k][dim+1]*this->fvGeom.boundaryFace[bfIdx].shapeValue[k];
                        FieldVector<Scalar,dim> grad(this->fvGeom.boundaryFace[bfIdx].grad[k]);
                        grad *= sol[k][dim];
                        grad *= varNData[k].density;
                        gradRhoX += grad;

                        grad = this->fvGeom.boundaryFace[bfIdx].grad[k];
                        grad *= sol[k][dim+1];
						gradP += grad;

						grad = this->fvGeom.boundaryFace[bfIdx].grad[k];

                        for (int comp = 0; comp < dim; comp++)
                        {
                            velocityValue[comp] += sol[k][comp]*this->fvGeom.boundaryFace[bfIdx].shapeValue[k];
                            rhoV[comp] += varNData[k].density*sol[k][comp]*this->fvGeom.boundaryFace[bfIdx].shapeValue[k];

                            FieldVector<Scalar,dim> gradVComp = grad;
                            gradVComp *= sol[k][comp];
                            velocityGradient[comp] += gradVComp;
                        }
                    }

                    Scalar xValue;
                    Scalar outward = velocityValue*this->fvGeom.boundaryFace[bfIdx].normal;
                    if (outward > 0)
                    	xValue = sol[node][dim];
                    else
                    	xValue = sol[node][dim];

                    for (int comp = 0; comp < dim; comp++)
                    	xV[comp] = rhoV[comp]*xValue;

                    FieldVector<Scalar,dim> DgradX(0);
                    D.mv(gradRhoX, DgradX);

//momentum balance
                    if (bctypeface[0] == BoundaryConditions::dirichlet)
                    {
/*
	                    FieldVector<Scalar,dim> gradVN(0);
        	            velocityGradient.umv(it->unitOuterNormal(faceLocal), gradVN);
                    	gradVN *= this->fvGeom.boundaryFace[bfIdx].area;
                    	//mutiply with mu !!!!!!!!!!!!! TODO
                    	for (int comp = 0; comp < dim; comp++)
                    	{
						//	FieldVector<Scalar,dim> gradVComp(0);
						//	for (int k = 0; k < this->fvGeom.numVertices; k++)
						//	{
						//		FieldVector<Scalar,dim> grad(this->fvGeom.boundaryFace[bfIdx].grad[k]);
						//		grad *= sol[k][comp];
						//		grad *= varNData[k].viscosity;
						//	}
						//
							FieldVector<Scalar,dim> pComp(0);
							pComp[comp] = -pressureValue;
							result[comp] = pComp*this->fvGeom.boundaryFace[bfIdx].normal;

                        	if (alpha == 0)
								result[comp] += gradVN[comp];// *this->fvGeom.boundaryFace[bfIdx].normal;
                    	}
   */
                    	for (int comp = 0; comp < dim; comp++)
						{
							FieldVector<Scalar,dim> gradVComp(0);
                    for (int k = 0; k < this->fvGeom.numVertices; k++)
                    {
                        FieldVector<Scalar,dim> grad(this->fvGeom.boundaryFace[bfIdx].grad[k]);
								grad *= sol[k][comp];
								grad *= varNData[k].viscosity;
								gradVComp = grad;
                    }

							FieldVector<Scalar,dim> pComp(0);
							pComp[comp] = -pressureValue;
							result[comp] = pComp*this->fvGeom.boundaryFace[bfIdx].normal;

							if (alpha == 0)
								result[comp] += gradVComp*this->fvGeom.boundaryFace[bfIdx].normal;
						}
					}
                    else
                    {
                        Scalar beaversJosephC = this->getImp().problem.beaversJosephC(this->fvGeom.boundaryFace[bfIdx].ipGlobal, element, it, this->fvGeom.boundaryFace[bfIdx].ipLocal);
                        if (beaversJosephC > 0) // realize Beavers-Joseph interface condition
                        {
                            FieldVector<Scalar,dim> tangentialV = velocityValue;

                            // normal n
                            FieldVector<Scalar,dim> normal = it->unitOuterNormal(faceLocal);

                            // v.n
                            Scalar normalComp = tangentialV*normal;

                            // (v.n)n
                            FieldVector<Scalar,dim> normalV = normal;
                            normalV *= normalComp;

                            // v_tau = v - (v.n)n
                            tangentialV -= normalV;

                            for (int comp = 0; comp < dim; comp++)
                            {
                            	result[comp] += 1.0/beaversJosephC*this->fvGeom.boundaryFace[bfIdx].area*tangentialV[comp];
                            }
                        }
                        else if (beaversJosephC < 0)
                        {
                            for (int comp = 0; comp < dim; comp++)
                            {
                                FieldVector<Scalar,dim>  pressVector(0);
                    			pressVector[comp] = -pressureValue;
                                result[comp] += pressVector*it->unitOuterNormal(faceLocal)*this->fvGeom.boundaryFace[bfIdx].area;
                            }
                        }
                    }
//transport
                    result[dim] = (xV - DgradX)*this->fvGeom.boundaryFace[bfIdx].normal;

//mass balance
                    FieldVector<Scalar,dim> massResidual = velocityValue;
					Scalar alphaH2 = alpha*this->fvGeom.subContVol[node].volume;
					SolutionVector source = problem.q(this->fvGeom.boundaryFace[bfIdx].ipGlobal, element, this->fvGeom.boundaryFace[bfIdx].ipLocal);
					FieldVector<Scalar,dim> dimSource;
					for (int comp = 0; comp < dim; comp++)
						dimSource[comp] = source[comp];
					dimSource *= alphaH2;
					massResidual += dimSource;
					result[dim+1] -= massResidual*it->unitOuterNormal(faceLocal)*this->fvGeom.boundaryFace[bfIdx].area;
                }
            }
        }
        return result;
    }


    void assembleBoundaryCondition(const Element& element)
    {
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
            //this->dirichletIndex[i] = 0;
        }

        // evaluate boundary conditions via intersection iterator
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
                            //this->getImp().problem.dirichletIndex(global, element, it, local, dirichletIdx); // eval bctype
                            if (bctypeface[equationNumber]!=BoundaryConditions::neumann)
                                break;
                            FieldVector<Scalar,numEq> J = this->getImp().problem.neumann(global, element, it, local);
//                            if (equationNumber < dim+1) {
                                J[equationNumber] *= this->fvGeom.boundaryFace[bfIdx].area;
                                this->b[nodeInElement][equationNumber] += J[equationNumber];
//                            }
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


            for (int equationNumber=0; equationNumber<dim+1; equationNumber++) {
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
                                //this->dirichletIndex[i][equationNumber] = dirichletIdx[equationNumber];

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
                                //this->dirichletIndex[i][equationNumber] = dirichletIdx[equationNumber];
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
    StokesTransportProblem<Grid,Scalar>& problem;
    double alpha;

    std::vector<VariableNodeData> varNData;
    std::vector<VariableNodeData> oldVarNData;

};
}
#endif
