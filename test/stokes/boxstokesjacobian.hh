#ifndef DUNE_BOXSTOKESJACOBIAN_HH
#define DUNE_BOXSTOKESJACOBIAN_HH

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
#include "dumux/stokes/stokesproblem.hh"

namespace Dune
{
//! A class for computing local jacobian matrices
/*! A class for computing local jacobian matrix for the
   Stokes equation

  div (-mu grad u + pI) = f;  div u = 0; in Omega

  u = g_D on Gamma_D; (-mu grad u + pI) n  = g_N on Gamma_N.

  Uses the box method.

  Template parameters are:

  - Grid     a DUNE grid type
  - Scalar    type used for return values
*/
template<class Grid, class Scalar, class BoxFunction = LeafP1Function<Grid, Scalar, Grid::dimension+1> >
class BoxStokesJacobian
    : public BoxJacobian<BoxStokesJacobian<Grid,Scalar,BoxFunction>,Grid,Scalar,Grid::dimension+1,BoxFunction>
{
    enum {dim=Grid::dimension};
    enum {numEq = dim+1};

    typedef typename Grid::Traits::template Codim<0>::Entity Element;
    typedef typename Grid::LeafGridView::IntersectionIterator IntersectionIterator;
    typedef typename Element::Geometry Geometry;
    typedef BoxStokesJacobian<Grid,Scalar,BoxFunction> ThisType;
    typedef typename LocalJacobian<ThisType,Grid,Scalar,numEq>::VBlockType SolutionVector;
    typedef BoxJacobian<ThisType,Grid,Scalar,numEq,BoxFunction> BoxJacobianType;
    typedef Dune::FVElementGeometry<Grid> FVElementGeometry;

public:

    //! Constructor
    BoxStokesJacobian (StokesProblem<Grid,Scalar>& params,
                       bool levelBoundaryAsDirichlet_, const Grid& grid,
                       BoxFunction& sol,
                       bool procBoundaryAsDirichlet_=true)
        : BoxJacobianType(levelBoundaryAsDirichlet_, grid, sol, procBoundaryAsDirichlet_),
          problem(params)
    {
        // (negative) parameter for stabilization of pressure
        alpha = -1;
        //parameter to take into account natural boundary conditions on Gamma_N
        beta = 1;
        this->analytic = false;
    }

    void clearVisited ()
    {
        return;
    }

    // calculate defect on the boundary for mass balance equation
    void localDefect(const Element& element, const SolutionVector* sol, bool withBC = true) {
        BoxJacobianType::localDefect(element, sol, withBC);

        Dune::GeometryType gt = element.geometry().type();
        const typename ReferenceElementContainer<Scalar,dim>::value_type& referenceElement = ReferenceElements<Scalar, dim>::general(gt);

        for (int vert=0; vert < this->fvGeom.numVertices; vert++) // begin loop over vertices / sub control volumes
            if (!this->fvGeom.subContVol[vert].inner)
                {
                    FieldVector<Scalar,dim> averagedNormal(0);
                    FieldMatrix<Scalar,dim,dim> velocityGradient(0);
                    int faces = 0;
                    IntersectionIterator endit = element.ileafend();
                    for (IntersectionIterator it = element.ileafbegin(); it!=endit; ++it)
                        {
                            if (it->boundary()) {
                                // get geometry type of face
                                GeometryType faceGT = it->geometryInInside().type();

                                // center in face's reference element
                                const FieldVector<Scalar,dim-1>& faceLocal = ReferenceElements<Scalar,dim-1>::general(faceGT).position(0,0);

                                // return the face number with respect to the generic reference element
                                int faceIdx = it->indexInInside();


                                int numVerticesOfFace = referenceElement.size(faceIdx, 1, dim);
                                for (int nodeInFace = 0; nodeInFace < numVerticesOfFace; nodeInFace++) {
                                    int nodeInElement = referenceElement.subEntity(faceIdx, 1, nodeInFace, dim);
                                    if (nodeInElement != vert)
                                        continue;
                                    int bfIdx = this->fvGeom.boundaryFaceIndex(faceIdx,    nodeInFace);

                                    FieldVector<Scalar,dim> local = this->fvGeom.boundaryFace[bfIdx].ipLocal;
                                    FieldVector<Scalar,dim> global = this->fvGeom.boundaryFace[bfIdx].ipGlobal;
                                    FieldVector<BoundaryConditions::Flags,numEq> bctypeface = this->getImp().problem.bctype(global, element, it, local);

                                    // Checks if boundary conditions for velocity are DIRICHLET
                                    bool onlyDirichlet = true;
                                    FieldVector<BoundaryConditions::Flags, numEq>  bctypeNode = this->getImp().problem.bctype(global, element, it, local);
                                    for (int i=0; i<dim; ++i)
                                    {
                                    	if (bctypeNode[i] != BoundaryConditions::dirichlet)
                                    		onlyDirichlet = false;
                                    }

                                    if (onlyDirichlet)
                                    {
                                    	FieldVector<Scalar,dim> normal = it->unitOuterNormal(faceLocal);
                                        normal *= this->fvGeom.boundaryFace[bfIdx].area;
                                        averagedNormal += normal;

                                        for (int k = 0; k < this->fvGeom.numVertices; k++)
                                        {
                                        	FieldVector<Scalar,dim> grad(this->fvGeom.boundaryFace[bfIdx].grad[k]);

                                        	for (int comp = 0; comp < dim; comp++)
                                        	{
                                        		FieldVector<Scalar,dim> gradVComp = grad;
                                        		gradVComp *= sol[k][comp];
                                        		velocityGradient[comp] += gradVComp;
                                        	}
                                        }
                                   }
                                    faces++;
                                }
                            }
                        }
                    Scalar defect = 0;
                    Scalar normalLength = averagedNormal.two_norm();
                    if (normalLength)
                        averagedNormal /= normalLength;

                    FieldVector<Scalar,dim>  tangent(0);
                    tangent[0] = averagedNormal[1];
                    tangent[1] = -averagedNormal[0];

                	// on the Dirichlet boundary Gamma_D we use momentum equation
                    // multiplied by the averaged normal vector
                    for (int k = 0; k < dim; k++)
                        defect += this->def[vert][k]*averagedNormal[k];

                    // in the corners, reduce the degree of the polynomials
                    // approximating pressure to 1 (remove term xy)
                    if (faces == 2 && this->fvGeom.numVertices == 4)
                    {
                    	this->def[vert][dim] = sol[0][dim] + sol[3][dim] - sol[1][dim] - sol[2][dim];
                    }
                    else
                       	if (this->bctype[vert][0] == BoundaryConditions::dirichlet)
                        {
                       		// from boundaryFlux we have p.n on the boundary
                       		// here we compute mu*gradVT*tangent instead of -mu*gradVN*n
                       		FieldVector<Scalar,dim>  gradVT(0);
                            velocityGradient.umv(tangent, gradVT);
                            gradVT *= normalLength;
                            gradVT *= elData.mu;
							defect -= gradVT*tangent;

                            this->def[vert][dim] = defect;
                        }
                    else // de-stabilize (remove alpha*grad p - alpha div f from computeFlux on the boundary)
                    	 // use mass balance equation without any stabilization on the boundary
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
                                        grad *= sol[k][dim];
                                        pressGradient += grad;
                                    }

                                    Scalar alphaH2 = 0.5*alpha*(this->fvGeom.subContVol[i].volume + this->fvGeom.subContVol[j].volume);
                                    pressGradient *= alphaH2;

                                    if (vert == i)
                                        this->def[vert][dim] += pressGradient*this->fvGeom.subContVolFace[face].normal;
                                    else
                                        this->def[vert][dim] -= pressGradient*this->fvGeom.subContVolFace[face].normal;
                                }

                            SolutionVector source = problem.q(this->fvGeom.subContVol[vert].global, element, this->fvGeom.subContVol[vert].local);
                            Scalar alphaH2 = alpha*this->fvGeom.subContVol[vert].volume;
                            this->def[vert][dim] -= alphaH2*source[dim]*this->fvGeom.subContVol[vert].volume;
                        }
                }

        return;
    }

    SolutionVector computeStorage (const Element& element, const SolutionVector* sol, int node, bool old = false)
    {
        SolutionVector result(0);

        return result;
    }

    SolutionVector computeSource (const Element& element, const SolutionVector* sol, const int& node)
    {
    	// calculate source term and take into account flux calculated on the boundary
    	// source term for mass balance equation: problem.q[dim] = div f
    	// should be multiplied with stabilization parameter alpha
    	SolutionVector result = problem.q(this->fvGeom.subContVol[node].global, element, this->fvGeom.subContVol[node].local);
        result *= -1.0;

        Scalar alphaH2 = alpha*this->fvGeom.subContVol[node].volume;
        result[dim] *= alphaH2;

        // used for Gamma_N to take into account natural boundary conditions
        // g_N*n = mu*gradVN*n - p*n
        // 0.5*beta*p, for the rest see boundaryFlux
        // factor 0.5 is used in the equation and RHS, and can be removed
        if (!this->fvGeom.subContVol[node].inner)
            result[dim] += beta*0.5*sol[node][dim]/this->fvGeom.subContVol[node].volume;

        SolutionVector flux = boundaryFlux(element, sol, node);

        flux /= this->fvGeom.subContVol[node].volume;
        result -= flux;

        return (result);
    }

    SolutionVector computeFlux (const Element& element, const SolutionVector* sol, int face)
    {
        SolutionVector flux(0);

        Scalar pressValue = 0;
        FieldVector<Scalar, dim> pressGradient(0);
        FieldVector<Scalar, dim> velocityValue(0);
        for (int k = 0; k < this->fvGeom.numVertices; k++) {
            pressValue += sol[k][dim]*this->fvGeom.subContVolFace[face].shapeValue[k];
            FieldVector<Scalar,dim> grad(this->fvGeom.subContVolFace[face].grad[k]);
            grad *= sol[k][dim];
            pressGradient += grad;
            for (int comp = 0; comp < dim; comp++)
                velocityValue[comp] += sol[k][comp]*this->fvGeom.subContVolFace[face].shapeValue[k];
        }

        // momentum balance:
        // (-mu grad u + p I)n
        for (int comp = 0; comp < dim; comp++)
            {
                FieldVector<Scalar, dim> gradVComp(0);
                for (int k = 0; k < this->fvGeom.numVertices; k++) {
                    FieldVector<Scalar,dim> grad(this->fvGeom.subContVolFace[face].grad[k]);
                    grad *= sol[k][comp];
                    gradVComp += grad;
                }
                gradVComp *= -elData.mu;

                FieldVector<Scalar, dim> pComp(0);
                pComp[comp] = pressValue;
                gradVComp += pComp;

                flux[comp] = gradVComp*this->fvGeom.subContVolFace[face].normal;
            }

        // mass balance:
        // (u + alpha grad p)n
        FieldVector<Scalar, dim> massResidual = velocityValue;
        int i = this->fvGeom.subContVolFace[face].i;
        int j = this->fvGeom.subContVolFace[face].j;
        Scalar alphaH2 = 0.5*alpha*(this->fvGeom.subContVol[i].volume + this->fvGeom.subContVol[j].volume);
        pressGradient *= alphaH2;
        massResidual += pressGradient;
        flux[dim] = massResidual*this->fvGeom.subContVolFace[face].normal;

        return flux;
    }

    void computeElementData (const Element& element)
    {
        // ASSUMING element-wise constant viscosity, evaluate mu at the cell center
        elData.mu = problem.mu(this->fvGeom.elementGlobal, element, this->fvGeom.elementLocal);
    };

    virtual void updateVariableData(const Element& element, const SolutionVector* sol, int i, bool old = false)
    {
        return;
    }

    void updateVariableData(const Element& element, const SolutionVector* sol, bool old = false)
    {
        return;
    }

    virtual void updateStaticData (const Element& element, const SolutionVector* sol)
    {
        return;
    }

    SolutionVector boundaryFlux(const Element& element, const SolutionVector* sol, int node)
    {
    	// calculate flux on the boundary for momentum balance equations
    	// and mass balance equation
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
                        FieldVector<Scalar,dim>  pressGradient(0);
                        FieldMatrix<Scalar,dim,dim> velocityGradient(0);
                        for (int vert = 0; vert < this->fvGeom.numVertices; vert++)
                        {
                            pressureValue += sol[vert][dim]*this->fvGeom.boundaryFace[bfIdx].shapeValue[vert];
                            FieldVector<Scalar,dim> grad(this->fvGeom.boundaryFace[bfIdx].grad[vert]);

                            for (int comp = 0; comp < dim; comp++)
                                {
                                    velocityValue[comp] += sol[vert][comp]*this->fvGeom.boundaryFace[bfIdx].shapeValue[vert];
                                    FieldVector<Scalar,dim> gradVComp = grad;
                                    gradVComp *= sol[vert][comp];
                                    velocityGradient[comp] += gradVComp;
                                }
                        }
                        FieldVector<Scalar, dim> massResidual = velocityValue;
                        result[dim] -= massResidual*it->unitOuterNormal(faceLocal)*this->fvGeom.boundaryFace[bfIdx].area;

                        FieldVector<Scalar,dim>  tangent(0);
                        tangent[0] = it->unitOuterNormal(faceLocal)[1];
                        tangent[1] = -(it->unitOuterNormal(faceLocal)[0]);
                        FieldVector<Scalar,dim>  gradVT(0);
                        velocityGradient.umv(tangent, gradVT);
                        gradVT *= this->fvGeom.boundaryFace[bfIdx].area;
                        gradVT *= elData.mu;

                        // for Dirichlet boundary Gamma_D we have: -pn
                        // the rest is done in localDefect: -mu*gradVT*tangent = mu*gradVN.n
                        if (bctypeface[0] == BoundaryConditions::dirichlet)
                            {
                                for (int comp = 0; comp < dim; comp++)
                                    {
                                        FieldVector<Scalar,dim>  pressVector(0);
                                        pressVector[comp] = -pressureValue;
                                        result[comp] += pressVector*it->unitOuterNormal(faceLocal)*this->fvGeom.boundaryFace[bfIdx].area;
                                    }
                            }
                        else // for the natural boundary Gamma_N we have
                        	 // g_N*n = mu*gradVN*n - p*n (this equation is added to mass balance with factor beta)
                        	 // -u*n - 0.5*beta*gradVT*tangent - 0.5*beta*g_N*n, the rest (0.5*beta*p) is in computeSource
                            {
                                Scalar beaversJosephC = this->getImp().problem.beaversJosephC(this->fvGeom.boundaryFace[bfIdx].ipGlobal, element, it, this->fvGeom.boundaryFace[bfIdx].ipLocal);
                                if (beaversJosephC == 0)
                                    {
                                        result[dim] -= beta*0.5*(gradVT*tangent)/this->fvGeom.boundaryFace[bfIdx].area;

                                        FieldVector<Scalar,numEq> neumann = this->getImp().problem.neumann(this->fvGeom.boundaryFace[bfIdx].ipGlobal, element, it, this->fvGeom.boundaryFace[bfIdx].ipLocal);

                                        for (int comp = 0; comp < dim; comp++)
                                            result[dim] -= beta*0.5*neumann[comp]*it->unitOuterNormal(faceLocal)[comp];
                                    }
                                else if (beaversJosephC > 0) // realize Beavers-Joseph interface condition
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
                                                result[comp] -= 1.0/beaversJosephC*this->fvGeom.boundaryFace[bfIdx].area*tangentialV[comp];
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
                    }
                }
            }

        return result;
    }


    void assembleBoundaryCondition(const Element& element, int k = 1)
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

                                if (bctypeface[equationNumber]!=BoundaryConditions::neumann)
                                    break;
                                FieldVector<Scalar,numEq> J = this->getImp().problem.neumann(global, element, it, local);

                                if (equationNumber < dim)
                                    {
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
                // (ii)  a processor boundary (i.e. neither boundary() nor neighbor() was true), or
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


                for (int equationNumber=0; equationNumber<dim; equationNumber++) {
                    for (int i=0; i<sfs.size(); i++) // loop over test function number
                        {
                            //this->dirichletIndex[i][equationNumber] = equationNumber;

                            if (sfs[i].codim()==0)
                                continue; // skip interior dof
                            if (sfs[i].codim()==1) // handle face dofs
                                {
                                    if (sfs[i].entity()==it->indexInInside()) {
                                        if (this->bctype[i][equationNumber] < bctypeface[equationNumber]) {
                                            this->bctype[i][equationNumber] = bctypeface[equationNumber];
                                            // this->dirichletIndex[i][equationNumber] = dirichletIdx[equationNumber];

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
                                            // this->dirichletIndex[i][equationNumber] = dirichletIdx[equationNumber];
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


    struct ElementData {
        Scalar mu;
    };

    ElementData elData;
    StokesProblem<Grid,Scalar>& problem;
    double alpha;
    double beta;
};
}
#endif
