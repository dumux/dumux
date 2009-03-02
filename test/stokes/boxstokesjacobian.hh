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
#include <dune/grid/utility/intersectiongetter.hh>

#include<dune/disc/shapefunctions/lagrangeshapefunctions.hh>
#include<dune/disc/operators/boundaryconditions.hh>

#include<dumux/operators/boxjacobian.hh>
#include "dumux/stokes/stokesproblem.hh"

namespace Dune
{
//! A class for computing local jacobian matrices
/*! A class for computing local jacobian matrix for the
  diffusion equation

  div j = q; j = -K grad u; in Omega

  u = g on Gamma1; j*n = J on Gamma2.

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
        alpha = -1.0e-3;
        this->analytic = false;
    }

    void clearVisited ()
    {
        return;
    }

    template<class TypeTag>
    void localDefect(const Element& element, const SolutionVector* sol, bool withBC = true) {
        BoxJacobianType::template localDefect<TypeTag>(element, sol, withBC);

        this->template assembleBC<TypeTag>(element);

        Dune::GeometryType gt = element.geometry().type();
        const typename ReferenceElementContainer<Scalar,dim>::value_type& referenceElement = ReferenceElements<Scalar, dim>::general(gt);

        for (int vert=0; vert < this->fvGeom.numVertices; vert++) // begin loop over vertices / sub control volumes
            if (!this->fvGeom.subContVol[vert].inner)
                {
                    typedef typename IntersectionIteratorGetter<Grid,TypeTag>::IntersectionIterator IntersectionIterator;

                    FieldVector<Scalar,dim> averagedNormal(0);
                    int faces = 0;
                    IntersectionIterator endit = IntersectionIteratorGetter<Grid, TypeTag>::end(element);
                    for (IntersectionIterator it = IntersectionIteratorGetter<Grid, TypeTag>::begin(element); it!=endit; ++it)
                        {
                            if (it->boundary()) {
                                // get geometry type of face
                                GeometryType faceGT = it->intersectionSelfLocal().type();

                                // center in face's reference element
                                const FieldVector<Scalar,dim-1>& faceLocal = ReferenceElements<Scalar,dim-1>::general(faceGT).position(0,0);

                                int faceIdx = it->numberInSelf();
                                int numVerticesOfFace = referenceElement.size(faceIdx, 1, dim);
                                for (int nodeInFace = 0; nodeInFace < numVerticesOfFace; nodeInFace++) {
                                    int nodeInElement = referenceElement.subEntity(faceIdx, 1, nodeInFace, dim);
                                    if (nodeInElement != vert)
                                        continue;
                                    int bfIdx = this->fvGeom.boundaryFaceIndex(faceIdx,    nodeInFace);

                                    FieldVector<Scalar,dim> local = this->fvGeom.boundaryFace[bfIdx].ipLocal;
                                    FieldVector<Scalar,dim> global = this->fvGeom.boundaryFace[bfIdx].ipGlobal;
                                    BoundaryConditions::Flags bctypeface = this->getImp().problem.bctype(global, element, it, local);

                                    if (bctypeface == BoundaryConditions::dirichlet)
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

                    //std::cout << this->fvGeom.subContVol[vert].global << ": N = " << averagedNormal << ", cond = " << this->bctype[vert][0] << std::endl;

                    if (faces == 2 && this->fvGeom.numVertices == 4)
                        {
                            //this->def[vert][1] = sol[0][1] + sol[3][1] - sol[1][1] - sol[2][1];
                            this->def[vert][dim] = sol[0][dim] + sol[3][dim] - sol[1][dim] - sol[2][dim];
                        }
                    else if (this->bctype[vert][0] == BoundaryConditions::dirichlet)
                        this->def[vert][dim] = defect;
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
                                        grad *= sol[k][dim];
                                        pressGradient += grad;
                                    }

                                    Scalar alphaH2 = 0.5*alpha*(this->fvGeom.subContVol[i].volume + this->fvGeom.subContVol[j].volume);
                                    pressGradient *= alphaH2;
                                    if (i == vert)
                                        this->def[vert][dim] += pressGradient*this->fvGeom.subContVolFace[face].normal;
                                    else
                                        this->def[vert][dim] -= pressGradient*this->fvGeom.subContVolFace[face].normal;
                                }
                        }
                }

        return;
    }

    SolutionVector computeM (const Element& element, const SolutionVector* sol, int node, bool old = false)
    {
        SolutionVector result(0);

        return result;
    }

    SolutionVector computeQ (const Element& element, const SolutionVector* sol, const int& node)
    {
        SolutionVector result = problem.q(this->fvGeom.subContVol[node].global, element, this->fvGeom.subContVol[node].local);
        result *= -1.0;

        Scalar alphaH2 = alpha*this->fvGeom.subContVol[node].volume;
        result[dim] *= alphaH2;

        SolutionVector flux = boundaryFlux(element, sol, node);

        flux /= this->fvGeom.subContVol[node].volume;
        result -= flux;

        return (result);
    }

    SolutionVector computeA (const Element& element, const SolutionVector* sol, int face)
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
        typedef typename IntersectionIteratorGetter<Grid,LeafTag>::IntersectionIterator IntersectionIterator;

        IntersectionIterator endit = IntersectionIteratorGetter<Grid, LeafTag>::end(element);
        for (IntersectionIterator it = IntersectionIteratorGetter<Grid, LeafTag>::begin(element); it!=endit; ++it)
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
                    int faceIdx = it->numberInSelf();

                    int numVerticesOfFace = referenceElement.size(faceIdx, 1, dim);
                    for (int nodeInFace = 0; nodeInFace < numVerticesOfFace; nodeInFace++) {
                        int nodeInElement = referenceElement.subEntity(faceIdx, 1, nodeInFace, dim);
                        if (node != nodeInElement)
                            continue;

                        int bfIdx = this->fvGeom.boundaryFaceIndex(faceIdx,    nodeInFace);

                        FieldVector<typename BoundaryConditions::Flags, numEq> bctypeface(BoundaryConditions::process);
                        bctypeface = this->getImp().problem.bctype(this->fvGeom.boundaryFace[bfIdx].ipGlobal, element, it, this->fvGeom.boundaryFace[bfIdx].ipLocal);

                        // get geometry type of face
                        GeometryType faceGT = it->intersectionSelfLocal().type();

                        // center in face's reference element
                        const FieldVector<Scalar,dim-1>& faceLocal = ReferenceElements<Scalar,dim-1>::general(faceGT).position(0,0);

                        Scalar pressureValue = 0;
                        FieldVector<Scalar,dim>  velocityValue(0);
                        FieldVector<Scalar,dim>  pressGradient(0);
                        FieldMatrix<Scalar,dim,dim> velocityGradient(0);
                        for (int vert = 0; vert < this->fvGeom.numVertices; vert++) {
                            pressureValue += sol[vert][dim]*this->fvGeom.boundaryFace[bfIdx].shapeValue[vert];
                            FieldVector<Scalar,dim> grad(this->fvGeom.boundaryFace[bfIdx].grad[vert]);
                            grad *= sol[vert][dim];
                            pressGradient += grad;

                            grad = this->fvGeom.boundaryFace[bfIdx].grad[vert];

                            for (int comp = 0; comp < dim; comp++)
                                {
                                    velocityValue[comp] += sol[vert][comp]*this->fvGeom.boundaryFace[bfIdx].shapeValue[vert];
                                    FieldVector<Scalar,dim> gradVComp = grad;
                                    gradVComp *= sol[vert][comp];
                                    velocityGradient[comp] += gradVComp;
                                }
                        }
                        FieldVector<Scalar, dim> massResidual = velocityValue;
                        Scalar alphaH2 = alpha*this->fvGeom.subContVol[node].volume;
                        SolutionVector source = problem.q(this->fvGeom.boundaryFace[bfIdx].ipGlobal, element, this->fvGeom.boundaryFace[bfIdx].ipLocal);
                        FieldVector<Scalar, dim> dimSource;
                        for (int comp = 0; comp < dim; comp++)
                            dimSource[comp] = source[comp];
                        dimSource *= alphaH2;
                        massResidual += dimSource;
                        result[dim] -= massResidual*it->unitOuterNormal(faceLocal)*this->fvGeom.boundaryFace[bfIdx].area;

                        if (bctypeface[0] == BoundaryConditions::dirichlet)
                            {
                                FieldVector<Scalar,dim>  gradVN(0);
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
                            }
                        else
                            {
                                Scalar beaversJosephC = this->getImp().problem.beaversJosephC(this->fvGeom.boundaryFace[bfIdx].ipGlobal, element, it, this->fvGeom.boundaryFace[bfIdx].ipLocal);
                                if (fabs(beaversJosephC) > 0) // realize Beavers-Joseph interface condition
                                    {
                                        FieldVector<Scalar,dim> tangentialV = velocityValue;
                                        //                            for (int comp = 0; comp < dim; comp++)
                                        //                            {
                                        //                                tangentialV[comp] = sol[node][comp];
                                        //                            }

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
                            }
                    }
                }
            }

        return result;
    }


    template<class TypeTag> void assembleBC(const Element& element) {
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
            this->dirichletIndex[i] = 0;
        }

        // evaluate boundary conditions via intersection iterator
        typedef typename IntersectionIteratorGetter<Grid,TypeTag>::IntersectionIterator IntersectionIterator;

        IntersectionIterator endit = IntersectionIteratorGetter<Grid, TypeTag>::end(element);
        for (IntersectionIterator it = IntersectionIteratorGetter<Grid, TypeTag>::begin(element); it!=endit; ++it)
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
                    int faceIdx = it->numberInSelf();
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
                                //                            std::cout << "faceIdx = " << faceIdx << ", nodeInElement = " << nodeInElement
                                //                                      << ", bfIdx = " << bfIdx << ", local = " << local << ", global = " << global
                                //                                      << ", bctypeface = " << bctypeface << std::endl;
                                if (bctypeface[equationNumber]!=BoundaryConditions::neumann)
                                    break;
                                FieldVector<Scalar,dim> J = this->getImp().problem.J(global, element, it, local);
                                //                            std::cout << "J = " << J << std::endl;
                                if (equationNumber < dim) {
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
                                    if (sfs[i].entity()==it->numberInSelf()) {
                                        if (this->bctype[i][equationNumber] < bctypeface[equationNumber]) {
                                            this->bctype[i][equationNumber] = bctypeface[equationNumber];
                                            this->dirichletIndex[i][equationNumber] = dirichletIdx[equationNumber];

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
                            for (int j=0; j<ReferenceElements<Scalar,dim>::general(gt).size(it->numberInSelf(), 1, sfs[i].codim()); j++)
                                if (sfs[i].entity()==ReferenceElements<Scalar,dim>::general(gt).subEntity(it->numberInSelf(), 1, j, sfs[i].codim()))
                                    {
                                        if (this->bctype[i][equationNumber] < bctypeface[equationNumber]) {
                                            this->bctype[i][equationNumber] = bctypeface[equationNumber];
                                            this->dirichletIndex[i][equationNumber] = dirichletIdx[equationNumber];
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
};
}
#endif
