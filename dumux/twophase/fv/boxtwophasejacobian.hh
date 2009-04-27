// $Id$

#ifndef DUNE_BOXTWOPHASEJACOBIAN_HH
#define DUNE_BOXTWOPHASEJACOBIAN_HH

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
#include "dumux/operators/localjacobian.hh"
#include "dumux/twophase/twophaseproblem_deprecated.hh"

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
  - Scalar    type used for return values
*/
template<class Grid, class Scalar, class BoxFunction = LeafP1Function<Grid, Scalar, 2> >
class BoxTwoPhaseLocalJacobian
    : public LocalJacobian<BoxTwoPhaseLocalJacobian<Grid,Scalar>,Grid,Scalar,2>
{
    typedef typename Grid::Traits::template Codim<0>::Entity Entity;
    typedef typename Entity::Geometry Geometry;
    typedef typename LocalJacobian<BoxTwoPhaseLocalJacobian<Grid,Scalar>,Grid,Scalar,2>::VBlockType SolutionVector;
    typedef typename LocalJacobian<BoxTwoPhaseLocalJacobian<Grid,Scalar>,Grid,Scalar,2>::MBlockType MBlockType;

public:
    // define the number of components of your system, this is used outside
    // to allocate the correct size of (dense) blocks with a FieldMatrix
    enum {dim=Grid::dimension};
    enum {numEq=2};
    enum {SIZE=LagrangeShapeFunctionSetContainer<Scalar,Scalar,dim>::maxsize};

    //! Constructor
    BoxTwoPhaseLocalJacobian (TwoPhaseProblem<Grid,Scalar>& params,
                              bool levelBoundaryAsDirichlet_, const Grid& grid,
                              BoxFunction& sol,
                              bool procBoundaryAsDirichlet_=true)
        : problem(params),levelBoundaryAsDirichlet(levelBoundaryAsDirichlet_),
          procBoundaryAsDirichlet(procBoundaryAsDirichlet_),
          currentSolution(sol), oldSolution(grid), dt(1)
    {}


    template<class TypeTag>
    void localDefect (const Entity& element, const SolutionVector* sol)
    {
        // extract some important parameters
        const Geometry& geometry = element.geometry();
        Dune::GeometryType gt = geometry.type();
        const typename Dune::LagrangeShapeFunctionSetContainer<Scalar,Scalar,dim>::value_type&
            sfs=Dune::LagrangeShapeFunctions<Scalar,Scalar,dim>::general(gt, 1);
        int size = sfs.size();
        this->setcurrentsize(size);

        // calculate secondary variables and set defect to zero
        Scalar saturationW[size];
        Scalar pC[size];
        SolutionVector mobility[size];
        Scalar mobN[size];
        SolutionVector density;
        for (int i=0; i < size; i++) {
            this->def[i] = 0;
            pC[i] = sol[i][1] - sol[i][0];
            saturationW[i] = problem.materialLaw().saturationW(pC[i]);
            mobility[i][0] = problem.materialLaw().mobW(saturationW[i]);
            mobility[i][1] = problem.materialLaw().mobN(1 - saturationW[i]);
        }
        density[0] = problem.materialLaw().wettingPhase.density();
        density[1] = problem.materialLaw().nonwettingPhase.density();

        // get cell volume
        Scalar cellVolume = geometry.volume();

        // cell center in reference element
        const Dune::FieldVector<Scalar,dim>
            elementLocal = Dune::ReferenceElements<Scalar,dim>::general(gt).position(0,0);

        // get global coordinate of cell center
        const Dune::FieldVector<Scalar,dim> elementGlobal = geometry.global(elementLocal);

        // ASSUMING element-wise constant permeability, evaluate K at the cell center
        const Dune::FieldMatrix<Scalar,dim,dim> K = problem.K(elementGlobal, element, elementLocal);

        // ASSUMING element-wise constant porosity, evaluate at the cell center
        double volumeFactor = problem.porosity(elementGlobal, element, elementLocal)*cellVolume/dt;
        for (int i=0; i < size; i++) // begin loop over vertices
        {
            // capillary pressure
            Scalar pCI = pC[i];
            Scalar pCOld = uold[i][1] - uold[i][0];

            // derivative of wetting phase saturation w.r.t. pC
            Scalar dSdP = problem.materialLaw().dSdP(pCI);

            // time derivative
            Scalar diffPC = dSdP*volumeFactor*(pCI - pCOld);
            this->def[i][0] += diffPC;
            this->def[i][1] -= diffPC;

            // local coordinate of vertex
            const FieldVector<Scalar,dim> vertexLocal = sfs[i].position();

            // get global coordinate of vertex
            const FieldVector<Scalar,dim> vertexGlobal = geometry.global(vertexLocal);

            // get source term
            FieldVector<Scalar, numEq> q = problem.q(vertexGlobal, element, vertexLocal);

            // add source to defect
            q *= cellVolume;
            this->def[i] -= q;

            for (int j=i+1; j < size; j++) // begin loop over neighboring vertices
            {
                // local coordinate of neighbor
                const FieldVector<Scalar,dim> neighborLocal = sfs[j].position();

                if (!gt.isSimplex()) {
                    // compute the local distance
                    Scalar distanceLocal = (vertexLocal - neighborLocal).two_norm();

                    // check whether the two vertices share a cell edge
                    if (distanceLocal > 1.01)
                        continue;
                }

                // get global coordinate of neighbor
                const FieldVector<Scalar,dim> neighborGlobal = geometry.global(neighborLocal);

                // compute the edge vector
                FieldVector<Scalar,dim>  edgeVector = neighborGlobal - vertexGlobal;

                // get distance between neighbors
                Scalar oneByDistanceGlobal = 1.0/edgeVector.two_norm();

                // normalize edge vector
                edgeVector *= oneByDistanceGlobal;

                // permeability in edge direction
                FieldVector<Scalar,dim> Kij(0);
                K.umv(edgeVector, Kij);

                // calculate pressure difference
                SolutionVector uDiff = sol[j] - sol[i];

                SolutionVector flux;
                for (int comp = 0; comp < numEq; comp++) {
                    // calculate pressure component gradient
                    FieldVector<Scalar, dim> uGrad(edgeVector);
                    uGrad *= oneByDistanceGlobal*uDiff[comp];

                    // adjust by gravity
                    FieldVector<Scalar, dim> gravity = problem.gravity();
                    gravity *= density[comp];
                    uGrad -= gravity;

                    // calculate the flux using upwind
                    Scalar outward = uGrad*Kij;
                    if (outward > 0)
                        flux[comp] = mobility[i][comp]*outward;
                    else
                        flux[comp] = mobility[j][comp]*outward;
                }

                // get the local edge center
                FieldVector<Scalar,dim> edgeLocal = vertexLocal + neighborLocal;
                edgeLocal *= 0.5;

                // get global coordinate of edge center
                const FieldVector<Scalar,dim> edgeGlobal = geometry.global(edgeLocal);

                // distance between cell center and edge center
                Scalar distanceEdgeCell = (elementGlobal - edgeGlobal).two_norm();

                ////////////////////////////////////////////////////////////
                // CAREFUL: only valid in 2D
                ////////////////////////////////////////////////////////////
                // obtain integrated Flux
                flux *= distanceEdgeCell;

                // add to defect
                this->def[i] -= flux;
                this->def[j] += flux;
            } // end loop over neighboring vertices
        } // end loop over vertices

        // adjust by density
        for (int i=0; i < size; i++) {
            this->def[i][0] *= density[0];
            this->def[i][1] *= density[1];
        }

        // assemble boundary conditions
        assembleBC<TypeTag> (element);

        // add to defect
        for (int i=0; i < size; i++)
            this->def[i] -= this->b[i];

        return;
    }

    void setLocalSolution (const Entity& element)
    {
        Dune::GeometryType gt = element.geometry().type();
        const typename Dune::LagrangeShapeFunctionSetContainer<Scalar,Scalar,dim>::value_type&
            sfs=Dune::LagrangeShapeFunctions<Scalar,Scalar,dim>::general(gt, 1);
        int size = sfs.size();
        this->setcurrentsize(size);

        for (int i = 0; i < size; i++)
            for (int comp = 0; comp < numEq; comp++) {
                this->u[i][comp] = currentSolution.evallocal(comp, element, sfs[i].position());
                uold[i][comp] = oldSolution.evallocal(comp, element, sfs[i].position());
            }

        return;
    }

    void setDt (double d)
    {
        dt = d;
    }

    void setOldSolution (BoxFunction& uOld)
    {
        *oldSolution = *uOld;
    }

private:
    template<class TypeTag>
    void assembleBC (const Entity& element)
    {
        // extract some important parameters
        Dune::GeometryType gt = element.geometry().type();
        const typename Dune::LagrangeShapeFunctionSetContainer<Scalar,Scalar,dim>::value_type&
            sfs=Dune::LagrangeShapeFunctions<Scalar,Scalar,dim>::general(gt,1);
        setcurrentsize(sfs.size());

        for (int i = 0; i < sfs.size(); i++)
            this->b[i] = 0;

        // determine quadrature order
        int p=2;
        // evaluate boundary conditions via intersection iterator
        typedef typename Grid::LeafGridView::IntersectionIterator
            IntersectionIterator;

        IntersectionIterator endit = element.ileafend();
        for (IntersectionIterator it = element.ileafbegin();
             it!=endit; ++it)
        {
            // if we have a neighbor then we assume there is no boundary (forget interior boundaries)
            // in level assemble treat non-level neighbors as boundary
            if (it.neighbor())
            {
                if (levelBoundaryAsDirichlet && it.outside()->level()==element.level())
                    continue;
                if (!levelBoundaryAsDirichlet)
                    continue;
            }

            // determine boundary condition type for this face, initialize with processor boundary
            FieldVector<typename BoundaryConditions::Flags, numEq> bctypeface(BoundaryConditions::process);

            // handle face on exterior boundary, this assumes there are no interior boundaries
            if (it.boundary())
            {
                Dune::GeometryType gtface = it.geometryInInside().type();
                for (size_t g=0; g<Dune::QuadratureRules<Scalar,dim-1>::rule(gtface,p).size(); ++g)
                {
                    const Dune::FieldVector<Scalar,dim-1>& facelocal = Dune::QuadratureRules<Scalar,dim-1>::rule(gtface,p)[g].position();
                    FieldVector<Scalar,dim> local = it.geometryInInside().global(facelocal);
                    FieldVector<Scalar,dim> global = it.geometry().global(facelocal);
                    bctypeface = problem.bctype(global,element,it,local); // eval bctype


                    if (bctypeface[0]!=BoundaryConditions::neumann) break;

                    SolutionVector J = problem.neumann(global,element,it,local);
                    if (J.two_norm() < 1e-10)
                        continue;
                    double weightface = Dune::QuadratureRules<Scalar,dim-1>::rule(gtface,p)[g].weight();
                    Scalar detjacface = it.geometry().integrationElement(facelocal);
                    for (int i=0; i<sfs.size(); i++) // loop over test function number
                        if (this->bctype[i][0]==BoundaryConditions::neumann)
                        {
                            //////////////////////////////////////////////////////////////////////////
                            // HACK: piecewise constants with respect to dual grid not implemented yet
                            // works only if exactly one quadrature point is located within each dual
                            // cell boundary (which should be the case for p = 2)
                            //////////////////////////////////////////////////////////////////////////
                            if (sfs[i].evaluateFunction(0,local) > 0.5) {
                                J *= weightface*detjacface;
                                this->b[i] -= J;
                            }
                        }
                }
                if (bctypeface[0]==BoundaryConditions::neumann) continue; // was a neumann face, go to next face
            }

            // If we are here, then it is
            // (i)   an exterior boundary face with Dirichlet condition, or
            // (ii)  a processor boundary (i.e. neither boundary() nor neighbor() was true), or
            // (iii) a level boundary in case of level-wise assemble
            // How processor boundaries are handled depends on the processor boundary mode
            if (bctypeface[0]==BoundaryConditions::process && procBoundaryAsDirichlet==false
                && levelBoundaryAsDirichlet==false)
                continue; // then it acts like homogeneous Neumann

            // now handle exterior or interior Dirichlet boundary
            for (int i=0; i<sfs.size(); i++) // loop over test function number
            {
                if (sfs[i].codim()==0) continue; // skip interior dof
                if (sfs[i].codim()==1) // handle face dofs
                {
                    if (sfs[i].entity()==it.indexInInside())
                    {
                        if (this->bctype[i][0]<bctypeface[0])
                        {
                            this->bctype[i].assign(bctypeface[0]);
                            if (bctypeface[0]==BoundaryConditions::process)
                                this->b[i] = 0;
                            if (bctypeface[0]==BoundaryConditions::dirichlet)
                            {
                                this->b[i] = 0;
                            }
                        }
                    }
                    continue;
                }
                // handle subentities of this face
                for (int j=0; j<ReferenceElements<Scalar,dim>::general(gt).size(it.indexInInside(),1,sfs[i].codim()); j++)
                    if (sfs[i].entity()==ReferenceElements<Scalar,dim>::general(gt).subEntity(it.indexInInside(),1,j,sfs[i].codim()))
                    {
                        if (this->bctype[i][0]<bctypeface[0])
                        {
                            this->bctype[i].assign(bctypeface[0]);
                            if (bctypeface[0]==BoundaryConditions::process)
                                this->b[i] = 0;
                            if (bctypeface[0]==BoundaryConditions::dirichlet)
                            {
                                this->b[i] = 0;
                            }
                        }
                    }
            }
        }
    }

    // parameters given in constructor
    TwoPhaseProblem<Grid,Scalar>& problem;
    bool levelBoundaryAsDirichlet;
    bool procBoundaryAsDirichlet;
    const BoxFunction& currentSolution;
    BoxFunction oldSolution;
    double dt;
    SolutionVector uold[SIZE];
};

/** @} */
}
#endif
