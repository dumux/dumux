// $Id$

#ifndef DUNE_BOXPARABOLICJACOBIAN_HH
#define DUNE_BOXPARABOLICJACOBIAN_HH

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
#include "dumux/operators/localjacobian.hh"
#include "dumux/parabolic/parabolicproblem.hh"

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
  - RT    type used for return values
*/
template<class G, class RT, class BoxFunction = LeafP1Function<G, RT> >
class BoxParabolicLocalJacobian
    : public LocalJacobian<BoxParabolicLocalJacobian<G,RT>,G,RT,1>
{
    typedef typename G::ctype DT;
    typedef typename G::Traits::template Codim<0>::Entity Entity;
    typedef typename Entity::Geometry Geometry;
    typedef typename LocalJacobian<BoxParabolicLocalJacobian<G,RT>,G,RT,1>::VBlockType VBlockType;
    typedef typename LocalJacobian<BoxParabolicLocalJacobian<G,RT>,G,RT,1>::MBlockType MBlockType;

public:
    // define the number of components of your system, this is used outside
    // to allocate the correct size of (dense) blocks with a FieldMatrix
    enum {n=G::dimension};
    enum {m=1};
    enum {SIZE=LagrangeShapeFunctionSetContainer<DT,RT,n>::maxsize};

    //! Constructor: do not change the first five arguments
    BoxParabolicLocalJacobian (ParabolicProblem<G,RT>& params,
                               bool levelBoundaryAsDirichlet_, const G& grid,
                               BoxFunction& sol,
                               bool procBoundaryAsDirichlet_=true)
        : problem(params),levelBoundaryAsDirichlet(levelBoundaryAsDirichlet_),
          procBoundaryAsDirichlet(procBoundaryAsDirichlet_),
          currentSolution(sol), oldSolution(grid), dt(1)
    {
        this->analytic = false;
    }


    template<class TypeTag>
    void localDefect (const Entity& e, const VBlockType* sol, bool withBC = true)
    {
        // extract some important parameters
        const Geometry& geometry = e.geometry();
        Dune::GeometryType gt = geometry.type();
        const typename Dune::LagrangeShapeFunctionSetContainer<DT,RT,n>::value_type&
            sfs=Dune::LagrangeShapeFunctions<DT,RT,n>::general(gt, 1);
        int size = sfs.size();
        this->setcurrentsize(size);

        // get cell volume
        DT cellVolume = geometry.volume();

        // cell center in reference element
        const Dune::FieldVector<DT,n>
            elementLocal = Dune::ReferenceElements<DT,n>::general(gt).position(0,0);

        // get global coordinate of cell center
        const Dune::FieldVector<DT,n> elementGlobal = geometry.global(elementLocal);

        // ASSUMING element-wise constant permeability, evaluate K at the cell center
        const Dune::FieldMatrix<DT,n,n> K = problem.K(elementGlobal, e, elementLocal);

        // set defect to zero
        for (int i=0; i < size; i++)
            this->def[i] = 0;

        double volumeFactor = cellVolume/dt;
        for (int i=0; i < size; i++) // begin loop over vertices
        {
            // time derivative
            this->def[i] += volumeFactor*(sol[i] - uold[i]);

            // local coordinate of vertex
            const FieldVector<DT,n> vertexLocal = sfs[i].position();

            // get global coordinate of vertex
            const FieldVector<DT,n> vertexGlobal = geometry.global(vertexLocal);

            // get source term
            RT q = problem.q(vertexGlobal, e, vertexLocal);

            // add source to defect
            this->def[i] -= cellVolume*q;

            for (int j=i+1; j < size; j++) // begin loop over neighboring vertices
            {
                // local coordinate of neighbor
                const FieldVector<DT,n> neighborLocal = sfs[j].position();

                if (!gt.isSimplex()) {
                    // compute the local distance
                    DT distanceLocal = (vertexLocal - neighborLocal).two_norm();

                    // check whether the two vertices share a cell edge
                    if (distanceLocal > 1.01)
                        continue;
                }

                // get global coordinate of neighbor
                const FieldVector<DT,n> neighborGlobal = geometry.global(neighborLocal);

                // compute the edge vector
                FieldVector<DT,n>  edgeVector = vertexGlobal - neighborGlobal;

                // get distance between neighbors
                DT oneByDistanceGlobal = 1.0/edgeVector.two_norm();

                // normalize edge vector
                edgeVector *= oneByDistanceGlobal;

                // calculate pressure difference
                DT uDiff = sol[j] - sol[i];

                // permeability in edge direction
                FieldVector<DT,n> Kij(0);
                K.umv(edgeVector, Kij);

                // calculate the flux
                RT flux = oneByDistanceGlobal*uDiff*(Kij*edgeVector);

                // get the local edge center
                FieldVector<DT,n> edgeLocal = vertexLocal + neighborLocal;
                edgeLocal *= 0.5;

                // get global coordinate of edge center
                const FieldVector<DT,n> edgeGlobal = geometry.global(edgeLocal);

                // distance between cell center and edge center
                DT distanceEdgeCell = (elementGlobal - edgeGlobal).two_norm();

                ////////////////////////////////////////////////////////////
                // CAREFUL: only valid in 2D
                ////////////////////////////////////////////////////////////
                // obtain integrated Flux
                RT integratedFlux = distanceEdgeCell*flux;

                // add to defect
                this->def[i] -= integratedFlux;
                this->def[j] += integratedFlux;
            } // end loop over neighboring vertices
        } // end loop over vertices

        // assemble boundary conditions
        if (withBC)
        {
            assembleBC<TypeTag> (e);

            // add to defect
            for (int i=0; i < size; i++)
                this->def[i] -= this->b[i];
        }

        return;
    }

    void setLocalSolution (const Entity& e)
    {
        Dune::GeometryType gt = e.geometry().type();
        const typename Dune::LagrangeShapeFunctionSetContainer<DT,RT,n>::value_type&
            sfs=Dune::LagrangeShapeFunctions<DT,RT,n>::general(gt, 1);
        int size = sfs.size();
        this->setcurrentsize(size);

        for (int i = 0; i < size; i++) {
            this->u[i] = currentSolution.evallocal(0, e, sfs[i].position());
            uold[i] = oldSolution.evallocal(0, e, sfs[i].position());
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

    template<class TypeTag>
    void assembleBC (const Entity& e)
    {
        // extract some important parameters
        Dune::GeometryType gt = e.geometry().type();
        const typename Dune::LagrangeShapeFunctionSetContainer<DT,RT,n>::value_type&
            sfs=Dune::LagrangeShapeFunctions<DT,RT,n>::general(gt,1);
        setcurrentsize(sfs.size());

        for (int i = 0; i < sfs.size(); i++)
            this->b[i] = 0;

        // determine quadrature order
        int p=2;
        // evaluate boundary conditions via intersection iterator
        typedef typename IntersectionIteratorGetter<G,TypeTag>::IntersectionIterator
            IntersectionIterator;

        IntersectionIterator endit = IntersectionIteratorGetter<G,TypeTag>::end(e);
        for (IntersectionIterator it = IntersectionIteratorGetter<G,TypeTag>::begin(e);
             it!=endit; ++it)
        {
            // if we have a neighbor then we assume there is no boundary (forget interior boundaries)
            // in level assemble treat non-level neighbors as boundary
            if (it->neighbor())
            {
                if (levelBoundaryAsDirichlet && it->outside()->level()==e.level())
                    continue;
                if (!levelBoundaryAsDirichlet)
                    continue;
            }

            // determine boundary condition type for this face, initialize with processor boundary
            typename BoundaryConditions::Flags bctypeface = BoundaryConditions::process;

            // handle face on exterior boundary, this assumes there are no interior boundaries
            if (it->boundary())
            {
                Dune::GeometryType gtface = it->intersectionSelfLocal().type();
                for (size_t g=0; g<Dune::QuadratureRules<DT,n-1>::rule(gtface,p).size(); ++g)
                {
                    const Dune::FieldVector<DT,n-1>& facelocal = Dune::QuadratureRules<DT,n-1>::rule(gtface,p)[g].position();
                    FieldVector<DT,n> local = it->intersectionSelfLocal().global(facelocal);
                    FieldVector<DT,n> global = it->intersectionGlobal().global(facelocal);
                    bctypeface = problem.bctype(global,e,local); // eval bctype


                    if (bctypeface!=BoundaryConditions::neumann) break;

                    RT J = problem.J(global,e,local);
                    if (fabs(J) < 1e-10)
                        continue;
                    double weightface = Dune::QuadratureRules<DT,n-1>::rule(gtface,p)[g].weight();
                    DT detjacface = it->intersectionGlobal().integrationElement(facelocal);
                    for (int i=0; i<sfs.size(); i++) // loop over test function number
                        if (this->bctype[i][0]==BoundaryConditions::neumann)
                        {
                            //////////////////////////////////////////////////////////////////////////
                            // HACK: piecewise constants with respect to dual grid not implemented yet
                            // works only if exactly one quadrature point is located within each dual
                            // cell boundary (which should be the case for p = 2)
                            //////////////////////////////////////////////////////////////////////////
                            if (sfs[i].evaluateFunction(0,local) > 0.5)
                                this->b[i] -= J*weightface*detjacface;
                        }
                }
                if (bctypeface==BoundaryConditions::neumann) continue; // was a neumann face, go to next face
            }

            // If we are here, then it is
            // (i)   an exterior boundary face with Dirichlet condition, or
            // (ii)  a processor boundary (i.e. neither boundary() nor neighbor() was true), or
            // (iii) a level boundary in case of level-wise assemble
            // How processor boundaries are handled depends on the processor boundary mode
            if (bctypeface==BoundaryConditions::process && procBoundaryAsDirichlet==false
                && levelBoundaryAsDirichlet==false)
                continue; // then it acts like homogeneous Neumann

            // now handle exterior or interior Dirichlet boundary
            for (int i=0; i<sfs.size(); i++) // loop over test function number
            {
                if (sfs[i].codim()==0) continue; // skip interior dof
                if (sfs[i].codim()==1) // handle face dofs
                {
                    if (sfs[i].entity()==it->numberInSelf())
                    {
                        if (this->bctype[i][0]<bctypeface)
                        {
                            this->bctype[i][0] = bctypeface;
                            if (bctypeface==BoundaryConditions::process)
                                this->b[i] = 0;
                            if (bctypeface==BoundaryConditions::dirichlet)
                            {
                                this->b[i] = 0;
                            }
                        }
                    }
                    continue;
                }
                // handle subentities of this face
                for (int j=0; j<ReferenceElements<DT,n>::general(gt).size(it->numberInSelf(),1,sfs[i].codim()); j++)
                    if (sfs[i].entity()==ReferenceElements<DT,n>::general(gt).subEntity(it->numberInSelf(),1,j,sfs[i].codim()))
                    {
                        if (this->bctype[i][0]<bctypeface)
                        {
                            this->bctype[i][0] = bctypeface;
                            if (bctypeface==BoundaryConditions::process)
                                this->b[i] = 0;
                            if (bctypeface==BoundaryConditions::dirichlet)
                            {
                                this->b[i] = 0;
                            }
                        }
                    }
            }
        }
    }

    // parameters given in constructor
    ParabolicProblem<G,RT>& problem;
    bool levelBoundaryAsDirichlet;
    bool procBoundaryAsDirichlet;
    const BoxFunction& currentSolution;
    BoxFunction oldSolution;
    double dt;
    VBlockType uold[SIZE];
};

/** @} */
}
#endif
