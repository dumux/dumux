// $Id$

#ifndef DUNE_P1PARABOLICJACOBIAN_HH
#define DUNE_P1PARABOLICJACOBIAN_HH

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
template<class G, class RT, class P1Function = LeafP1Function<G, RT> >
class P1ParabolicLocalJacobian
    : public LocalJacobian<P1ParabolicLocalJacobian<G,RT>,G,RT,1>
{
    typedef typename G::ctype DT;
    typedef typename G::Traits::template Codim<0>::Entity Entity;
    typedef typename LocalJacobian<P1ParabolicLocalJacobian<G,RT>,G,RT,1>::VBlockType VBlockType;
    typedef typename LocalJacobian<P1ParabolicLocalJacobian<G,RT>,G,RT,1>::MBlockType MBlockType;

public:
    // define the number of components of your system, this is used outside
    // to allocate the correct size of (dense) blocks with a FieldMatrix
    enum {n=G::dimension};
    enum {m=1};
    enum {SIZE=LagrangeShapeFunctionSetContainer<DT,RT,n>::maxsize};

    //! Constructor
    P1ParabolicLocalJacobian (ParabolicProblem<G,RT>& params,
                              bool levelBoundaryAsDirichlet_, const G& grid,
                              P1Function& sol,
                              bool procBoundaryAsDirichlet_=true)
        : problem(params),levelBoundaryAsDirichlet(levelBoundaryAsDirichlet_),
          procBoundaryAsDirichlet(procBoundaryAsDirichlet_),
          currentSolution(sol), oldSolution(grid), dt(1)
    {}


    template<class TypeTag>
    void localDefect (const Entity& e, const VBlockType* sol)
    {
        // extract some important parameters
        Dune::GeometryType gt = e.geometry().type();
        const typename Dune::LagrangeShapeFunctionSetContainer<DT,RT,n>::value_type&
            sfs=Dune::LagrangeShapeFunctions<DT,RT,n>::general(gt, 1);
        int size = sfs.size();
        this->setcurrentsize(size);

        MBlockType locstiff[size][size];
        MBlockType locmass[size][size];
        for (int i = 0; i < size; i++) {
            this->bctype[i][0] = BoundaryConditions::neumann;
            this->b[i] = 0;
            for (int j = 0; j < size; j++) {
                locstiff[i][j] = 0;
                locmass[i][j] = 0;
            }
        }

        // Loop over all quadrature points and assemble matrix and right hand side
        int p=2;
        //if (gt.isSimplex()) p=1;
        for (size_t g=0; g<Dune::QuadratureRules<DT,n>::rule(gt,p).size(); ++g) // run through all quadrature points
        {
            const Dune::FieldVector<DT,n>&
                local = Dune::QuadratureRules<DT,n>::rule(gt,p)[g].position(); // pos of integration point
            Dune::FieldVector<DT,n> global = e.geometry().global(local);     // ip in global coordinates
            const Dune::FieldMatrix<DT,n,n>
                jac = e.geometry().jacobianInverseTransposed(local);           // eval jacobian inverse
            Dune::FieldMatrix<DT,n,n> K = problem.K(global,e,local);   // eval diffusion tensor
            double weight = Dune::QuadratureRules<DT,n>::rule(gt,p)[g].weight();// weight of quadrature point
            DT detjac = e.geometry().integrationElement(local);              // determinant of jacobian
            RT q = problem.q(global,e,local);                                // source term
            RT factor = weight*detjac;

            // evaluate gradients at Gauss points
            Dune::FieldVector<DT,n> grad[SIZE], temp, gv;
            RT sfsVal[size];
            for (int i=0; i<size; i++)
            {
                for (int l=0; l<n; l++)
                    temp[l] = sfs[i].evaluateDerivative(0,l,local);
                grad[i] = 0;
                jac.umv(temp,grad[i]); // transform gradient to global coordinates
                sfsVal[i] = sfs[i].evaluateFunction(0,local);
            }

            for (int i=0; i<size; i++) // loop over test function number
            {
                // rhs
                this->b[i] += q*sfs[i].evaluateFunction(0,local)*factor;

                // matrix
                gv = 0;    K.umv(grad[i],gv); // multiply with diffusion tensor
                locstiff[i][i] += (grad[i]*gv)*factor;
                locmass[i][i] += sfsVal[i]*sfsVal[i]*factor/dt;
                for (int j=0; j<i; j++)
                {
                    RT entry = (grad[j]*gv)*factor;
                    locstiff[i][j] += entry;
                    locstiff[j][i] += entry;
                    entry = sfsVal[i]*sfsVal[j]*factor/dt;
                    locmass[i][j] += entry;
                    locmass[j][i] += entry;
                }
            }
        }

        assembleBC<TypeTag> (e);

        for (int i=0; i<size; i++) {
            this->def[i] = -this->b[i];
            for (int j=0; j<size; j++)
                this->def[i] += (locmass[i][j] + locstiff[i][j])*sol[j]- locmass[i][j]*uold[j];
        }


        //std::cout << locmass[0][0] << ", " << locmass[0][1] << ", " << locmass[0][2] << ", " << locmass[0][3] << std::endl;
        //std::cout << locmass[1][0] << ", " << locmass[1][1] << ", " << locmass[1][2] << ", " << locmass[1][3] << std::endl;
        //std::cout << locmass[2][0] << ", " << locmass[2][1] << ", " << locmass[2][2] << ", " << locmass[2][3] << std::endl;
        //std::cout << locmass[3][0] << ", " << locmass[3][1] << ", " << locmass[3][2] << ", " << locmass[3][3] << std::endl << std::endl;
        //std::cout << "sol = " << sol[0] << ", " << sol[1] << ", " << sol[2] << ", " << sol[3] << std::endl;
        //std::cout << "uold = " << uold[0] << ", " << uold[1] << ", " << uold[2] << ", " << uold[3] << std::endl;
        //std::cout << "defect = " << this->def[0] << ", " << this->def[1] << ", " << this->def[2] << ", " << this->def[3] << std::endl;
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

    void setOldSolution (P1Function& uOld)
    {
        *oldSolution = *uOld;
    }

private:
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
        if (gt.isSimplex()) p=1;
        // evaluate boundary conditions via intersection iterator
        typedef typename IntersectionIteratorGetter<G,TypeTag>::IntersectionIterator
            IntersectionIterator;

        IntersectionIterator endit = IntersectionIteratorGetter<G,TypeTag>::end(e);
        for (IntersectionIterator it = IntersectionIteratorGetter<G,TypeTag>::begin(e);
             it!=endit; ++it)
        {
            // if we have a neighbor then we assume there is no boundary (forget interior boundaries)
            // in level assemble treat non-level neighbors as boundary
            if (it.neighbor())
            {
                if (levelBoundaryAsDirichlet && it.outside()->level()==e.level())
                    continue;
                if (!levelBoundaryAsDirichlet)
                    continue;
            }

            // determine boundary condition type for this face, initialize with processor boundary
            typename BoundaryConditions::Flags bctypeface = BoundaryConditions::process;

            // handle face on exterior boundary, this assumes there are no interior boundaries
            if (it.boundary())
            {
                Dune::GeometryType gtface = it.intersectionSelfLocal().type();
                for (size_t g=0; g<Dune::QuadratureRules<DT,n-1>::rule(gtface,p).size(); ++g)
                {
                    const Dune::FieldVector<DT,n-1>& facelocal = Dune::QuadratureRules<DT,n-1>::rule(gtface,p)[g].position();
                    FieldVector<DT,n> local = it.intersectionSelfLocal().global(facelocal);
                    FieldVector<DT,n> global = it.intersectionGlobal().global(facelocal);
                    bctypeface = problem.bctype(global,e,local); // eval bctype


                    if (bctypeface!=BoundaryConditions::neumann) break;

                    RT J = problem.J(global,e,local);
                    double weightface = Dune::QuadratureRules<DT,n-1>::rule(gtface,p)[g].weight();
                    DT detjacface = it.intersectionGlobal().integrationElement(facelocal);
                    for (int i=0; i<sfs.size(); i++) // loop over test function number
                        if (this->bctype[i][0]==BoundaryConditions::neumann)
                        {
                            this->b[i] -= J*sfs[i].evaluateFunction(0,local)*weightface*detjacface;
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
                    if (sfs[i].entity()==it.numberInSelf())
                    {
                        if (this->bctype[i][0]<bctypeface)
                        {
                            this->bctype[i][0] = bctypeface;
                            if (bctypeface==BoundaryConditions::process)
                                this->b[i] = 0;
                            if (bctypeface==BoundaryConditions::dirichlet)
                            {
                                Dune::FieldVector<DT,n> global = e.geometry().global(sfs[i].position());
                                this->b[i] = 0;
                            }
                        }
                    }
                    continue;
                }
                // handle subentities of this face
                for (int j=0; j<ReferenceElements<DT,n>::general(gt).size(it.numberInSelf(),1,sfs[i].codim()); j++)
                    if (sfs[i].entity()==ReferenceElements<DT,n>::general(gt).subEntity(it.numberInSelf(),1,j,sfs[i].codim()))
                    {
                        if (this->bctype[i][0]<bctypeface)
                        {
                            this->bctype[i][0] = bctypeface;
                            if (bctypeface==BoundaryConditions::process)
                                this->b[i] = 0;
                            if (bctypeface==BoundaryConditions::dirichlet)
                            {
                                Dune::FieldVector<DT,n> global = e.geometry().global(sfs[i].position());
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
    const P1Function& currentSolution;
    P1Function oldSolution;
    double dt;
    VBlockType uold[SIZE];
};

/** @} */
}
#endif
