// $Id$

#ifndef DUNE_P1GROUNDWATER_HH
#define DUNE_P1GROUNDWATER_HH

#include <map>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <vector>
#include <sstream>

#include <dune/common/exceptions.hh>
#include <dune/grid/common/grid.hh>
#include <dune/grid/common/referenceelements.hh>
#include <dune/common/geometrytype.hh>
#include <dune/grid/common/quadraturerules.hh>

#include <dune/disc/shapefunctions/lagrangeshapefunctions.hh>
#include <dune/disc/operators/boundaryconditions.hh>
#include <dune/disc/operators/localstiffness.hh>
#include "dumux/fractionalflow/variableclass.hh"
#include "dumux/fractionalflow/fractionalflowproblem.hh"
#include "dumux/diffusion/diffusionproblem.hh"

/**
 * @file
 * @brief  compute local stiffness matrix for conforming finite elements for diffusion equation
 * @author Peter Bastian
 */


namespace Dune
{
/** @addtogroup DISC_Disc
 *
 * @{
 */
/**
 * @brief compute local stiffness matrix for conforming finite elements for diffusion equation
 *
 */


//! A class for computing local stiffness matrices
/*! A class for computing local stiffness matrix for the
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
template<class G, class RT, class Problem = FractionalFlowProblem<G, RT, VariableClass<G,RT> > >
class GroundwaterEquationLocalStiffness
    : public LinearLocalStiffness<typename G::LevelGridView,RT,1>
{
    template<int dim>
    struct ElementLayout
    {
        bool contains (Dune::GeometryType gt)
        {
            return gt.dim() == dim;
        }
    };

    // grid types
    typedef typename G::LevelGridView GV;
    typedef typename G::ctype DT;
    typedef typename G::Traits::template Codim<0>::Entity Entity;
    typedef typename G::Traits::LevelIndexSet IS;
    typedef Dune::MultipleCodimMultipleGeomTypeMapper<GV,ElementLayout> EM;
    typedef Dune::VariableClass<G, RT> VC;

public:
    // define the number of components of your system, this is used outside
    // to allocate the correct size of (dense) blocks with a FieldMatrix
    enum {n=G::dimension};
    enum {m=1};
    enum {SIZE=LagrangeShapeFunctionSetContainer<DT,RT,n>::maxsize};

    //! Constructor
    GroundwaterEquationLocalStiffness (Problem& params,
                                       bool levelBoundaryAsDirichlet_, const G& grid,
                                       int level = 0,
                                       bool procBoundaryAsDirichlet_=true)
        : problem(params),levelBoundaryAsDirichlet(levelBoundaryAsDirichlet_),
          procBoundaryAsDirichlet(procBoundaryAsDirichlet_), elementmapper(grid.levelView(level))
    {}


    //! assemble local stiffness matrix for given element and order
    /*! On exit the following things have been done:
      - The stiffness matrix for the given entity and polynomial degree has been assembled and is
      accessible with the mat() method.
      - The boundary conditions have been evaluated and are accessible with the bc() method
      - The right hand side has been assembled. It contains either the value of the essential boundary
      condition or the assembled source term and neumann boundary condition. It is accessible via the rhs() method.
      @param[in]  e    a codim 0 entity reference
      @param[in]  k    order of Lagrange basis
    */
    void assemble (const Entity& e, int k=1)
    {
        // extract some important parameters
        Dune::GeometryType gt = e.geometry().type();
        const typename Dune::LagrangeShapeFunctionSetContainer<DT,RT,n>::value_type&
            sfs=Dune::LagrangeShapeFunctions<DT,RT,n>::general(gt,k);
        setcurrentsize(sfs.size());

        // clear assemble data
        for (int i=0; i<sfs.size(); i++)
        {
            this->b[i] = 0;
            this->bctype[i][0] = BoundaryConditions::neumann;
            for (int j=0; j<sfs.size(); j++)
                this->A[i][j] = 0;
        }

        assembleV(e,k);
        assembleBC(e,k);
    }

    //! assemble only boundary conditions for given element
    /*! On exit the following things have been done:
      - The boundary conditions have been evaluated and are accessible with the bc() method
      - The right hand side contains either the value of the essential boundary
      condition or the assembled neumann boundary condition. It is accessible via the rhs() method.
      @param[in]  e    a codim 0 entity reference
      @param[in]  k    order of Lagrange basis
    */
    void assembleBoundaryCondition (const Entity& e, int k=1)
    {
        // extract some important parameters
        Dune::GeometryType gt = e.geometry().type();
        const typename Dune::LagrangeShapeFunctionSetContainer<DT,RT,n>::value_type&
            sfs=Dune::LagrangeShapeFunctions<DT,RT,n>::general(gt,k);
        setcurrentsize(sfs.size());

        // clear assemble data
        for (int i=0; i<sfs.size(); i++)
        {
            this->b[i] = 0;
            this->bctype[i][0] = BoundaryConditions::neumann;
        }

        this->assembleBC(e,k);
    }

private:

    void assembleV (const Entity& e, int k=1)
    {
        // extract some important parameters
        Dune::GeometryType gt = e.geometry().type();
        const typename Dune::LagrangeShapeFunctionSetContainer<DT,RT,n>::value_type&
            sfs=Dune::LagrangeShapeFunctions<DT,RT,n>::general(gt,k);
        setcurrentsize(sfs.size());

        // Loop over all quadrature points and assemble matrix and right hand side
        int p=2;
        if (gt.isSimplex()) p=1;
        if (k>1) p=2*(k-1);
        int elemId = elementmapper.map(e);
        for (size_t g=0; g<Dune::QuadratureRules<DT,n>::rule(gt,p).size(); ++g) // run through all quadrature points
        {
            const Dune::FieldVector<DT,n>&
                local = Dune::QuadratureRules<DT,n>::rule(gt,p)[g].position(); // pos of integration point
            Dune::FieldVector<DT,n> global = e.geometry().global(local);     // ip in global coordinates
            const Dune::FieldMatrix<DT,n,n>
                jac = e.geometry().jacobianInverseTransposed(local);           // eval jacobian inverse
            Dune::FieldMatrix<DT,n,n> K = problem.K(global,e,local);   // eval diffusion tensor
            K *= problem.materialLaw_.mobTotal(
                    problem.variables.sat(global,e,local), 
                    global, 
                    e, 
                    local);
            double weight = Dune::QuadratureRules<DT,n>::rule(gt,p)[g].weight();// weight of quadrature point
            DT detjac = e.geometry().integrationElement(local);              // determinant of jacobian
            RT q = problem.source(global,e,local);                                // source term
            RT factor = weight*detjac;

            // evaluate gradients at Gauss points
            Dune::FieldVector<DT,n> grad[SIZE], temp, gv;
            for (int i=0; i<sfs.size(); i++)
            {
                for (int l=0; l<n; l++)
                    temp[l] = sfs[i].evaluateDerivative(0,l,local);
                grad[i] = 0;
                jac.umv(temp,grad[i]); // transform gradient to global ooordinates
            }

            for (int i=0; i<sfs.size(); i++) // loop over test function number
            {
                // rhs
                this->b[i] += q*sfs[i].evaluateFunction(0,local)*factor;

                // matrix
                gv = 0;    K.umv(grad[i],gv); // multiply with diffusion tensor
                this->A[i][i] += (grad[i]*gv)*factor;
                for (int j=0; j<i; j++)
                {
                    RT t = (grad[j]*gv)*factor;
                    this->A[i][j] += t;
                    this->A[j][i] += t;
                }
            }
        }
    }


    void assembleBC (const Entity& e, int k=1)
    {
        // extract some important parameters
        Dune::GeometryType gt = e.geometry().type();
        const typename Dune::LagrangeShapeFunctionSetContainer<DT,RT,n>::value_type&
            sfs=Dune::LagrangeShapeFunctions<DT,RT,n>::general(gt,k);
        setcurrentsize(sfs.size());

        // determine quadrature order
        int p=2;
        if (gt.isSimplex()) p=1;
        if (k>1) p=2*(k-1);

        // evaluate boundary conditions via intersection iterator
        typedef typename G::LeafGridView::IntersectionIterator
            IntersectionIterator;

        IntersectionIterator endit = e.ileafend();
        for (IntersectionIterator it = e.ileafbegin();
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
                Dune::GeometryType gtface = it->geometryInInside().type();
                for (size_t g=0; g<Dune::QuadratureRules<DT,n-1>::rule(gtface,p).size(); ++g)
                {
                    const Dune::FieldVector<DT,n-1>& facelocal = Dune::QuadratureRules<DT,n-1>::rule(gtface,p)[g].position();
                    FieldVector<DT,n> local = it->geometryInInside().global(facelocal);
                    FieldVector<DT,n> global = it->geometry().global(facelocal);
                    bctypeface = problem.bctype(global,e,local); // eval bctype


                    if (bctypeface!=BoundaryConditions::neumann) break;

                    RT J = problem.neumannPress(global,e,local);
                    double weightface = Dune::QuadratureRules<DT,n-1>::rule(gtface,p)[g].weight();
                    DT detjacface = it->geometry().integrationElement(facelocal);
                    for (int i=0; i<sfs.size(); i++) // loop over test function number
                        if (this->bctype[i][0]==BoundaryConditions::neumann)
                        {
                            this->b[i] += J*sfs[i].evaluateFunction(0,local)*weightface*detjacface;
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
                    if (sfs[i].entity()==it->indexInInside())
                    {
                        if (this->bctype[i][0]<bctypeface)
                        {
                            this->bctype[i][0] = bctypeface;
                            if (bctypeface==BoundaryConditions::process)
                                this->b[i] = 0;
                            if (bctypeface==BoundaryConditions::dirichlet)
                            {
                                Dune::FieldVector<DT,n> global = e.geometry().global(sfs[i].position());
                                this->b[i] = problem.dirichletPress(global,e,sfs[i].position());
                            }
                        }
                    }
                    continue;
                }
                // handle subentities of this face
                for (int j=0; j<ReferenceElements<DT,n>::general(gt).size(it->indexInInside(),1,sfs[i].codim()); j++)
                    if (sfs[i].entity()==ReferenceElements<DT,n>::general(gt).subEntity(it->indexInInside(),1,j,sfs[i].codim()))
                    {
                        if (this->bctype[i][0]<bctypeface)
                        {
                            this->bctype[i][0] = bctypeface;
                            if (bctypeface==BoundaryConditions::process)
                                this->b[i] = 0;
                            if (bctypeface==BoundaryConditions::dirichlet)
                            {
                                Dune::FieldVector<DT,n> global = e.geometry().global(sfs[i].position());
                                this->b[i] = problem.dirichletPress(global,e,sfs[i].position());
                            }
                        }
                    }
            }
        }
    }

    // parameters given in constructor
    Problem& problem;
    bool levelBoundaryAsDirichlet;
    bool procBoundaryAsDirichlet;
    EM elementmapper;
};

/** @} */
}
#endif
