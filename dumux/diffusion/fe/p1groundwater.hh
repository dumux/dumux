// $Id$
/*****************************************************************************
 *   Copyright (C) 2007-2009 by Peter Bastian, Bernd Flemisch                *
 *   Institute of Hydraulic Engineering                                      *
 *   University of Stuttgart, Germany                                        *
 *   email: <givenname>.<name>@iws.uni-stuttgart.de                          *
 *                                                                           *
 *   This program is free software; you can redistribute it and/or modify    *
 *   it under the terms of the GNU General Public License as published by    *
 *   the Free Software Foundation; either version 2 of the License, or       *
 *   (at your option) any later version, as long as this copyright notice    *
 *   is included in its original form.                                       *
 *                                                                           *
 *   This program is distributed WITHOUT ANY WARRANTY.                       *
 *****************************************************************************/
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
#include "dumux/fractionalflow/variableclass2p.hh"
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
  - Scalar    type used for return values
*/
template<class GridView, class Scalar, class Problem = FractionalFlowProblem<GridView, Scalar, VariableClass<GridView,Scalar> > >
class GroundwaterEquationLocalStiffness
    : public LinearLocalStiffness<GridView,Scalar,1>
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
    typedef typename GridView::Traits::template Codim<0>::Entity Element;
    typedef Dune::MultipleCodimMultipleGeomTypeMapper<GridView,ElementLayout> EM;
    typedef Dune::VariableClass<GridView, Scalar> VC;
    typedef typename GridView::IntersectionIterator IntersectionIterator;

public:
    // define the number of components of your system, this is used outside
    // to allocate the correct size of (dense) blocks with a FieldMatrix
    enum {dim=GridView::dimension};
    enum {numEq=1};
    enum {SIZE=LagrangeShapeFunctionSetContainer<Scalar,Scalar,dim>::maxsize};

    //! Constructor
    GroundwaterEquationLocalStiffness (Problem& params,
                                       bool levelBoundaryAsDirichlet_, const GridView& gridView,
                                       int level = 0,
                                       bool procBoundaryAsDirichlet_=true)
        : problem(params),levelBoundaryAsDirichlet(levelBoundaryAsDirichlet_),
          procBoundaryAsDirichlet(procBoundaryAsDirichlet_), elementmapper(gridView), gridView_(gridView)
    {}


    //! assemble local stiffness matrix for given element and order
    /*! On exit the following things have been done:
      - The stiffness matrix for the given entity and polynomial degree has been assembled and is
      accessible with the mat() method.
      - The boundary conditions have been evaluated and are accessible with the bc() method
      - The right hand side has been assembled. It contains either the value of the essential boundary
      condition or the assembled source term and neumann boundary condition. It is accessible via the rhs() method.
      @param[in]  element    a codim 0 entity reference
      @param[in]  k    order of Lagrange basis
    */
    void assemble (const Element& element, int k=1)
    {
        // extract some important parameters
        Dune::GeometryType gt = element.geometry().type();
        const typename Dune::LagrangeShapeFunctionSetContainer<Scalar,Scalar,dim>::value_type&
            sfs=Dune::LagrangeShapeFunctions<Scalar,Scalar,dim>::general(gt,k);
        setcurrentsize(sfs.size());

        // clear assemble data
        for (int i=0; i<sfs.size(); i++)
        {
            this->b[i] = 0;
            this->bctype[i][0] = BoundaryConditions::neumann;
            for (int j=0; j<sfs.size(); j++)
                this->A[i][j] = 0;
        }

        assembleV(element,k);
        assembleBC(element,k);
    }

    //! assemble only boundary conditions for given element
    /*! On exit the following things have been done:
      - The boundary conditions have been evaluated and are accessible with the bc() method
      - The right hand side contains either the value of the essential boundary
      condition or the assembled neumann boundary condition. It is accessible via the rhs() method.
      @param[in]  element    a codim 0 entity reference
      @param[in]  k    order of Lagrange basis
    */
    void assembleBoundaryCondition (const Element& element, int k=1)
    {
        // extract some important parameters
        Dune::GeometryType gt = element.geometry().type();
        const typename Dune::LagrangeShapeFunctionSetContainer<Scalar,Scalar,dim>::value_type&
            sfs=Dune::LagrangeShapeFunctions<Scalar,Scalar,dim>::general(gt,k);
        setcurrentsize(sfs.size());

        // clear assemble data
        for (int i=0; i<sfs.size(); i++)
        {
            this->b[i] = 0;
            this->bctype[i][0] = BoundaryConditions::neumann;
        }

        this->assembleBC(element,k);
    }

private:

    void assembleV (const Element& element, int k=1)
    {
        // extract some important parameters
        Dune::GeometryType gt = element.geometry().type();
        const typename Dune::LagrangeShapeFunctionSetContainer<Scalar,Scalar,dim>::value_type&
            sfs=Dune::LagrangeShapeFunctions<Scalar,Scalar,dim>::general(gt,k);
        setcurrentsize(sfs.size());

        // Loop over all quadrature points and assemble matrix and right hand side
        int p=2;
        if (gt.isSimplex()) p=1;
        if (k>1) p=2*(k-1);
//        int elemId = elementmapper.map(element);
        for (size_t g=0; g<Dune::QuadratureRules<Scalar,dim>::rule(gt,p).size(); ++g) // run through all quadrature points
        {
            const Dune::FieldVector<Scalar,dim>&
                local = Dune::QuadratureRules<Scalar,dim>::rule(gt,p)[g].position(); // pos of integration point
            Dune::FieldVector<Scalar,dim> global = element.geometry().global(local);     // ip in global coordinates
            const Dune::FieldMatrix<Scalar,dim,dim>
                jac = element.geometry().jacobianInverseTransposed(local);           // eval jacobian inverse
            Dune::FieldMatrix<Scalar,dim,dim> K = problem.soil().K(global,element,local);   // eval diffusion tensor
            K *= problem.materialLaw().mobTotal(
                    problem.variables().satElement(element),
                    global,
                    element,
                    local);
            double weight = Dune::QuadratureRules<Scalar,dim>::rule(gt,p)[g].weight();// weight of quadrature point
            Scalar detjac = element.geometry().integrationElement(local);              // determinant of jacobian
            Scalar q = problem.sourcePress(global,element,local);                                // source term
            Scalar factor = weight*detjac;

            // evaluate gradients at Gauss points
            Dune::FieldVector<Scalar,dim> grad[SIZE], temp, gv;
            for (int i=0; i<sfs.size(); i++)
            {
                for (int l=0; l<dim; l++)
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
                    Scalar t = (grad[j]*gv)*factor;
                    this->A[i][j] += t;
                    this->A[j][i] += t;
                }
            }
        }
    }


    void assembleBC (const Element& element, int k=1)
    {
        // extract some important parameters
        Dune::GeometryType gt = element.geometry().type();
        const typename Dune::LagrangeShapeFunctionSetContainer<Scalar,Scalar,dim>::value_type&
            sfs=Dune::LagrangeShapeFunctions<Scalar,Scalar,dim>::general(gt,k);
        setcurrentsize(sfs.size());

        // determine quadrature order
        int p=2;
        if (gt.isSimplex()) p=1;
        if (k>1) p=2*(k-1);

        // evaluate boundary conditions via intersection iterator
        IntersectionIterator endit = gridView_.template iend(element);
        for (IntersectionIterator it = gridView_.template ibegin(element); it!=endit; ++it)
        {
            // if we have a neighbor then we assume there is no boundary (forget interior boundaries)
            // in level assemble treat non-level neighbors as boundary
            if (it->neighbor())
            {
                if (levelBoundaryAsDirichlet && it->outside()->level()==element.level())
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
                for (size_t g=0; g<Dune::QuadratureRules<Scalar,dim-1>::rule(gtface,p).size(); ++g)
                {
                    const Dune::FieldVector<Scalar,dim-1>& facelocal = Dune::QuadratureRules<Scalar,dim-1>::rule(gtface,p)[g].position();
                    FieldVector<Scalar,dim> local = it->geometryInInside().global(facelocal);
                    FieldVector<Scalar,dim> global = it->geometry().global(facelocal);
                    bctypeface = problem.bctypePress(global,element,local); // eval bctype


                    if (bctypeface!=BoundaryConditions::neumann) break;

                    Scalar J = problem.neumannPress(global,element,local);
                    double weightface = Dune::QuadratureRules<Scalar,dim-1>::rule(gtface,p)[g].weight();
                    Scalar detjacface = it->geometry().integrationElement(facelocal);
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
            // (ii)  a processor boundary (i.element. neither boundary() nor neighbor() was true), or
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
                                Dune::FieldVector<Scalar,dim> global = element.geometry().global(sfs[i].position());
                                this->b[i] = problem.dirichletPress(global,element,sfs[i].position());
                            }
                        }
                    }
                    continue;
                }
                // handle subentities of this face
                for (int j=0; j<ReferenceElements<Scalar,dim>::general(gt).size(it->indexInInside(),1,sfs[i].codim()); j++)
                    if (sfs[i].entity()==ReferenceElements<Scalar,dim>::general(gt).subEntity(it->indexInInside(),1,j,sfs[i].codim()))
                    {
                        if (this->bctype[i][0]<bctypeface)
                        {
                            this->bctype[i][0] = bctypeface;
                            if (bctypeface==BoundaryConditions::process)
                                this->b[i] = 0;
                            if (bctypeface==BoundaryConditions::dirichlet)
                            {
                                Dune::FieldVector<Scalar,dim> global = element.geometry().global(sfs[i].position());
                                this->b[i] = problem.dirichletPress(global,element,sfs[i].position());
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
    const GridView& gridView_;
};

/** @} */
}
#endif
