/*****************************************************************************
 *   Copyright (C) 2007 by Peter Bastian                                     *
 *   Institute of Parallel and Distributed System                            *
 *   Department Simulation of Large Systems                                  *
 *   University of Stuttgart, Germany                                        *
 *                                                                           *
 *   Copyright (C) 2008 by Andreas Lauser, Bernd Flemisch                    *
 *   Institute of Hydraulic Engineering                                      *
 *   University of Stuttgart, Germany                                        *
 *   email: and _at_ poware.org                                              *
 *                                                                           *
 *   This program is free software; you can redistribute it and/or modify    *
 *   it under the terms of the GNU General Public License as published by    *
 *   the Free Software Foundation; either version 2 of the License, or       *
 *   (at your option) any later version, as long as this copyright notice    *
 *   is included in its original form.                                       *
 *                                                                           *
 *   This program is distributed WITHOUT ANY WARRANTY.                       *
 *****************************************************************************/
/**
 * @file
 * @brief  compute local jacobian matrix for conforming finite elements for diffusion equation
 * @author Peter Bastian
 */
#ifndef DUMUX_BOX_JACOBIAN_HH
#define DUMUX_BOX_JACOBIAN_HH

#include <dune/disc/operators/localstiffness.hh>

#include <dumux/fvgeometry/fvelementgeometry.hh>

#include <boost/format.hpp>

#ifdef HAVE_VALGRIND
#include <valgrind/memcheck.h>
#elif !defined VALGRIND_CHECK_MEM_IS_DEFINED
#define VALGRIND_CHECK_MEM_IS_DEFINED(addr, s)
#define VALGRIND_MAKE_MEM_UNDEFINED(addr, s)
#endif // HAVE_VALGRIND


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

  - Problem        The specification of the problem (grid, boundary conditions, ...)
  - JacobianImp    The Box model specific part of the Jacobian
*/
template<class Problem, 
         class BoxTraitsT, 
         class JacobianImp,
         class VertexData>
class BoxJacobian : public Dune::LocalStiffness<typename Problem::DomainTraits::Grid::LeafGridView,
                                                typename Problem::DomainTraits::Scalar,
                                                BoxTraitsT::numEq>

{
private:
    typedef typename Problem::DomainTraits DomainTraits;
    typedef BoxTraitsT                     BoxTraits;

    enum {
        numEq = BoxTraits::numEq,

        dim     = DomainTraits::dim,
        dimWorld    = DomainTraits::dimWorld
    };
    typedef BoxJacobian ThisType;

    typedef typename Problem::DomainTraits::Grid::LeafGridView    LeafGridView;
    typedef typename DomainTraits::Grid                           Grid;
    typedef typename DomainTraits::Scalar                         Scalar;
    typedef typename DomainTraits::CoordScalar                    CoordScalar;
    typedef typename DomainTraits::LocalPosition                  LocalPosition;
    typedef typename DomainTraits::GlobalPosition                 GlobalPosition;

    typedef typename LeafGridView::template Codim<0>::Entity      Element;
    typedef typename Element::EntityPointer                       ElementPointer;
    typedef typename DomainTraits::ReferenceElement               ReferenceElement;
    typedef typename DomainTraits::IntersectionIterator           IntersectionIterator;

    typedef typename Element::Geometry                     Geometry;

    typedef typename BoxTraits::SpatialFunction            SpatialFunction;
    typedef typename BoxTraits::SolutionVector             SolutionVector;
    typedef typename BoxTraits::BoundaryTypeVector         BoundaryTypeVector;
    typedef typename BoxTraits::FVElementGeometry          FVElementGeometry;
    typedef typename BoxTraits::LocalFunction              LocalFunction;

    typedef typename BoxTraits::ShapeFunctionSetContainer  ShapeFunctionSetContainer;
    typedef typename ShapeFunctionSetContainer::value_type ShapeFunctionSet;

    typedef typename std::vector<VertexData> VertexDataArray;

public:
    BoxJacobian(Problem &problem)
        : problem_(problem),
          curElementPtr_(problem.elementEnd()),
          
          curElemDat_(BoxTraits::ShapeFunctionSetContainer::maxsize),
          prevElemDat_(BoxTraits::ShapeFunctionSetContainer::maxsize),

          curSolution_(NULL),
          oldSolution_(NULL)
    {
    }

    /*!
     * \brief Set the global solution of the current time step.
     */
    void setCurSolution(SpatialFunction *uCur)
    {
        curSolution_ = uCur;
    }

    /*!
     * \brief Set the global solution of the previous time step.
     *
     * This is required for implicit Euler time integration.
     *
     * TODO: If we would use a different time integration scheme
     *       we might need more old solutions, so it would be nice
     *       to have a non-hacky way to accomplish this.
     */
    void setOldSolution(SpatialFunction *uOld)
    {
        oldSolution_ = uOld;
    }

    /*!
     * \brief Assemble the linear system of equations for the
     *        verts of a element, given a local solution 'localU'.
     */
    void assemble(const Element &element, const LocalFunction& localU, int orderOfShapeFns = 1)
    {
        // set the current grid element
        asImp_().setCurrentElement(element);

#if HAVE_VALGRIND
        for (size_t i = 0; i < localU.size(); ++i)
            VALGRIND_CHECK_MEM_IS_DEFINED(&localU[i],
                                          sizeof(SolutionVector));
#endif // HAVE_VALGRIND

        LocalFunction mutableLocalU(localU);
        assemble_(element, mutableLocalU, orderOfShapeFns);
    }

    /*!
     * \brief Assemble the linear system of equations for the
     *        verts of a element.
     */
    void assemble(const Element &element, int orderOfShapeFns = 1)
    {
        // set the current grid element
        asImp_().setCurrentElement(element);

        int numVertices = curElementGeom_.numVertices;
        LocalFunction localU(numVertices);
        restrictToElement(localU, *curSolution_);
        assemble_(element, localU, orderOfShapeFns);
    }
    
    /*!
     * \brief Express the boundary conditions for a element in terms
     *        of a linear equation system.
     */
    void assembleBoundaryCondition(const Element &element, int orderOfShapeFns=1)
    {
        // set the current grid element
        asImp_().setCurrentElement(element);

        // reset the right hand side and the boundary
        // condition type vector
        resetRhs_();

        // set the boundary types
        updateBoundaryTypes_();

        // apply the neumann conditions
        applyBoundaryCondition_();
    };

    void updateBoundaryTypes(const Element &element)
    {
        // set the current grid element
        asImp_().setCurrentElement(element);
        updateBoundaryTypes_();
    }



    void applyBoundaryCondition_()
    {
        Dune::GeometryType      geoType = curElement_().geometry().type();
        const ReferenceElement &refElem = DomainTraits::referenceElement(geoType);

        // evaluate boundary conditions for all intersections of
        // the current element
        IntersectionIterator isIt = curElement_().ileafbegin();
        const IntersectionIterator &endIt = curElement_().ileafend();
        for (; isIt != endIt; ++isIt)
        {
            // handle only faces on boundaries.
            if (!isIt->boundary())
                continue;

            // Assemble the boundary for all verts of the
            // current face
            int faceIdx = isIt->indexInInside();
            int numFaceVerts = refElem.size(faceIdx, 1, dim);
            for (int faceVertIdx = 0;
                 faceVertIdx < numFaceVerts;
                 ++faceVertIdx)
            {
                int elemVertIdx = refElem.subEntity(faceIdx,
                                                    1,
                                                    faceVertIdx,
                                                    dim);
                int boundaryFaceIdx =
                    curElementGeom_.boundaryFaceIndex(faceIdx,
                                                      faceVertIdx);

                // handle boundary conditions for a single
                // sub-control volume face
                applyBoundaryCondition_(isIt,
                                        elemVertIdx,
                                        boundaryFaceIdx);
            }
        }
    }

    // handle boundary conditions for a single
    // sub-control volume face
    void applyBoundaryCondition_(const IntersectionIterator &isIt,
                                 int scvIdx,
                                 int boundaryFaceIdx)
    {
        // temporary vector to store the neumann boundaries
        SolutionVector values(0.0);
        bool wasEvaluated = false;

        // loop over all primary variables to deal with mixed
        // boundary conditions
        for (int eqIdx = 0; eqIdx < numEq; ++eqIdx)
        {
            if (this->bctype[scvIdx][eqIdx]
                != BoundaryConditions::neumann)
            {
                // dirichlet boundary conditions are treated as
                // zeros on the right hand side
                this->b[scvIdx][eqIdx] = 0;
                continue;
            }

            // if we are here we've got a neumann boundary condition

            // make sure that we evaluate call the problem's
            // neumann() method exactly once if
            if (!wasEvaluated)
            {
                // make sure that we only evaluate
                // the neumann fluxes once
                wasEvaluated = true;

                problem_.neumann(values,
                                 curElement_(),
                                 curElementGeom_,
                                 isIt,
                                 scvIdx,
                                 boundaryFaceIdx);
                // TODO (?): multiple integration
                // points
                values *= curElementGeom_.boundaryFace[boundaryFaceIdx].area;
            }
            VALGRIND_CHECK_MEM_IS_DEFINED(&values[eqIdx], 
                                          sizeof(Scalar));

            this->b[scvIdx][eqIdx] += values[eqIdx];
            VALGRIND_CHECK_MEM_IS_DEFINED(&this->b[scvIdx][eqIdx], 
                                          sizeof(Scalar));
        }
    }

    /*!
     * \brief Compute the local residual, i.e. the right hand side
     *        of an equation we would like to have zero.
     */
    void evalLocalResidual(LocalFunction &residual,
                           bool withBoundary = true)
    {
        // reset residual
        for (int i = 0; i < curElementGeom_.numVertices; i++) {
            residual[i] = 0;
        }

        // evaluate the local rate
        for (int i=0; i < curElementGeom_.numVertices; i++)
        {
            SolutionVector massContrib(0), tmp(0);

            // mass balance within the element. this is the
            // $\frac{m}{\partial t}$ term if using implicit
            // euler as time discretization.
            //
            // TODO (?): we might need a more explicit way for
            // doing the time discretization...
            this->asImp_().computeStorage(massContrib, i, false);
            this->asImp_().computeStorage(tmp, i, true);

            massContrib -= tmp;
            massContrib *= curElementGeom_.subContVol[i].volume/problem_.timeStepSize();
            residual[i] += massContrib;

            // subtract the source term from the local rate
            SolutionVector source;
            this->asImp_().computeSource(source, i);
            source *= curElementGeom_.subContVol[i].volume;
            residual[i] -= source;
            VALGRIND_CHECK_MEM_IS_DEFINED(&residual[i],
                                          sizeof(SolutionVector));
        }

        // calculate the mass flux over the faces and subtract
        // it from the local rates
        for (int k = 0; k < curElementGeom_.numEdges; k++)
        {
            int i = curElementGeom_.subContVolFace[k].i;
            int j = curElementGeom_.subContVolFace[k].j;

            SolutionVector flux;
            VALGRIND_MAKE_MEM_UNDEFINED(&flux, sizeof(SolutionVector));
            this->asImp_().computeFlux(flux, k);
            VALGRIND_CHECK_MEM_IS_DEFINED(&flux, sizeof(SolutionVector));

            // subtract fluxes from the local mass rates of
            // the respective sub control volume adjacent to
            // the face.
            residual[i] -= flux;
            residual[j] += flux;
        }

        if (withBoundary) {
            assembleBoundaryCondition(this->curElement_());
            for (int i = 0; i < curElementGeom_.numVertices; i++) {
                residual[i] += this->b[i];
            }
        }

#if HAVE_VALGRIND
        for (int i=0; i < curElementGeom_.numVertices; i++)
            VALGRIND_CHECK_MEM_IS_DEFINED(&residual[i],
                                          sizeof(SolutionVector));
#endif // HAVE_VALGRIND
    }

    /*!
     * \brief Set the current grid element.
     */
    void setCurrentElement(const Element &element)
    {
        if (curElementPtr_ != element) {
            // update the FV element geometry if the current element
            // is to be changed
            curElementPtr_ = element;
            curElementGeom_.update(curElement_());

            // tell LocalStiffness (-> parent class) the current
            // number of degrees of freedom
            this->setcurrentsize(curElementGeom_.numVertices);
        }
    };


    /*!
     * \brief Restrict the global function 'globalFn' to the vertices
     *        of the current element, save the result to 'dest'.
     */
    void restrictToElement(LocalFunction &dest,
                           const SpatialFunction &globalFn) const
    {
        // we assert that the i-th shape function is
        // associated to the i-th vert of the element.
        int n = curElement_().template count<dim>();
        for (int i = 0; i < n; i++) {
            dest[i] = (*globalFn)[problem_.vertexIdx(curElement_(), i)];
        }
    }

    /*!
     * \brief Initialize the static data of all elements. The current
     *        solution is not yet valid when this method is
     *        called. (updateStaticData() is called as soon the
     *        current solution is valid.)
     *
     * This method should be overwritten by the child class if
     * necessary.
     */
    void initStaticData()
    { };

    /*!
     * \brief Update the static data of all elements with the current solution
     *
     * This method should be overwritten by the child class if
     * necessary.
     */
    void updateStaticData(SpatialFunction &curSol, SpatialFunction &oldSol)
    { };

    /*!
     * \brief Update the model specific vertex data of a whole
     *        element.
     */
    void updateElementData_(VertexDataArray &dest, const LocalFunction &sol, bool isOldSol)
    {
#ifdef ENABLE_VALGRIND
        for (int i = 0; i < elemDat.size(); ++i)
            VALGRIND_MAKE_MEM_UNDEFINED(&elemDat[i], 
                                        sizeof(VertexData));
#endif // ENABLE_VALGRIND

        int numVertices = this->curElement_().template count<dim>();
        for (int vertIdx = 0; vertIdx < numVertices; vertIdx++) {
            dest[vertIdx].update(sol[vertIdx],
                                 this->curElement_(),
                                 vertIdx,
                                 isOldSol,
                                 asImp_());
        }
    }


    /*!
     * \brief This returns the finite volume geometric information of the current
     *        element. <b>Use this method with great care!</br>
     */
    const FVElementGeometry &curFvElementGeometry() const
    { return curElementGeom_; }

    /*!
     * \brief Set current local solution
     */
    void setCurrentSolution(const LocalFunction &sol)
    {
        asImp_().updateElementData_(curElemDat_, sol, false);
    }

    /*!
     * \brief Set local solution of the last time step
     */
    void setPreviousSolution(const LocalFunction &sol)
    {
        asImp_().updateElementData_(prevElemDat_, sol, true);
    }

    /*!
     * \brief Vary a single component of a single vert of the
     *        local solution for the current element.
     *
     * This method is a optimization, since if varying a single
     * component at a degree of freedom not the whole element cache
     * needs to be recalculated. (Updating the element cache is very
     * expensive since material laws need to be evaluated.)
     */
    void deflectCurrentSolution(LocalFunction &curSol, 
                                int vertIdx,
                                int eqIdx,
                                Scalar value)
    {
        // stash away the orignial vertex data so that we do not need
        // to re calculate them if restoreCurSolution() is called.
        curVertexDataStash_ = curElemDat_[vertIdx];

        // recalculate the vertex data for the box which should be
        // changed
        curSol[vertIdx][eqIdx] = value;
        
        VALGRIND_MAKE_MEM_UNDEFINED(&curElemDat_[vertIdx], 
                                    sizeof(VertexData));

        curElemDat_[vertIdx].template update<JacobianImp>(curSol[vertIdx], 
                                                          this->curElement_(),
                                                          vertIdx, 
                                                          false,
                                                          asImp_());
        VALGRIND_CHECK_MEM_IS_DEFINED(&curElemDat_[vertIdx], 
                                      sizeof(VertexData));
    }

    /*!
     * \brief Restore the local jacobian to the state before
     *        deflectCurSolution() was called.
     *
     * This only works if deflectSolution was only called with
     * (vert, component) as arguments.
     */
    void restoreCurrentSolution(LocalFunction &curSol, 
                                int vertIdx, 
                                int eqIdx,
                                Scalar origValue)
    {
        curSol[vertIdx][eqIdx] = origValue;
        curElemDat_[vertIdx] = curVertexDataStash_;
    };

    /*!
     * \brief Returns a reference to the problem.
     */
    Problem &problem()
    { return problem_; };

    /*!
     * \brief Returns a reference to the problem.
     */
    const Problem &problem() const
    { return problem_; };

protected:
    const Element &curElement_() const
    { return *curElementPtr_; }

    // The problem we would like to solve
    Problem &problem_;

    ElementPointer    curElementPtr_;
    FVElementGeometry curElementGeom_;

    // current and previous element data. (this is model specific.)
    VertexDataArray  curElemDat_;
    VertexDataArray  prevElemDat_;
    
    // temporary variable to store the variable vertex data
    VertexData   curVertexDataStash_;

    // global solution of the current and of the last timestep
    const SpatialFunction *curSolution_;
    const SpatialFunction *oldSolution_;

private:
    void updateBoundaryTypes_()
    {
        Dune::GeometryType      geoType = curElement_().geometry().type();
        const ReferenceElement &refElem = DomainTraits::referenceElement(geoType);

        int numVerts = curElement_().template count<dim>();
        for (int i = 0; i < numVerts; ++i)
            for (int comp=0; comp < numEq; comp++)
                this->bctype[i][comp] = BoundaryConditions::neumann;

        // evaluate boundary conditions
        IntersectionIterator isIt = curElement_().ileafbegin();
        const IntersectionIterator &endIt = curElement_().ileafend();
        for (; isIt != endIt; ++isIt)
        {
            // Ignore non- boundary faces.
            if (!isIt->boundary())
                continue;

            // Set the bctype for all vertices of the face
            int faceIdx = isIt->indexInInside();
            int numFaceVerts = refElem.size(faceIdx, 1, dim);
            for (int faceVertIdx = 0;
                 faceVertIdx < numFaceVerts;
                 faceVertIdx++)
            {
                int elemVertIdx = refElem.subEntity(faceIdx,
                                                    1,
                                                    faceVertIdx,
                                                    dim);
                int boundaryFaceIdx =
                    curElementGeom_.boundaryFaceIndex(faceIdx,
                                                      faceVertIdx);

                // set the boundary types
                BoundaryTypeVector tmp;
                problem_.boundaryTypes(tmp,
                                       curElement_(),
                                       curFvElementGeometry(),
                                       isIt,
                                       elemVertIdx,
                                       boundaryFaceIdx);

                // copy boundary type to the bctype array.
                for (int eqIdx = 0; eqIdx < numEq; ++eqIdx) {
                    // make sure that dirichlet boundaries have
                    // priority over neumann ones
                    if (this->bctype[elemVertIdx][eqIdx]
                        == BoundaryConditions::dirichlet)
                    {
                        continue;
                    }

                    this->bctype[elemVertIdx][eqIdx] = tmp[eqIdx];
                    VALGRIND_CHECK_MEM_IS_DEFINED(&this->bctype[elemVertIdx][eqIdx], 
                                                  sizeof(Dune::BoundaryConditions));
                }
            }
        }
    };

    void assemble_(const Element &element, LocalFunction& localU, int orderOfShapeFns = 1)
    {
        // set the current grid element
        asImp_().setCurrentElement(element);

        // reset the right hand side and the bctype array
        resetRhs_();

        int numVertices = curElementGeom_.numVertices;

        // restrict the previous global solution to the current element
        LocalFunction localOldU(numVertices);
        restrictToElement(localOldU, *oldSolution_);

        this->asImp_().setCurrentSolution(localU);
        this->asImp_().setPreviousSolution(localOldU);

        // approximate the local stiffness matrix numerically
        // TODO: do this analytically if possible
        LocalFunction residUPlusEps(numVertices);
        LocalFunction residUMinusEps(numVertices);
        for (int j = 0; j < numVertices; j++)
        {
            for (int eqIdx = 0; eqIdx < numEq; eqIdx++)
            {
                Scalar eps = std::max(fabs(1e-5*localU[j][eqIdx]), 1e-5);
                Scalar uJ = localU[j][eqIdx];

                // vary the eqIdx-th equation at the element's j-th
                // vertex and calculate the residual, don't include
                // the boundary conditions
                asImp_().deflectCurrentSolution(localU, j, eqIdx, uJ + eps);
                evalLocalResidual(residUPlusEps, false);
                asImp_().restoreCurrentSolution(localU, j, eqIdx, uJ);

                asImp_().deflectCurrentSolution(localU, j, eqIdx, uJ - eps);
                evalLocalResidual(residUMinusEps, false);
                asImp_().restoreCurrentSolution(localU, j, eqIdx, uJ);

 
                // calculate the gradient when varying the
                // eqIdx-th primary variable at the j-th vert
                // of the element and fill the required fields of
                // the LocalStiffness base class.
                residUPlusEps -= residUMinusEps;

                residUPlusEps /= 2*eps;
                updateLocalStiffness_(j,
                                      eqIdx,
                                      residUPlusEps);
            }
        }

        // calculate the right hand side
        LocalFunction residU(numVertices);
        evalLocalResidual(residU, true); // include boundary conditions for the residual

        for (int i=0; i < numVertices; i++) {
            for (int eqIdx=0; eqIdx < numEq; eqIdx++) {
                // TODO: in most cases this is not really a boundary
                // vertex but an interior vertex, so
                // Dune::BoundaryConditions::neumann is misleading...
                if (this->bctype[i][eqIdx] == Dune::BoundaryConditions::neumann) {
                    this->b[i][eqIdx] = residU[i][eqIdx];
                }
            }
        }
    };

    void updateLocalStiffness_(int j,
                               int varJIdx,
                               const LocalFunction &stiffness)
    {
        for (int i = 0; i < curElementGeom_.numVertices; i++) {
            for (int varIIdx = 0; varIIdx < numEq; varIIdx++) {
                // A[i][j][varIIdx][eqIdx] is the
                // approximate rate of change of the
                // unknown 'varIIdx' at vertex 'i' if the
                // unknown 'varJIdx' at vertex 'j' is
                // slightly deflected from the current
                // local solution
                this->A[i][j][varIIdx][varJIdx] = stiffness[i][varIIdx];
            }
        }
    }

    // reset the right hand side
    void resetRhs_()
    {
        int numVertices = this->currentsize();
        for (int i=0; i < numVertices; i++) {
            this->bctype[i].assign(Dune::BoundaryConditions::neumann);
            this->b[i] = 0;
        }
    };

    JacobianImp &asImp_()
    { return *static_cast<JacobianImp*>(this); }

    const JacobianImp &asImp_() const
    { return *static_cast<const JacobianImp*>(this); }
};
}

#endif
