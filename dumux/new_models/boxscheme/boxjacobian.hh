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
    template<class Problem, class BoxTraitsT, class JacobianImp>
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
        typedef BoxJacobian<Problem, BoxTraits, JacobianImp> ThisType;

        typedef typename Problem::DomainTraits::Grid::LeafGridView    LeafGridView;
        typedef typename DomainTraits::Grid                           Grid;
        typedef typename DomainTraits::Scalar                         Scalar;
        typedef typename DomainTraits::CoordScalar                    CoordScalar;
        typedef typename DomainTraits::LocalPosition                  LocalPosition;
        typedef typename DomainTraits::GlobalPosition                     GlobalPosition;

        typedef typename LeafGridView::template Codim<0>::Entity      Element;
        typedef typename Element::EntityPointer                       ElementPointer;
        typedef typename DomainTraits::ReferenceElement               ReferenceElement;
        typedef typename DomainTraits::IntersectionIterator           IntersectionIterator;

        typedef typename Element::Geometry                        Geometry;

        typedef typename BoxTraits::SpatialFunction            SpatialFunction;
        typedef typename BoxTraits::SolutionVector             SolutionVector;
        typedef typename BoxTraits::BoundaryTypeVector         BoundaryTypeVector;
        typedef typename BoxTraits::FVElementGeometry          FVElementGeometry;
        typedef typename BoxTraits::LocalFunction              LocalFunction;

        typedef typename BoxTraits::ShapeFunctionSetContainer  ShapeFunctionSetContainer;
        typedef typename ShapeFunctionSetContainer::value_type ShapeFunctionSet;

    public:
        BoxJacobian(Problem &problem)
            : problem_(problem),
              curElementPtr_(problem.elementEnd()),

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
                asImp_()->setCurrentElement(element);

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
                asImp_()->setCurrentElement(element);

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
                asImp_()->setCurrentElement(element);

                // reset the right hand side and the boundary
                // condition type vector
                resetRhs_();

                Dune::GeometryType          geoType = curElement_().geometry().type();
                const ReferenceElement &refElem = DomainTraits::referenceElement(geoType);

                // temporary vector to store the neumann boundaries
                SolutionVector fluxes(0.0);

                // evaluate boundary conditions
                IntersectionIterator isIt = curElement_().ileafbegin();
                const IntersectionIterator &endIt = curElement_().ileafend();
                for (; isIt != endIt; ++isIt)
                {
                    // handle only faces on exterior boundaries. This
                    // assumes there are no interior boundaries.
                    if (!isIt->boundary())
                        continue;

                    // Assemble the boundary for all verts of the
                    // current face
                    int faceIdx = isIt->numberInSelf();
                    int numVerticesOfFace = refElem.size(faceIdx, 1, dim);
                    for (int vertInFace = 0;
                         vertInFace < numVerticesOfFace;
                         vertInFace++)
                    {
                        int vertInElement = refElem.subEntity(faceIdx,
                                                              1,
                                                              vertInFace,
                                                              dim);
                        int bfIdx = curElementGeom_.boundaryFaceIndex(faceIdx, vertInFace);
                        const LocalPosition &local = curElementGeom_.boundaryFace[bfIdx].ipLocal;
                        const GlobalPosition &global = curElementGeom_.boundaryFace[bfIdx].ipGlobal;

                        // set the boundary types
                        // TODO: better parameters: element, FVElementGeometry, bfIdx
                        BoundaryTypeVector tmp;
                        problem_.boundaryTypes(tmp,
                                               curElement_(),
                                               isIt,
                                               global,
                                               local);

                        // handle boundary conditions
                        bool neumannEvaluated = false;
                        for (int bcIdx = 0; bcIdx < numEq; ++bcIdx) {
                            // set the bctype of the LocalStiffness
                            // base class
                            this->bctype[vertInElement][bcIdx] = tmp[bcIdx];

                            // neumann boundaries
                            if (tmp[bcIdx] == BoundaryConditions::neumann) {
                                if (!neumannEvaluated) {
                                    // make sure that we only evaluate
                                    // the neumann fluxes once
                                    neumannEvaluated = true;
                                    // TODO: better Parameters: element, FVElementGeometry, bfIdx
                                    problem_.neumann(fluxes,
                                                     curElement_(),
                                                     isIt,
                                                     global,
                                                     local);
                                    fluxes *= curElementGeom_.boundaryFace[bfIdx].area;
                                }
                                this->b[vertInElement][bcIdx] += fluxes[bcIdx];
                            }
                            // dirichlet and processor boundaries
                            else {
                                // set right hand side to 0
                                this->b[vertInElement][bcIdx] = Scalar(0);
                            }
                        }
                    }
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
                    this->asImp_()->computeStorage(massContrib, i, false);
                    this->asImp_()->computeStorage(tmp, i, true);

                    massContrib -= tmp;
                    massContrib *= curElementGeom_.subContVol[i].volume/problem_.timeStepSize();
                    residual[i] += massContrib;

                    // subtract the source term from the local rate
                    SolutionVector source;
                    this->asImp_()->computeSource(source, i);
                    source *= curElementGeom_.subContVol[i].volume;
                    residual[i] -= source;
                }

                // calculate the mass flux over the faces and subtract
                // it from the local rates
                for (int k = 0; k < curElementGeom_.numEdges; k++)
                {
                    int i = curElementGeom_.subContVolFace[k].i;
                    int j = curElementGeom_.subContVolFace[k].j;

                    SolutionVector flux;
                    this->asImp_()->computeFlux(flux, k);

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
            }


        /*!
         * \brief Restrict the global function 'globalFn' to the verts
         *        of the current element, save the result to 'dest'.
         */
        void restrictToElement(LocalFunction &dest,
                            const SpatialFunction &globalFn) const
            {
                // we assert that the i-th shape function is
                // associated to the i-th vert of the element.
                int n = curElement_().template count<dim>();
                for (int i = 0; i < n; i++) {
                    dest[i] = (*globalFn)[problem_.vertIdx(curElement_(), i)];
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

        // clear the visited flag of all verts, required to make the
        // old-style models work in conjunction with the new newton
        // method.  HACK: remove
        void clearVisited() DUNE_DEPRECATED {};


    protected:

        // set the element which is currently considered as the local
        // stiffness matrix
        bool setCurrentElement_(const Element &e)
        {
            if (curElementPtr_ != e) {
                curElementPtr_ = e;
                curElementGeom_.update(curElement_());
                // tell the LocalStiffness (-> parent class) the current
                // number of degrees of freedom
                setcurrentsize(curElementGeom_.numVertices);
                return true;
            }
            return false;
        }

        const Element &curElement_() const
            { return *curElementPtr_; }

        // The problem we would like to solve
        Problem &problem_;

        ElementPointer       curElementPtr_;
        FVElementGeometry curElementGeom_;

        // global solution of the current and of the last timestep
        const SpatialFunction *curSolution_;
        const SpatialFunction *oldSolution_;

    private:
        void assemble_(const Element &element, LocalFunction& localU, int orderOfShapeFns = 1)
            {
                // set the current grid element
                asImp_()->setCurrentElement(element);

                // reset the right hand side and the bctype array
                resetRhs_();

                int numVertices = curElementGeom_.numVertices;

                // restrict the previous global solution to the current element
                LocalFunction localOldU(numVertices);
                restrictToElement(localOldU, *oldSolution_);

                this->asImp_()->setParams(element, localU, localOldU);

                // approximate the local stiffness matrix numerically
                // TODO: do this analytically if possible
                LocalFunction residUPlusEps(numVertices);
                LocalFunction residUMinusEps(numVertices);
                for (int j = 0; j < numVertices; j++)
                {
                    for (int comp = 0; comp < numEq; comp++)
                    {
                        Scalar eps = std::max(fabs(1e-5*localU[j][comp]), 1e-5);
                        Scalar uJ = localU[j][comp];

                        // vary the comp-th component at the element's j-th vert and
                        // calculate the residual, don't include the boundary
                        // conditions
                        this->asImp_()->deflectCurSolution(j, comp, uJ + eps);
                        evalLocalResidual(residUPlusEps, false);

                        this->asImp_()->deflectCurSolution(j, comp, uJ - eps);
                        evalLocalResidual(residUMinusEps, false);

                        // restore the current local solution to the state before
                        // varyCurSolution() has been called
                        this->asImp_()->restoreCurSolution(j, comp);

                        // calculate the gradient when varying the
                        // comp-th primary variable at the j-th vert
                        // of the element and fill the required fields of
                        // the LocalStiffness base class.
                        residUPlusEps -= residUMinusEps;

                        residUPlusEps /= 2*eps;
                        updateLocalStiffness_(j,
                                              comp,
                                              residUPlusEps);
                    }
                }

                // calculate the right hand side
                LocalFunction residU(numVertices);
                evalLocalResidual(residU, true); // include boundary conditions for the residual

                for (int i=0; i < numVertices; i++) {
                    for (int comp=0; comp < numEq; comp++) {
                        // TODO: in most cases this is not really a
                        // boundary element but an interior element, so
                        // Dune::BoundaryConditions::neumann is
                        // misleading...
                        if (this->bctype[i][comp] == Dune::BoundaryConditions::neumann) {
                            this->b[i][comp] = residU[i][comp];
                        }
                    }
                }
            };

        void updateLocalStiffness_(int j,
                                   int comp,
                                   const LocalFunction &stiffness)
            {
                for (int i = 0; i < curElementGeom_.numVertices; i++) {
                    for (int compi = 0; compi < numEq; compi++) {
                        // A[i][j][compi][comp] is the
                        // approximate rate of change of the
                        // unknown 'compi' at vert 'i' if the
                        // component 'comp' at vert 'j' is
                        // slightly deflected from the current
                        // local solution
                        this->A[i][j][compi][comp] = stiffness[i][compi];
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

    private:
        JacobianImp *asImp_()
            { return static_cast<JacobianImp*>(this); }

        const JacobianImp *asImp_() const
            { return static_cast<const JacobianImp*>(this); }
    };
}

#endif
