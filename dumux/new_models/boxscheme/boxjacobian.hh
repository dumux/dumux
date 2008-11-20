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
                                                    BoxTraitsT::PrimaryVariables>
    {
    private:
        typedef typename Problem::DomainTraits DomainTraits;
        typedef BoxTraitsT                     BoxTraits;

        enum {
            PrimaryVariables = BoxTraits::PrimaryVariables,

            GridDim     = DomainTraits::GridDim,
            WorldDim    = DomainTraits::WorldDim
        };
        typedef BoxJacobian<Problem, BoxTraits, JacobianImp> ThisType;

        typedef typename Problem::DomainTraits::Grid::LeafGridView    LeafGridView;
        typedef typename DomainTraits::Grid                           Grid;
        typedef typename DomainTraits::Scalar                         Scalar;
        typedef typename DomainTraits::CoordScalar                    CoordScalar;
        typedef typename DomainTraits::LocalCoord                     LocalCoord;
        typedef typename DomainTraits::WorldCoord                     WorldCoord;

        typedef typename LeafGridView::template Codim<0>::Entity      Cell;
        typedef typename Cell::EntityPointer                          CellPointer;
        typedef typename DomainTraits::ReferenceElement               ReferenceElement;
        typedef typename DomainTraits::IntersectionIterator           IntersectionIterator;
        typedef typename DomainTraits::IntersectionIteratorGetter     IntersectionIteratorGetter;

        typedef typename Cell::Geometry                        Geometry;

        typedef typename BoxTraits::SpatialFunction            SpatialFunction;
        typedef typename BoxTraits::UnknownsVector             UnknownsVector;
        typedef typename BoxTraits::BoundaryTypeVector         BoundaryTypeVector;
        typedef typename BoxTraits::FVElementGeometry          FVElementGeometry;
        typedef typename BoxTraits::LocalFunction              LocalFunction;

        typedef typename BoxTraits::ShapeFunctionSetContainer  ShapeFunctionSetContainer;
        typedef typename ShapeFunctionSetContainer::value_type ShapeFunctionSet;

    public:
        BoxJacobian(Problem &problem)
            : problem_(problem),
              curCellPtr_(problem.cellEnd()),

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
         *        nodes of a cell, given a local solution 'localU'.
         */
        void assemble(const Cell &cell, const LocalFunction& localU, int orderOfShapeFns = 1)
            {
                // set the current grid cell
                asImp_()->setCurrentCell(cell);

                LocalFunction mutableLocalU(localU);
                assemble_(cell, mutableLocalU, orderOfShapeFns);
            }

        /*!
         * \brief Assemble the linear system of equations for the
         *        nodes of a cell.
         */
        void assemble(const Cell &cell, int orderOfShapeFns = 1)
            {
                // set the current grid cell
                asImp_()->setCurrentCell(cell);

                int numVertices = curCellGeom_.nNodes;
                LocalFunction localU(numVertices);
                restrictToCell(localU, *curSolution_); 
                assemble_(cell, localU, orderOfShapeFns);
            }
        
        /*!
         * \brief Express the boundary conditions for a cell in terms
         *        of a linear equation system.
         */
        void assembleBoundaryCondition(const Cell &cell, int orderOfShapeFns=1)
            {
                // set the current grid cell
                asImp_()->setCurrentCell(cell);

                // reset the right hand side and the boundary
                // condition type vector
                resetRhs_();
                
                Dune::GeometryType          geoType = curCell_().geometry().type();
                const ReferenceElement &refElem = DomainTraits::referenceElement(geoType);
                
                // evaluate boundary conditions
                IntersectionIterator endIt = IntersectionIteratorGetter::end(curCell_());
                IntersectionIterator faceIt = IntersectionIteratorGetter::begin(curCell_());
                for (; faceIt != endIt; ++faceIt)
                {
                    // handle only faces on exterior boundaries. This
                    // assumes there are no interior boundaries.
                    if (!faceIt.boundary())
                        continue;
                    
                    // Assemble the boundary for all nodes of the
                    // current face
                    int faceIdx = faceIt.numberInSelf();
                    int nNodesOfFace = refElem.size(faceIdx, 1, GridDim);
                    for (int nodeInFace = 0; 
                         nodeInFace < nNodesOfFace;
                         nodeInFace++)
                    {
                        int nodeInElement = refElem.subEntity(faceIdx,
                                                              1, 
                                                              nodeInFace, 
                                                              GridDim);
                        int bfIdx = curCellGeom_.boundaryFaceIndex(faceIdx, nodeInFace);
                        const LocalCoord &local = curCellGeom_.boundaryFace[bfIdx].ipLocal;
                        const WorldCoord &global = curCellGeom_.boundaryFace[bfIdx].ipGlobal;
                        
                        // set the boundary types
                        // TODO: better parameters: cell, FVElementGeometry, bfIndex
                        BoundaryTypeVector tmp;
                        problem_.boundaryTypes(tmp,
                                               curCell_(),
                                               faceIt,
                                               global,
                                               local);
                        // TODO: mixed boundary conditions
                        this->bctype[nodeInElement].assign(tmp[0]);
                        
                        // handle neumann boundary conditions
                        if (tmp[0] == BoundaryConditions::neumann) {
                            // TODO: better Parameters: cell, FVElementGeometry, bfIndex
                            UnknownsVector J;
                            problem_.neumann(J,
                                             curCell_(),
                                             faceIt, 
                                             global,
                                             local);
                            J *= curCellGeom_.boundaryFace[bfIdx].area;

                            this->b[nodeInElement] += J;
                        }
                        // handle dirichlet and process boundaries
                        else if (tmp[0] == BoundaryConditions::dirichlet || 
                                 tmp[0] == BoundaryConditions::process)
                        {
                            // right hand side
                            this->b[nodeInElement] = Scalar(0);
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
                for (int i = 0; i < curCellGeom_.nNodes; i++) {
                    residual[i] = 0;
                }

                // evaluate the local rate
                for (int i=0; i < curCellGeom_.nNodes; i++)
                {
                    UnknownsVector massContrib(0), tmp(0);

                    // mass balance within the cell. this is the
                    // $\frac{m}{\partial t}$ term if using implicit
                    // euler as time discretization. 
                    //
                    // TODO (?): we might need a more explicit way for
                    // doing the time discretization...
                    this->asImp_()->localRate(massContrib, i, false);
                    this->asImp_()->localRate(tmp, i, true);

                    massContrib -= tmp;
                    massContrib *= curCellGeom_.subContVol[i].volume/problem_.timeStepSize();
                    residual[i] += massContrib;
                    
                    // subtract the source term from the local rate
                    UnknownsVector q;
                    problem_.sourceTerm(q,
                                        curCell_(),
                                        curCellGeom_,
                                        i);
                    q *= curCellGeom_.subContVol[i].volume;
                    residual[i] -= q;
                }

                // calculate the mass flux over the faces and subtract
                // it from the local rates
                for (int k = 0; k < curCellGeom_.nEdges; k++)
                {
                    int i = curCellGeom_.subContVolFace[k].i;
                    int j = curCellGeom_.subContVolFace[k].j;

                    UnknownsVector flux;
                    this->asImp_()->fluxRate(flux, k);

                    // subtract fluxes from the local mass rates of
                    // the respective sub control volume adjacent to
                    // the face.
                    residual[i] -= flux;
                    residual[j] += flux;
                }
                
                if (withBoundary) {
                    assembleBoundaryCondition(this->curCell_());
                    for (int i = 0; i < curCellGeom_.nNodes; i++) {
                        residual[i] += this->b[i];
                    }
                }
            }


        /*!
         * \brief Restrict the global function 'globalFn' to the nodes
         *        of the current cell, save the result to 'dest'.
         */
        void restrictToCell(LocalFunction &dest,
                            const SpatialFunction &globalFn) const
            {
                // we assert that the i-th shape function is
                // associated to the i-th node of the cell.
                int n = curCell_().template count<GridDim>();
                for (int i = 0; i < n; i++) {
                    dest[i] = (*globalFn)[problem_.nodeIndex(curCell_(), i)];
                }
            }

        /*!
         * \brief Initialize the static data of all cells. The current
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
         * \brief Update the static data of all cells with the current solution
         *
         * This method should be overwritten by the child class if
         * necessary.
         */
        void updateStaticData(SpatialFunction &curSol, SpatialFunction &oldSol)
            { };

        // clear the visited flag of all nodes, required to make the
        // old-style models work in conjunction with the new newton
        // method.  HACK: remove
        void clearVisited() DUNE_DEPRECATED {};


    protected:
                  
        // set the cell which is currently considered as the local
        // stiffness matrix
        bool setCurrentCell_(const Cell &e)
        {
            if (curCellPtr_ != e) {
                curCellPtr_ = e;
                curCellGeom_.update(curCell_());
                // tell the LocalStiffness (-> parent class) the current
                // number of degrees of freedom
                setcurrentsize(curCellGeom_.nNodes);
                return true;
            }
            return false;
        }

        const Cell &curCell_() const
            { return *curCellPtr_; }

        // The problem we would like to solve
        Problem &problem_;

        CellPointer       curCellPtr_;
        FVElementGeometry curCellGeom_;

        // global solution of the current and of the last timestep
        const SpatialFunction *curSolution_;
        const SpatialFunction *oldSolution_;

    private:
        void assemble_(const Cell &cell, LocalFunction& localU, int orderOfShapeFns = 1)
            {
                // set the current grid cell
                asImp_()->setCurrentCell(cell);

                // reset the right hand side and the bctype array
                resetRhs_();

                int numVertices = curCellGeom_.nNodes;

                // restrict the previous global solution to the current cell
                LocalFunction localOldU(numVertices);
                restrictToCell(localOldU, *oldSolution_);

                this->asImp_()->setParams(cell, localU, localOldU);

                // approximate the local stiffness matrix numerically
                // TODO: do this analytically if possible
                LocalFunction residUPlusEps(numVertices);
                LocalFunction residUMinusEps(numVertices);
                for (int j = 0; j < numVertices; j++)
                {
                    for (int comp = 0; comp < PrimaryVariables; comp++)
                    {
                        Scalar eps = std::max(fabs(1e-5*localU[j][comp]), 1e-5);
                        Scalar uJ = localU[j][comp];
                        
                        // vary the comp-th component at the cell's j-th node and
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
                        // comp-th primary variable at the j-th node
                        // of the cell and fill the required fields of
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
                    for (int comp=0; comp < PrimaryVariables; comp++) {
                        // TODO: in most cases this is not really a
                        // boundary cell but an interior cell, so
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
                for (int i = 0; i < curCellGeom_.nNodes; i++) {
                    for (int compi = 0; compi < PrimaryVariables; compi++) {
                        // A[i][j][compi][comp] is the
                        // approximate rate of change of the
                        // unknown 'compi' at node 'i' if the
                        // component 'comp' at node 'j' is
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
