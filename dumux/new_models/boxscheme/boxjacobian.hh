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
        typedef typename DomainTraits::CellReferenceElements          CellReferenceElements;
        typedef typename DomainTraits::CellReferenceElement           CellReferenceElement;
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
            : _problem(problem),
              _curCellPtr(problem.cellEnd()),

              _curSolution(NULL),
              _oldSolution(NULL)
            {
            }

        /*!
         * \brief Set the global solution of the current time step.
         */
        void setCurSolution(SpatialFunction *uCur)
            {
                _curSolution = uCur;
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
                _oldSolution = uOld;
            }

        /*!
         * \brief Assemble the linear system of equations for the
         *        nodes of a cell, given a local solution 'localU'.
         */
        void assemble(const Cell &cell, const LocalFunction& localU, int orderOfShapeFns = 1)
            {
                // set the current grid cell
                _asImp()->setCurrentCell(cell);

                LocalFunction mutableLocalU(localU);
                _assemble(cell, mutableLocalU, orderOfShapeFns);
            }

        /*!
         * \brief Assemble the linear system of equations for the
         *        nodes of a cell.
         */
        void assemble(const Cell &cell, int orderOfShapeFns = 1)
            {
                // set the current grid cell
                _asImp()->setCurrentCell(cell);

                int numVertices = _curCellGeom.nNodes;
                LocalFunction localU(numVertices);
                restrictToCell(localU, *_curSolution); 
                _assemble(cell, localU, orderOfShapeFns);
            }
        
        /*!
         * \brief Express the boundary conditions for a cell in terms
         *        of a linear equation system.
         */
        void assembleBoundaryCondition(const Cell &cell, int orderOfShapeFns=1)
            {
                // set the current grid cell
                _asImp()->setCurrentCell(cell);

                // reset the right hand side and the boundary
                // condition type vector
                _resetRhs();
                
                Dune::GeometryType          geoType = _curCell().geometry().type();
                const CellReferenceElement &refElem = CellReferenceElements::general(geoType);
                
                // evaluate boundary conditions
                IntersectionIterator endIt = IntersectionIteratorGetter::end(_curCell());
                IntersectionIterator faceIt = IntersectionIteratorGetter::begin(_curCell());
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
                        int bfIdx = _curCellGeom.boundaryFaceIndex(faceIdx, nodeInFace);
                        const LocalCoord &local = _curCellGeom.boundaryFace[bfIdx].ipLocal;
                        const WorldCoord &global = _curCellGeom.boundaryFace[bfIdx].ipGlobal;
                        
                        // set the boundary types
                        // TODO: better parameters: cell, FVElementGeometry, bfIndex
                        BoundaryTypeVector tmp;
                        _problem.boundaryTypes(tmp,
                                               _curCell(),
                                               faceIt,
                                               global,
                                               local);
                        // TODO: mixed boundary conditions
                        this->bctype[nodeInElement].assign(tmp[0]);
                        
                        // handle neumann boundary conditions
                        if (tmp[0] == BoundaryConditions::neumann) {
                            // TODO: better Parameters: cell, FVElementGeometry, bfIndex
                            UnknownsVector J;
                            _problem.neumann(J,
                                             _curCell(),
                                             faceIt, 
                                             global,
                                             local);
                            J *= _curCellGeom.boundaryFace[bfIdx].area;

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
                for (int i = 0; i < _curCellGeom.nNodes; i++) {
                    residual[i] = 0;
                }

                // evaluate the local rate
                for (int i=0; i < _curCellGeom.nNodes; i++)
                {
                    UnknownsVector massContrib(0), tmp(0);

                    // mass balance within the cell. this is the
                    // $\frac{m}{\partial t}$ term if using implicit
                    // euler as time discretization. 
                    //
                    // TODO (?): we might need a more explicit way for
                    // doing the time discretization...
                    this->_asImp()->localRate(massContrib, i, false);
                    this->_asImp()->localRate(tmp, i, true);

                    massContrib -= tmp;
                    massContrib *= _curCellGeom.subContVol[i].volume/_problem.timeStepSize();
                    residual[i] += massContrib;
                    
                    // subtract the source term from the local rate
                    UnknownsVector q;
                    _problem.sourceTerm(q,
                                        _curCell(),
                                        _curCellGeom,
                                        i);
                    q *= _curCellGeom.subContVol[i].volume;
                    residual[i] -= q;
                }

                // calculate the mass flux over the faces and subtract
                // it from the local rates
                for (int k = 0; k < _curCellGeom.nEdges; k++)
                {
                    int i = _curCellGeom.subContVolFace[k].i;
                    int j = _curCellGeom.subContVolFace[k].j;

                    UnknownsVector flux;
                    this->_asImp()->fluxRate(flux, k);

                    // subtract fluxes from the local mass rates of
                    // the respective sub control volume adjacent to
                    // the face.
                    residual[i] -= flux;
                    residual[j] += flux;
                }
                
                if (withBoundary) {
                    assembleBoundaryCondition(this->_curCell());
                    for (int i = 0; i < _curCellGeom.nNodes; i++) {
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
                int n = _curCell().template count<GridDim>();
                for (int i = 0; i < n; i++) {
                    dest[i] = (*globalFn)[_problem.nodeIndex(_curCell(), i)];
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
        bool _setCurrentCell(const Cell &e)
        {
            if (_curCellPtr != e) {
                _curCellPtr = e;
                _curCellGeom.update(_curCell());
                // tell the LocalStiffness (-> parent class) the current
                // number of degrees of freedom
                setcurrentsize(_curCellGeom.nNodes);
                return true;
            }
            return false;
        }

        const Cell &_curCell() const
            { return *_curCellPtr; }

        // The problem we would like to solve
        Problem &_problem;

        CellPointer       _curCellPtr;
        FVElementGeometry _curCellGeom;

        // global solution of the current and of the last timestep
        const SpatialFunction *_curSolution;
        const SpatialFunction *_oldSolution;

    private:
        void _assemble(const Cell &cell, LocalFunction& localU, int orderOfShapeFns = 1)
            {
                // set the current grid cell
                _asImp()->setCurrentCell(cell);

                // reset the right hand side and the bctype array
                _resetRhs();

                int numVertices = _curCellGeom.nNodes;

                // restrict the previous global solution to the current cell
                LocalFunction localOldU(numVertices);
                restrictToCell(localOldU, *_oldSolution);

                this->_asImp()->setParams(cell, localU, localOldU);

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
                        this->_asImp()->deflectCurSolution(j, comp, uJ + eps);
                        evalLocalResidual(residUPlusEps, false);
         
                        this->_asImp()->deflectCurSolution(j, comp, uJ - eps);
                        evalLocalResidual(residUMinusEps, false);

                        // restore the current local solution to the state before
                        // varyCurSolution() has been called 
                        this->_asImp()->restoreCurSolution(j, comp);

                        // calculate the gradient when varying the
                        // comp-th primary variable at the j-th node
                        // of the cell and fill the required fields of
                        // the LocalStiffness base class.
                        residUPlusEps -= residUMinusEps;
                        residUPlusEps /= 2*eps;
                        _updateLocalStiffness(j, 
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

        void _updateLocalStiffness(int j, 
                                   int comp,
                                   const LocalFunction &stiffness)
            {
                for (int i = 0; i < _curCellGeom.nNodes; i++) {
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
        void _resetRhs()
            {
                int numVertices = this->currentsize();
                for (int i=0; i < numVertices; i++) {
                    this->bctype[i].assign(Dune::BoundaryConditions::neumann);
                    this->b[i] = 0;
                }
            };

    private:
        JacobianImp *_asImp()
            { return static_cast<JacobianImp*>(this); }

        const JacobianImp *_asImp() const
            { return static_cast<const JacobianImp*>(this); }
    };
}

#endif
