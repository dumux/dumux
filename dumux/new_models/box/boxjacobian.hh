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
        BoxJacobian(Problem &problem,
                    bool levelBoundaryAsDirichlet,
                    bool procBoundaryAsDirichlet=true)
            : _problem(problem),
              _curCellPtr(problem.cellEnd()),

              _curSolution(NULL),
              _oldSolution(NULL),

              _levelBoundaryAsDirichlet(levelBoundaryAsDirichlet),
              _processBoundaryAsDirichlet(procBoundaryAsDirichlet)
            {
            }

        /*!
         * \brief Set the global solution of the current time step.
         */
        void setCurrentSolution(SpatialFunction *uCur)
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

//        void clearVisited()
//            { }

        /*!
         * \brief Assemble the linear system of equations for the
         *        nodes of a cell, given a local solution 'localU'.
         */
        void assemble(const Cell &cell, const LocalFunction& localU, int orderOfShapeFns = 1)
            {
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
                
                if (!cell.hasBoundaryIntersections())
                    return;

                Dune::GeometryType      geoType = _curCell().geometry().type();
                const CellReferenceElement &refElem = CellReferenceElements::general(geoType);
                
                // evaluate boundary conditions via intersection iterator
                IntersectionIterator endIt = IntersectionIteratorGetter::end(_curCell());
                IntersectionIterator faceIt = IntersectionIteratorGetter::begin(_curCell());
                for (; faceIt != endIt; ++faceIt)
                {

                    // if we have a neighbor then we assume there is
                    // no boundary (forget interior boundaries) in
                    // level assemble treat non-level neighbors as
                    // boundary
                    if (faceIt.neighbor())
                    {
                        if (!_levelBoundaryAsDirichlet)
                            // we are not supposed to handle
                            // boundaries within a level.
                            continue;
                        else if (faceIt.outside()->level() == _curCell().level())
                            // we are supposed to handle boundaries
                            // within a level, but we are not at the
                            // level boundary of the cell
                            continue;
                    }


                    // determine boundary condition type for this
                    // face, initialize as process boundary
                    BoundaryTypeVector faceBCType(Dune::BoundaryConditions::process);

                    // handle face on exterior boundary, this assumes
                    // there are no interior boundaries
                    if (faceIt.boundary())
                    {
                        _assembleExteriorBC(faceBCType, faceIt, refElem);
                        
                        if (faceBCType[0]==Dune::BoundaryConditions::neumann)
                            continue; // was a neumann face, go to next face
                    }

                    // If we are here, then faceIt is
                    // (i)   an exterior boundary face with Dirichlet condition, or
                    // (ii)  a process boundary (i.e. neither boundary() nor neighbor() was true), or
                    // (iii) a level boundary in case of level-wise assemble
                    // How process boundaries are handled depends on the process boundary mode
                    if (faceBCType[0]==Dune::BoundaryConditions::process
                        && !_processBoundaryAsDirichlet
                        && !_levelBoundaryAsDirichlet)
                        continue; // then faceIt acts like homogeneous Neumann


                    // now handle exterior or interior Dirichlet boundary
                    _assembleDirichletBC(faceBCType,
                                         faceIt,
                                         refElem);
                }
            }

        /*!
         * \brief Compute the local residual, i.e. the right hand side
         *        of an equation we would like to have zero.
         */
        void evalLocalResidual(LocalFunction &residual) const
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
            }


        /*!
         * \brief Restrict the global function 'globalFn' to the nodes
         *        of the current cell, save the result to 'dest'.
         */
        void restrictToCell(LocalFunction &dest,
                            const SpatialFunction &globalFn)
            {
#if 1
                // we assert that the i-th shape function is
                // associated to the i-th node of the cell.
                int n = _curCell().template count<GridDim>();
                for (int i = 0; i < n; i++) {
                    dest[i] = (*globalFn)[_problem.nodeIndex(_curCell(), i)];
                }
#else
                // TODO: globalFn is probably not required to be
                // evaluated at the center position of a shape
                // function since all other shape functions are 0 at
                // this location. thus we just copy the values at the
                // nodes of the cell, but if higher order shape
                // functions are involved, they might not be centered
                // at the cell's nodes, so this code might actually be
                // necessarry...
                Dune::GeometryType geoType = _curCell().geometry().type();
                const ShapeFunctionSet &shapeFns = BoxTraits::shapeFunctions()(geoType, 1);
                //int size = ;
                //dest.resize(size);

                for (int i = 0; i < shapeFns.size(); i++) {
                    for (int comp = 0; comp < PrimaryVariables; comp++) {
                        dest[i][comp] = globalFn.evallocal(comp,
                                                           _curCell(),
                                                           shapeFns[i].position());
                    }
                }
#endif
            }


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

        bool         _levelBoundaryAsDirichlet;
        bool         _processBoundaryAsDirichlet;

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

                LocalFunction residUPlusEps(numVertices);
                LocalFunction residUMinusEps(numVertices);
                for (int j = 0; j < numVertices; j++)
                {
                    for (int comp = 0; comp < PrimaryVariables; comp++)
                    {
                        Scalar eps = std::max(fabs(1e-3 * localU[j][comp]), 1e-3);
                        Scalar uJ = localU[j][comp];
                        
                        // vary the comp-th component at the cell's j-th node and
                        // calculate the residual.
                        this->_asImp()->deflectCurSolution(j, comp, uJ + eps);
                        evalLocalResidual(residUPlusEps);
                        
                        this->_asImp()->deflectCurSolution(j, comp, uJ - eps);
                        evalLocalResidual(residUMinusEps);

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

                // assemble boundary conditions
                assembleBoundaryCondition(cell);
                
                // calculate the right hand side
                LocalFunction residU(numVertices);
                evalLocalResidual(residU);
                for (int i=0; i < numVertices; i++) {
                    for (int comp=0; comp < PrimaryVariables; comp++) {
                        if (this->bctype[i][comp]==Dune::BoundaryConditions::neumann) {
                            this->b[i][comp] += residU[i][comp];
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


       void _assembleExteriorBC(BoundaryTypeVector &faceBCType,
                                const IntersectionIterator &faceIt,
                                const CellReferenceElement &refElem)
            {
                int faceIdx = faceIt.numberInSelf();
                int nVerticesOfFace = refElem.size(faceIdx, 1, GridDim);
                for (int nodeInFace = 0; nodeInFace < nVerticesOfFace; nodeInFace++)
                {
                    int nodeInElement = refElem.subEntity(faceIdx, 1, nodeInFace, GridDim);
                    
                    // TODO: "mixed" boundary conditions are not yet
                    // possible. either it's neumann for all primary
                    // variables or it's dirchlet.
                    if (this->bctype[nodeInElement][0] == Dune::BoundaryConditions::neumann) 
                    {
                        int bfIdx = _curCellGeom.boundaryFaceIndex(faceIdx, nodeInFace);
                        const LocalCoord &local = _curCellGeom.boundaryFace[bfIdx].ipLocal;
                        WorldCoord global = _curCellGeom.boundaryFace[bfIdx].ipGlobal;
                        
                        _problem.boundaryTypes(faceBCType,
                                               _curCell(),
                                               faceIt,
                                               global,
                                               local);
                        
                        
                        if (faceBCType[0]!=Dune::BoundaryConditions::neumann)
                            break;
                        
                        UnknownsVector J;
                        _problem.neumann(J,
                                         _curCell(),
                                         faceIt,
                                         global,
                                         local);
                        
                        J *= _curCellGeom.boundaryFace[bfIdx].area;
                        this->b[nodeInElement] += J;
                    }
                }
            }
        

        // assemble a dirichlet boundary condition for a face which
        // may be either exterior or within the grid (if level or
        // process boundaries are required/used)
        void _assembleDirichletBC(const BoundaryTypeVector &faceBCType,
                                  const IntersectionIterator &faceIt,
                                  const CellReferenceElement &refElem)
            {
#if 1
                int nFaceVerts = refElem.size(faceIt.numberInSelf(),
                                              1,
                                              GridDim);
                // loop over all vertices of the face
                for (int i=0; i < nFaceVerts; ++ i) {
                    // find the local node index of the cell using
                    // the local node index of the face
                    int j = refElem.subEntity(faceIt.numberInSelf(),
                                              1,
                                              i,
                                              GridDim);
                    
                    // the associated node of the j-th shape
                    // function is a node of the current
                    // face.
                    this->bctype[j] = faceBCType[0];
                    this->b[j] = 0;
                }
#else
                Dune::GeometryType      geoType = _curCell().geometry().type();
                const ShapeFunctionSet &shapeFnSet = BoxTraits::shapeFunctions()(geoType, 1);

                // loop over shape functions
                for (int i=0; i<shapeFnSet.size(); i++)
                {
                    if (shapeFnSet[i].codim() != GridDim)
                        continue; // skip interior and face degrees of freedom
                    
                    // loop over all vertices of the face in order to
                    // make sure the entity associated with the
                    // current shape function is also a node of the
                    // current face.
                    //
                    // TODO: there must be a better way to do this??
                    int nFaceVerts = refElem.size(faceIt.numberInSelf(),
                                                  1,
                                                  shapeFnSet[i].codim());
                    for (int j=0; j < nFaceVerts; j++)
                    {
                        if (shapeFnSet[i].entity() ==
                            refElem.subEntity(faceIt.numberInSelf(),
                                              1,
                                              j,
                                              shapeFnSet[i].codim()))
                        {
                            // the associated node of the i-th shape
                            // function is a node of the current
                            // face.
                            this->bctype[i] = faceBCType[0];
                            this->b[i] = 0;
                            break;
                        }
                    }
                }
#endif
            }

    private:
        JacobianImp *_asImp()
            { return static_cast<JacobianImp*>(this); }

        const JacobianImp *_asImp() const
            { return static_cast<const JacobianImp*>(this); }
    };
}

#endif
