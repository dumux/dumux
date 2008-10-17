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

      - _Problem        The specification of the problem (grid, boundary conditions, ...)
      - _JacobianImp    The Box model specific part of the Jacobian
    */
    template<class Problem, class BoxTraitsT, class JacobianImp>
    class BoxJacobian : public Dune::LocalStiffness<typename Problem::DomainTraits::Grid::LeafGridView,
                                                    typename Problem::DomainTraits::Scalar,
                                                    BoxTraitsT::NumUnknowns>
    {
    private:
        typedef typename Problem::DomainTraits DomainTraits;
        typedef BoxTraitsT                     BoxTraits;

        enum {
            NumUnknowns = BoxTraits::NumUnknowns,

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
        typedef typename DomainTraits::CellReferenceElements          CellReferenceElements;
        typedef typename DomainTraits::CellReferenceElement           CellReferenceElement;
        typedef typename DomainTraits::IntersectionIterator           IntersectionIterator;
        typedef typename DomainTraits::IntersectionIteratorGetter     IntersectionIteratorGetter;

        typedef typename Cell::Geometry                  Geometry;

        typedef typename BoxTraits::SpatialFunction      SpatialFunction;
        typedef typename BoxTraits::UnknownsVector       UnknownsVector;
        typedef typename BoxTraits::BoundaryTypeVector   BoundaryTypeVector;
        typedef typename BoxTraits::FVElementGeometry    FVElementGeometry;
        typedef typename BoxTraits::ShapeFnSets          ShapeFnSets;
        typedef typename BoxTraits::ShapeFnSet           ShapeFnSet;
        typedef typename BoxTraits::LocalFunction        LocalFunction;
        typedef typename BoxTraits::CachedCellData       CachedCellData;
        typedef typename BoxTraits::CachedSubContVolData CachedSubContVolData;

    public:
        BoxJacobian(Problem &problem,
                    bool levelBoundaryAsDirichlet,
                    bool procBoundaryAsDirichlet=true)
            : _problem(problem),

              _curCellPtr(NULL),

              _curSolution(NULL),
              _oldSolution(NULL),

              _levelBoundaryAsDirichlet(levelBoundaryAsDirichlet),
              _processBoundaryAsDirichlet(procBoundaryAsDirichlet),

              _dt(1)
            {
            }

        // TODO/HACK: keeping track of the time step length is already
        // job of the problem, so we should actually call it back
        void setDt (double d)
            { _dt = d; }
        double getDt () const
            { return _dt; }
        // END TODO/HACK

        void setCurrentSolution(SpatialFunction *uCur)
            {
                _curSolution = uCur;
            }

        void setOldSolution(SpatialFunction *uOld)
            {
                _oldSolution = uOld;
            }

        void clearVisited()
            { }

        // TODO/FIXME: this is only valid for linear problems where
        // the local stiffness matrix is independend of the current
        // solution. We need to implement this properly, but this
        // should at least make the thing compile...
        typedef Dune::FieldVector<Scalar, BoxTraitsT::NumUnknowns> VBlockType;
        void assemble(const Cell &cell, const Dune::BlockVector<VBlockType>& localSolution, int orderOfShapeFns = 1)
            {
                assemble(cell, orderOfShapeFns);
            }

        void assemble(const Cell &cell, int orderOfShapeFns = 1)
            {
                // we assert that TypeTag == Dune::LeafTag, else this
                // doesn't make any sense.

                // set the current cell to the one we are supposed
                // to assemble the local stiffnes matrix.
                setCurrentCell(cell);

                // reset the right hand side and the bctype array
                _resetRhs();

                LocalFunction localU, localOldU;
                CachedCellData uCache;
                evalLocal(localU, *_curSolution);
                evalLocal(localOldU, *_oldSolution);

                JacobianImp::updateCellCache(uCache,
                                              _problem,
                                              _curCell(),
                                              _curCellGeom,
                                              localU);
                
                LocalFunction  uPlusEps,      uMinusEps;
                CachedCellData uPlusMinusEpsCache;
                LocalFunction defUPlusEps, defUMinusEps;
                uPlusMinusEpsCache = uCache;
                int numVertices = _curCellGeom.nNodes;
                for (int j = 0; j < numVertices; j++)
                {
                    for (int comp = 0; comp < NumUnknowns; comp++)
                    {
                        Scalar eps = std::max(fabs(1e-3*localU.atSubContVol[j][comp]), 1e-3);

                        uPlusEps = localU;
                        uMinusEps = localU;
                        uPlusEps.atSubContVol[j][comp]  += eps;
                        uMinusEps.atSubContVol[j][comp] -= eps;

                        // copy uCache into uPlusMinusEpsCache, but only evaluate
                        // component j fully using the local solution uPlusEps
                        JacobianImp::partialUpdateCellCache(uPlusMinusEpsCache,
                                                             uCache,
                                                             _problem,
                                                             _curCell(),
                                                             _curCellGeom,
                                                             uPlusEps,
                                                             j);
                        evalLocalDefect(defUPlusEps,
                                        uPlusEps, uPlusMinusEpsCache,
                                        localOldU);


                        // copy uCache into uPlusMinusEpsCache, but only evaluate
                        // component j fully using the local solution uMinusEps
                        JacobianImp::partialUpdateCellCache(uPlusMinusEpsCache,
                                                             uCache,
                                                             _problem,
                                                             _curCell(),
                                                             _curCellGeom,
                                                             uMinusEps,
                                                             j);
                        evalLocalDefect(defUMinusEps,
                                        uMinusEps, uPlusMinusEpsCache,
                                        localOldU);
                        
                        
/*                        if (HACKY_HACK) {
                            printf("varying parameter %d at vertex %d of cell %d\n",
                                   comp, 
                                   j,
                                   cellIdx);
                            printf("node 14: pw:%f/sn:%f vs pw:%f/sn:%f ///// node 1:pw:%f/sn%f vs pw:%f/sn:%f\n", 
                                   defUPlusEps.atSubContVol[0][0], 
                                   defUMinusEps.atSubContVol[0][1],
                                   
                                   defUPlusEps.atSubContVol[0][0], 
                                   defUMinusEps.atSubContVol[0][1],
                                   
                                   defUPlusEps.atSubContVol[1][0], 
                                   defUMinusEps.atSubContVol[1][1],
                                   
                                   defUPlusEps.atSubContVol[1][0], 
                                   defUMinusEps.atSubContVol[1][1]);
                        }
*/


                        Scalar deltaX = 2*eps;
                        for (int i = 0; i < numVertices; i++) {
                            for (int compi = 0; compi < NumUnknowns; compi++) {
                                Scalar deltaY = defUPlusEps.atSubContVol[i][compi] -
                                                 defUMinusEps.atSubContVol[i][compi];

                                // A[i][j][compi][comp] is the
                                // approximate rate of change of the
                                // unknown 'compi' at vertex 'i' if the
                                // component 'comp' at vertex 'j' is
                                // slightly deflected from the current
                                // local solution
                                this->A[i][j][compi][comp] = deltaY/deltaX;
                            }
                        }
                    }
                }

                // assemble boundary conditions
                assembleBoundaryCondition(_curCell());
                
                // calculate the right hand side
                LocalFunction defU;
                evalLocalDefect(defU,
                                localU, uCache,
                                localOldU);
                for (int i=0; i < numVertices; i++) {
                    for (int comp=0; comp < NumUnknowns; comp++) {
                        if (this->bctype[i][comp]==Dune::BoundaryConditions::neumann) {
                            this->b[i][comp] += defU.atSubContVol[i][comp];
                        }
                    }
                }
            }

        void assembleBoundaryCondition(const Cell &cell, int orderOfShapeFns=1)
            {
                setCurrentCell(cell);
                // reset the right hand side and the boundary
                // condition type vector
                _resetRhs();
                
                if (!cell.hasBoundaryIntersections())
                    return;

                Dune::GeometryType      geoType = _curCell().geometry().type();
                const ShapeFnSet       &shapeFnSet = ShapeFnSets::general(geoType, 1);
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
                                         refElem,
                                         shapeFnSet);
                }
            }

        // set the cell which is currently considered as the local
        // stiffness matrix
        void setCurrentCell(const Cell &e)
        {
            // TODO: This function is called whenever assemble() or
            // assembleBC() is called. Since both are often called
            // consecutively, and I suspect that
            // FVElementGeometry.update() is a pretty expensive
            // operation, it would be nice if there was a way to find
            // out whether FVElementGeometry really needs an update...
            _curCellPtr = &e;
            _curCellGeom.update(_curCell());
            // tell the LocalStiffness (-> parent class) the current
            // number of degrees of freedom
            setcurrentsize(_curCellGeom.nNodes);
        }

        // compute the local defect between 'solNew' and 'solOld' at
        // the vertices of the current cell. solOld must be the local
        // solution of the last time step.
        void evalLocalDefect(LocalFunction &defect,
                             const LocalFunction &solNew,
                             const CachedCellData &solNewCache,
                             const LocalFunction &solOld) const
            {
                // reset defect
                for (int i = 0; i < _curCellGeom.nNodes; i++) {
                    defect.atSubContVol[i] = 0;
                }

                // evaluate the defect of the mass balance
                for (int i=0; i < _curCellGeom.nNodes; i++)
                {
                    UnknownsVector massContrib, tmp;

                    // mass balance within the cell
                    JacobianImp::evalMassBalance(massContrib,
                                                  _problem,
                                                  _curCell(),
                                                  _curCellGeom,
                                                  solNew,
                                                  i);
                    JacobianImp::evalMassBalance(tmp,
                                                  _problem,
                                                  _curCell(),
                                                  _curCellGeom,
                                                  solOld,
                                                  i);
                    
                    massContrib -= tmp;
                    massContrib *= _curCellGeom.subContVol[i].volume/getDt();
                    defect.atSubContVol[i] += massContrib;

                    // calculate the mass which was injected since the last timestep
                    // and subtract it from the defect due to the mass balance.
                    UnknownsVector q;
                    _problem.sourceTerm(q,
                                        _curCell(),
                                        _curCellGeom,
                                        i);
                    q *= _curCellGeom.subContVol[i].volume;
                    defect.atSubContVol[i] -= q;

//                    if (HACKY_HACK) {
//                        printf("q: %f/%f\n", q[0], q[1]);
//                    }
                }

                // calculate defect of mass flux over the faces
                for (int k = 0; k < _curCellGeom.nEdges; k++)
                {
                    int i = _curCellGeom.subContVolFace[k].i;
                    int j = _curCellGeom.subContVolFace[k].j;

                    UnknownsVector flux;
                    JacobianImp::evalMassFlux(flux,
                                               _problem,
                                               _curCell(),
                                               _curCellGeom,
                                               solNew,
                                               solNewCache,
                                               k);

/*                    if (HACKY_HACK) {
                        std::cout << boost::format("flux: %g/%g\n")%flux[0]%flux[1];
                    }
*/

                    // add to defect
                    defect.atSubContVol[i] -= flux;
                    defect.atSubContVol[j] += flux;
                }
            }

        //! same as above but less efficient if the update of the
        //! entire cell cache is a slow operation. Overloaded for
        //! convenience.
        void evalLocalDefect(LocalFunction &defect,
                             const LocalFunction &solNew,
                             const LocalFunction &solOld) const
            {
                CachedCellData solNewCache;
                JacobianImp::updateCellCache(solNewCache,
                                              _problem,
                                              _curCell(),
                                              _curCellGeom,
                                              solNew);
                evalLocalDefect(defect, solNew, solNewCache, solOld);
            }

        // evaluate 'globalFn' for the current cell, save the result
        // to 'dest'
        void evalLocal(LocalFunction &dest,
                       const SpatialFunction &globalFn)
            {
                Dune::GeometryType geoType = _curCell().geometry().type();
                const ShapeFnSet &shapeFns = ShapeFnSets::general(geoType, 1);

                int size = shapeFns.size();
//                dest.resize(size);
                for (int i = 0; i < size; i++)
                    for (int comp = 0; comp < NumUnknowns; comp++)
                        dest.atSubContVol[i][comp] = globalFn.evallocal(comp,
                                                                        _curCell(),
                                                                        shapeFns[i].position());
            }


    private:
        // reset the right hand side
        void _resetRhs()
            {
                int numVertices = this->currentsize();
                for (int i=0; i < numVertices; i++) {
                    this->bctype[i].assign(Dune::BoundaryConditions::neumann);
                    this->b[i] = 0;
                }
            };

        const Cell &_curCell() const
            { return *_curCellPtr; }


       void _assembleExteriorBC(BoundaryTypeVector &faceBCType,
                                const IntersectionIterator &faceIt,
                                const CellReferenceElement &refElem)
            {
                int faceIdx = faceIt.numberInSelf();
                int nVerticesOfFace = refElem.size(faceIdx, 1, GridDim);
                for (int vertexInFace = 0; vertexInFace < nVerticesOfFace; vertexInFace++)
                {
                    int vertexInElement = refElem.subEntity(faceIdx, 1, vertexInFace, GridDim);
                    
                    // TODO: "mixed" boundary conditions are not yet
                    // possible. either it's neumann for all primary
                    // variables or it's dirchlet.
                    if (this->bctype[vertexInElement][0] == Dune::BoundaryConditions::neumann) {
                        int bfIdx = _curCellGeom.boundaryFaceIndex(faceIdx, vertexInFace);
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
                        this->b[vertexInElement] += J;
                    }
                }
            }
        

        // assemble a dirichlet boundary condition for a face which
        // may be either exterior or within the grid (if level or
        // process boundaries are required/used)
        void _assembleDirichletBC(const BoundaryTypeVector &faceBCType,
                                  const IntersectionIterator &faceIt,
                                  const CellReferenceElement &refElem,
                                  const ShapeFnSet &shapeFnSet)
            {
                // loop over shape functions
                for (int i=0; i<shapeFnSet.size(); i++)
                {
                    if (shapeFnSet[i].codim() != GridDim)
                        continue; // skip interior and face degrees of freedom
                    
                    // loop over all vertices of the face in order to
                    // make sure the entity associated with the
                    // current shape function is also a vertex of the
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
                            // the associated vertex of the i-th shape
                            // function is a vertex of the current
                            // face.
                            this->bctype[i] = faceBCType[0];
                            this->b[i] = 0;
                            break;
                        }
                    }
                }
            }

        // The problem we would like to solve
        Problem &_problem;

        const Cell             *_curCellPtr;
        FVElementGeometry       _curCellGeom;

        // global solution of the current and of the last timestep
        const SpatialFunction *_curSolution;
        const SpatialFunction *_oldSolution;

        bool         _levelBoundaryAsDirichlet;
        bool         _processBoundaryAsDirichlet;

        double      _dt;
    };
}

#endif
