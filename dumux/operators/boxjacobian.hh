// $Id$

#ifndef DUNE_BOXJACOBIAN_HH
#define DUNE_BOXJACOBIAN_HH

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
#include<dune/grid/common/mcmgmapper.hh>

#include<dune/disc/shapefunctions/lagrangeshapefunctions.hh>
#include<dune/disc/operators/boundaryconditions.hh>
#include<dune/disc/functions/p1function.hh>

#include"localjacobian.hh"

/**
 * @file
 * @brief  compute local jacobian matrix for conforming finite elements for diffusion equation
 * @author Peter Bastian
 */

namespace Dune {
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
  - Scalar    type used for return values
*/
template<class Imp, class Grid, class Scalar, int numEq,
         class BoxFunction = LeafP1Function<Grid, Scalar, numEq> >
class BoxJacobian :
        public LocalJacobian<Imp,Grid,Scalar,numEq> {
    // mapper: one data element per vertex
    template<int dim> struct P1Layout
    {
        bool contains(Dune::GeometryType gt) {
            return gt.dim() == 0;
        }
    };

    typedef typename Grid::Traits::template Codim<0>::Entity Element;
    typedef typename Grid::template Codim<0>::LeafIntersectionIterator IntersectionIterator;
    typedef typename Element::Geometry Geometry;
    typedef typename LocalJacobian<Imp,Grid,Scalar,numEq>::VBlockType VBlockType;
    typedef typename LocalJacobian<Imp,Grid,Scalar,numEq>::MBlockType MBlockType;
    typedef MultipleCodimMultipleGeomTypeMapper<typename Grid::LeafGridView, P1Layout> VertexMapper;

public:
    // define the number of components of your system, this is used outside
    // to allocate the correct size of (dense) blocks with a FieldMatrix
    enum {dim=Grid::dimension};
    enum {SIZE=LagrangeShapeFunctionSetContainer<Scalar,Scalar,dim>::maxsize};

    //! Constructor
    BoxJacobian(bool levelBoundaryAsDirichlet_, const Grid& grid,
                BoxFunction& sol, bool procBoundaryAsDirichlet_=true) :
        vertexMapper(grid.leafView()),
        levelBoundaryAsDirichlet(levelBoundaryAsDirichlet_),
        procBoundaryAsDirichlet(procBoundaryAsDirichlet_),
        currentSolution(sol), oldSolution(grid), dt(1)
        {    }

    //**********************************************************
    //*                                                                            *
    //*    Computation of the local defect                                *
    //*                                                                            *
    //**********************************************************
    void localDefect(const Element& e, const VBlockType* sol, bool withBC = true) {
        for (int i=0; i < this->fvGeom.numVertices; i++) // begin loop over vertices / sub control volumes
        {
            // implicit Euler
            bool old = true;
            VBlockType massContrib = computeM(e, this->uold, i, old);
            massContrib *= -1.0;
            this->def[i] = massContrib;
        }

        //updateVariableData(e, sol);
        for (int i=0; i < this->fvGeom.numVertices; i++) // begin loop over vertices / sub control volumes
        {
            VBlockType massContrib = computeM(e, sol, i);
            this->def[i] += massContrib;
            this->def[i] *= this->fvGeom.subContVol[i].volume/dt;

            // get source term
            VBlockType q = computeQ(e, sol, i);
            q *= this->fvGeom.subContVol[i].volume;
            this->def[i] -= q;
        } // end loop over vertices / sub control volumes

        for (int k = 0; k < this->fvGeom.numEdges; k++) // begin loop over edges / sub control volume faces
        {
            int i = this->fvGeom.subContVolFace[k].i;
            int j = this->fvGeom.subContVolFace[k].j;

            VBlockType flux = computeA(e, sol, k);

            // add to defect
            this->def[i] -= flux;
            this->def[j] += flux;
            //                std::cout << "i = " << i << ", j = " << j << ", flux = " << flux << std::endl;
        } // end loop over edges / sub control volume faces

        if (withBC) {
            // assemble boundary conditions
            this->getImp().assembleBoundaryCondition(e);
            //                for (int i = 0; i < this->fvGeom.numVertices; i++)
            //                    std::cout << this->fvGeom.subContVol[i].global << ": cond = " << this->bctype[i][0] << std::endl;

            //                assembleBC(e);

            // add to defect
            for (int i=0; i < this->fvGeom.numVertices; i++) {
                this->def[i] += this->b[i];
            }
        }

        return;
    }

    void setLocalSolution(const Element& e)
    {
        Dune::GeometryType gt = e.geometry().type();
        const typename Dune::LagrangeShapeFunctionSetContainer<Scalar,Scalar,dim>::value_type
            &sfs=Dune::LagrangeShapeFunctions<Scalar, Scalar, dim>::general(gt, 1);
        int size = sfs.size();
        this->setcurrentsize(size);

        for (int i = 0; i < size; i++)
            for (int comp = 0; comp < numEq; comp++) {
                this->u[i][comp] = currentSolution.evallocal(comp, e, sfs[i].position());
                this->uold[i][comp] = oldSolution.evallocal(comp, e, sfs[i].position());
            }

        return;
    }

    void localToGlobal(const Element& e, const VBlockType* sol)
    {
        int size = e.template count<dim>();
        for (int i = 0; i < size; i++) {
            int globalIdx = vertexMapper.template map<dim>(e, i);
            (*currentSolution)[globalIdx] = this->u[i];
        }
    }


    void setDt(double d) {
        dt = d;

        return;
    }

    double getDt() {
        return dt;
    }

    void setOldSolution(BoxFunction& uOld) {
        *oldSolution = *uOld;
    }

    VBlockType computeM(const Element& e, const VBlockType* sol, int node, bool old = false) {
        return this->getImp().computeM(e, sol, node, old);
    }

    VBlockType computeQ(const Element& e, const VBlockType* sol, int node) {
        return this->getImp().computeQ(e, sol, node);
    }

    VBlockType computeA(const Element& e, const VBlockType* sol, int face) {
        return this->getImp().computeA(e, sol, face);
    }

    // analog to EvalStaticData in MUFTE
    virtual void updateStaticData(const Element& e, VBlockType* sol) {
        return this->getImp().updateStaticData(e, sol);
    }

    virtual void assembleBoundaryCondition(const Element& e, int k = 1)
    {
        Dune::GeometryType gt = e.geometry().type();
        const typename Dune::LagrangeShapeFunctionSetContainer<Scalar,Scalar,dim>::value_type
            &sfs=Dune::LagrangeShapeFunctions<Scalar, Scalar, dim>::general(gt, 1);
        setcurrentsize(sfs.size());
        this->fvGeom.update(e);

        const typename ReferenceElementContainer<Scalar,dim>::value_type
            &referenceElement = ReferenceElements<Scalar, dim>::general(gt);

        for (int i = 0; i < sfs.size(); i++) {
            this->bctype[i].assign(BoundaryConditions::neumann);
            this->b[i] = 0;
	    //            this->dirichletIndex[i] = 0;
        }

        // evaluate boundary conditions via intersection iterator
        IntersectionIterator endit = e.ileafend();
        for (IntersectionIterator it = e.ileafbegin(); it!=endit; ++it)
        {
            // if we have a neighbor then we assume there is no boundary (forget interior boundaries)
            // in level assemble treat non-level neighbors as boundary
            if (it->neighbor()) {
                if (levelBoundaryAsDirichlet && it->outside()->level()==e.level())
                    continue;
                if (!levelBoundaryAsDirichlet)
                    continue;
            }

            // determine boundary condition type for this face, initialize with processor boundary
            FieldVector<typename BoundaryConditions::Flags, numEq> bctypeface(BoundaryConditions::process);
            // FieldVector<int,numEq> dirichletIdx(0);

            // handle face on exterior boundary, this assumes there are no interior boundaries
            if (it->boundary()) {
                int faceIdx = it->numberInInside();
                //                 std::cout << "faceIdx = " << faceIdx << ", beginning: " << std::endl;
                //                 for (int i = 0; i < 4; i++)
                //                   std::cout << "bctype[" << i << "] = " << this->bctype[i] << std::endl;

                int numVerticesOfFace = referenceElement.size(faceIdx, 1, dim);
                for (int nodeInFace = 0; nodeInFace < numVerticesOfFace; nodeInFace++) {
                    int nodeInElement = referenceElement.subEntity(faceIdx, 1, nodeInFace, dim);
                    for (int equationNumber = 0; equationNumber < numEq; equationNumber++) {
                        if (this->bctype[nodeInElement][equationNumber] == BoundaryConditions::neumann) {
                            int bfIdx = this->fvGeom.boundaryFaceIndex(faceIdx,    nodeInFace);
                            FieldVector<Scalar,dim> local = this->fvGeom.boundaryFace[bfIdx].ipLocal;
                            FieldVector<Scalar,dim> global = this->fvGeom.boundaryFace[bfIdx].ipGlobal;
                            bctypeface = this->getImp().problem.bctype(global, e, it, local); // eval bctype
                            // this->getImp().problem.dirichletIndex(global, e, it, local, dirichletIdx); // eval bctype
                            //                                                     std::cout << "faceIdx = " << faceIdx << ", nodeInElement = " << nodeInElement
                            //                                                           << ", bfIdx = " << bfIdx << ", local = " << local << ", global = " << global
                            //                                                           << ", bctypeface = " << bctypeface << std::endl;
                            if (bctypeface[equationNumber]!=BoundaryConditions::neumann)
                                break;
                            VBlockType J = this->getImp().problem.J(global, e, it, local);
                            J[equationNumber] *= this->fvGeom.boundaryFace[bfIdx].area;
                            this->b[nodeInElement][equationNumber] += J[equationNumber];
                        }
                    }
                }

                bool nface(true); // check if face is a neumann face
                for(int i=0; i<numEq; i++)
                {
                    if(bctypeface[i] != BoundaryConditions::neumann)
                        nface = false; // was not a neumann face
                }
                if(nface == true)
                    continue; // was a neumann face, go to next face
            }

            // If we are here, then it is
            // (i)   an exterior boundary face with Dirichlet condition, or
            // (ii)  a processor boundary (i.e. neither boundary() nor neighbor() was true), or
            // (iii) a level boundary in case of level-wise assemble
            // How processor boundaries are handled depends on the processor boundary mode

            bool pface(false);  // check if face is a process boundary
            for(int i=0; i<numEq; i++)
            {
                if (bctypeface[i]==BoundaryConditions::process
                    && procBoundaryAsDirichlet==false
                    && levelBoundaryAsDirichlet==false)
                {
                    pface = true;
                    break;
                }
            }
            if(pface == true)
                continue;   // if face is a process boundary it acts like homogeneous Neumann


            for (int equationNumber=0; equationNumber<numEq; equationNumber++) {
                for (int i=0; i<sfs.size(); i++) // loop over test function number
                {
                    //this->dirichletIndex[i][equationNumber] = equationNumber;

                    //std::cout<<"i = "<<i<<std::endl;
                    if (sfs[i].codim()==0)
                        continue; // skip interior dof
                    if (sfs[i].codim()==1) // handle face dofs
                    {
                        if (sfs[i].entity()==it->numberInInside()) {
                            if (this->bctype[i][equationNumber] < bctypeface[equationNumber]) {
                                this->bctype[i][equationNumber] = bctypeface[equationNumber];
                                // this->dirichletIndex[i][equationNumber] = dirichletIdx[equationNumber];

                                if (bctypeface[equationNumber] == BoundaryConditions::process)
                                    this->b[i][equationNumber] = 0;
                                if (bctypeface[equationNumber] == BoundaryConditions::dirichlet) {
                                    this->b[i][equationNumber] = 0;
                                }
                            }
                        }
                        continue;
                    }
                    // handle subentities of this face
                    for (int j=0; j<ReferenceElements<Scalar,dim>::general(gt).size(it->numberInInside(), 1, sfs[i].codim()); j++)
                        if (sfs[i].entity()==ReferenceElements<Scalar,dim>::general(gt).subEntity(it->numberInInside(), 1, j, sfs[i].codim()))
                        {
                            if (this->bctype[i][equationNumber] < bctypeface[equationNumber]) {
                                this->bctype[i][equationNumber] = bctypeface[equationNumber];
                                // this->dirichletIndex[i][equationNumber] = dirichletIdx[equationNumber];
                                if (bctypeface[equationNumber] == BoundaryConditions::process)
                                    this->b[i][equationNumber] = 0;
                                if (bctypeface[equationNumber] == BoundaryConditions::dirichlet) {
                                    this->b[i][equationNumber] = 0;
                                }
                            }
                        }
                }
            }
        }

    }

    // parameters given in constructor
    VertexMapper vertexMapper;
    bool levelBoundaryAsDirichlet;
    bool procBoundaryAsDirichlet;
    BoxFunction& currentSolution;
    BoxFunction oldSolution;

public:
    double dt;
};
}
#endif
