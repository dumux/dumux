#ifndef DUNE_BOXTRANSPORTJACOBIAN_HH
#define DUNE_BOXTRANSPORTJACOBIAN_HH

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
#include <dune/grid/utility/intersectiongetter.hh>

#include<dune/disc/shapefunctions/lagrangeshapefunctions.hh>
#include<dune/disc/operators/boundaryconditions.hh>

#include<dumux/operators/boxjacobian.hh>
#include<dumux/transport/transportproblem.hh>

namespace Dune
{
  //! A class for computing local jacobian matrices
  /*! A class for computing local jacobian matrix for the
    diffusion equation

        div j = q; j = -K grad u; in Omega

        u = g on Gamma1; j*n = J on Gamma2.

    Uses the box method.

    Template parameters are:

    - Grid     a DUNE grid type
    - Scalar    type used for return values
  */
  template<class Grid, class Scalar, class VC, class BoxFunction = LeafP1Function<Grid, Scalar, 1> >
  class BoxTransportJacobian
    : public BoxJacobian<BoxTransportJacobian<Grid,Scalar,VC,BoxFunction>,Grid,Scalar,1,BoxFunction>
  {
    typedef typename Grid::Traits::template Codim<0>::Entity Entity;
    enum {dim=Grid::dimension};
    enum {numEq = 1};
    enum {SIZE=LagrangeShapeFunctionSetContainer<Scalar,Scalar,dim>::maxsize};
    typedef typename Entity::Geometry Geometry;
    typedef BoxTransportJacobian<Grid,Scalar,VC,BoxFunction> ThisType;
    typedef typename LocalJacobian<ThisType,Grid,Scalar,1>::VBlockType VBlockType;
    typedef Dune::FVElementGeometry<Grid> FVElementGeometry;
    typedef FieldMatrix<Scalar,dim,dim> FMatrix;
    typedef FieldVector<Scalar,dim> FVector;

  public:

    //! Constructor
    BoxTransportJacobian (TransportProblem<Grid,Scalar,VC>& params,
                  bool levelBoundaryAsDirichlet_, const Grid& grid,
                  BoxFunction& sol,
                  bool procBoundaryAsDirichlet_=true)
    : BoxJacobian<ThisType,Grid,Scalar,1,BoxFunction>(levelBoundaryAsDirichlet_, grid, sol, procBoundaryAsDirichlet_),
      problem(params), varNData(SIZE), oldVarNData(SIZE)
    {
      this->analytic = false;
    }

    void clearVisited ()
    {
        return;
    }

    VBlockType computeM (const Entity& element, const VBlockType* sol, int node, bool old = false)
    {
        return (elData.porosity*sol[node]);
    }

    VBlockType computeQ (const Entity& element, const VBlockType* sol, const int& node)
    {
        return 0;
    }

    VBlockType computeA (const Entity& element, const VBlockType* sol, int face)
    {
        FieldVector<Scalar, dim> gradS(0);
        Scalar lambdaBarFace = 0;
        Scalar dPdSFace = 0;
        for (int k = 0; k < this->fvGeom.numVertices; k++) {
            FieldVector<Scalar,dim> grad(this->fvGeom.subContVolFace[face].grad[k]);
            grad *= sol[k];
            gradS += grad;
            lambdaBarFace += varNData[k].lambdaBar*this->fvGeom.subContVolFace[face].shapeValue[k];
            dPdSFace += varNData[k].dPdS*this->fvGeom.subContVolFace[face].shapeValue[k];
        }

        gradS *= lambdaBarFace*dPdSFace;
        FieldVector<Scalar,dim> KGradS(0);
        elData.K.umv(gradS, KGradS);

        VBlockType flux = KGradS*this->fvGeom.subContVolFace[face].normal;

        int i = this->fvGeom.subContVolFace[face].i;
        int j = this->fvGeom.subContVolFace[face].j;
        Scalar outward = elData.velocity*this->fvGeom.subContVolFace[face].normal;
        if (outward > 0)
            flux -= varNData[i].fractionalW*outward;
        else
            flux -= varNData[j].fractionalW*outward;


        return flux;
    }

    void computeElementData (const Entity& element)
    {
        // ASSUMING element-wise constant permeability, evaluate K at the cell center
        elData.K = problem.soil().K(this->fvGeom.elementGlobal, element, this->fvGeom.elementLocal);
        elData.porosity = problem.soil().porosity(this->fvGeom.elementGlobal, element, this->fvGeom.elementLocal);

        switch (dim)
        {
        case 1:
            elData.velocity[0] = 0.5*(problem.variables.vTotal(element, 0)[0] + problem.variables.vTotal(element, 1)[0]);
            break;
        case 2:
            elData.velocity[0] = 0.5*(problem.variables.vTotal(element, 0)[0] + problem.variables.vTotal(element, 1)[0]);
            elData.velocity[1] = 0.5*(problem.variables.vTotal(element, 2)[1] + problem.variables.vTotal(element, 3)[1]);
            break;
        }
    };

    // the members of the struct are defined here
    struct VariableNodeData
    {
        VBlockType lambdaBar;
        VBlockType dPdS;
        VBlockType fractionalW;
    };

    void updateVariableData(const Entity& element, const VBlockType* sol,
            int i, std::vector<VariableNodeData>& varData)
    {
        FVector& global = this->fvGeom.subContVol[i].global;
        FVector& local = this->fvGeom.subContVol[i].local;

        Scalar mobilityW = problem.materialLaw().mobW(sol[i], global, element, local);
        Scalar mobilityN = problem.materialLaw().mobN(1.0 - sol[i], global, element, local);

        varData[i].dPdS = problem.materialLaw().dPdS(sol[i], global, element, local);
        varData[i].fractionalW = problem.materialLaw().fractionalW(sol[i], global, element, local);
        varData[i].lambdaBar = mobilityW*mobilityN/(mobilityW + mobilityN);

              return;
    }

     virtual void updateVariableData(const Entity& element, const VBlockType* sol, int i, bool old = false)
    {
         if (old) {
             updateVariableData(element, sol, i, oldVarNData);
         }
         else
             updateVariableData(element, sol, i, varNData);

         return;
    }

    void updateVariableData(const Entity& element, const VBlockType* sol, bool old = false)
    {
        int size = this->fvGeom.numVertices;

        for (int i = 0; i < size; i++)
                updateVariableData(element, sol, i, old);

        return;
    }

    virtual void updateStaticData (const Entity& element, const VBlockType* sol)
    {
        return;
    }

    template<class TypeTag> void assembleBC(const Entity& e) {
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
            this->dirichletIndex[i] = 0;
        }

        // evaluate boundary conditions via intersection iterator
        typedef typename IntersectionIteratorGetter<Grid,TypeTag>::IntersectionIterator IntersectionIterator;

        IntersectionIterator endit = IntersectionIteratorGetter<Grid, TypeTag>::end(e);
        for (IntersectionIterator it = IntersectionIteratorGetter<Grid, TypeTag>::begin(e); it!=endit; ++it)
        {
            // if we have a neighbor then we assume there is no boundary (forget interior boundaries)
            // in level assemble treat non-level neighbors as boundary
            if (it->neighbor()) {
                if (this->levelBoundaryAsDirichlet && it->outside()->level()==e.level())
                    continue;
                if (!this->levelBoundaryAsDirichlet)
                    continue;
            }

            // determine boundary condition type for this face, initialize with processor boundary
            FieldVector<typename BoundaryConditions::Flags, numEq> bctypeface(BoundaryConditions::process);
            FieldVector<int,numEq> dirichletIdx(0);

            // handle face on exterior boundary, this assumes there are no interior boundaries
            if (it->boundary()) {
                int faceIdx = it->numberInSelf();
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
                            bctypeface = this->getImp().problem.bctypeSat(global, e, local); // eval bctype

                            if (bctypeface[equationNumber]!=BoundaryConditions::neumann)
                                break;
                            VBlockType J = this->getImp().problem.neumannSat(global, e, local, 0);
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
                        && this->procBoundaryAsDirichlet==false
                        && this->levelBoundaryAsDirichlet==false)
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
                        if (sfs[i].entity()==it->numberInSelf()) {
                            if (this->bctype[i][equationNumber] < bctypeface[equationNumber]) {
                                this->bctype[i][equationNumber] = bctypeface[equationNumber];
                                this->dirichletIndex[i][equationNumber] = dirichletIdx[equationNumber];

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
                    for (int j=0; j<ReferenceElements<Scalar,dim>::general(gt).size(it->numberInSelf(), 1, sfs[i].codim()); j++)
                        if (sfs[i].entity()==ReferenceElements<Scalar,dim>::general(gt).subEntity(it->numberInSelf(), 1, j, sfs[i].codim()))
                        {
                            if (this->bctype[i][equationNumber] < bctypeface[equationNumber]) {
                                this->bctype[i][equationNumber] = bctypeface[equationNumber];
                                this->dirichletIndex[i][equationNumber] = dirichletIdx[equationNumber];
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

    struct ElementData {
        FieldMatrix<Scalar,dim,dim> K;
        Scalar porosity;
        FVector velocity;
       };

       ElementData elData;
       std::vector<VariableNodeData> varNData;
       std::vector<VariableNodeData> oldVarNData;
    TransportProblem<Grid,Scalar,VC>& problem;
  };
}
#endif
