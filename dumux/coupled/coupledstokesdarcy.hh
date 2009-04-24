// $Id$

#ifndef DUNE_COUPLEDSTOKESDARCY_HH
#define DUNE_COUPLEDSTOKESDARCY_HH

#include "coupledmodel.hh"
#include "dumux/fvgeometry/fvelementgeometry.hh"

namespace Dune
{
/** \todo Please doc me! */
template<class StokesModel, class DarcyModel>
class CoupledStokesDarcy : public CoupledModel<StokesModel, DarcyModel, CoupledStokesDarcy<StokesModel, DarcyModel> > {
public:
    template<int dim>
    struct NodeLayout
    {
        bool contains(GeometryType gt) {
            return gt.dim() == 0;
        }
    };
    template<int dim>
    struct ElementLayout
    {
        bool contains(GeometryType gt) {
            return gt.dim() == dim;
        }
    };

    // typedefs depending on base class
    typedef CoupledModel<StokesModel, DarcyModel, CoupledStokesDarcy<StokesModel, DarcyModel> > BaseType;
    typedef typename BaseType::FirstGrid StokesGrid;
    typedef typename BaseType::SecondGrid DarcyGrid;

    // new typedefs
    enum {dim = StokesGrid::dimension};
    enum {vOrder = StokesModel::v_ordr};
    typedef typename StokesGrid::LeafGridView StokesGV;
    typedef typename StokesGV::IndexSet StokesIS;
    typedef typename StokesGV::template Codim<0>::Iterator StokesIterator;
    typedef typename StokesGrid::template Codim<dim>::EntityPointer StokesVPointer;
    typedef typename StokesGV::IntersectionIterator StokesIntersectionIterator;
    typedef MultipleCodimMultipleGeomTypeMapper<StokesGV,ElementLayout> StokesEM;
    typedef typename StokesGrid::HostGridType StokesHostGrid;
    typedef typename StokesHostGrid::template Codim<0>::EntityPointer StokesHostPointer;

    typedef typename DarcyGrid::LeafGridView DarcyGV;
    typedef typename DarcyGV::IndexSet DarcyIS;
    typedef typename DarcyGV::template Codim<0>::Iterator DarcyIterator;
    typedef typename DarcyGV::template Codim<dim>::Iterator DarcyVIterator;
    typedef typename DarcyGV::IntersectionIterator DarcyIntersectionIterator;
    typedef MultipleCodimMultipleGeomTypeMapper<DarcyGV,NodeLayout> DarcyVM;

    typedef typename DarcyGrid::HostGridType HostGrid;
    typedef typename HostGrid::template Codim<0>::EntityPointer HostPointer;
    typedef typename HostGrid::template Codim<dim>::EntityPointer HostVPointer;
    typedef MonomialShapeFunctionSet<double,double,dim> ShapeFunctionSet;

    template <class ASDType, class ADSType>
    void assembleCoupling(ASDType& A_SD, ADSType& A_DS)
    {
        // intialize rowsizes with 0
        for (typename ASDType::size_type i = 0; i < A_SD.N(); i++)
            A_SD.setrowsize(i, 0);
        for (typename ADSType::size_type i = 0; i < A_DS.N(); i++)
            A_DS.setrowsize(i, 0);

        // loop 1 over all elements of the Stokes grid to set the rowsizes of the coupling matrices
        StokesIterator endIt = (this->stokesGrid_).template leafend<0>();
        for (StokesIterator stokesIt = (this->stokesGrid_).template leafbegin<0>(); stokesIt != endIt; ++stokesIt)
        {
            GeometryType gt = stokesIt->geometry().type();
            const typename ReferenceElementContainer<double,dim>::value_type& referenceElement = ReferenceElements<double, dim>::general(gt);

            StokesIntersectionIterator endIsIt = stokesIt->ileafend();
            for (StokesIntersectionIterator stokesIsIt = stokesIt->ileafbegin(); stokesIsIt != endIsIt; ++ stokesIsIt)
            {
                if (stokesIsIt->neighbor())
                    continue; // only work on border entities

                // In the following, it is determined whether the intersection really is on the interface.
                // Every node of the corresponding face is checked whether it also belongs to the Darcy grid.
                int faceIdx = stokesIsIt->indexInInside();
                int numVerticesOfFace = referenceElement.size(faceIdx, 1, dim);
                int numVerticesInDarcyGrid = 0;
                int darcyIds[numVerticesOfFace];
                for (int nodeInFace = 0; nodeInFace < numVerticesOfFace; nodeInFace++)
                {
                    int nodeInElement = referenceElement.subEntity(faceIdx, 1, nodeInFace, dim);

                    // get the node pointer on the Stokes grid
                    const StokesVPointer& stokesVPointer = (*stokesIt).template subEntity<dim>(nodeInElement);

                    // get the node pointer on the host grid
                    const HostVPointer& hostVPointer = (this->stokesGrid_).template getHostEntity<dim>(*stokesVPointer);

                    // check if the node is also part of the darcy grid
                    if (!((this->darcyGrid_).template contains<dim>(hostVPointer)))
                        break; // otherwise this face is not at the interface

                    // get the index of the node with respect to the Darcy matrix
                    darcyIds[numVerticesInDarcyGrid] = darcyVM_.map(*stokesVPointer);

                    numVerticesInDarcyGrid++;
                }

                if (numVerticesInDarcyGrid < numVerticesOfFace) // then this face is not at the interface
                    continue;

                // get the index of the element with respect to the Stokes matrix
                int stokesId = stokesEM_.map(*stokesIt);

                // set rowsize for coupling Stokes <- Darcy
                A_SD.setrowsize(stokesId, numVerticesOfFace);

                // set rowsize for coupling Darcy <- Stokes
                for (int nodeInFace = 0; nodeInFace < numVerticesOfFace; nodeInFace++)
                    A_DS.incrementrowsize(darcyIds[nodeInFace]);
            }
        } // end loop 1 over all elements of the Stokes grid
        A_SD.endrowsizes();
        A_DS.endrowsizes();


        // loop 2 over all elements of the Stokes grid to set the indices of the nonzero entries of the coupling matrices
        for (StokesIterator stokesIt = (this->stokesGrid_).template leafbegin<0>(); stokesIt != endIt; ++stokesIt)
        {
            GeometryType gt = stokesIt->geometry().type();
            const typename ReferenceElementContainer<double,dim>::value_type& referenceElement = ReferenceElements<double, dim>::general(gt);

            StokesIntersectionIterator endIsIt = stokesIt->ileafend();
            for (StokesIntersectionIterator stokesIsIt = stokesIt->ileafbegin(); stokesIsIt != endIsIt; ++ stokesIsIt)
            {
                if (stokesIsIt->neighbor())
                    continue; // only work on border entities

                // In the following, it is determined whether the intersection really is on the interface.
                // Every node of the corresponding face is checked whether it also belongs to the Darcy grid.
                int faceIdx = stokesIsIt->indexInInside();
                int numVerticesOfFace = referenceElement.size(faceIdx, 1, dim);
                int numVerticesInDarcyGrid = 0;
                int darcyIds[numVerticesOfFace];
                for (int nodeInFace = 0; nodeInFace < numVerticesOfFace; nodeInFace++)
                {
                    int nodeInElement = referenceElement.subEntity(faceIdx, 1, nodeInFace, dim);

                    // get the node pointer on the Stokes grid
                    const StokesVPointer& stokesVPointer = (*stokesIt).template subEntity<dim>(nodeInElement);

                    // get the node pointer on the host grid
                    const HostVPointer& hostVPointer = (this->stokesGrid_).template getHostEntity<dim>(*stokesVPointer);

                    // check if the node is also part of the darcy grid
                    if (!((this->darcyGrid_).template contains<dim>(hostVPointer)))
                        break; // otherwise this face is not at the interface

                    // get the index of the node with respect to the Darcy matrix
                    darcyIds[numVerticesInDarcyGrid] = darcyVM_.map(*stokesVPointer);

                    numVerticesInDarcyGrid++;
                }

                if (numVerticesInDarcyGrid < numVerticesOfFace) // then this face is not at the interface
                    continue;

                // get the index of the element with respect to the Stokes matrix
                int stokesId = stokesEM_.map(*stokesIt);

                // set indices of nonzero entries
                for (int nodeInFace = 0; nodeInFace < numVerticesOfFace; nodeInFace++) {
                    A_SD.addindex(stokesId, darcyIds[nodeInFace]);
                    A_DS.addindex(darcyIds[nodeInFace], stokesId);
                }
            }
        } // end loop 2 over all elements of the Stokes grid
        A_SD.endindices();
        A_DS.endindices();
        A_SD = 0;
        A_DS = 0;

        // loop 3 over all elements of the Stokes grid to actually fill the coupling matrices
        for (StokesIterator stokesIt = (this->stokesGrid_).template leafbegin<0>(); stokesIt != endIt; ++stokesIt)
        {
            GeometryType gt = stokesIt->geometry().type();
            const typename ReferenceElementContainer<double,dim>::value_type& referenceElement = ReferenceElements<double, dim>::general(gt);

            FVElementGeometry<StokesGrid> fvGeom;
            fvGeom.update(*stokesIt);

            StokesIntersectionIterator endIsIt = stokesIt->ileafend();
            for (StokesIntersectionIterator stokesIsIt = stokesIt->ileafbegin(); stokesIsIt != endIsIt; ++ stokesIsIt)
            {
                if (stokesIsIt->neighbor())
                    continue; // only work on border entities

                // In the following, it is determined whether the intersection really is on the interface.
                // Every node of the corresponding face is checked whether it also belongs to the Darcy grid.
                int faceIdx = stokesIsIt->indexInInside();
                int numVerticesOfFace = referenceElement.size(faceIdx, 1, dim);
                int numVerticesInDarcyGrid = 0;
                int darcyIds[numVerticesOfFace];
                int nodeInElement[numVerticesOfFace];
                for (int nodeInFace = 0; nodeInFace < numVerticesOfFace; nodeInFace++) {
                    int nInEl = referenceElement.subEntity(faceIdx, 1, nodeInFace, dim);

                    // get the node pointer on the Stokes grid
                    const StokesVPointer& stokesVPointer = (*stokesIt).template subEntity<dim>(nInEl);

                    // get the node pointer on the host grid
                    const HostVPointer& hostVPointer = (this->stokesGrid_).template getHostEntity<dim>(*stokesVPointer);

                    // check if the node is also part of the darcy grid
                    if (!((this->darcyGrid_).template contains<dim>(hostVPointer)))
                        break; // otherwise this face is not at the interface

                    // get the index of the node with respect to the Darcy matrix
                    darcyIds[numVerticesInDarcyGrid] = darcyVM_.map(*stokesVPointer);

                    // save the local index of the node inside the element
                    nodeInElement[numVerticesInDarcyGrid] = nInEl;

                    numVerticesInDarcyGrid++;
                }

                if (numVerticesInDarcyGrid < numVerticesOfFace) // then this face is not at the interface
                    continue;

                // get the index of the element with respect to the Stokes matrix
                int stokesId = stokesEM_.map(*stokesIt);


                // get the geometry type of the face
                GeometryType geomTypeBoundary = stokesIsIt->geometryInInside().type();

                // The coupling Stokes <- Darcy is realized in the FE way.
                // The unknown Darcy pressure is the piecewise linear FE interpolant.
                int qOrder = vOrder + 1;
                ShapeFunctionSet velShapeFuncSet(vOrder);
                for(unsigned int qNode = 0; qNode < QuadratureRules<double,dim-1>::rule(geomTypeBoundary,qOrder).size(); ++qNode)
                {
                    const FieldVector<double,dim-1>& qLocalDimM1 = QuadratureRules<double,dim-1>::rule(geomTypeBoundary,qOrder)[qNode].position();
                    FieldVector<double,dim> qLocal = stokesIsIt->geometryInInside().global(qLocalDimM1);
                    double qWeight = QuadratureRules<double,dim-1>::rule(geomTypeBoundary,qOrder)[qNode].weight();
                    double qDetJac = stokesIsIt->geometry().integrationElement(qLocalDimM1);
                    FieldVector<double,dim> normal = stokesIsIt->unitOuterNormal(qLocalDimM1);
                    const typename LagrangeShapeFunctionSetContainer<double,double,dim>::value_type&
                        pressShapeFuncSet = LagrangeShapeFunctions<double,double,dim>::general(gt,1);

                    for(int comp = 0; comp < dim; ++comp) // loop over the velocity components
                    {
                        for (int i = 0; i < velShapeFuncSet.size(); ++i) // loop over the scalar velocity basis functions
                        {
                            int ii = comp*velShapeFuncSet.size() + i;
                            double velShapeValue = velShapeFuncSet[i].evaluateFunction(0,qLocal);
                            for (int j = 0; j < numVerticesOfFace; ++j) // loop over interface nodes
                            {
                                int jInElement = nodeInElement[j];
                                double pressShapeValue = pressShapeFuncSet[jInElement].evaluateFunction(0,qLocal);
                                double entrySD = (pressShapeValue*(velShapeValue*normal[comp]))* qDetJac * qWeight;
                                int jInMatrix = darcyIds[j];
                                (A_SD[stokesId][jInMatrix])[ii] -= entrySD;
                            }
                        }
                    }
                }

                // The coupling Darcy <- Stokes is realized in the FV way.
                // The unknown Stokes normal velocity is evaluated in the center of the subcontrolvolume face.
                const FieldVector<double,dim-1>& faceLocalDimM1 = ReferenceElements<double,dim-1>::general(geomTypeBoundary).position(0,0);
                FieldVector<double,dim> normal = stokesIsIt->unitOuterNormal(faceLocalDimM1);
                for(int comp = 0; comp < dim; ++comp) // loop over the velocity components
                {
                    for (int i = 0; i < velShapeFuncSet.size(); ++i) // loop over the scalar velocity basis functions
                    {
                        int ii = comp*velShapeFuncSet.size() + i;
                        for (int j = 0; j < numVerticesOfFace; ++j) // loop over interface nodes
                        {
                            int bfIdx = fvGeom.boundaryFaceIndex(faceIdx, j);
                            FieldVector<double,dim> qLocal = fvGeom.boundaryFace[bfIdx].ipLocal;
                            double velShapeValue = velShapeFuncSet[i].evaluateFunction(0,qLocal);
                            double entryDS = (velShapeValue*normal[comp])*fvGeom.boundaryFace[bfIdx].area;
                            int jInMatrix = darcyIds[j];
                            (A_DS[jInMatrix][stokesId])[0][ii] += entryDS;
                        }
                    }
                }
            }
        } // end loop 3 over all elements of the Stokes grid
    }

    virtual void solve()
    {
        typedef typename BaseType::GlobalMatrixType GlobalMatrixType;
        typedef typename BaseType::GlobalVectorType GlobalVectorType;

        typedef MatrixAdapter<GlobalMatrixType,GlobalVectorType,GlobalVectorType> Operator;
        Operator op(this->A);
        double red=1E-14;
        SeqPardiso<GlobalMatrixType,GlobalVectorType,GlobalVectorType> pardiso(this->A);
        LoopSolver<GlobalVectorType> solver(op,pardiso,red,10000,1);
        //      SeqILU0<GlobalMatrixType,GlobalVectorType,GlobalVectorType> ilu(this->A, 0.9);
        //      BiCGSTABSolver<GlobalVectorType> solver(op,ilu,red,10000,1);
        InverseOperatorResult r;
        solver.apply(this->u, this->f, r);

        int rowsInBlock1 = BaseType::FirstMatrixType::block_type::rows;
        int nOfBlockRows1 = (this->firstModel_).matrix().N();
        DarcyIterator dummyIT2((this->darcyGrid_).template leafbegin<0>());
        DarcyIntersectionIterator dummyIS2(dummyIT2->ileafbegin());
        DarcyVIterator endItV2 = (this->darcyGrid_).template leafend<dim>();
        for (DarcyVIterator it = (this->darcyGrid_).template leafbegin<dim>(); it != endItV2; ++it)
        {
            FieldVector<double,dim> globalCoord = (*it).geometry().corner(0);
            BoundaryConditions::Flags bctype = (this->secondModel_).problem.bctype(globalCoord, *dummyIT2, dummyIS2, globalCoord);
            int darcyId = darcyVM_.map(*it);
            this->u[rowsInBlock1*nOfBlockRows1 + darcyId] = -this->u[rowsInBlock1*nOfBlockRows1 + darcyId];

            if (bctype == BoundaryConditions::dirichlet)
            {
                double dirichletBC = (this->secondModel_).problem.g(globalCoord, *dummyIT2, dummyIS2, globalCoord);
                this->u[rowsInBlock1*nOfBlockRows1 + darcyId] += dirichletBC;
            }
        }

        // transfer to local solution vectors
        const typename BaseType::FirstMatrixType::block_type::size_type colsInBlock1 = BaseType::FirstMatrixType::block_type::cols;
        const typename BaseType::SecondMatrixType::block_type::size_type colsInBlock2 = BaseType::SecondMatrixType::block_type::cols;
        for (int i = 0; i < this->firstModel_.sol().size(); i++)
            for (typename BaseType::FirstMatrixType::block_type::size_type k = 0; k < colsInBlock1; k++)
                this->firstModel_.sol()[i][k] = this->u[i*colsInBlock1 + k];
        for (int i = 0; i < this->secondModel_.sol().size(); i++)
            for (typename BaseType::SecondMatrixType::block_type::size_type k = 0; k < colsInBlock2; k++)
                this->secondModel_.sol()[i][k] = this->u[colsInBlock1*this->firstModel_.sol().size() + i*colsInBlock2 + k];
    }

    virtual void vtkout (const char* name, int k)
    {
        char localName[128];
        sprintf(localName, "%s_stokes", name);
        this->firstModel_.vtkout(localName, k);
        sprintf(localName, "%s_darcy", name);
        this->secondModel_.vtkout(localName, k);
    }

    CoupledStokesDarcy(const StokesGrid& stokesGrid, StokesModel& stokesModel,
                       const DarcyGrid& darcyGrid, DarcyModel& darcyModel,
                       bool assembleGlobalSystem)
        : BaseType(stokesGrid, stokesModel, darcyGrid, darcyModel, assembleGlobalSystem),
          stokesGrid_(this->firstGrid_), darcyGrid_(this->secondGrid_),
          stokesEM_(stokesGrid_.leafView()), darcyVM_(darcyGrid_.leafView())
    {}

private:
    const StokesGrid& stokesGrid_;
    const DarcyGrid& darcyGrid_;
    StokesEM stokesEM_;
    DarcyVM darcyVM_;
};

}

#endif
