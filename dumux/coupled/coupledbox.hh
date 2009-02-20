// $Id: coupledfirstsecond.hh 1128 2009-02-06 09:07:10Z anneli $

#ifndef DUNE_COUPLEDBOX_HH
#define DUNE_COUPLEDBOX_HH

#include "coupledmodel.hh"
#include "dumux/fvgeometry/fvelementgeometry.hh"

namespace Dune
{
/** \todo Please doc me! */
template<class FirstModel, class SecondModel, class Imp>
class CoupledBox : public CoupledModel<FirstModel, SecondModel, Imp> {
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
    typedef CoupledModel<FirstModel, SecondModel, Imp> BaseType;
    typedef typename BaseType::FirstGrid FirstGrid;
    typedef typename BaseType::SecondGrid SecondGrid;
    typedef typename FirstModel::MatrixType FirstMatrixType;
    typedef typename SecondModel::MatrixType SecondMatrixType;
    enum{firstNumEq = FirstMatrixType::block_type::rows};
    enum{secondNumEq = SecondMatrixType::block_type::rows};

    // new typedefs
    enum {dim = FirstGrid::dimension};
    typedef typename FirstGrid::LeafGridView FirstGV;
    typedef typename FirstGV::IndexSet FirstIS;
    typedef typename FirstGV::template Codim<0>::Iterator FirstIterator;
    typedef typename FirstGrid::Traits::template Codim<0>::Entity FirstElement;
    typedef typename FirstGrid::template Codim<dim>::EntityPointer FirstVPointer;
    typedef typename IntersectionIteratorGetter<FirstGrid,LeafTag>::IntersectionIterator FirstIntersectionIterator;
    typedef MultipleCodimMultipleGeomTypeMapper<FirstGrid,FirstIS,NodeLayout> FirstVertexMapper;
    typedef typename FirstGrid::HostGridType FirstHostGrid;
    typedef typename FirstHostGrid::template Codim<0>::EntityPointer FirstHostPointer;

    typedef typename SecondGrid::LeafGridView SecondGV;
    typedef typename SecondGV::IndexSet SecondIS;
    typedef typename SecondGV::template Codim<0>::Iterator SecondIterator;
    typedef typename SecondGrid::Traits::template Codim<0>::Entity SecondElement;
    typedef typename SecondGV::template Codim<dim>::Iterator SecondVIterator;
    typedef typename IntersectionIteratorGetter<SecondGrid,LeafTag>::IntersectionIterator SecondIntersectionIterator;
    typedef MultipleCodimMultipleGeomTypeMapper<SecondGrid,SecondIS,NodeLayout> SecondVertexMapper;

    typedef typename SecondGrid::HostGridType HostGrid;
    typedef typename HostGrid::template Codim<0>::EntityPointer HostPointer;
    typedef typename HostGrid::template Codim<dim>::EntityPointer HostVPointer;

    template <class FirstFV, class SecondFV>
    void localCoupling12(FirstFV& firstSol, SecondFV& secondSol, int firstIndex, int secondIndex,
            FieldVector<double, dim> qGlobal, FirstFV& result)
    {
        this->getImp().template localCoupling12<FirstFV,SecondFV>(firstSol, secondSol,
                firstIndex, secondIndex, qGlobal, result);
    }

    template <class FirstFV, class SecondFV>
    void localCoupling21(FirstFV& firstSol, SecondFV& secondSol, int firstIndex, int secondIndex,
            FieldVector<double, dim> qGlobal, SecondFV& result)
    {
        this->getImp().template localCoupling21<FirstFV,SecondFV>(firstSol, secondSol,
                firstIndex, secondIndex, qGlobal, result);
    }

    template <class A12Type, class A21Type>
    void assembleCoupling(A12Type& A_12, A21Type& A_21)
    {
        // intialize rowsizes with 0
        for (typename A12Type::size_type i = 0; i < A_12.N(); i++)
            A_12.setrowsize(i, 0);
        for (typename A21Type::size_type i = 0; i < A_21.N(); i++)
            A_21.setrowsize(i, 0);

        std::vector<bool> visitedFirst(A_12.N());
        for (unsigned int i = 0; i < A_12.N(); i++)
            visitedFirst[i] = false;
        std::vector<bool> visitedSecond(A_21.N());
        for (unsigned int i = 0; i < A_21.N(); i++)
            visitedSecond[i] = false;
        // loop 1 over all elements of the First grid to set the rowsizes of the coupling matrices
        FirstIterator endIt = (this->firstGrid()).template leafend<0>();
        for (FirstIterator firstIt = (this->firstGrid()).template leafbegin<0>(); firstIt != endIt; ++firstIt)
        {
            GeometryType gt = firstIt->geometry().type();
            const typename ReferenceElementContainer<double,dim>::value_type& referenceElement = ReferenceElements<double, dim>::general(gt);

            FirstIntersectionIterator endIsIt = firstIt->ileafend();
            for (FirstIntersectionIterator firstIsIt = firstIt->ileafbegin(); firstIsIt != endIsIt; ++ firstIsIt)
            {
                if (firstIsIt->neighbor())
                    continue; // only work on border entities

                // In the following, it is determined whether the intersection really is on the interface.
                // Every node of the corresponding face is checked whether it also belongs to the Second grid.
                int faceIdx = firstIsIt->numberInSelf();
                int numVerticesOfFace = referenceElement.size(faceIdx, 1, dim);
                int numVerticesInSecondGrid = 0;
                int secondIds[numVerticesOfFace];
                int firstIds[numVerticesOfFace];
                for (int nodeInFace = 0; nodeInFace < numVerticesOfFace; nodeInFace++)
                {
                    int nodeInElement = referenceElement.subEntity(faceIdx, 1, nodeInFace, dim);

                    // get the node pointer on the First grid
                    const FirstVPointer& firstVPointer = (*firstIt).template entity<dim>(nodeInElement);

                    // get the node pointer on the host grid
                    const HostVPointer& hostVPointer = (this->firstGrid()).template getHostEntity<dim>(*firstVPointer);

                    // check if the node is also part of the second grid
                    if (!((this->secondGrid()).template contains<dim>(hostVPointer)))
                        break; // otherwise this face is not at the interface

                    // get the index of the node with respect to the first matrix
                    firstIds[numVerticesInSecondGrid] = firstVertexMapper_.map(*firstVPointer);

                    // get the index of the node with respect to the Second matrix
                    secondIds[numVerticesInSecondGrid] = secondVertexMapper_.map(*firstVPointer);

                    numVerticesInSecondGrid++;
                }

                if (numVerticesInSecondGrid < numVerticesOfFace) // then this face is not at the interface
                    continue;

                for (int nodeInFace = 0; nodeInFace < numVerticesOfFace; nodeInFace++)
                {
                    int firstIndex = firstIds[nodeInFace];

                    if (!visitedFirst[firstIndex])
                    {
                        A_12.incrementrowsize(firstIndex);
                        visitedFirst[firstIndex] = true;
                    }

                    A_12.incrementrowsize(firstIndex, numVerticesOfFace - 1);

                    int secondIndex = secondIds[nodeInFace];

                    if (!visitedSecond[secondIndex])
                    {
                        A_21.incrementrowsize(secondIndex);
                        visitedSecond[secondIndex] = true;
                    }

                    A_21.incrementrowsize(secondIndex, numVerticesOfFace - 1);
                }
            }
        } // end loop 1 over all elements of the First grid
        A_12.endrowsizes();
        A_21.endrowsizes();


        // loop 2 over all elements of the First grid to set the indices of the nonzero entries of the coupling matrices
        for (FirstIterator firstIt = (this->firstGrid()).template leafbegin<0>(); firstIt != endIt; ++firstIt)
        {
            GeometryType gt = firstIt->geometry().type();
            const typename ReferenceElementContainer<double,dim>::value_type& referenceElement = ReferenceElements<double, dim>::general(gt);

            FirstIntersectionIterator endIsIt = firstIt->ileafend();
            for (FirstIntersectionIterator firstIsIt = firstIt->ileafbegin(); firstIsIt != endIsIt; ++ firstIsIt)
            {
                if (firstIsIt->neighbor())
                    continue; // only work on border entities

                // In the following, it is determined whether the intersection really is on the interface.
                // Every node of the corresponding face is checked whether it also belongs to the Second grid.
                int faceIdx = firstIsIt->numberInSelf();
                int numVerticesOfFace = referenceElement.size(faceIdx, 1, dim);
                int numVerticesInSecondGrid = 0;
                int secondIds[numVerticesOfFace];
                int firstIds[numVerticesOfFace];
                for (int nodeInFace = 0; nodeInFace < numVerticesOfFace; nodeInFace++)
                {
                    int nodeInElement = referenceElement.subEntity(faceIdx, 1, nodeInFace, dim);

                    // get the node pointer on the First grid
                    const FirstVPointer& firstVPointer = (*firstIt).template entity<dim>(nodeInElement);

                    // get the node pointer on the host grid
                    const HostVPointer& hostVPointer = (this->firstGrid()).template getHostEntity<dim>(*firstVPointer);

                    // check if the node is also part of the second grid
                    if (!((this->secondGrid()).template contains<dim>(hostVPointer)))
                        break; // otherwise this face is not at the interface

                    // get the index of the node with respect to the first matrix
                    firstIds[numVerticesInSecondGrid] = firstVertexMapper_.map(*firstVPointer);

                    // get the index of the node with respect to the Second matrix
                    secondIds[numVerticesInSecondGrid] = secondVertexMapper_.map(*firstVPointer);

                    numVerticesInSecondGrid++;
                }

                if (numVerticesInSecondGrid < numVerticesOfFace) // then this face is not at the interface
                    continue;

                for (int nodeInFaceFirst = 0; nodeInFaceFirst < numVerticesOfFace; nodeInFaceFirst++)
                {
                    int firstIndex = firstIds[nodeInFaceFirst];

                    for (int nodeInFaceSecond = 0; nodeInFaceSecond < numVerticesOfFace; nodeInFaceSecond++)
                    {
                        int secondIndex = secondIds[nodeInFaceSecond];

                        A_12.addindex(firstIndex, secondIndex);
                        A_21.addindex(secondIndex, firstIndex);
                    }
                }
            }
        } // end loop 2 over all elements of the First grid
        A_12.endindices();
        A_21.endindices();
        A_12 = 0;
        A_21 = 0;

        // loop 3 over all elements of the First grid to actually fill the coupling matrices
        for (FirstIterator firstIt = (this->firstGrid()).template leafbegin<0>(); firstIt != endIt; ++firstIt)
        {
            GeometryType gt = firstIt->geometry().type();
            const typename ReferenceElementContainer<double,dim>::value_type& referenceElement = ReferenceElements<double, dim>::general(gt);

            FirstIntersectionIterator endIsIt = firstIt->ileafend();
            for (FirstIntersectionIterator firstIsIt = firstIt->ileafbegin(); firstIsIt != endIsIt; ++ firstIsIt)
            {
                if (firstIsIt->neighbor())
                    continue; // only work on border entities

                // In the following, it is determined whether the intersection really is on the interface.
                // Every node of the corresponding face is checked whether it also belongs to the Second grid.
                int faceIdx = firstIsIt->numberInSelf();
                int numVerticesOfFace = referenceElement.size(faceIdx, 1, dim);
                int numVerticesInSecondGrid = 0;
                int secondIds[numVerticesOfFace];
                int firstIds[numVerticesOfFace];
                for (int nodeInFace = 0; nodeInFace < numVerticesOfFace; nodeInFace++)
                {
                    int nodeInElement = referenceElement.subEntity(faceIdx, 1, nodeInFace, dim);

                    // get the node pointer on the First grid
                    const FirstVPointer& firstVPointer = (*firstIt).template entity<dim>(nodeInElement);

                    // get the node pointer on the host grid
                    const HostVPointer& hostVPointer = (this->firstGrid()).template getHostEntity<dim>(*firstVPointer);

                    // check if the node is also part of the second grid
                    if (!((this->secondGrid()).template contains<dim>(hostVPointer)))
                        break; // otherwise this face is not at the interface

                    // get the index of the node with respect to the first matrix
                    firstIds[numVerticesInSecondGrid] = firstVertexMapper_.map(*firstVPointer);

                    // get the index of the node with respect to the Second matrix
                    secondIds[numVerticesInSecondGrid] = secondVertexMapper_.map(*firstVPointer);

                    numVerticesInSecondGrid++;
                }

                if (numVerticesInSecondGrid < numVerticesOfFace) // then this face is not at the interface
                    continue;

                FVElementGeometry<FirstGrid> firstFVGeom;
                firstFVGeom.update(*firstIt);

                // get the geometry type of the face
                GeometryType geomTypeBoundary = firstIsIt->intersectionSelfLocal().type();

                // The coupling is realized in the FV way.
                // The unknown First normal velocity is evaluated in the center of the subcontrolvolume face.
                const FieldVector<double,dim-1>& faceLocalDimM1 = ReferenceElements<double,dim-1>::general(geomTypeBoundary).position(0,0);
                FieldVector<double,dim> normal = firstIsIt->unitOuterNormal(faceLocalDimM1);

                for (int nodeInFace = 0; nodeInFace < numVerticesOfFace; ++nodeInFace) // loop over interface nodes
                {
                    int bfIdx = firstFVGeom.boundaryFaceIndex(faceIdx, nodeInFace);
                    FieldVector<double,dim> firstQLocal = firstFVGeom.boundaryFace[bfIdx].ipLocal;
                    FieldVector<double,dim> qGlobal = firstFVGeom.boundaryFace[bfIdx].ipGlobal;

                    int firstIndex = firstIds[nodeInFace];
                    int secondIndex = secondIds[nodeInFace];

                    FieldVector<double, firstNumEq> firstSol = (this->firstModel_.sol())[firstIndex];
                    FieldVector<double, secondNumEq> secondSol = (this->secondModel_.sol())[secondIndex];

                    // calculate the entries of A_12
                    for (int j = 0; j < secondNumEq; j++)
                    {
                        double eps = std::max(fabs(1e-5*secondSol[j]), 1e-5);

                        secondSol[j] += eps;
                        FieldVector<double, firstNumEq> c12PlusEps(0);
                        localCoupling12(firstSol, secondSol, firstIndex, secondIndex, qGlobal, c12PlusEps);

                        secondSol[j] -= 2.0*eps;
                        FieldVector<double, firstNumEq> c12MinusEps(0);
                        localCoupling12(firstSol, secondSol, firstIndex, secondIndex, qGlobal, c12MinusEps);

                        for (int i = 0; i < firstNumEq; i++)
                            A_12[firstIndex][secondIndex][i][j] = 0.5/eps*(c12PlusEps[i] - c12MinusEps[i]);

                        secondSol[j] += eps;
                    }

                    // calculate the entries of A_21
                    for (int j = 0; j < firstNumEq; j++)
                    {
                        double eps = std::max(fabs(1e-5*firstSol[j]), 1e-5);

                        firstSol[j] += eps;
                        FieldVector<double, secondNumEq> c21PlusEps(0);
                        localCoupling21(firstSol, secondSol, firstIndex, secondIndex, qGlobal, c21PlusEps);

                        firstSol[j] -= 2.0*eps;
                        FieldVector<double, secondNumEq> c21MinusEps(0);
                        localCoupling21(firstSol, secondSol, firstIndex, secondIndex, qGlobal, c21MinusEps);

                        for (int i = 0; i < secondNumEq; i++)
                            A_21[secondIndex][firstIndex][i][j] = 0.5/eps*(c21PlusEps[i] - c21MinusEps[i]);

                        firstSol[j] += eps;
                    }
                }
            }
        } // end loop 3 over all elements of the First grid
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

        //        int rowsInBlock1 = BaseType::FirstMatrixType::block_type::rows;
        //        int nOfBlockRows1 = (this->firstModel_).matrix().N();
        //        SecondIterator dummyIT2((this->secondGrid()).template leafbegin<0>());
        //        SecondIntersectionIterator dummyIS2(IntersectionIteratorGetter<SecondGrid,LeafTag>::begin(*dummyIT2));
        //        SecondVIterator endItV2 = (this->secondGrid()).template leafend<dim>();
        //        for (SecondVIterator it = (this->secondGrid()).template leafbegin<dim>(); it != endItV2; ++it)
        //        {
        //            FieldVector<double,dim> globalCoord = (*it).geometry().corner(0);
        //            BoundaryConditions::Flags bctype = (this->secondModel_).problem.bctype(globalCoord, *dummyIT2, dummyIS2, globalCoord);
        //            int secondId = secondVertexMapper_.map(*it);
        //            this->u[rowsInBlock1*nOfBlockRows1 + secondId] = -this->u[rowsInBlock1*nOfBlockRows1 + secondId];
        //
        //            if (bctype == BoundaryConditions::dirichlet)
        //            {
        //                double dirichletBC = (this->secondModel_).problem.g(globalCoord, *dummyIT2, dummyIS2, globalCoord);
        //                this->u[rowsInBlock1*nOfBlockRows1 + secondId] += dirichletBC;
        //            }
        //        }

        // transfer to local solution vectors
        const typename BaseType::FirstMatrixType::block_type::size_type colsInBlock1 = BaseType::FirstMatrixType::block_type::cols;
        const typename BaseType::SecondMatrixType::block_type::size_type colsInBlock2 = BaseType::SecondMatrixType::block_type::cols;
        for (unsigned int i = 0; i < this->firstModel_.sol().size(); i++)
            for (typename BaseType::FirstMatrixType::block_type::size_type k = 0; k < colsInBlock1; k++)
                this->firstModel_.sol()[i][k] -= this->u[i*colsInBlock1 + k];
        for (unsigned int i = 0; i < this->secondModel_.sol().size(); i++)
            for (typename BaseType::SecondMatrixType::block_type::size_type k = 0; k < colsInBlock2; k++)
                this->secondModel_.sol()[i][k] -= this->u[colsInBlock1*this->firstModel_.sol().size() + i*colsInBlock2 + k];
    }

    virtual void vtkout (const char* name, int k)
    {
        char localName[128];
        sprintf(localName, "%s_first", name);
        this->firstModel_.vtkout(localName, k);
        sprintf(localName, "%s_second", name);
        this->secondModel_.vtkout(localName, k);
    }

    CoupledBox(const FirstGrid& firstGrid, FirstModel& firstModel,
            const SecondGrid& secondGrid, SecondModel& secondModel,
            bool assembleGlobalSystem)
    : BaseType(firstGrid, firstModel, secondGrid, secondModel, assembleGlobalSystem),
    firstVertexMapper_(firstGrid, firstGrid.leafIndexSet()),
    secondVertexMapper_(secondGrid, secondGrid.leafIndexSet())
    {}

private:
    FirstVertexMapper firstVertexMapper_;
    SecondVertexMapper secondVertexMapper_;
};

}

#endif
