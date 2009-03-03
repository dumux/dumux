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
    typedef typename SecondGrid::template Codim<0>::EntityPointer SecondElementPointer;
    typedef typename SecondGrid::template Codim<dim>::EntityPointer SecondVPointer;
    typedef typename SecondGV::template Codim<dim>::Iterator SecondVIterator;
    typedef typename IntersectionIteratorGetter<SecondGrid,LeafTag>::IntersectionIterator SecondIntersectionIterator;
    typedef MultipleCodimMultipleGeomTypeMapper<SecondGrid,SecondIS,NodeLayout> SecondVertexMapper;

    typedef typename SecondGrid::HostGridType HostGrid;
    typedef typename HostGrid::template Codim<0>::Entity HostElement;
    typedef typename HostGrid::template Codim<0>::EntityPointer HostPointer;
    typedef typename HostGrid::template Codim<dim>::EntityPointer HostVPointer;
    typedef typename IntersectionIteratorGetter<HostGrid,LeafTag>::IntersectionIterator HostIntersectionIterator;

    template <class FirstFV, class SecondFV>
    void localCoupling12(FirstFV& firstSol, SecondFV& secondSol, int firstIndex, int secondIndex,
                         FieldVector<double, dim> qGlobal, FieldVector<double, dim> normal, FirstFV& result)
    {
        this->getImp().template localCoupling12<FirstFV,SecondFV>(firstSol, secondSol,
                                                                  firstIndex, secondIndex, qGlobal, normal, result);
    }

    template <class FV, class FVGrad, class ElementT>
    void localBoundaryDefect1(FV& firstSol, FVGrad& firstSolGrad, int firstIndex,
                              FieldVector<double, dim> qGlobal, const ElementT& element, FieldVector<double, dim> qLocal,
                              FieldVector<double, dim> normal, FV& result)
    {
        this->getImp().template localBoundaryDefect1<FV,FVGrad>(firstSol, firstSolGrad,
                                                                firstIndex, qGlobal, element, qLocal, normal, result);
    }

    template <class FirstFV, class SecondFV>
    void localCoupling21(FirstFV& firstSol, SecondFV& secondSol, int firstIndex, int secondIndex,
                         FieldVector<double, dim> qGlobal, FieldVector<double, dim> normal, SecondFV& result)
    {
        this->getImp().template localCoupling21<FirstFV,SecondFV>(firstSol, secondSol,
                                                                  firstIndex, secondIndex, qGlobal, normal, result);
    }

    template <class FirstFVVec, class FirstFV>
    void calculateFirstSolAtQ(const FirstFVVec& firstSol, int subCVF, FirstFV& firstSolAtQ)
    {
        if (subCVF == 0)
        {
            firstSolAtQ = firstSol[0];
            //            firstSolAtQ *= 3.0;
            //            firstSolAtQ += firstSol[1];
            //            firstSolAtQ *= 0.25;
        }
        else
        {
            firstSolAtQ = firstSol[1];
            //            firstSolAtQ *= 3.0;
            //            firstSolAtQ += firstSol[0];
            //            firstSolAtQ *= 0.25;
        }
    }

    template <class Element, class FVGeom, class Gradient>
    void calculateFirstSolGradientAtQ(const Element& element, const FVGeom& fvGeom, int bfIdx, Gradient& result)
    {
        result = 0;

        for (int vert = 0; vert < fvGeom.numVertices; vert++)
        {
            FieldVector<double, dim> shapeGrad = fvGeom.boundaryFace[bfIdx].grad[vert];

            const FirstVPointer& vertexPointer = element.template entity<dim>(vert);

            int globalId = firstVertexMapper_.map(*vertexPointer);

            FieldVector<double, firstNumEq> vertexSol = this->firstModel().sol()[globalId];

            //std::cout << "vert " << vert << ": sol " << vertexSol << ", shapeGrad = " << shapeGrad << std::endl;

            for (int equation = 0; equation < firstNumEq; equation++)
            {
                FieldVector<double, dim> gradient = shapeGrad;
                gradient *= vertexSol[equation];

                result[equation] += gradient;
            }
        }

        return;
    }

    template <class SecondFVVec, class SecondFV>
    void calculateSecondSolAtQ(const SecondFVVec& secondSol, int subCVF, SecondFV& secondSolAtQ)
    {
        if (subCVF == 0)
        {
            secondSolAtQ = secondSol[0];
            //            secondSolAtQ *= 3.0;
            //            secondSolAtQ += secondSol[1];
            //            secondSolAtQ *= 0.25;
        }
        else
        {
            secondSolAtQ = secondSol[1];
            //            secondSolAtQ *= 3.0;
            //            secondSolAtQ += secondSol[0];
            //            secondSolAtQ *= 0.25;
        }
    }

    template <class A12Type, class A21Type>
    void assembleCoupling(A12Type& A_12, A21Type& A_21)
    {
        count ++;
        //std::cout << "count = " << count << std::endl;
        //printvector(std::cout, this->firstModel().sol(), "Stokes solution", "row", 3, 1, 3);
        //printvector(std::cout, this->secondModel().sol(), "Darcy solution", "row", 100, 1, 3);

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

                std::vector<FieldVector<double, firstNumEq> > firstSol(numVerticesOfFace);
                std::vector<FieldVector<double, secondNumEq> > secondSol(numVerticesOfFace);
                for (int nodeInFace = 0; nodeInFace < numVerticesOfFace; ++nodeInFace) // loop over interface nodes
                {
                    firstSol[nodeInFace] = this->firstModel().sol()[firstIds[nodeInFace]];
                    secondSol[nodeInFace] = this->secondModel().sol()[secondIds[nodeInFace]];
                }

                for (int subCVF = 0; subCVF < numVerticesOfFace; ++subCVF) // loop over interface nodes
                {
                    int bfIdx = firstFVGeom.boundaryFaceIndex(faceIdx, subCVF);
                    FieldVector<double,dim> firstQLocal = firstFVGeom.boundaryFace[bfIdx].ipLocal;
                    FieldVector<double,dim> qGlobal = firstFVGeom.boundaryFace[bfIdx].ipGlobal;
                    FieldVector<double,dim> normal = firstFVGeom.boundaryFace[bfIdx].normal;

                    int nodeInElement = referenceElement.subEntity(faceIdx, 1, subCVF, dim);
                    FieldVector<double, dim> nodeGlobal = firstIt->geometry().corner(nodeInElement);
                    FieldVector<double, dim> nodeLocal = firstIt->geometry().local(nodeGlobal);

                    // WARNING: ASSUMES scalar BoundaryConditions flags for firstModel
                    typename BoundaryConditions::Flags bctypeNode = this->firstModel().problem.bctype(nodeGlobal, *firstIt, firstIsIt, nodeLocal);

                    // WARNING: ASSUMES that Dirichlet nodes are always on both sides
                    if (bctypeNode == BoundaryConditions::neumann)
                    {
                        // add coupling defect to rhs:
                        int firstIndex = firstIds[subCVF];
                        int secondIndex = secondIds[subCVF];
                        FieldVector<double, firstNumEq>  firstSolAtQ(0);
                        FieldMatrix<double, firstNumEq, dim>  firstSolGradientAtQ(0);
                        FieldVector<double, secondNumEq>  secondSolAtQ(0);
                        calculateFirstSolAtQ(firstSol, subCVF, firstSolAtQ);
                        calculateFirstSolGradientAtQ(*firstIt, firstFVGeom, bfIdx, firstSolGradientAtQ);
                        calculateSecondSolAtQ(secondSol, subCVF, secondSolAtQ);
                        FieldVector<double, firstNumEq> boundaryDefect1(0);
                        //                        if (count < 3)
                        //                            localBoundaryDefect1(firstSolAtQ, firstSolGradientAtQ, firstIndex, qGlobal, *firstIt, firstQLocal, normal, boundaryDefect1);
                        //                        else
                        localCoupling12(firstSolAtQ, secondSolAtQ, firstIndex, secondIndex, qGlobal, normal, boundaryDefect1);
                        (this->firstModel().rhs())[firstIndex] += boundaryDefect1;
                        FieldVector<double, secondNumEq> boundaryDefect2(0);
                        localCoupling21(firstSolAtQ, secondSolAtQ, firstIndex, secondIndex, qGlobal, normal, boundaryDefect2);
                        (this->secondModel().rhs())[secondIndex] += boundaryDefect2;

                        // calculate the entries of A_12
                        firstIndex = firstIds[subCVF];
                        for (int secondNode = 0; secondNode < numVerticesOfFace; ++secondNode)
                        {
                            secondIndex = secondIds[secondNode];

                            for (int j = 0; j < secondNumEq; j++)
                            {
                                double eps = std::max(fabs(1e-5*secondSol[secondNode][j]), 1e-5);

                                secondSol[secondNode][j] += eps;
                                calculateSecondSolAtQ(secondSol, subCVF, secondSolAtQ);

                                FieldVector<double, firstNumEq> coupling12PlusEps(0);
                                localCoupling12(firstSolAtQ, secondSolAtQ, firstIndex, secondIndex, qGlobal, normal, coupling12PlusEps);

                                secondSol[secondNode][j] -= 2.0*eps;
                                calculateSecondSolAtQ(secondSol, subCVF, secondSolAtQ);

                                FieldVector<double, firstNumEq> coupling12MinusEps(0);
                                localCoupling12(firstSolAtQ, secondSolAtQ, firstIndex, secondIndex, qGlobal, normal, coupling12MinusEps);

                                for (int i = 0; i < firstNumEq; i++)
                                    A_12[firstIndex][secondIndex][i][j] += 0.5/eps*(coupling12PlusEps[i] - coupling12MinusEps[i]);

                                secondSol[secondNode][j] += eps;
                            }
                        }

                        // calculate the entries of A_21
                        secondIndex = secondIds[subCVF];
                        firstSolAtQ = 0;
                        secondSolAtQ = 0;
                        calculateSecondSolAtQ(secondSol, subCVF, secondSolAtQ);
                        for (int firstNode = 0; firstNode < numVerticesOfFace; ++firstNode)
                        {
                            firstIndex = firstIds[firstNode];

                            for (int j = 0; j < firstNumEq; j++)
                            {
                                double eps = std::max(fabs(1e-5*firstSol[firstNode][j]), 1e-5);

                                firstSol[firstNode][j] += eps;
                                calculateFirstSolAtQ(firstSol, subCVF, firstSolAtQ);

                                FieldVector<double, secondNumEq> coupling21PlusEps(0);
                                localCoupling21(firstSolAtQ, secondSolAtQ, firstIndex, secondIndex, qGlobal, normal, coupling21PlusEps);

                                firstSol[firstNode][j] -= 2.0*eps;
                                calculateFirstSolAtQ(firstSol, subCVF, firstSolAtQ);

                                FieldVector<double, secondNumEq> coupling21MinusEps(0);
                                localCoupling21(firstSolAtQ, secondSolAtQ, firstIndex, secondIndex, qGlobal, normal, coupling21MinusEps);

                                for (int i = 0; i < secondNumEq; i++)
                                    A_21[secondIndex][firstIndex][i][j] += 0.5/eps*(coupling21PlusEps[i] - coupling21MinusEps[i]);

                                firstSol[firstNode][j] += eps;
                            }
                        }
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
        InverseOperatorResult r;
        solver.apply(this->u, this->f, r);

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
          secondVertexMapper_(secondGrid, secondGrid.leafIndexSet()), count(0)
    {}

private:
    FirstVertexMapper firstVertexMapper_;
    SecondVertexMapper secondVertexMapper_;
    int count;
};

}

#endif
