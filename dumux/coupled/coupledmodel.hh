// $Id$

#ifndef DUNE_COUPLEDMODEL_HH
#define DUNE_COUPLEDMODEL_HH

#include "dumux/pardiso/pardiso.hh"

namespace Dune
{
/** \todo Please doc me! */

template<class FirstModel, class SecondModel, class Imp>
class CoupledModel {
public:
    typedef typename FirstModel::GridType FirstGrid;
    typedef typename FirstModel::MatrixType FirstMatrixType;
    typedef typename FirstMatrixType::RowIterator FirstRowIterator;
    typedef typename FirstMatrixType::ColIterator FirstColIterator;

    typedef typename SecondModel::GridType SecondGrid;
    typedef typename SecondModel::MatrixType SecondMatrixType;
    typedef typename SecondMatrixType::RowIterator SecondRowIterator;
    typedef typename SecondMatrixType::ColIterator SecondColIterator;

    typedef FieldMatrix<double,1,1> GlobalMBlockType;
    typedef BCRSMatrix<GlobalMBlockType> GlobalMatrixType;
    typedef FieldVector<double,1> GlobalVBlockType;
    typedef BlockVector<GlobalVBlockType> GlobalVectorType;

    virtual void initial()
    {
        firstModel_.initial();
        secondModel_.initial();

        if (assembleGlobalSystem_)
        {
            //            firstModel_.assemble();
            //            secondModel_.assemble();

            const typename FirstMatrixType::block_type::size_type rowsInBlock1 = FirstMatrixType::block_type::rows;
            const typename FirstMatrixType::block_type::size_type colsInBlock1 = FirstMatrixType::block_type::cols;
            const typename SecondMatrixType::block_type::size_type rowsInBlock2 = SecondMatrixType::block_type::rows;
            const typename SecondMatrixType::block_type::size_type colsInBlock2 = SecondMatrixType::block_type::cols;
            typedef FieldMatrix<double,rowsInBlock1,colsInBlock2> BlockType12;
            typedef BCRSMatrix<BlockType12> A12Type;
            typedef FieldMatrix<double,rowsInBlock2,colsInBlock1> BlockType21;
            typedef BCRSMatrix<BlockType21> A21Type;
            typename FirstMatrixType::size_type nOfBlockCols1 = firstModel_.matrix().M();
            typename FirstMatrixType::size_type nOfBlockRows1 = firstModel_.matrix().N();
            typename SecondMatrixType::size_type nOfBlockCols2 = secondModel_.matrix().M();
            typename SecondMatrixType::size_type nOfBlockRows2 = secondModel_.matrix().N();
            A12Type A12(nOfBlockRows1, nOfBlockCols2, A12Type::random);
            A21Type A21(nOfBlockRows2, nOfBlockCols1, A21Type::random);
            assembleCoupling(A12, A21);

            A.setBuildMode(GlobalMatrixType::random);
            typename GlobalMatrixType::size_type colsGlobal, rowsGlobal;
            colsGlobal = colsInBlock1*nOfBlockCols1 + colsInBlock2*nOfBlockCols2;
            rowsGlobal = rowsInBlock1*nOfBlockRows1 + rowsInBlock2*nOfBlockRows2;
            A.setSize(rowsGlobal, colsGlobal);
            u.resize(colsGlobal);
            f.resize(rowsGlobal);

            FirstRowIterator endIBlock1 = A11.end();
            for (FirstRowIterator iBlock = A11.begin(); iBlock != endIBlock1; ++iBlock)
            {
                for (typename FirstMatrixType::block_type::size_type iLocal = 0; iLocal < rowsInBlock1; iLocal++)
                {
                    typename GlobalMatrixType::size_type rowSize = colsInBlock1*iBlock->size() + colsInBlock2*A12.getrowsize(iBlock.index());
                    A.setrowsize(iBlock.index()*rowsInBlock1 + iLocal, rowSize);
                }
            }
            SecondRowIterator endIBlock2 = A22.end();
            for (SecondRowIterator iBlock = A22.begin(); iBlock != endIBlock2; ++iBlock)
            {
                for (typename SecondMatrixType::block_type::size_type iLocal = 0; iLocal < rowsInBlock2; iLocal++)
                {
                    typename GlobalMatrixType::size_type rowSize = colsInBlock2*iBlock->size() + colsInBlock1*A21.getrowsize(iBlock.index());
                    A.setrowsize(rowsInBlock1*nOfBlockRows1 + iBlock.index()*rowsInBlock2 + iLocal, rowSize);
                }
            }
            A.endrowsizes();

            endIBlock1 = A11.end();
            for (FirstRowIterator iBlock = A11.begin(); iBlock != endIBlock1; ++iBlock)
            {
                FirstColIterator endJBlock = iBlock->end();
                for (FirstColIterator jBlock = iBlock->begin(); jBlock != endJBlock; ++jBlock)
                {
                    for (typename FirstMatrixType::block_type::size_type iLocal = 0; iLocal < rowsInBlock1; iLocal++)
                    {
                        typename GlobalMatrixType::size_type globalRow = iBlock.index()*rowsInBlock1 + iLocal;
                        for (typename FirstMatrixType::block_type::size_type jLocal = 0; jLocal < colsInBlock1; jLocal++)
                        {
                            typename GlobalMatrixType::size_type globalCol = jBlock.index()*colsInBlock1 + jLocal;
                            A.addindex(globalRow, globalCol);
                            u[globalCol] = (firstModel_.sol())[jBlock.index()][jLocal];
                        }
                    }
                }
            }

            typename A12Type::RowIterator endIBlock12 = A12.end();
            for (typename A12Type::RowIterator iBlock = A12.begin(); iBlock != endIBlock12; ++iBlock)
            {
                typename A12Type::ColIterator endJBlock = iBlock->end();
                for (typename A12Type::ColIterator jBlock = iBlock->begin(); jBlock != endJBlock; ++jBlock)
                {
                    for (typename A12Type::block_type::size_type iLocal = 0; iLocal < rowsInBlock1; iLocal++)
                    {
                        typename GlobalMatrixType::size_type globalRow = iBlock.index()*rowsInBlock1 + iLocal;
                        for (typename A12Type::block_type::size_type jLocal = 0; jLocal < colsInBlock2; jLocal++)
                        {
                            typename GlobalMatrixType::size_type globalCol = colsInBlock1*nOfBlockCols1 + jBlock.index()*colsInBlock2 + jLocal;
                            A.addindex(globalRow, globalCol);
                        }
                    }
                }
            }
            endIBlock2 = A22.end();
            for (SecondRowIterator iBlock = A22.begin(); iBlock != endIBlock2; ++iBlock)
            {
                SecondColIterator endJBlock = iBlock->end();
                for (SecondColIterator jBlock = iBlock->begin(); jBlock != endJBlock; ++jBlock)
                {
                    for (typename SecondMatrixType::block_type::size_type iLocal = 0; iLocal < rowsInBlock2; iLocal++)
                    {
                        typename GlobalMatrixType::size_type globalRow = rowsInBlock1*nOfBlockRows1 + iBlock.index()*rowsInBlock2 + iLocal;
                        for (typename SecondMatrixType::block_type::size_type jLocal = 0; jLocal < colsInBlock2; jLocal++)
                        {
                            typename GlobalMatrixType::size_type globalCol = colsInBlock1*nOfBlockCols1 + jBlock.index()*colsInBlock2 + jLocal;
                            A.addindex(globalRow, globalCol);
                            u[globalCol] = (secondModel_.sol())[jBlock.index()][jLocal];
                        }
                    }
                }
            }
            typename A21Type::RowIterator endIBlock21 = A21.end();
            for (typename A21Type::RowIterator iBlock = A21.begin(); iBlock != endIBlock21; ++iBlock)
            {
                typename A21Type::ColIterator endJBlock = iBlock->end();
                for (typename A21Type::ColIterator jBlock = iBlock->begin(); jBlock != endJBlock; ++jBlock)
                {
                    for (typename A21Type::block_type::size_type iLocal = 0; iLocal < rowsInBlock2; iLocal++)
                    {
                        typename GlobalMatrixType::size_type globalRow = rowsInBlock1*nOfBlockRows1 + iBlock.index()*rowsInBlock2 + iLocal;
                        for (typename A21Type::block_type::size_type jLocal = 0; jLocal < colsInBlock1; jLocal++)
                        {
                            typename GlobalMatrixType::size_type globalCol = jBlock.index()*colsInBlock1 + jLocal;
                            A.addindex(globalRow, globalCol);
                        }
                    }
                }
            }
            A.endindices();

        }

        uOld = u;
    }

    template <class A12Type, class A21Type>
    void assembleCoupling(A12Type& A12, A21Type& A21)
    {
        getImp().template assembleCoupling<A12Type,A21Type>(A12, A21);
    }

    virtual void assemble()
    {
        firstModel_.solOldNewtonStep() = firstModel_.sol();
        secondModel_.solOldNewtonStep() = secondModel_.sol();
        firstModel_.assemble();
        secondModel_.assemble();
        firstModel_.sol() = firstModel_.solOldNewtonStep();
        secondModel_.sol() = secondModel_.solOldNewtonStep();

        const typename FirstMatrixType::block_type::size_type rowsInBlock1 = FirstMatrixType::block_type::rows;
        const typename FirstMatrixType::block_type::size_type colsInBlock1 = FirstMatrixType::block_type::cols;
        const typename SecondMatrixType::block_type::size_type rowsInBlock2 = SecondMatrixType::block_type::rows;
        const typename SecondMatrixType::block_type::size_type colsInBlock2 = SecondMatrixType::block_type::cols;
        typedef FieldMatrix<double,rowsInBlock1,colsInBlock2> BlockType12;
        typedef BCRSMatrix<BlockType12> A12Type;
        typedef FieldMatrix<double,rowsInBlock2,colsInBlock1> BlockType21;
        typedef BCRSMatrix<BlockType21> A21Type;

        typename FirstMatrixType::size_type nOfBlockCols1 = firstModel_.matrix().M();
        typename FirstMatrixType::size_type nOfBlockRows1 = firstModel_.matrix().N();
        typename SecondMatrixType::size_type nOfBlockCols2 = secondModel_.matrix().M();
        typename SecondMatrixType::size_type nOfBlockRows2 = secondModel_.matrix().N();
        A12Type A12(nOfBlockRows1, nOfBlockCols2, A12Type::random);
        A21Type A21(nOfBlockRows2, nOfBlockCols1, A21Type::random);
        assembleCoupling(A12, A21);

        if (assembleGlobalSystem_)
        {
            FirstRowIterator endIBlock1 = A11.end();
            for (FirstRowIterator iBlock = A11.begin(); iBlock != endIBlock1; ++iBlock)
            {
                FirstColIterator endJBlock = iBlock->end();
                for (FirstColIterator jBlock = iBlock->begin(); jBlock != endJBlock; ++jBlock)
                {
                    for (typename FirstMatrixType::block_type::size_type iLocal = 0; iLocal < rowsInBlock1; iLocal++)
                    {
                        typename GlobalMatrixType::size_type globalRow = iBlock.index()*rowsInBlock1 + iLocal;
                        f[globalRow] = (firstModel_.rhs())[iBlock.index()][iLocal];
                        for (typename FirstMatrixType::block_type::size_type jLocal = 0; jLocal < colsInBlock1; jLocal++)
                        {
                            typename GlobalMatrixType::size_type globalCol = jBlock.index()*colsInBlock1 + jLocal;
                            A[globalRow][globalCol] = (*jBlock)[iLocal][jLocal];
                            u[globalCol] = (firstModel_.sol())[jBlock.index()][jLocal];
                        }
                    }
                }
            }
            typename A12Type::RowIterator endIBlock12 = A12.end();
            for (typename A12Type::RowIterator iBlock = A12.begin(); iBlock != endIBlock12; ++iBlock)
            {
                typename A12Type::ColIterator endJBlock = iBlock->end();
                for (typename A12Type::ColIterator jBlock = iBlock->begin(); jBlock != endJBlock; ++jBlock)
                {
                    for (typename A12Type::block_type::size_type iLocal = 0; iLocal < rowsInBlock1; iLocal++)
                    {
                        typename GlobalMatrixType::size_type globalRow = iBlock.index()*rowsInBlock1 + iLocal;
                        for (typename A12Type::block_type::size_type jLocal = 0; jLocal < colsInBlock2; jLocal++)
                        {
                            typename GlobalMatrixType::size_type globalCol = colsInBlock1*nOfBlockCols1 + jBlock.index()*colsInBlock2 + jLocal;
                            A[globalRow][globalCol] = (*jBlock)[iLocal][jLocal];
                        }
                    }
                }
            }
            SecondRowIterator endIBlock2 = A22.end();
            for (SecondRowIterator iBlock = A22.begin(); iBlock != endIBlock2; ++iBlock)
            {
                SecondColIterator endJBlock = iBlock->end();
                for (SecondColIterator jBlock = iBlock->begin(); jBlock != endJBlock; ++jBlock)
                {
                    for (typename SecondMatrixType::block_type::size_type iLocal = 0; iLocal < rowsInBlock2; iLocal++)
                    {
                        typename GlobalMatrixType::size_type globalRow = rowsInBlock1*nOfBlockRows1 + iBlock.index()*rowsInBlock2 + iLocal;
                        f[globalRow] = (secondModel_.rhs())[iBlock.index()][iLocal];
                        for (typename SecondMatrixType::block_type::size_type jLocal = 0; jLocal < colsInBlock2; jLocal++)
                        {
                            typename GlobalMatrixType::size_type globalCol = colsInBlock1*nOfBlockCols1 + jBlock.index()*colsInBlock2 + jLocal;
                            A[globalRow][globalCol] = (*jBlock)[iLocal][jLocal];
                            u[globalCol] = (secondModel_.sol())[jBlock.index()][jLocal];
                        }
                    }
                }
            }
            typename A21Type::RowIterator endIBlock21 = A21.end();
            for (typename A21Type::RowIterator iBlock = A21.begin(); iBlock != endIBlock21; ++iBlock)
            {
                typename A21Type::ColIterator endJBlock = iBlock->end();
                for (typename A21Type::ColIterator jBlock = iBlock->begin(); jBlock != endJBlock; ++jBlock)
                {
                    for (typename A21Type::block_type::size_type iLocal = 0; iLocal < rowsInBlock2; iLocal++)
                    {
                        typename GlobalMatrixType::size_type globalRow = rowsInBlock1*nOfBlockRows1 + iBlock.index()*rowsInBlock2 + iLocal;
                        for (typename A21Type::block_type::size_type jLocal = 0; jLocal < colsInBlock1; jLocal++)
                        {
                            typename GlobalMatrixType::size_type globalCol = jBlock.index()*colsInBlock1 + jLocal;
                            A[globalRow][globalCol] = (*jBlock)[iLocal][jLocal];
                        }
                    }
                }
            }
        }
    }

    virtual void solve()
    {
        typedef MatrixAdapter<GlobalMatrixType,GlobalVectorType,GlobalVectorType> Operator;
        Operator op(A);
        double red=1E-14;
        SeqILU0<GlobalMatrixType,GlobalVectorType,GlobalVectorType> ilu0(A,1.0);
        BiCGSTABSolver<GlobalVectorType> solver(op,ilu0,red,10000,1);
        //         SeqPardiso<GlobalMatrixType,GlobalVectorType,GlobalVectorType> pardiso(A);
        //         LoopSolver<GlobalVectorType> solver(op,pardiso,red,10000,1);
        InverseOperatorResult r;
        solver.apply(u, f, r);

        // transfer to local solution vectors
        const typename FirstMatrixType::block_type::size_type colsInBlock1 = FirstMatrixType::block_type::cols;
        const typename SecondMatrixType::block_type::size_type colsInBlock2 = SecondMatrixType::block_type::cols;
        for (unsigned int i = 0; i < firstModel_.sol().size(); i++)
            for (typename FirstMatrixType::block_type::size_type k = 0; k < colsInBlock1; k++)
                firstModel_.sol()[i][k] = u[i*colsInBlock1 + k];
        for (unsigned int i = 0; i < secondModel_.sol().size(); i++)
            for (typename SecondMatrixType::block_type::size_type k = 0; k < colsInBlock2; k++)
                secondModel_.sol()[i][k] = u[colsInBlock1*firstModel_.sol().size() + i*colsInBlock2 + k];
    }

    const FirstGrid& firstGrid()
    {
        return firstGrid_;
    }

    const SecondGrid& secondGrid()
    {
        return secondGrid_;
    }

    FirstModel& firstModel()
    {
        return firstModel_;
    }

    SecondModel& secondModel()
    {
        return secondModel_;
    }

    GlobalMatrixType& matrix()
    {
        return A;
    }

    GlobalVectorType& rhs()
    {
        return f;
    }

    GlobalVectorType& sol()
    {
        return u;
    }

    GlobalVectorType& solOld()
    {
        return uOld;
    }

    virtual void vtkout (const char* name, int k)
    {
        char localName[128];
        sprintf(localName, "%s_firstModel", name);
        firstModel_.vtkout(localName, k);
        sprintf(localName, "%s_secondModel", name);
        secondModel_.vtkout(localName, k);
    }

    CoupledModel(const FirstGrid& firstGrid, FirstModel& firstModel,
                 const SecondGrid& secondGrid, SecondModel& secondModel,
                 bool assembleGlobalSystem)
        : firstGrid_(firstGrid), secondGrid_(secondGrid), firstModel_(firstModel), secondModel_(secondModel),
          A11(firstModel_.matrix()), A22(secondModel_.matrix()), assembleGlobalSystem_(assembleGlobalSystem)
    {}

protected:
    const FirstGrid& firstGrid_;
    const SecondGrid& secondGrid_;
    FirstModel& firstModel_;
    SecondModel& secondModel_;
    FirstMatrixType& A11;
    SecondMatrixType& A22;
    bool assembleGlobalSystem_;
    GlobalMatrixType A;
    GlobalVectorType u;
    GlobalVectorType f;
    GlobalVectorType uOld;

    Imp& getImp()
    {
        return static_cast<Imp&>(*this);
    }
    const Imp& getImp() const
    {
        return static_cast<const Imp&>(*this);
    }
};

}

#endif
