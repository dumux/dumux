// $Id$

#ifndef DUNE_COUPLEDMODEL_HH
#define DUNE_COUPLEDMODEL_HH

namespace Dune
{

template<class FirstModel, class SecondModel>
class CoupledModel {
public:
	typedef typename FirstModel::GridType FirstGrid;
	typedef typename FirstModel::MatrixType FirstMatrixType;
    typedef typename FirstMatrixType::RowIterator FirstRowIterator;
    typedef typename FirstMatrixType::ColIterator FirstColIterator;
	typedef typename FirstGrid::LeafGridView FirstGV;
    typedef typename FirstGV::IndexSet FirstIS;
	typedef typename FirstGV::template Codim<0>::Iterator FirstIterator;

	typedef typename SecondModel::GridType SecondGrid;
	typedef typename SecondModel::MatrixType SecondMatrixType;
    typedef typename SecondMatrixType::RowIterator SecondRowIterator;
    typedef typename SecondMatrixType::ColIterator SecondColIterator;
	typedef typename SecondGrid::LeafGridView SecondGV;
    typedef typename SecondGV::IndexSet SecondIS;
	typedef typename SecondGV::template Codim<0>::Iterator SecondIterator;

	typedef FieldMatrix<double,1,1> GlobalMBlockType;
	typedef BCRSMatrix<GlobalMBlockType> GlobalMatrixType;
	typedef FieldVector<double,1> GlobalVBlockType;
	typedef BlockVector<GlobalVBlockType> GlobalVectorType;

	template <class A12Type>
	void assembleA12(A12Type& A12)
	{
		for (typename A12Type::size_type i = 0; i < A12.M(); i++)
			A12.setrowsize(i, 0);
		A12.endrowsizes();
		A12.endindices();
	}

	template <class A21Type>
	void assembleA21(A21Type& A21)
	{
		for (typename A21Type::size_type i = 0; i < A21.M(); i++)
			A21.setrowsize(i, 0);
		A21.endrowsizes();
		A21.endindices();
	}

	void assemblePostProcess()
	{}

	void initial()
	{
		firstModel_.initial();
		secondModel_.initial();
	}

	void assemble()
	{
		firstModel_.assemble();
		secondModel_.assemble();

		const typename FirstMatrixType::block_type::size_type mFirst = FirstMatrixType::block_type::rows;
		const typename FirstMatrixType::block_type::size_type nFirst = FirstMatrixType::block_type::cols;
		const typename SecondMatrixType::block_type::size_type mSecond = SecondMatrixType::block_type::rows;
		const typename SecondMatrixType::block_type::size_type nSecond = SecondMatrixType::block_type::cols;
		std::cout << "nFirst = " << nFirst << ", mFirst = " << mFirst << ", nSecond = " << nSecond << ", mSecond = " << mSecond << std::endl;
		typedef FieldMatrix<double,mFirst,nSecond> BlockType12;
		typedef BCRSMatrix<BlockType12> A12Type;
		typedef FieldMatrix<double,mSecond,nFirst> BlockType21;
		typedef BCRSMatrix<BlockType21> A21Type;

		typename FirstMatrixType::size_type firstNBlock = firstModel_.matrix().N();
		typename FirstMatrixType::size_type firstMBlock = firstModel_.matrix().M();
		typename SecondMatrixType::size_type secondNBlock = secondModel_.matrix().N();
		typename SecondMatrixType::size_type secondMBlock = secondModel_.matrix().M();
		A12Type A12(firstMBlock, secondNBlock, A12Type::random);
		A21Type A21(secondMBlock, firstNBlock, A21Type::random);
		assembleA12(A12);
		assembleA21(A21);

		assemblePostProcess();

		if (assembleGlobalSystem_)
		{
			A.setBuildMode(GlobalMatrixType::random);
			typename GlobalMatrixType::size_type nGlobal, mGlobal;
			nGlobal = nFirst*firstNBlock + nSecond*secondNBlock;
			mGlobal = mFirst*firstMBlock + mSecond*secondMBlock;
			std::cout << "nGlobal = " << nGlobal << ", mGlobal = " << mGlobal << std::endl;
			A.setSize(mGlobal, nGlobal);
			u.resize(nGlobal);
			f.resize(mGlobal);

	    	FirstRowIterator endIBlock1 = A11.end();
	    	for (FirstRowIterator iBlock = A11.begin(); iBlock != endIBlock1; ++iBlock)
	    	{
				for (typename FirstMatrixType::block_type::size_type iLocal = 0; iLocal < mFirst; iLocal++)
				{
					typename GlobalMatrixType::size_type rowSize = nFirst*iBlock->size() + nSecond*A12.getrowsize(iBlock.index());
					A.setrowsize(iBlock.index()*mFirst + iLocal, rowSize);
					std::cout << "1: row " << iBlock.index()*mFirst + iLocal << ", size " << rowSize << std::endl;
				}
	    	}
	    	SecondRowIterator endIBlock2 = A22.end();
	    	for (SecondRowIterator iBlock = A22.begin(); iBlock != endIBlock2; ++iBlock)
	    	{
				for (typename SecondMatrixType::block_type::size_type iLocal = 0; iLocal < mSecond; iLocal++)
				{
					typename GlobalMatrixType::size_type rowSize = nSecond*iBlock->size() + nFirst*A21.getrowsize(iBlock.index());
					A.setrowsize(mFirst*firstMBlock + iBlock.index()*mSecond + iLocal, rowSize);
					std::cout << "2: row " << mFirst*firstMBlock + iBlock.index()*mSecond + iLocal << ", size " << rowSize << std::endl;
				}
	    	}
	    	A.endrowsizes();


	    	endIBlock1 = A11.end();
	    	for (FirstRowIterator iBlock = A11.begin(); iBlock != endIBlock1; ++iBlock)
	    	{
	    		FirstColIterator endJBlock = iBlock->end();
	    		for (FirstColIterator jBlock = iBlock->begin(); jBlock != endJBlock; ++jBlock)
	    		{
					for (typename FirstMatrixType::block_type::size_type iLocal = 0; iLocal < mFirst; iLocal++)
					{
						typename GlobalMatrixType::size_type globalRow = iBlock.index()*mFirst + iLocal;
						for (typename FirstMatrixType::block_type::size_type jLocal = 0; jLocal < nFirst; jLocal++)
						{
							typename GlobalMatrixType::size_type globalCol = jBlock.index()*nFirst + jLocal;
							A.addindex(globalRow, globalCol);
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
					for (typename A12Type::block_type::size_type iLocal = 0; iLocal < mFirst; iLocal++)
					{
						typename GlobalMatrixType::size_type globalRow = nFirst*firstNBlock + iBlock.index()*mFirst + iLocal;
						for (typename A12Type::block_type::size_type jLocal = 0; jLocal < nSecond; jLocal++)
						{
							typename GlobalMatrixType::size_type globalCol = jBlock.index()*nSecond + jLocal;
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
					for (typename SecondMatrixType::block_type::size_type iLocal = 0; iLocal < mSecond; iLocal++)
					{
						typename GlobalMatrixType::size_type globalRow = mFirst*firstMBlock + iBlock.index()*mSecond + iLocal;
						for (typename SecondMatrixType::block_type::size_type jLocal = 0; jLocal < nSecond; jLocal++)
						{
							typename GlobalMatrixType::size_type globalCol = nFirst*firstNBlock + jBlock.index()*nSecond + jLocal;
							A.addindex(globalRow, globalCol);
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
					for (typename A21Type::block_type::size_type iLocal = 0; iLocal < mSecond; iLocal++)
					{
						typename GlobalMatrixType::size_type globalRow = mFirst*firstMBlock + iBlock.index()*mSecond + iLocal;
						for (typename A21Type::block_type::size_type jLocal = 0; jLocal < nFirst; jLocal++)
						{
							typename GlobalMatrixType::size_type globalCol = jBlock.index()*nFirst + jLocal;
							A.addindex(globalRow, globalCol);
						}
					}
	    		}
			}
	    	A.endindices();

	    	endIBlock1 = A11.end();
	    	for (FirstRowIterator iBlock = A11.begin(); iBlock != endIBlock1; ++iBlock)
	    	{
	    		FirstColIterator endJBlock = iBlock->end();
	    		for (FirstColIterator jBlock = iBlock->begin(); jBlock != endJBlock; ++jBlock)
	    		{
					for (typename FirstMatrixType::block_type::size_type iLocal = 0; iLocal < mFirst; iLocal++)
					{
						typename GlobalMatrixType::size_type globalRow = iBlock.index()*mFirst + iLocal;
						f[globalRow] = (firstModel_.rhs())[iBlock.index()][iLocal];
						for (typename FirstMatrixType::block_type::size_type jLocal = 0; jLocal < nFirst; jLocal++)
						{
							typename GlobalMatrixType::size_type globalCol = jBlock.index()*nFirst + jLocal;
							A[globalRow][globalCol] = (*jBlock)[iLocal][jLocal];
							u[globalCol] = (firstModel_.sol())[jBlock.index()][jLocal];
						}
					}
	    		}
			}
	    	endIBlock12 = A12.end();
	    	for (typename A12Type::RowIterator iBlock = A12.begin(); iBlock != endIBlock12; ++iBlock)
	    	{
	    		typename A12Type::ColIterator endJBlock = iBlock->end();
	    		for (typename A12Type::ColIterator jBlock = iBlock->begin(); jBlock != endJBlock; ++jBlock)
	    		{
					for (typename A12Type::block_type::size_type iLocal = 0; iLocal < mFirst; iLocal++)
					{
						typename GlobalMatrixType::size_type globalRow = nFirst*firstNBlock + iBlock.index()*mFirst + iLocal;
						for (typename A12Type::block_type::size_type jLocal = 0; jLocal < nSecond; jLocal++)
						{
							typename GlobalMatrixType::size_type globalCol = jBlock.index()*nSecond + jLocal;
							A[globalRow][globalCol] = (*jBlock)[iLocal][jLocal];
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
					for (typename SecondMatrixType::block_type::size_type iLocal = 0; iLocal < mSecond; iLocal++)
					{
						typename GlobalMatrixType::size_type globalRow = mFirst*firstMBlock + iBlock.index()*mSecond + iLocal;
						f[globalRow] = (secondModel_.rhs())[iBlock.index()][iLocal];
						for (typename SecondMatrixType::block_type::size_type jLocal = 0; jLocal < nSecond; jLocal++)
						{
							typename GlobalMatrixType::size_type globalCol = nFirst*firstNBlock + jBlock.index()*nSecond + jLocal;
							A[globalRow][globalCol] = (*jBlock)[iLocal][jLocal];
							u[globalCol] = (secondModel_.sol())[jBlock.index()][jLocal];
						}
					}
	    		}
			}
	    	endIBlock21 = A21.end();
	    	for (typename A21Type::RowIterator iBlock = A21.begin(); iBlock != endIBlock21; ++iBlock)
	    	{
	    		typename A21Type::ColIterator endJBlock = iBlock->end();
	    		for (typename A21Type::ColIterator jBlock = iBlock->begin(); jBlock != endJBlock; ++jBlock)
	    		{
					for (typename A21Type::block_type::size_type iLocal = 0; iLocal < mSecond; iLocal++)
					{
						typename GlobalMatrixType::size_type globalRow = mFirst*firstMBlock + iBlock.index()*mSecond + iLocal;
						for (typename A21Type::block_type::size_type jLocal = 0; jLocal < nFirst; jLocal++)
						{
							typename GlobalMatrixType::size_type globalCol = jBlock.index()*nFirst + jLocal;
							A[globalRow][globalCol] = (*jBlock)[iLocal][jLocal];
						}
					}
	    		}
			}


		}
	}

	FirstModel& firstModel()
	{
		return firstModel_;
	}

	SecondModel& secondModel()
	{
		return secondModel_;
	}

	const GlobalMatrixType& matrix()
	{
		return A;
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
};

}

#endif
