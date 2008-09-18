// $Id$ 

#ifndef DUNE_COUPLEDMODEL_HH
#define DUNE_COUPLEDMODEL_HH

namespace Dune
{

template<class FirstModel, class SecondModel>
class CoupledModel {
public:
	typedef typename FirstModel::ProblemType FirstProblem;
	typedef typename FirstModel::GridType FirstGrid;
	typedef typename FirstModel::MatrixType FirstMatrixType;
    typedef typename FirstMatrixType::RowIterator FirstRowIterator;
    typedef typename FirstMatrixType::ColIterator FirstColIterator;
	typedef typename FirstGrid::LeafGridView FirstGV;
    typedef typename FirstGV::IndexSet FirstIS;
	typedef typename FirstGV::template Codim<0>::Iterator FirstIterator;

	typedef typename SecondModel::ProblemType SecondProblem;
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
		for (typename A12Type::size_type i = 0; i < A12.N(); i++)
			A12.setrowsize(i, 0);
		A12.endrowsizes();
		A12.endindices();
	}
	
	template <class A21Type>
	void assembleA21(A21Type& A21) 
	{
		for (typename A21Type::size_type i = 0; i < A21.N(); i++)
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
		
		const typename FirstMatrixType::block_type::size_type nFirst = FirstMatrixType::block_type::rows; 
		const typename FirstMatrixType::block_type::size_type mFirst = FirstMatrixType::block_type::cols; 
		const typename SecondMatrixType::block_type::size_type nSecond = SecondMatrixType::block_type::rows; 
		const typename SecondMatrixType::block_type::size_type mSecond = SecondMatrixType::block_type::cols; 
		std::cout << "nFirst = " << nFirst << ", mFirst = " << mFirst << ", nSecond = " << nSecond << ", mSecond = " << mSecond << std::endl;
		typedef FieldMatrix<double,nFirst,mSecond> BlockType12;
		typedef BCRSMatrix<BlockType12> A12Type;
		typedef FieldMatrix<double,nSecond,mFirst> BlockType21;
		typedef BCRSMatrix<BlockType21> A21Type;
		
		typename FirstMatrixType::size_type firstNBlock = (*(firstModel_.A)).N();
		typename FirstMatrixType::size_type firstMBlock = (*(firstModel_.A)).M();
		typename SecondMatrixType::size_type secondNBlock = (*(secondModel_.A)).N();
		typename SecondMatrixType::size_type secondMBlock = (*(secondModel_.A)).M();
		A12Type A12(firstNBlock, secondMBlock, A12Type::random);
		A21Type A21(secondNBlock, firstMBlock, A21Type::random);
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
			A.setSize(nGlobal, mGlobal);
			u.resize(mGlobal);
			f.resize(nGlobal);
			
	    	FirstRowIterator endIBlock1 = A11.end();
	    	for (FirstRowIterator iBlock = A11.begin(); iBlock != endIBlock1; ++iBlock) 
	    	{
				for (typename FirstMatrixType::block_type::size_type iLocal = 0; iLocal < nFirst; iLocal++)
				{
					typename GlobalMatrixType::size_type rowSize = mFirst*iBlock->size() + mSecond*A12.getrowsize(iBlock.index());
					A.setrowsize(iBlock.index()*nFirst + iLocal, rowSize);
				}
	    	}
	    	SecondRowIterator endIBlock2 = A22.end();
	    	for (SecondRowIterator iBlock = A22.begin(); iBlock != endIBlock2; ++iBlock) 
	    	{
				for (typename SecondMatrixType::block_type::size_type iLocal = 0; iLocal < nSecond; iLocal++)
				{
					typename GlobalMatrixType::size_type rowSize = mSecond*iBlock->size() + mFirst*A21.getrowsize(iBlock.index());
					A.setrowsize(nFirst*firstNBlock + iBlock.index()*nSecond + iLocal, rowSize);
				}
	    	}
	    	A.endrowsizes();
	    	
	    	endIBlock1 = A11.end();
	    	for (FirstRowIterator iBlock = A11.begin(); iBlock != endIBlock1; ++iBlock) 
	    	{
	    		FirstColIterator endJBlock = iBlock->end();
	    		for (FirstColIterator jBlock = iBlock->begin(); jBlock != endJBlock; ++jBlock) 
	    		{
					for (typename FirstMatrixType::block_type::size_type iLocal = 0; iLocal < nFirst; iLocal++)
					{
						typename GlobalMatrixType::size_type globalRow = iBlock.index()*nFirst + iLocal;
						f[globalRow] = (*(firstModel_.f))[iBlock.index()][iLocal];
						for (typename FirstMatrixType::block_type::size_type jLocal = 0; jLocal < mFirst; jLocal++)
						{
							typename GlobalMatrixType::size_type globalCol = jBlock.index()*mFirst + jLocal;
							A.addindex(globalRow, globalCol);
						}	    
					}
	    		}
			}
	    	endIBlock1 = A12.end();
	    	for (FirstRowIterator iBlock = A12.begin(); iBlock != endIBlock1; ++iBlock) 
	    	{
	    		SecondColIterator endJBlock = iBlock->end();
	    		for (SecondColIterator jBlock = iBlock->begin(); jBlock != endJBlock; ++jBlock) 
	    		{
					for (typename A12Type::block_type::size_type iLocal = 0; iLocal < nFirst; iLocal++)
					{
						typename GlobalMatrixType::size_type globalRow = mFirst*firstMBlock + iBlock.index()*nFirst + iLocal;
						for (typename A12Type::block_type::size_type jLocal = 0; jLocal < mSecond; jLocal++)
						{
							typename GlobalMatrixType::size_type globalCol = jBlock.index()*mSecond + jLocal;
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
					for (typename SecondMatrixType::block_type::size_type iLocal = 0; iLocal < nSecond; iLocal++)
					{
						typename GlobalMatrixType::size_type globalRow = nFirst*firstNBlock + iBlock.index()*nSecond + iLocal;
						f[globalRow] = (*(secondModel_.f))[iBlock.index()][iLocal];
						for (typename SecondMatrixType::block_type::size_type jLocal = 0; jLocal < mSecond; jLocal++)
						{
							typename GlobalMatrixType::size_type globalCol = mFirst*firstMBlock + jBlock.index()*mSecond + jLocal;
							A.addindex(globalRow, globalCol);
						}	    
					}
	    		}
			}
	    	endIBlock2 = A21.end();
	    	for (SecondRowIterator iBlock = A21.begin(); iBlock != endIBlock2; ++iBlock) 
	    	{
	    		FirstColIterator endJBlock = iBlock->end();
	    		for (FirstColIterator jBlock = iBlock->begin(); jBlock != endJBlock; ++jBlock) 
	    		{
					for (typename A21Type::block_type::size_type iLocal = 0; iLocal < nSecond; iLocal++)
					{
						typename GlobalMatrixType::size_type globalRow = nFirst*firstNBlock + iBlock.index()*nSecond + iLocal;
						for (typename A21Type::block_type::size_type jLocal = 0; jLocal < mFirst; jLocal++)
						{
							typename GlobalMatrixType::size_type globalCol = jBlock.index()*mFirst + jLocal;
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
					for (typename FirstMatrixType::block_type::size_type iLocal = 0; iLocal < nFirst; iLocal++)
					{
						typename GlobalMatrixType::size_type globalRow = iBlock.index()*nFirst + iLocal;
						f[globalRow] = (*(firstModel_.f))[iBlock.index()][iLocal];
						for (typename FirstMatrixType::block_type::size_type jLocal = 0; jLocal < mFirst; jLocal++)
						{
							typename GlobalMatrixType::size_type globalCol = jBlock.index()*mFirst + jLocal;
							A[globalRow][globalCol] = (*jBlock)[iLocal][jLocal];
							u[globalCol] = (*(firstModel_.u))[jBlock.index()][jLocal];
						}	    
					}
	    		}
			}
	    	endIBlock1 = A12.end();
	    	for (FirstRowIterator iBlock = A12.begin(); iBlock != endIBlock1; ++iBlock) 
	    	{
	    		SecondColIterator endJBlock = iBlock->end();
	    		for (SecondColIterator jBlock = iBlock->begin(); jBlock != endJBlock; ++jBlock) 
	    		{
					for (typename A12Type::block_type::size_type iLocal = 0; iLocal < nFirst; iLocal++)
					{
						typename GlobalMatrixType::size_type globalRow = mFirst*firstMBlock + iBlock.index()*nFirst + iLocal;
						for (typename A12Type::block_type::size_type jLocal = 0; jLocal < mSecond; jLocal++)
						{
							typename GlobalMatrixType::size_type globalCol = jBlock.index()*mSecond + jLocal;
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
					for (typename SecondMatrixType::block_type::size_type iLocal = 0; iLocal < nSecond; iLocal++)
					{
						typename GlobalMatrixType::size_type globalRow = nFirst*firstNBlock + iBlock.index()*nSecond + iLocal;
						f[globalRow] = (*(secondModel_.f))[iBlock.index()][iLocal];
						for (typename SecondMatrixType::block_type::size_type jLocal = 0; jLocal < mSecond; jLocal++)
						{
							typename GlobalMatrixType::size_type globalCol = mFirst*firstMBlock + jBlock.index()*mSecond + jLocal;
							A[globalRow][globalCol] = (*jBlock)[iLocal][jLocal];
							u[globalCol] = (*(secondModel_.u))[jBlock.index()][jLocal];
						}	    
					}
	    		}
			}
	    	endIBlock2 = A21.end();
	    	for (SecondRowIterator iBlock = A21.begin(); iBlock != endIBlock2; ++iBlock) 
	    	{
	    		FirstColIterator endJBlock = iBlock->end();
	    		for (FirstColIterator jBlock = iBlock->begin(); jBlock != endJBlock; ++jBlock) 
	    		{
					for (typename A21Type::block_type::size_type iLocal = 0; iLocal < nSecond; iLocal++)
					{
						typename GlobalMatrixType::size_type globalRow = nFirst*firstNBlock + iBlock.index()*nSecond + iLocal;
						for (typename A21Type::block_type::size_type jLocal = 0; jLocal < mFirst; jLocal++)
						{
							typename GlobalMatrixType::size_type globalCol = jBlock.index()*mFirst + jLocal;
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
	A11(*(firstModel_.A)), A22(*(secondModel_.A)), assembleGlobalSystem_(assembleGlobalSystem)
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
