// $Id:$
#ifndef DUNE_VARIABLEMATRIX_HH
#define DUNE_VARIABLEMATRIX_HH

#include<vector>

namespace Dune
{

template<class Scalar, class EntryType> class VariableMatrix
{
public:
    typedef Dune::VariableMatrix<Scalar, EntryType> ThisType;
    typedef EntryType field_type;

    VariableMatrix(Scalar initVel = 0.0, int initRowNum = 0, int initColNum = 0) :
        rowSize_(initRowNum)
    {
        variableMatrix_.resize(rowSize_);
        for (int i = 0; i < rowSize_; i++)
        {
            variableMatrix_[i].resize(initColNum);
            for (int j = 0; j < initColNum; j++)
            {
                variableMatrix_[i][j] = initVel;
            }
        }
    }

    std::vector<EntryType>& operator[](int rowNum)
    {
        if (rowNum >= rowSize_)
        {
            checkRowSizeValue();

            if (rowNum >= rowSize_)
            {
                DUNE_THROW(RangeError, "rowNumber exceeds matrix size!");
            }
        }
        return variableMatrix_[rowNum];
    }

    const std::vector<EntryType>& operator[](int rowNum) const
    {
        if (rowNum >= rowSize_)
        {
            int testRealSize = variableMatrix_.size();
            if (rowNum >= testRealSize)
            {
                DUNE_THROW(RangeError, "rowNumber exceeds matrix size!");
            }
        }
        return variableMatrix_[rowNum];
    }

    void operator +=(ThisType& summand)
    {
        checkRowSizeValue();

        if (rowSize_ != summand.rowSize())
        {
            DUNE_THROW(RangeError, "not the same rowNumbers!");
        }
        ThisType sum;
        sum.resizeRows(rowSize_);
        for (int i = 0; i < rowSize_; i++)
        {
            int columnSize = std::max(variableMatrix_[i].size(),
                    summand[i].size());
            std::cout<<"columnSize = "<<columnSize<<std::endl;
            std::cout << "columnSizeA = " << variableMatrix_[i].size()
                    << std::endl;
            std::cout << "columnSizeB = " << summand[i].size() << std::endl;
            sum[i].resize(columnSize);
        }
        sum = *this + summand;
        *this = sum;
        std::cout << "test" << std::endl;
        return;
    }

    ThisType operator+(const ThisType& summand) const
    {
        std::cout << "I am here 2" << std::endl;
        ThisType sum;
        sum.resizeRows(rowSize_);
        for (int i = 0; i < rowSize_; i++)
        {
            int columnSize = std::max(variableMatrix_[i].size(),
                    summand[i].size());
            std::cout << "columnSizeA = " << variableMatrix_[i].size()
                    << std::endl;
            std::cout << "columnSizeB = " << summand[i].size() << std::endl;
            std::cout << "columnSize = " << columnSize << std::endl;
            sum[i].resize(columnSize);
        }
        for (int i = 0; i < rowSize_; i++)
        {
            std::cout << "I am here 3" << std::endl;
            int columnSizeA = variableMatrix_[i].size();
            int columnSizeB = summand[i].size();
            int columnSize = std::max(columnSizeA, columnSizeB);
            for (int j = 0; j < columnSize; j++)
            {
                std::cout << "I am here 4" << std::endl;
                if (columnSizeA < columnSize)
                {
                    if (j < columnSizeA)
                    {
                        sum[i][j] = variableMatrix_[i][j] + summand[i][j];
                    }
                    else
                    {
                        sum[i][j] = summand[i][j];
                    }
                }
                if (columnSizeB < columnSize)
                {
                    if (j < columnSizeB)
                    {
                        sum[i][j] = variableMatrix_[i][j] + summand[i][j];
                    }
                    else
                    {
                        sum[i][j] = variableMatrix_[i][j];
                    }
                }
                if (columnSizeA == columnSizeB)
                {
                    sum[i][j] = variableMatrix_[i][j] + summand[i][j];
                }
            }
        }
        std::cout << "test" << std::endl;
        return sum;
    }

    void addEntry(EntryType newEntry, int rowNum, int columnNum)
    {
//        std::cout<<"rowNum = "<<rowNum<<"columnNum = "<<columnNum<<std::endl;
        checkRowSizeValue();
        while (rowSize_ < rowNum + 1)
        {
            EntryType in(0);
            std::vector<EntryType> addEntry(0, in);
            variableMatrix_.push_back(addEntry);
            rowSize_++;
        }
        int columnSize = variableMatrix_[rowNum].size();
//        std::cout<<"columnSizeAdd = "<<columnSize<<std::endl;
        while (columnSize < columnNum + 1)
        {
            EntryType addEntry(0);
            variableMatrix_[rowNum].push_back(addEntry);
            columnSize++;
        }
        EntryType oldEntry = variableMatrix_[rowNum][columnNum];
        variableMatrix_[rowNum][columnNum] = oldEntry + newEntry;

        return;
    }

    int rowSize()
    {
        return rowSize_;
    }

    void resizeRows(int newSize)
    {
        variableMatrix_.resize(newSize);
        rowSize_ = newSize;
        return;
    }
    void checkRowSizeValue()
    {
        int testRealSize = variableMatrix_.size();
        if (testRealSize != rowSize_)
        {
            rowSize_ = testRealSize;
        }
        return;
    }

    void resizeColumnX(int newSize, int rowNum)
    {
        variableMatrix_[rowNum].resize(newSize);
        return;
    }

    std::vector<std::vector<EntryType> >& writeMatrix()
    {
        return variableMatrix_;
    }

    const std::vector<std::vector<EntryType> >& getMatrix()
    {
        return variableMatrix_;
    }

private:
    int rowSize_;
    std::vector<std::vector<EntryType> > variableMatrix_;
};
}
#endif
