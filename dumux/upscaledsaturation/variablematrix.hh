#ifndef DUNE_VARIABLEMATRIX_HH
#define DUNE_VARIABLEMATRIX_HH

#include<vector>

namespace Dune
{

template<class Scalar, class EntryType> class VariableMatrix
{
public:

    VariableMatrix(Scalar initVel = 0.0, int initRowNum = 1, int initColNum = 1) :
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
            int testRealSize = variableMatrix_.size();
            if (testRealSize != rowSize_)
            {
                rowSize_ = testRealSize;
            }
            if (rowNum >= rowSize_)
            {
                DUNE_THROW(RangeError, "rowNumber exceeds matrix size!");
            }
        }
        return variableMatrix_[rowNum];
    }

    std::vector<std::vector<EntryType> > operator+(
            std::vector<std::vector<EntryType> >& summand) const
    {
        std::vector<std::vector<EntryType> > sum;
        sum.resize(rowSize_);
        for (int i = 0; i < rowSize_; i++)
        {
            int columnSize = variableMatrix_[i].size();
            sum[i].resize(columnSize);
        }
        for (int i = 0; i < rowSize_; i++)
        {
            int columnSize = variableMatrix_[i].size();
            for (int j = 0; j < columnSize; j++)
            {
                sum[i][j] = variableMatrix_[i][j] + summand[i][j];
            }
        }
        return sum;
    }

    void addEntry(EntryType newEntry, int rowNum, int columnNum)
    {
        if (rowSize_ < rowNum + 1)
        {
            std::vector<EntryType> addEntry(0);
            addEntry = 0;
            variableMatrix_.push_back(addEntry);
        }
        int columnSize = variableMatrix_[rowNum].size();
        if (columnSize < columnNum + 1)
        {
            variableMatrix_[rowNum].push_back(newEntry);
        }
        else
        {
            EntryType oldEntry = variableMatrix_[rowNum][columnNum];
            variableMatrix_[rowNum][columnNum] = oldEntry + newEntry;
        }
        rowSize_ = std::max(rowSize_, rowNum + 1);
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

    void resizeColumnX(int newSize, int rowNum)
    {
        variableMatrix_[rowNum].resize(newSize);
        return;
    }

    std::vector<std::vector<EntryType> >& writeMatrix()
    {
        return variableMatrix_;
    }

    std::vector<std::vector<EntryType> >& getMatrix() const
    {
        return variableMatrix_;
    }

private:
    int rowSize_;
    std::vector<std::vector<EntryType> > variableMatrix_;
};
}
#endif
