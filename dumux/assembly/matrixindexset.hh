// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUMUX_SIMPLEMATRIXINDEXSET_HH
#define DUMUX_SIMPLEMATRIXINDEXSET_HH

#include <vector>
#include <set>

namespace Dune {

//TODO: figure out how to handle this copy from Dune

  /** \brief Stores the nonzero entries in a sparse matrix */
  class SimpleMatrixIndexSet
  {

  public:
    typedef std::size_t size_type;

    /** \brief Default constructor */
    SimpleMatrixIndexSet() : rows_(0), cols_(0)
    {}

    /** \brief Constructor setting the matrix size */
    SimpleMatrixIndexSet(size_type rows, size_type cols) : rows_(rows), cols_(cols) {
      indices_.resize(rows_);
    }

    /** \brief Reset the size of an index set */
    void resize(size_type rows, size_type cols) {
      rows_ = rows;
      cols_ = cols;
      indices_.resize(rows_);
    }

    /** \brief Add an index to the index set */
    void add(size_type i, size_type j) {
      indices_[i].insert(j);
    }

    /** \brief Return the number of entries */
    size_type size() const {
      size_type entries = 0;
      for (size_type i=0; i<rows_; i++)
        entries += indices_[i].size();

      return entries;
    }

    /** \brief Return the number of rows */
    size_type rows() const {return rows_;}


    /** \brief Return the number of entries in a given row */
    size_type rowsize(size_type row) const {return indices_[row].size();}

    /** \brief Import all nonzero entries of a sparse matrix into the index set
        \tparam MatrixType Needs to be BCRSMatrix<...>
        \param m reference to the MatrixType object
        \param rowOffset don't write to rows<rowOffset
        \param colOffset don't write to cols<colOffset
     */
    template <class MatrixType>
    void import(const MatrixType& m, size_type rowOffset=0, size_type colOffset=0) {

      typedef typename MatrixType::row_type RowType;
      typedef typename RowType::ConstIterator ColumnIterator;

      for (size_type rowIdx=0; rowIdx<m.N(); rowIdx++) {

        const RowType& row = m[rowIdx];

        ColumnIterator cIt    = row.begin();
        ColumnIterator cEndIt = row.end();

        for(; cIt!=cEndIt; ++cIt)
          add(rowIdx+rowOffset, cIt.index()+colOffset);

      }

    }

    void getIndices(std::vector<std::set<size_type> >& indicesParam) const {
        indicesParam = indices_;
    }

    void setIndices(const std::vector<std::set<size_type> >& indicesParam){
        indices_ = indicesParam;
    }

    std::vector<std::set<size_type> >* getPtrToIndices(){
        return &indices_;
    }

    /** \brief Initializes a BCRSMatrix with the indices contained
        in this MatrixIndexSet
        \tparam MatrixType Needs to be BCRSMatrix<...>
        \param matrix reference to the MatrixType object
     */
    template <class MatrixType>
    void exportIdx(MatrixType& matrix) const {

      matrix.setSize(rows_, cols_);
      matrix.setBuildMode(MatrixType::random);

      for (size_type i=0; i<rows_; i++)
        matrix.setrowsize(i, indices_[i].size());

      matrix.endrowsizes();

      for (size_type i=0; i<rows_; i++) {

        typename std::set<size_type>::iterator it = indices_[i].begin();
        for (; it!=indices_[i].end(); ++it)
          matrix.addindex(i, *it);

      }

      matrix.endindices();

    }

  private:

    std::vector<std::set<size_type> > indices_;

    size_type rows_, cols_;

  };


} // end namespace Dune

#endif
