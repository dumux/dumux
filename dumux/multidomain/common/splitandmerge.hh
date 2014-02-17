// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*****************************************************************************
 *   See the file COPYING for full copying permissions.                      *
 *                                                                           *
 *   This program is free software: you can redistribute it and/or modify    *
 *   it under the terms of the GNU General Public License as published by    *
 *   the Free Software Foundation, either version 2 of the License, or       *
 *   (at your option) any later version.                                     *
 *                                                                           *
 *   This program is distributed in the hope that it will be useful,         *
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of          *
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the            *
 *   GNU General Public License for more details.                            *
 *                                                                           *
 *   You should have received a copy of the GNU General Public License       *
 *   along with this program.  If not, see <http://www.gnu.org/licenses/>.   *
 *****************************************************************************/
/**
 * \file
 * \brief Some methods required by several classes of the coupling model.
 */
#ifndef DUMUX_SPLIT_AND_MERGE_HH
#define DUMUX_SPLIT_AND_MERGE_HH

#include "multidomainproperties.hh"
#include <dumux/common/valgrind.hh>

/*
* \brief docme
*/

namespace Dumux
{
/*!
 * \addtogroup ModelCoupling
 */
// \{

/*!
 * \brief Some methods required by several classes of the coupling model.
 */
template<class TypeTag>
class SplitAndMerge
{
    typedef typename GET_PROP_TYPE(TypeTag, SubDomain1TypeTag) SubTypeTag1;
    typedef typename GET_PROP_TYPE(TypeTag, SubDomain2TypeTag) SubTypeTag2;

    typedef typename GET_PROP_TYPE(TypeTag, SolutionVector) SolutionVector;
    typedef typename GET_PROP_TYPE(SubTypeTag1, SolutionVector) SolutionVector1;
    typedef typename GET_PROP_TYPE(SubTypeTag2, SolutionVector) SolutionVector2;

    typedef typename GET_PROP_TYPE(TypeTag, JacobianMatrix) JacobianMatrix;
    typedef typename GET_PROP_TYPE(SubTypeTag1, JacobianMatrix) JacobianMatrix1;
    typedef typename GET_PROP_TYPE(SubTypeTag2, JacobianMatrix) JacobianMatrix2;


    enum {
            numEq1 = GET_PROP_VALUE(SubTypeTag1, NumEq),
            numEq2 = GET_PROP_VALUE(SubTypeTag2, NumEq)
    };
public:
    /*!
     * \brief Merge two solution vectors of the sub models into a
     *        global vector: nonoverlapping case.
     *
     * \param vec1 docme
     * \param vec2 docme
     * \param dest docme
     *
     */
    static void mergeSolVectors(const SolutionVector1 &vec1,
                                 const SolutionVector2 &vec2,
                                 SolutionVector &dest)
    {
//        printvector(std::cout, vec1, "vec1", "row", 200, 1, 3);
//        printvector(std::cout, vec2, "vec2", "row", 200, 1, 3);

        int nDofs1 = vec1.size();
        int nDofs2 = vec2.size();

        // merge the global vector into two individual ones
        for (int i = 0; i < nDofs1; ++i)
            for (int j = 0; j < numEq1; j++)
            {
				dest[i*numEq1 + j][0] = vec1[i][j];

                Valgrind::CheckDefined(dest[i*numEq1 + j][0]);
            }

        for (int i = 0; i < nDofs2; ++i)
            for (int j = 0; j < numEq2; j++)
            {
				dest[nDofs1*numEq1 + i*numEq2 + j][0] = vec2[i][j];

                Valgrind::CheckDefined(dest[nDofs1*numEq1 + i*numEq2 + j][0]);
            }

//        printvector(std::cout, dest, "dest", "row", 200, 1, 3);

    }

    /*!
     * \brief Split a global solution vector into two solution vectors
     *        of the sub models: nonoverlapping case.
     *
     * \param vec docme
     * \param dest1 docme
     * \param dest2 docme
     *
     */
    static void splitSolVector(const SolutionVector &vec,
                               SolutionVector1 &dest1,
                               SolutionVector2 &dest2)
    {
//        printvector(std::cout, vec, "vec", "row", 200, 1, 3);

        int nDofs1 = dest1.size();
        int nDofs2 = dest2.size();

        for (int i = 0; i < nDofs1; ++i)
            for (int j = 0; j < numEq1; j++)
				dest1[i][j] = vec[i*numEq1 + j][0];

        for (int i = 0; i < nDofs2; ++i)
            for (int j = 0; j < numEq2; j++)
				dest2[i][j] = vec[nDofs1*numEq1 + i*numEq2 + j][0];

//        printvector(std::cout, dest1, "dest1", "row", 200, 1, 3);
//        printvector(std::cout, dest2, "dest2", "row", 200, 1, 3);

    }

    /*!
     * \brief Merge two solution vectors of the sub models into a
     *        global vector: more general case.
     *
     * \param vec1 docme
     * \param vec2 docme
     * \param dest docme
     * \param subDOFToCoupledDOF docme
     *
     */
    static void mergeSolVectors(const SolutionVector1 &vec1,
                                 const SolutionVector2 &vec2,
                                 SolutionVector &dest,
                                 const std::vector<int>& subDOFToCoupledDOF)
    {
        int nDofs1 = vec1.size();
        int nDofs2 = vec2.size();
//        printvector(std::cout, vec1, "vec1", "row", 200, 1, 3);
//        printvector(std::cout, vec2, "vec2", "row", 200, 1, 3);

        for (int i = 0; i < nDofs1; ++i)
            for (int j = 0; j < numEq1; j++)
                dest[i*numEq1 + j] = vec1[i][j];

        for (int i = 0; i < nDofs2; ++i)
        {
            for (int j = numEq1; j < numEq2; j++)
                dest[nDofs1*numEq1 + i*(numEq2-numEq1) + (j - numEq1)] = vec2[i][j];
        }
//        printvector(std::cout, dest, "dest", "row", 200, 1, 3);
    }

    /*!
     * \brief Split a global solution vector into two solution vectors
     *        of the sub models: more general case.
     *
     * \param vec docme
     * \param dest1 docme
     * \param dest2 docme
     * \param subDOFToCoupledDOF docme
     *
     */
    static void splitSolVector(const SolutionVector &vec,
                               SolutionVector1 &dest1,
                               SolutionVector2 &dest2,
                               const std::vector<int>& subDOFToCoupledDOF)
    {
        int nDofs1 = dest1.size();
        int nDofs2 = dest2.size();

//        printvector(std::cout, vec, "vec", "row", 200, 1, 3);
        for (int i = 0; i < nDofs1; ++i)
            for (int j = 0; j < numEq1; j++)
                dest1[i][j] = vec[i*numEq1 + j];

        for (int i = 0; i < nDofs2; ++i)
        {
            int blockIdxCoupled = subDOFToCoupledDOF[i];
            for (int j = 0; j < numEq1; j++)
            {
                dest2[i][j] = vec[blockIdxCoupled*numEq1 + j];
            }

            for (int j = numEq1; j < numEq2; j++)
                dest2[i][j] = vec[nDofs1*numEq1 + i*(numEq2-numEq1) + (j - numEq1)];
        }
//        printvector(std::cout, dest1, "dest1", "row", 200, 1, 3);
//        printvector(std::cout, dest2, "dest2", "row", 200, 1, 3);

    }


    /*!
     * \brief Merge individual jacobian matrices of the sub models
     *        into a global jacobian matrix.
     *
     * \param M1 docme
     * \param M2 docme
     * \param M docme
     *
     */
    static void mergeMatrices(const JacobianMatrix1 &M1,
                               const JacobianMatrix2 &M2,
                               JacobianMatrix &M)
    {
    	DUNE_THROW(Dune::NotImplemented, "mergeMatrices in coupled common");
    }

    /*!
     * \brief Copy a sub matrix into into the main diagonal of the global matrix.
     *
     * \param Msub docme
     * \param M docme
     * \param offset docme?
     *
     */
    template <class SubMatrix>
    static void copyMatrix(const SubMatrix &Msub,
                           JacobianMatrix &M,
                           size_t offset)
    {
        // loop over all rows of the submatrix
        typedef typename SubMatrix::ConstRowIterator RowIterator;
        typedef typename SubMatrix::ConstColIterator ColIterator;
        RowIterator endRow = Msub.end();
        for (RowIterator row = Msub.begin(); row != endRow; ++row) {
            // loop over columns of the current row
            ColIterator endCol = row->end();
            for (ColIterator col = row->begin(); col != endCol; ++ col) {
                // copy entry in the global matrix
                M[row.index() + offset][col.index() + offset]
                    = *col;
            }
        }
    }

};
// \}

} // namespace Dumux

#endif
