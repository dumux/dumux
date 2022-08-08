// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*****************************************************************************
 *   See the file COPYING for full copying permissions.                      *
 *                                                                           *
 *   This program is free software: you can redistribute it and/or modify    *
 *   it under the terms of the GNU General Public License as published by    *
 *   the Free Software Foundation, either version 3 of the License, or       *
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
/*!
 * \file
 * \ingroup Linear
 * \brief A parallel version of a linear operator
 */
#ifndef DUMUX_LINEAR_PARALLEL_LINEAR_OPERATOR_HH
#define DUMUX_LINEAR_PARALLEL_LINEAR_OPERATOR_HH

#include <dune/common/hybridutilities.hh>
#include <dune/istl/operators.hh>
#include <dune/istl/bcrsmatrix.hh>

#include <dumux/parallel/parallel_for.hh>

namespace Dumux {

/*!
 * \brief Adapter to turn a matrix into a linear operator.
 * Adapts a matrix to the assembled linear operator interface
 */
template<class M, class X, class Y>
class ThreadParallelMatrixAdapter : public Dune::MatrixAdapter<M,X,Y>
{
    using ParentType = Dune::MatrixAdapter<M,X,Y>;
public:
    //! export types
    typedef M matrix_type;
    typedef X domain_type;
    typedef Y range_type;
    typedef typename X::field_type field_type;

    //! constructor: just store a reference to a matrix
    explicit ThreadParallelMatrixAdapter (const M& A)
    : ParentType(A) {}

    //! constructor: store an std::shared_ptr to a matrix
    explicit ThreadParallelMatrixAdapter (std::shared_ptr<const M> A)
    : ParentType(A) {}

    //! apply operator to x:  \f$ y = A(x) \f$
    void apply (const X& x, Y& y) const override
    {
        y = 0.0;

        auto& A = this->getmat();
        using namespace Dune::Hybrid;
        forEach(integralRange(Dune::Hybrid::size(y)), [&](auto&& i)
        {
            forEach(integralRange(Dune::Hybrid::size(x)), [&](auto&& j)
            {
                const auto& mat = A[i][j];
                const auto& xx = x[j];
                auto& yy = y[i];

                Dumux::parallelFor(mat.N(), [&](const std::size_t ii)
                {
                    const auto& row = mat[ii];
                    const auto endj = row.end();
                    for (auto jj=row.begin(); jj!=endj; ++jj)
                    {
                        const auto& xj = Dune::Impl::asVector(xx[jj.index()]);
                        auto&& yi = Dune::Impl::asVector(yy[ii]);
                        Dune::Impl::asMatrix(*jj).umv(xj, yi);
                    }
                });
            });
        });
    }

    //! apply operator to x, scale and add:  \f$ y = y + \alpha A(x) \f$
    void applyscaleadd (field_type alpha, const X& x, Y& y) const override
    {
        auto& A = this->getmat();
        using namespace Dune::Hybrid;
        forEach(integralRange(Dune::Hybrid::size(y)), [&](auto&& i)
        {
            forEach(integralRange(Dune::Hybrid::size(x)), [&](auto&& j)
            {
                const auto& mat = A[i][j];
                const auto& xx = x[j];
                auto& yy = y[i];

                Dumux::parallelFor(mat.N(), [&](const std::size_t ii)
                {
                    const auto& row = mat[ii];
                    const auto endj = row.end();
                    for (auto jj=row.begin(); jj!=endj; ++jj)
                    {
                        const auto& xj = Dune::Impl::asVector(xx[jj.index()]);
                        auto&& yi = Dune::Impl::asVector(yy[ii]);
                        Dune::Impl::asMatrix(*jj).usmv(alpha, xj, yi);
                    }
                });
            });
        });
    }

    //! Category of the scalar product (see SolverCategory::Category)
    Dune::SolverCategory::Category category() const override
    {
      return Dune::SolverCategory::sequential;
    }
};

/**
* @brief A nonoverlapping operator with communication object.
*/
template<class M, class X, class Y, class C>
class ParallelStokesOperator : public Dune::MatrixAdapter<M,X,Y>
{
    static constexpr std::size_t N = 2;

    using ParentType = Dune::MatrixAdapter<M,X,Y>;

    using A = std::decay_t<decltype(std::declval<M>()[Dune::Indices::_0][Dune::Indices::_0])>;
    using U = std::decay_t<decltype(std::declval<X>()[Dune::Indices::_0])>;

    using P = std::decay_t<decltype(std::declval<M>()[Dune::Indices::_1][Dune::Indices::_1])>;
    using V = std::decay_t<decltype(std::declval<X>()[Dune::Indices::_1])>;

public:
    //! \brief The type of the matrix we operate on.
    typedef M matrix_type;
    //! \brief The type of the domain.
    typedef X domain_type;
    //! \brief The type of the range.
    typedef Y range_type;
    //! \brief The field type of the range
    typedef typename X::field_type field_type;
    //! \brief The type of the communication object
    typedef C communication_type;

    /**
     * @brief constructor: just store a reference to a matrix.
     *
     * @param A The assembled matrix.
     * @param com The communication object for syncing owner and copy
     * data points. (E.~g. OwnerOverlapCommunication )
     */
    ParallelStokesOperator (const matrix_type& A_, const std::array<std::shared_ptr<const communication_type>, 2>& comms)
    : ParentType(Dune::stackobject_to_shared_ptr(A_))
    , buildComm_({ true, true })
    , comms_(comms)
    {}

   /*
    * \brief Apply operator to x:  \f$ y = A(x) \f$
    * \param x vector to apply the operator to
    * \param y vector to store the result
    *
    * x is assumed to be in a consistent representation
    * that is all dofs _stored_ on this process have the correct value of the global vector
    * Blatt and Bastian (2009) Definition (2.3) https://doi.org/10.1504/IJCSE.2008.021112
    * At the exit of this function, y has to be in a unique additive representation (Def 2.5)
    * This is achieved by masking out all indices which are not owned by this process.
    * Applying the (local) operator doesn't require communication if the operator
    * pattern and entries are extended using the ParallelISTLHelper for non-overlapping decompositions.
    */
    void apply (const X& x, Y& y) const override
    {
        using namespace Dune::Indices;

        y = 0;

        applyNonoverlapping_(1.0, x, y[_0], _0, _0);
        this->getmat()[_0][_1].umv(x[_1], y[_0]);
        // make A*x consistent
        comms_[0]->addOwnerCopyToOwnerCopy(y[_0], y[_0]);

        this->getmat()[_1][_1].umv(x[_1], y[_1]);
        this->getmat()[_1][_0].umv(x[_0], y[_1]);
        comms_[1]->project(y[_1]);
    }

    //! apply operator to x, scale and add:  \f$ y = y + \alpha A(x) \f$
    void applyscaleadd (field_type alpha, const X& x, Y& y) const override
    {
        using namespace Dune::Indices;
        using namespace Dune::Hybrid;

        // y is already consistent so compute alpha*A*x separately
        auto y1(y[_0]);
        y[_0] = 0.0;
        applyNonoverlapping_(alpha, x, y[_0], _0, _0);
        this->getmat()[_0][_1].usmv(alpha, x[_1], y[_0]);
        // make alpha*A*x consistent
        comms_[0]->addOwnerCopyToOwnerCopy(y[_0], y[_0]);
        y[_0] += y1;

        this->getmat()[_1][_1].usmv(alpha, x[_1], y[_1]);
        this->getmat()[_1][_0].usmv(alpha, x[_0], y[_1]);
        comms_[1]->project(y[_1]);
    }

    //! Category of the linear operator (see SolverCategory::Category)
    Dune::SolverCategory::Category category() const override
    {
        return Dune::SolverCategory::overlapping;
    }

private:

    template<std::size_t rId, std::size_t cId>
    void applyNonoverlapping_(field_type alpha, const X& x, U& y, Dune::index_constant<rId> rowID, Dune::index_constant<cId> colID) const
    {
        static_assert(rId == 0 && cId == 0, "Only implemented for first row and column");

        // this is from Dune::NonoverlappingSchwarzOperator but rewritten/reformated

        // at the beginning make a multimap "borderContribution_".
        // process has i and j as border dofs but is not the owner
        // => only contribute to Ax if i,j is in borderContribution_
        if (buildComm_[colID])
        {
            // create masks that contain the attribute of the dof on the local process
            auto& mask = std::get<colID>(mask_);
            const auto& parallelIndexSet = std::get<colID>(comms_)->indexSet();
            mask.assign(x[colID].size(), 1); // initialize owner
            for (auto it = parallelIndexSet.begin(); it != parallelIndexSet.end(); ++it)
            {
                if (it->local().attribute() == Dune::OwnerOverlapCopyAttributeSet::copy)
                    mask[it->local().local()] = 0; // copy (not owner)
                else if (it->local().attribute() == Dune::OwnerOverlapCopyAttributeSet::overlap)
                    mask[it->local().local()] = 2; // overlap (not owner)
            }

            std::get<colID>(borderContribution_).clear();

            // temporary data structures
            // key: local index i; value: process that owns i
            std::map<int, int> owner;
            // for each local index make the multimap remoteIndexMap:
            // key: local index i, data: pair of process that knows i and pointer to remote index entry
            using RemoteIndexListIterator = typename C::RI::RemoteIndexList::const_iterator;
            std::multimap<int, std::pair<int, RemoteIndexListIterator>> remoteIndexMap;

            auto& subA = this->getmat()[colID][colID];
            for (auto rowIt = subA.begin(); rowIt != subA.end(); ++rowIt)
            {
                // skip owner and overlap
                if (std::get<colID>(mask_)[rowIt.index()] == 0)
                {
                    const auto& remoteIndices = std::get<colID>(comms_)->remoteIndices();
                    for (auto rIt = remoteIndices.begin(); rIt != remoteIndices.end(); ++rIt)
                    {
                        auto& remoteIndexList = *(rIt->second.first); // ??
                        for (auto rIdxIt = remoteIndexList.begin(); rIdxIt != remoteIndexList.end(); ++rIdxIt)
                        {
                            if (rIdxIt->attribute() != Dune::OwnerOverlapCopyAttributeSet::overlap)
                            {
                                if (rIdxIt->localIndexPair().local().local() == rowIt.index())
                                {
                                    remoteIndexMap.insert(std::make_pair(
                                        rowIt.index(), std::make_pair(rIt->first, rIdxIt)
                                    ));

                                    if (rIdxIt->attribute() == Dune::OwnerOverlapCopyAttributeSet::owner)
                                        owner.insert(
                                            std::make_pair(rowIt.index(), rIt->first)
                                        );
                                }
                            }
                        }
                    }
                }
            }

            for (auto rowIt = subA.begin(); rowIt != subA.end(); ++rowIt)
            {
                // handle copy dofs of the row (i.e. dofs we don't own)
                if (std::get<rowID>(mask_)[rowIt.index()] == 0)
                {
                    const int ownerRank = owner.find(rowIt.index())->second;
                    const auto [beginRemoteRowIt, endRemoteRowIt] = remoteIndexMap.equal_range(rowIt.index());

                    for (auto colIt = rowIt->begin(); colIt != rowIt->end(); ++colIt)
                    {
                        // handle copy dofs in the column (i.e. dofs we don't own)
                        if (std::get<colID>(mask_)[colIt.index()] == 0)
                        {
                            bool addBorderContribution = true;
                            const auto [beginRemoteColIt, endRemoteColIt] = remoteIndexMap.equal_range(colIt.index());

                            for (auto remoteRowIt = beginRemoteRowIt; remoteRowIt != endRemoteRowIt; ++remoteRowIt)
                            {
                                for (auto remoteColIt = beginRemoteColIt; remoteColIt != endRemoteColIt; ++remoteColIt)
                                {
                                    // if the ranks are the same
                                    if (remoteColIt->second.first == remoteRowIt->second.first)
                                    {
                                        if (remoteColIt->second.second->attribute() == Dune::OwnerOverlapCopyAttributeSet::owner
                                            || remoteColIt->second.first == ownerRank
                                            || remoteColIt->second.first < std::get<colID>(comms_)->communicator().rank())
                                        {
                                            addBorderContribution = false;
                                            break;
                                        }
                                    }
                                }
                            }

                            // don't contribute to A[row]*x if
                            // 1. the owner of j has i as interior/border dof
                            // 2. ownerRank has j as interior/border dof
                            // 3. there is another process with smaller rank that has i and j as interior/border dofs
                            //
                            // if the owner of j does not have i as interior/border dof,
                            // it will not be taken into account
                            if (addBorderContribution)
                                std::get<colID>(borderContribution_).insert(std::pair<int, int>{ rowIt.index(), colIt.index() });
                        }
                    }
                }
            }

            buildComm_[colID] = false;
        }

        // Compute alpha*A[rowID]*x for the non-overlapping case
        auto& subA = this->getmat()[rowID][colID];
        for (auto rowIt = subA.begin(); rowIt != subA.end(); ++rowIt)
        {
            // the case of a copy dof
            // dof is on border because copies only exist on the border
            // and it means that another process is the owner
            if (std::get<rowID>(mask_)[rowIt.index()] == 0)
            {
                for (auto colIt = rowIt->begin(); colIt != rowIt->end(); ++colIt)
                {
                    // the column dof is owned by us
                    if (std::get<colID>(mask_)[colIt.index()] == 1) // j is owner => then sum entries
                        Dune::Impl::asMatrix(*colIt).usmv(alpha, x[colID][colIt.index()], y[rowIt.index()]);

                    // the column dof is a copy: check if we have a border contribution, and if yes, sum entries
                    else if (std::get<colID>(mask_)[colIt.index()] == 0)
                    {
                        const auto [beginColBCIt, endColBCIt] = std::get<colID>(borderContribution_).equal_range(rowIt.index());
                        for (auto bcColIt = beginColBCIt; bcColIt != endColBCIt; ++bcColIt)
                            if (bcColIt->second == int(colIt.index()))
                                Dune::Impl::asMatrix(*colIt).usmv(alpha, x[colID][colIt.index()], y[rowIt.index()]);
                    }
                }
            }

            // the case of an owned dof
            else if (std::get<rowID>(mask_)[rowIt.index()] == 1)
            {
                for (auto colIt = rowIt->begin(); colIt != rowIt->end(); ++colIt)
                {
                    if (std::get<colID>(mask_)[colIt.index()] != 2) // anything but overlap
                        Dune::Impl::asMatrix(*colIt).usmv(alpha, x[colID][colIt.index()], y[rowIt.index()]);
                }
            }
        }
    }

    // for now only one row
    mutable std::array<bool, N> buildComm_;
    mutable std::array<std::vector<int>, N> mask_;
    mutable std::array<std::multimap<int, int>, N> borderContribution_;

    std::array<std::shared_ptr<const communication_type>, N> comms_;
};

/**
 * \brief Scalar product for overlapping Schwarz methods.
 *
 * Consistent vectors in interior and border are assumed.
 * \tparam  X The type of the sequential vector to use for the left hand side,
 * e.g. BlockVector or another type fulfilling the ISTL
 * vector interface.
 * \tparam C The type of the communication object.
 * This must either be OwnerOverlapCopyCommunication or a type
 * implementing the same interface.
 */
template<class X, class C, std::size_t N>
class ParallelScalarProduct : public Dune::ScalarProduct<X>
{
public:
    //! \brief The type of the vector to compute the scalar product on.
    //!
    //! E.g. BlockVector or another type fulfilling the ISTL
    //! vector interface.
    typedef X domain_type;
    //!  \brief The field type used by the vector type domain_type.
    typedef typename X::field_type field_type;
    typedef typename Dune::FieldTraits<field_type>::real_type real_type;
    //! \brief The type of the communication object.
    //!
    //! This must either be OwnerOverlapCopyCommunication or a type
    //! implementing the same interface.
    typedef C communication_type;

    /*!
     * \param com The communication object for syncing overlap and copy
     * data points.
     * \param cat parallel solver category (nonoverlapping or overlapping)
     */
    ParallelScalarProduct (const std::array<std::shared_ptr<const communication_type>, N>& comms)
    : comms_(comms)
    {}

    /*! \brief Dot product of two vectors.
     * It is assumed that the vectors are consistent on the interior+border
     * partition. According to Blatt and Bastian (2009)
     * https://doi.org/10.1504/IJCSE.2008.021112 they only have to be in a
     ' valid representation (i.e. all dofs owned by the process have the same value as the global vector)
     */
    field_type dot (const X& x, const X& y) const override
    {
        field_type result = 0.0;

        using namespace Dune::Hybrid;
        forEach(integralRange(Dune::Hybrid::size(y)), [&](auto&& i)
        {
            field_type subResult = 0.0;
            comms_[i()]->dot(x[i], y[i], subResult);
            result += subResult;
        });

        return result;
    }

    /*! \brief Norm of a right-hand side vector.
     */
    real_type norm (const X& x) const override
    {
        return std::sqrt(dot(x, x));
    }

    //! Category of the scalar product (see SolverCategory::Category)
    Dune::SolverCategory::Category category() const override
    {
        return Dune::SolverCategory::overlapping;
    }

private:
    std::array<std::shared_ptr<const communication_type>, N> comms_;
};

} // end namespace Dumux

#endif
