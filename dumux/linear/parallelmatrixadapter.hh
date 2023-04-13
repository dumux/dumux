// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup Linear
 * \brief A parallel version of a linear operator
 */
#ifndef DUMUX_LINEAR_PARALLEL_MATRIX_ADAPTER_HH
#define DUMUX_LINEAR_PARALLEL_MATRIX_ADAPTER_HH

#include <dune/common/hybridutilities.hh>
#include <dune/istl/operators.hh>
#include <dune/istl/bcrsmatrix.hh>
#include <dumux/parallel/parallel_for.hh>

namespace Dumux {

/*!
 * \brief Adapter to turn a multi-type matrix into a thread-parallel linear operator.
 * Adapts a matrix to the assembled linear operator interface
 */
template<class M, class X, class Y>
class ParallelMultiTypeMatrixAdapter : public Dune::MatrixAdapter<M,X,Y>
{
    using ParentType = Dune::MatrixAdapter<M,X,Y>;
public:
    //! export types
    typedef M matrix_type;
    typedef X domain_type;
    typedef Y range_type;
    typedef typename X::field_type field_type;

    //! constructor: just store a reference to a matrix
    explicit ParallelMultiTypeMatrixAdapter (const M& A)
    : ParentType(A) {}

    //! constructor: store an std::shared_ptr to a matrix
    explicit ParallelMultiTypeMatrixAdapter (std::shared_ptr<const M> A)
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
};

} // end namespace Dumux

#endif
