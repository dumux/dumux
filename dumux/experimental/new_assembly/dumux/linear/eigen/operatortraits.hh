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
 * \brief Operator traits specialization for eigen
 */
#ifndef DUMUX_LINEAR_EIGEN_OPERATOR_TRAITS_HH
#define DUMUX_LINEAR_EIGEN_OPERATOR_TRAITS_HH

#include <type_traits>

#include <dumux/experimental/new_assembly/dumux/linear/operator.hh>
#include <dumux/experimental/new_assembly/dumux/linear/eigen/detail.hh>

namespace Dumux {
namespace Eigen_ {

template<LinearSystem::Detail::EigenMatrix M>
class EigenMatrixOperator
{
public:
    template<typename _M> requires(
        std::is_lvalue_reference_v<_M> and
        std::constructible_from<const M&, _M>)
    EigenMatrixOperator(_M&& m)
    : m_(m)
    {}

    operator const M&() const
    { return m_; }

private:
    const M& m_;
};

} // namespace Eigen_

namespace Linear::Traits {

template<LinearSystem::Detail::EigenMatrix M,
         LinearSystem::Detail::EigenVector V>
struct MatrixOperator<M, V>
: public std::type_identity<Eigen_::EigenMatrixOperator<M>>
{};

} // namespace Linear::Traits
} // namespace Dumux

#endif
