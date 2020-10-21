// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/****************************************************************************
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
 * \ingroup Nonlinear
 * \brief Backends for the Newton for different variable types
 */
#ifndef DUMUX_NEWTON_VARIABLES_BACKEND_HH
#define DUMUX_NEWTON_VARIABLES_BACKEND_HH

#include <array>
#include <utility>
#include <type_traits>

#include <dune/common/indices.hh>
#include <dune/common/typetraits.hh>
#include <dune/common/hybridutilities.hh>
#include <dune/common/std/type_traits.hh>
#include <dune/istl/bvector.hh>

// forward declaration
namespace Dune {

template<class... Args>
class MultiTypeBlockVector;

} // end namespace Dune

namespace Dumux {

/*!
 * \ingroup Nonlinear
 * \brief Class providing operations with primary variable vectors
 */
template<class DofVector, class Enable = void>
class NewtonDofBackend;

/*!
 * \ingroup Nonlinear
 * \brief Specialization providing Newton operations for scalar/number types
 */
template<class Scalar>
class NewtonDofBackend<Scalar, std::enable_if_t<Dune::IsNumber<Scalar>::value, Scalar>>
{
public:
    using DofVector = Scalar; //!< the type of the dofs parametrizing the variables object

    static std::size_t size(const DofVector& d)
    { return 1; }

    static DofVector makeDofVector(const DofVector& d)
    { return d; }

    static DofVector makeZeroDofVector(std::size_t size)
    { return 0.0; }
};

/*!
 * \ingroup Nonlinear
 * \brief Specialization providing Newton operations for block vectors
 */
template<class BT>
class NewtonDofBackend<Dune::BlockVector<BT>>
{

public:
    using DofVector = Dune::BlockVector<BT>; //!< the type of the dofs parametrizing the variables object

    static std::size_t size(const DofVector& d)
    { return d.size(); }

    static DofVector makeDofVector(const DofVector& d)
    { return d; }

    static DofVector makeZeroDofVector(std::size_t size)
    { DofVector d; d.resize(size); return d; }
};

/*!
 * \ingroup Nonlinear
 * \brief Specialization providing Newton operations for multitype block vectors
 */
template<class... Blocks>
class NewtonDofBackend<Dune::MultiTypeBlockVector<Blocks...>>
{
    using DV = Dune::MultiTypeBlockVector<Blocks...>;
    static constexpr auto numBlocks = DV::size();

    using VectorSizeInfo = std::array<std::size_t, numBlocks>;

public:
    using DofVector = DV; //!< the type of the dofs parametrizing the variables object

    static VectorSizeInfo size(const DofVector& d)
    {
        VectorSizeInfo result;
        using namespace Dune::Hybrid;
        forEach(std::make_index_sequence<numBlocks>{}, [&](auto i) {
            result[i] = d[Dune::index_constant<i>{}].size();
        });
        return result;
    }

    static DofVector makeDofVector(const DofVector& d)
    { return d; }

    static DofVector makeZeroDofVector(const VectorSizeInfo& size)
    {
        DofVector result;
        using namespace Dune::Hybrid;
        forEach(std::make_index_sequence<numBlocks>{}, [&](auto i) {
            result[Dune::index_constant<i>{}].resize(size[i]);
        });
        return result;
    }
};

namespace Impl {

template<class Vars>
using SolutionVectorType = typename Vars::SolutionVector;

template<class Vars, bool varsExportSolution>
class NewtonVariablesBackend;

/*!
 * \ingroup Nonlinear
 * \brief Class providing Newton operations for scalar/number types
 * \note We assume the variables being simply a dof vector if we
 *       do not find the variables class to export `SolutionVector`.
 */
template<class Vars>
class NewtonVariablesBackend<Vars, false>
: public NewtonDofBackend<Vars>
{
    using ParentType = NewtonDofBackend<Vars>;

public:
    using Variables = Vars;
    using typename ParentType::DofVector;

    static void update(Variables& v, const DofVector& dofs)
    { v = dofs; }

    //! operations on variables
    static const DofVector& getDofVector(const Variables& v)
    { return v; }
};

/*!
 * \file
 * \ingroup Nonlinear
 * \brief Class providing Newton operations for generic variable
 *        types, possibly containing also secondary variables.
 */
template<class Vars>
class NewtonVariablesBackend<Vars, true>
: public NewtonDofBackend<typename Vars::SolutionVector>
{
public:
    using DofVector = typename Vars::SolutionVector;
    using Variables = Vars; //!< the type of the variables object

    static void update(Variables& v, const DofVector& dofs)
    { v.update(dofs); }

    //! operations on variables
    static const DofVector& getDofVector(const Variables& v)
    { return v.dofs(); }
};
} // end namespace Impl

/*!
 * \ingroup Nonlinear
 * \brief Class providing Newton operations for generic variable classes
 *        that represent the state of a numerical solution, possibly
 *        consisting of primary/secondary variables and information on
 *        the time level.
 */
template<class Vars>
using NewtonVariablesBackend = Impl::NewtonVariablesBackend<Vars, Dune::Std::is_detected_v<Impl::SolutionVectorType, Vars>>;

} // end namespace Dumux

#endif
