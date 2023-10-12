// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup Core
 * \brief Backends for operations on different solution vector types
 *        or more generic variable classes to be used in places where
 *        several different types/layouts should be supported.
 */
#ifndef DUMUX_COMMON_VARIABLES_BACKEND_HH
#define DUMUX_COMMON_VARIABLES_BACKEND_HH

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
 * \ingroup Core
 * \brief Class providing operations with primary variable vectors
 */
template<class DofVector, bool isScalar = Dune::IsNumber<DofVector>::value>
class DofBackend;

/*!
 * \ingroup Core
 * \brief Specialization providing operations for scalar/number types
 */
template<class Scalar>
class DofBackend<Scalar, true>
{
public:
    using DofVector = Scalar; //!< the type of the dofs parametrizing the variables object
    using SizeType = std::size_t;

    //! Return the number of entries in the dof vector
    static SizeType size(const DofVector& d)
    { return 1; }

    //! Make a zero-initialized dof vector instance
    static DofVector zeros(SizeType size)
    { return 0.0; }

    //! Perform axpy operation (y += a * x)
    template<class OtherDofVector>
    static void axpy(Scalar a, const OtherDofVector& x, DofVector& y)
    { y += a*x; }
};

/*!
 * \ingroup Core
 * \brief Specialization providing operations for block vectors
 * \tparam Vector a type that is
 *   - default-constructible
 *   - has size() member
 *   - has resize(0) member
 *   - has axpy(a, x) member
 */
template<class Vector>
class DofBackend<Vector, false>
{
public:
    using DofVector = Vector; //!< the type of the dofs parametrizing the variables object
    using SizeType = std::size_t;

    //! Return the number of entries in the dof vector
    static SizeType size(const DofVector& d)
    { return d.size(); }

    //! Make a zero-initialized dof vector instance
    static DofVector zeros(SizeType size)
    { DofVector d; d.resize(size); return d; }

    //! Perform axpy operation (y += a * x)
    template<class OtherDofVector>
    static void axpy(typename DofVector::field_type a, const OtherDofVector& x, DofVector& y)
    {
        for (typename DofVector::size_type i = 0; i < y.size(); ++i)
            y[i].axpy(a, x[i]);
    }
};

/*!
 * \ingroup Core
 * \brief Specialization providing operations for multitype block vectors
 */
template<class... Blocks>
class DofBackend<Dune::MultiTypeBlockVector<Blocks...>, false>
{
    using DV = Dune::MultiTypeBlockVector<Blocks...>;
    static constexpr auto numBlocks = DV::size();

    using VectorSizeInfo = std::array<std::size_t, numBlocks>;

public:
    using DofVector = DV; //!< the type of the dofs parametrizing the variables object
    using SizeType = VectorSizeInfo;

    //! Return the number of entries in the sub-dof-vectors
    static SizeType size(const DofVector& d)
    {
        VectorSizeInfo result;
        using namespace Dune::Hybrid;
        forEach(std::make_index_sequence<numBlocks>{}, [&](auto i) {
            result[i] = d[Dune::index_constant<i>{}].size();
        });
        return result;
    }

    //! Make a zero-initialized dof vector instance
    static DofVector zeros(const SizeType& size)
    {
        DofVector result;
        using namespace Dune::Hybrid;
        forEach(std::make_index_sequence<numBlocks>{}, [&](auto i) {
            result[Dune::index_constant<i>{}].resize(size[i]);
        });
        return result;
    }

    //! Perform axpy operation (y += a * x)
    template<class Scalar, class OtherDofVector, std::enable_if_t< Dune::IsNumber<Scalar>::value, int> = 0>
    static void axpy(Scalar a, const OtherDofVector& x, DofVector& y)
    {
        using namespace Dune::Hybrid;
        forEach(std::make_index_sequence<numBlocks>{}, [&](auto i) {
            DofBackend<std::decay_t<decltype(y[Dune::index_constant<i>{}])>>::axpy(
                a, x[Dune::index_constant<i>{}], y[Dune::index_constant<i>{}]
            );
        });
    }
};

namespace Detail {

template<class Vars>
using SolutionVectorType = typename Vars::SolutionVector;

template<class Vars, bool varsExportSolution>
class VariablesBackend;

/*!
 * \ingroup Core
 * \brief Class providing operations for primary variable vector/scalar types
 * \note We assume the variables being simply a dof vector if we
 *       do not find the variables class to export `SolutionVector`.
 */
template<class Vars>
class VariablesBackend<Vars, false>
: public DofBackend<Vars>
{
    using ParentType = DofBackend<Vars>;

public:
    using Variables = Vars;
    using typename ParentType::DofVector;

    //! update to new solution vector
    static void update(Variables& v, const DofVector& dofs)
    { v = dofs; }

    //! return const reference to dof vector
    static const DofVector& dofs(const Variables& v)
    { return v; }

    //! return reference to dof vector
    static DofVector& dofs(Variables& v)
    { return v; }
};

/*!
 * \ingroup Core
 * \brief Class providing operations for generic variable classes,
 *        containing primary and possibly also secondary variables.
 */
template<class Vars>
class VariablesBackend<Vars, true>
: public DofBackend<typename Vars::SolutionVector>
{
public:
    using DofVector = typename Vars::SolutionVector;
    using Variables = Vars; //!< the type of the variables object

    //! update to new solution vector
    static void update(Variables& v, const DofVector& dofs)
    { v.update(dofs); }

    //! return const reference to dof vector
    static const DofVector& dofs(const Variables& v)
    { return v.dofs(); }

    //! return reference to dof vector
    static DofVector& dofs(Variables& v)
    { return v.dofs(); }
};
} // end namespace Detail

/*!
 * \ingroup Core
 * \brief Class providing operations for generic variable classes
 *        that represent the state of a numerical solution, possibly
 *        consisting of primary/secondary variables and information on
 *        the time level.
 */
template<class Vars>
using VariablesBackend = Detail::VariablesBackend<Vars, Dune::Std::is_detected_v<Detail::SolutionVectorType, Vars>>;

} // end namespace Dumux

#endif
