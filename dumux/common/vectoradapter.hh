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
 * \ingroup Common
 * \brief docme
 */
#ifndef DUMUX_COMMON_VECTOR_ADAPTER_HH
#define DUMUX_COMMON_VECTOR_ADAPTER_HH

#include <memory>
#include <utility>
#include <type_traits>
#include <initializer_list>

#include <dune/common/typetraits.hh>

namespace Dumux {
namespace VectorAdaptation {

template<typename V>
class VectorStorageImpl
{
    static_assert(!std::is_pointer_v<V>, "This expects values or references");

public:
    static constexpr bool isConst = std::is_const_v<std::remove_reference_t<V>>;
    static constexpr bool isReference = std::is_lvalue_reference_v<V>;

    using StoredType = V;
    using Vector = std::decay_t<V>;

    template<bool d = std::is_default_constructible_v<Vector>, std::enable_if_t<d, bool> = true>
    VectorStorageImpl() {}

    template<bool isRef = isReference, std::enable_if_t<isRef, bool> = true>
    VectorStorageImpl(StoredType v) : v_(v) {}

    template<bool isRef = isReference, std::enable_if_t<!isRef, bool> = true>
    VectorStorageImpl(const Vector& v) : v_(v) {}

    template<bool isRef = isReference, std::enable_if_t<!isRef, bool> = true>
    VectorStorageImpl(Vector&& v) : v_(std::move(v)) {}

    template<bool c = isConst, std::enable_if_t<!c, bool> = true>
    Vector& get() { return v_; }
    const Vector& get() const { return v_; }

private:
    StoredType v_;
};

template<typename Vector>
using VectorStorage = VectorStorageImpl<Vector>;

template<typename Vector>
using VectorViewStorage = VectorStorageImpl<const Vector&>;

template<typename Vector>
using VectorProxyStorage = VectorStorageImpl<Vector&>;

//! Generic vector adapter class
template<typename Storage,
         typename Factory,
         typename Access,
         typename Operators>
class Adapter
{
    static constexpr bool isConstNative = Storage::isConst;

public:
    static constexpr bool hasStaticSize = Access::hasStaticSize;
    static constexpr bool isReference = Storage::isReference;

    using Native = typename Storage::Vector;

    //! Construction using the factory
    template<typename... Args,
             bool isRef = isReference, std::enable_if_t<!isRef, bool> = true>
    Adapter(Args&&... args)
    : storage_(Factory::make(std::forward<Args>(args)...))
    {}

    //! Construction from reference
    template<bool isRef = isReference, std::enable_if_t<isRef, bool> = true>
    Adapter(typename Storage::StoredType v)
    : storage_(v)
    {}

    /*!
     * \name Size queries and modifications
     */
    // \{

    template<bool isStatic = hasStaticSize, std::enable_if_t<isStatic, bool> = true>
    static constexpr auto size() { return Access::size(); }

    template<bool isStatic = hasStaticSize, std::enable_if_t<!isStatic, bool> = true>
    auto size() const { return Access::size(native(*this)); }

    template<bool isDynamic = !hasStaticSize, std::enable_if_t<isDynamic, bool> = true>
    void resize(std::size_t n) { return Access::resize(native(*this), n); }

    // \}

    /*!
     * \name Vector access & iterators
     */
    // \{

    template<typename Index>
    decltype(auto) operator[] (const Index& i)
    { return Access::get(native(*this), i); }

    template<typename Index>
    decltype(auto) operator[] (const Index& i) const
    { return Access::get(native(*this), i); }

    decltype(auto) begin() { return Access::begin(native(*this)); }
    decltype(auto) begin() const { return Access::begin(native(*this)); }

    decltype(auto) end() { return Access::end(native(*this)); }
    decltype(auto) end() const { return Access::end(native(*this)); }

    // \}

    /*!
     * \name Assignment operations
     */
    // \{

    template<typename OtherVector>
    Adapter& operator=(const OtherVector& other)
    { Factory::assign(native(*this), other); return *this; }

    template<typename... Args>
    Adapter& operator=(const Adapter<Args...>& other)
    { Factory::assign(native(*this), native(other)); return *this; }

    // \}


    /*!
     * \name Operators interface
     */
    // \{

    //! Multiplication with a scalar
    template<typename Scalar, std::enable_if_t<Dune::IsNumber<Scalar>::value, bool> = true>
    Adapter& operator*=(const Scalar& value)
    { Operators::multiplyWithScalar(native(*this), value); return *this; }

    //! Division by a scalar
    template<typename Scalar, std::enable_if_t<Dune::IsNumber<Scalar>::value, bool> = true>
    Adapter& operator/=(const Scalar& value)
    { Operators::divideByScalar(native(*this), value); return *this; }

    //! add another native vector
    // (TODO: What about e.g. BlockVector with different BlockType?)
    Adapter& operator+=(const Native& other)
    { Operators::addVector(native(*this), other); return *this; }

    //! add another adapted vector
    template<typename... Args>
    Adapter& operator+=(const Adapter<Args...>& other)
    { Operators::addVector(native(*this), native(other)); return *this; }

    //! subtract another native vector
    // (TODO: What about e.g. BlockVector with different BlockType?)
    Adapter& operator-=(const Native& other)
    { Operators::subtractVector(native(*this), other); return *this; }

    //! subtract another adapted vector
    template<typename... Args>
    Adapter& operator-=(const Adapter<Args...>& other)
    { Operators::subtractVector(native(*this), native(other)); return *this; }

    // \}

    //! Return reference to the native vector type
    template<bool isConst = isConstNative, std::enable_if_t<!isConst, bool> = true>
    friend Native& native(Adapter& a)
    { return a.storage_.get(); }

    //! Return const reference to the native vector type
    friend const Native& native(const Adapter& a)
    { return a.storage_.get(); }

private:
    Storage storage_;
};

// Traits to be specialized by different vector types
namespace Traits {

template<typename Vector> struct VectorFactory;
template<typename Vector> struct VectorAccess;
template<typename Vector> struct VectorOperators;

} // end namespace Traits
} // end namespace VectorAdaptation

//! Vector adapter using specialized traits
template<typename Vector>
using VectorAdapter = VectorAdaptation::Adapter<
    typename VectorAdaptation::VectorStorage<Vector>,
    typename VectorAdaptation::Traits::VectorFactory<Vector>::type,
    typename VectorAdaptation::Traits::VectorAccess<Vector>::type,
    typename VectorAdaptation::Traits::VectorOperators<Vector>::type
>;

//! Vector view adapter class using specialized traits
template<typename Vector>
using VectorViewAdapter = VectorAdaptation::Adapter<
    typename VectorAdaptation::VectorViewStorage<Vector>,
    typename VectorAdaptation::Traits::VectorFactory<Vector>::type,
    typename VectorAdaptation::Traits::VectorAccess<Vector>::type,
    typename VectorAdaptation::Traits::VectorOperators<Vector>::type
>;

//! Vector proxy adapter class using specialized traits
template<typename Vector>
using VectorProxyAdapter = VectorAdaptation::Adapter<
    typename VectorAdaptation::VectorProxyStorage<Vector>,
    typename VectorAdaptation::Traits::VectorFactory<Vector>::type,
    typename VectorAdaptation::Traits::VectorAccess<Vector>::type,
    typename VectorAdaptation::Traits::VectorOperators<Vector>::type
>;

// default implementations of the vector adapter helpers
namespace VectorAdaptation {

template<typename Vector>
class DefaultStaticVectorFactory
{
public:
    static Vector make() { return Vector{}; }

    template<typename T>
    static Vector make(const std::initializer_list<T>& v)
    { return {v}; }

    template<typename OtherVector>
    static void assign(Vector& v, const OtherVector& other)
    { v = other; }
};

template<typename Vector>
class DefaultDynamicVectorFactory
: public DefaultStaticVectorFactory<Vector>
{
public:
    using DefaultStaticVectorFactory<Vector>::make;

    static Vector make(std::size_t n)
    { return {n}; }
};


namespace Detail {

template<typename Vector>
class DefaultRange
{
public:
    static auto begin(Vector& v) { return v.begin(); }
    static auto begin(const Vector& v) { return v.begin(); }

    static auto end(Vector& v) { return v.end(); }
    static auto end(const Vector& v) { return v.end(); }
};

template<typename Vector>
class DefaultAccess
{
public:
    template<typename Index>
    static decltype(auto) get(Vector& v, Index i)
    { return v[i]; }

    template<typename Index>
    static decltype(auto) get(const Vector& v, Index i)
    { return v[i]; }
};

} // end namespace Detail

template<typename Vector>
class DefaultDynamicVectorAccess
: public Detail::DefaultRange<Vector>
, public Detail::DefaultAccess<Vector>
{
public:
    static constexpr bool hasStaticSize = false;

    static auto size(const Vector& v)
    { return v.size(); }

    static auto resize(Vector& v, std::size_t n)
    { v.resize(n); }
};

template<typename Vector>
class DefaultStaticVectorAccess
: public Detail::DefaultRange<Vector>
, public Detail::DefaultAccess<Vector>
{
public:
    static constexpr bool hasStaticSize = true;

    static auto size(const Vector& v) { return Vector::size(); }
    static constexpr auto size() { return Vector::size(); }
};

template<typename Vector>
class DefaultVectorOperators
{
public:
    template<typename Scalar>
    static void multiplyWithScalar(Vector& v, const Scalar& value)
    { v *= value; }

    template<typename Scalar>
    static void divideByScalar(Vector& v, const Scalar& value)
    { v /= value; }

    template<typename OtherVector>
    static void addVector(Vector& v, const OtherVector& other)
    { v += other; }

    template<typename OtherVector>
    static void subtractVector(Vector& v, const OtherVector& other)
    { v -= other; }
};

} // end namespace VectorAdaptation
} // end namespace Dumux


// specialize the traits for Dune vectors
namespace Dune {

template<typename K, int size> class FieldVector;
template<typename K, typename A> class DynamicVector;
template<typename BT, typename A> class BlockVector;
template<typename... Rows> class MultiTypeBlockVector;

} // end namespace Dune

namespace Dumux::VectorAdaptation::Traits {

template<typename BT>
struct VectorFactory<Dune::BlockVector<BT>>
{ using type = DefaultDynamicVectorFactory<Dune::BlockVector<BT>>; };

template<typename FT>
struct VectorFactory<Dune::DynamicVector<FT>>
{ using type = DefaultDynamicVectorFactory<Dune::DynamicVector<FT>>; };

template<typename FT, int size>
struct VectorFactory<Dune::FieldVector<FT, size>>
{ using type = DefaultStaticVectorFactory<Dune::FieldVector<FT, size>>; };

template<typename... Rows>
struct VectorFactory<Dune::MultiTypeBlockVector<Rows...>>
{ using type = DefaultStaticVectorFactory<Dune::MultiTypeBlockVector<Rows...>>; };

template<typename BT>
struct VectorAccess<Dune::BlockVector<BT>>
{ using type = DefaultDynamicVectorAccess<Dune::BlockVector<BT>>; };

template<typename FT>
struct VectorAccess<Dune::DynamicVector<FT>>
{ using type = DefaultDynamicVectorAccess<Dune::DynamicVector<FT>>; };

template<typename FT, int size>
struct VectorAccess<Dune::FieldVector<FT, size>>
{ using type = DefaultStaticVectorAccess<Dune::FieldVector<FT, size>>; };

template<typename... Rows>
struct VectorAccess<Dune::MultiTypeBlockVector<Rows...>>
{ using type = DefaultStaticVectorAccess<Dune::MultiTypeBlockVector<Rows...>>; };

template<typename Vector>
struct VectorOperators
{ using type = DefaultVectorOperators<Vector>; };

} // end namespace Dumux::VectorAdaptation::Traits

#endif
