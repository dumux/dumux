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
#ifndef DUMUX_COMMON_VECTOR_WITH_STATE_HH
#define DUMUX_COMMON_VECTOR_WITH_STATE_HH

#include <vector>
#include <utility>
#include <type_traits>
#include <initializer_list>

#include <dumux/common/vectoradapter.hh>

namespace Dumux {

namespace Detail {

template<typename Vector, typename State, bool isPrimaryVariablesVector>
class ProxyStorageImpl;

template<typename T>
using IndexedType = decltype(std::declval<T>()[0]);

template<typename T>
inline constexpr bool isPrimaryVariablesVector
    = Dune::IsNumber<std::decay_t<IndexedType<T>>>::value;

template<typename Vector, typename State>
using ProxyStorage = ProxyStorageImpl<
    Vector,
    State,
    isPrimaryVariablesVector<Vector> // TODO: requires operator[]
>;

} // end namespace Detail

//! vectors with state implementation
template<typename VectorAdapter,
         typename State,
         typename ProxyStorage>
class VectorWithStateImpl
{
    static constexpr bool useProxies = ProxyStorage::useProxies;
    static constexpr bool indexableState = Dune::IsIndexable<State>::value;

    using RawState = std::decay_t<State>;

public:
    static constexpr bool hasStaticSize = VectorAdapter::hasStaticSize;
    static constexpr bool isReference = VectorAdapter::isReference;

    using Native = typename VectorAdapter::Native;

    template<typename TheState, typename... Args,
             std::enable_if_t<std::is_convertible_v<TheState, State>, bool> = true>
    VectorWithStateImpl(TheState&& state, Args&&... args)
    : vectorAdapter_(std::forward<Args>(args)...)
    , state_(std::forward<State>(state))
    { updateProxies_(); }

    template<typename... Args, typename S = State,
             std::enable_if_t<
                (std::is_default_constructible_v<S> && std::is_constructible_v<VectorAdapter, Args...>)
             , bool> = true>
    VectorWithStateImpl(Args&&... args)
    : vectorAdapter_(std::forward<Args>(args)...)
    {
        if constexpr (indexableState)
            state_.resize(size()); // TODO: Outsource as well? This requires dynamic state vectors
        updateProxies_();
    }

    /*!
     * \name Size queries and modifications
     */
    // \{

    template<bool isStatic = hasStaticSize, std::enable_if_t<isStatic, bool> = true>
    static constexpr auto size() { return VectorAdapter::size(); }

    template<bool isStatic = hasStaticSize, std::enable_if_t<!isStatic, bool> = true>
    auto size() const { return vectorAdapter_.size(); }

    template<bool isStatic = hasStaticSize, std::enable_if_t<!isStatic, bool> = true>
    void resize(std::size_t n)
    {
        vectorAdapter_.resize(n);
        if constexpr (indexableState)
            state_.resize(size()); // TODO: Outsource as well? This requires dynamic state vectors
        updateProxies_();
    }

    // \}

    /*!
     * \name Vector access & iterators
     */
    // \{

    template<typename Index>
    decltype(auto) operator[] (const Index& i)
    {
        if constexpr (useProxies)
            return proxyStorage_.get(i);
        else
            return vectorAdapter_[i];
    }

    template<typename Index>
    decltype(auto) operator[] (const Index& i) const
    {
        if constexpr (useProxies)
            return proxyStorage_.get(i);
        else
            return vectorAdapter_[i];
    }

    decltype(auto) begin()
    {
        if constexpr (useProxies)
            return proxyStorage_.begin();
        else
            return vectorAdapter_.begin();
    }

    decltype(auto) begin() const
    {
        if constexpr (useProxies)
            return proxyStorage_.constBegin();
        else
            return vectorAdapter_.begin();
    }

    decltype(auto) end()
    {
        if constexpr (useProxies)
            return proxyStorage_.end();
        else
            return vectorAdapter_.end();
    }

    decltype(auto) end() const
    {
        if constexpr (useProxies)
            return proxyStorage_.constEnd();
        else
            return vectorAdapter_.end();
    }

    // \}

    /*!
     * \name Assignment operations
     */
    // \{

    // TODO: This should not be permitted as it does not modify the states?
    //       But we probably do this somewhere in dumux also with stateful vectors
    template<typename OtherVector>
    VectorWithStateImpl& operator=(const OtherVector& other)
    {
        vectorAdapter_ = other;
        updateProxies_();
        return *this;
    }

    template<typename V, typename S, typename B>
    VectorWithStateImpl& operator=(const VectorWithStateImpl<V, S, B>& other)
    {
        vectorAdapter_ = native(other);
        state_ = other.state_;
        updateProxies_();
        return *this;
    }

    // \}


    /*!
     * \name Operators interface
     */
    // \{

    //! Multiplication with a scalar
    template<typename Scalar, std::enable_if_t<Dune::IsNumber<Scalar>::value, bool> = true>
    VectorWithStateImpl& operator*=(const Scalar& value) { vectorAdapter_ *= value; return *this; }

    //! Division by a scalar
    template<typename Scalar, std::enable_if_t<Dune::IsNumber<Scalar>::value, bool> = true>
    VectorWithStateImpl& operator/=(const Scalar& value) { vectorAdapter_ /= value; return *this; }

    //! add another native vector
    // (TODO: What about e.g. BlockVector with different BlockType?)
    VectorWithStateImpl& operator+=(const Native& other) { vectorAdapter_ += other; return *this; }

    //! add another vector with state (TODO: this sould be prohibited? Adding states??)
    template<typename V, typename S, typename B>
    VectorWithStateImpl& operator+=(const VectorWithStateImpl<V, S, B>& other)
    { vectorAdapter_ += native(other); return *this; }

    //! subtract another native vector
    // (TODO: What about e.g. BlockVector with different BlockType?)
    VectorWithStateImpl& operator-=(const Native& other)
    { vectorAdapter_ -= other; return *this; }

    //! subtract another vector with state (TODO: this sould be prohibited? Subtracting states??)
    template<typename V, typename S, typename B>
    VectorWithStateImpl& operator-=(const VectorWithStateImpl<V, S, B>& other)
    { vectorAdapter_ -= native(other); return *this; }

    // \}

    RawState& state() { return state_; }
    const RawState& state() const { return state_; }
    void setState(const RawState& s) { state_ = s; }

    friend Native& native(VectorWithStateImpl& v)
    { return native(v.vectorAdapter_); }

    friend const Native& native(const VectorWithStateImpl& v)
    { return native(v.vectorAdapter_); }

private:
    void updateProxies_()
    {
        if constexpr (useProxies)
            proxyStorage_.update(native(vectorAdapter_), state_);
    }

    VectorAdapter vectorAdapter_;
    ProxyStorage proxyStorage_;
    State state_;
};

//! Convenience alias for vectors with state
template<typename Vector, typename State, typename Adapter = VectorAdapter<Vector>>
using VectorWithState = VectorWithStateImpl<
    Adapter, State, typename Detail::ProxyStorage<Vector, State>
>;


namespace Detail {

template<typename Vector, typename State>
class ProxyStorageImpl<Vector, State, /*isPrimaryVariablesVector*/true>
{
public:
    static constexpr bool useProxies = false;
    int get(std::size_t i) { return throw_(); }
    int get(std::size_t i) const { return throw_(); }

private:

    int throw_()
    {
        DUNE_THROW(Dune::InvalidStateException,
                   "Empty proxy storage shouldn't be called");
    }

};

template<typename Vector, typename State>
class ProxyStorageImpl<Vector, State, /*isPrimaryVariablesVector*/false>
{
    template<typename T>
    using ConstRef = std::add_lvalue_reference_t<
        std::add_const_t<std::remove_reference_t<T>>
    >;

    static_assert(Dune::IsIndexable<Vector>::value);
    static_assert(Dune::IsIndexable<State>::value);

    using VectorReturn = IndexedType<Vector>;
    using StateReturn = IndexedType<State>;
    using ConstVectorReturn = ConstRef<VectorReturn>;
    using ConstStateReturn = ConstRef<StateReturn>;

    static_assert(std::is_lvalue_reference_v<VectorReturn>);
    static_assert(std::is_lvalue_reference_v<StateReturn>);
    static_assert(std::is_lvalue_reference_v<ConstVectorReturn>);
    static_assert(std::is_lvalue_reference_v<ConstStateReturn>);

    using VectorReturnDecayType = std::decay_t<VectorReturn>;
    using StateReturnDecayType = std::decay_t<StateReturn>;

    // TODO: This now hardcodes the adapter using the traits...
    using Proxy = VectorWithStateImpl<
        VectorProxyAdapter<VectorReturnDecayType>,
        StateReturn,
        ProxyStorage<VectorReturnDecayType, StateReturnDecayType>
    >;

    using ConstProxy = VectorWithStateImpl<
        VectorViewAdapter<VectorReturnDecayType>,
        ConstStateReturn,
        ProxyStorage<VectorReturnDecayType, StateReturnDecayType>
    >;

public:
    static constexpr bool useProxies = true;

    // TODO: state/vector access??
    void update(Vector& v, State& s)
    {
        proxy_.clear();
        proxy_.reserve(v.size());
        for (std::size_t i = 0; i < v.size(); ++i)
            proxy_.emplace_back(s[i], v[i]);
        updateConst_(v, s);
    }

    Proxy& get(std::size_t i)
    { return proxy_[i]; }

    const ConstProxy& get(std::size_t i) const
    { return constProxy_[i]; }

    auto begin() { return proxy_.begin(); }
    auto end() { return proxy_.end(); }

    auto constBegin() const { return constProxy_.begin(); }
    auto constEnd() const { return constProxy_.end(); }

private:
    // TODO: state/vector access??
    void updateConst_(const Vector& v, const State& s)
    {
        constProxy_.clear();
        constProxy_.reserve(v.size());
        for (std::size_t i = 0; i < v.size(); ++i)
            constProxy_.emplace_back(s[i], v[i]);
    }

    std::vector<Proxy> proxy_;
    std::vector<ConstProxy> constProxy_;
};

} // end namespace Detail
} // end namespace Dumux

#endif
