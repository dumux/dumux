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
#ifndef DUMUX_COMMON_SOLUTIONVECTOR_HH
#define DUMUX_COMMON_SOLUTIONVECTOR_HH

#include <dumux/common/numeqvector.hh>
#include <dune/istl/bvector.hh>
#include <dune/common/std/type_traits.hh>

namespace Dumux {

// forward declare
template<class PV, class S>
class SwitchablePrimaryVariables;

template<class PV>
class NumEqVectorTraits;

} // end namespace Dumux

namespace Dumux::Istl {

// TODO: can probably be replaced with class BlockVector : public BlockVectorType
// This does not really add any value. Maybe delete completely and just use the dune type directly.
template<class BlockVectorType>
class BlockVector
{
public:
    using field_type = typename BlockVectorType::field_type;
    using block_type = typename BlockVectorType::block_type;
    using allocator_type = typename BlockVectorType::allocator_type;
    using size_type = typename BlockVectorType::size_type;
    using Iterator = typename BlockVectorType::Iterator;
    using ConstIterator = typename BlockVectorType::ConstIterator;

    BlockVector() = default;

    //! make vector with _n components
    explicit BlockVector (size_type n) : blockVector_(n)
    {}

    /** \brief Construct from a std::initializer_list */
    BlockVector (const std::initializer_list<block_type>& l) : blockVector_(l)
    {}

    template<typename S>
    BlockVector (size_type n, S capacity) : blockVector_(n, capacity)
    {}

    block_type& operator [](int i)
    { return blockVector_[i]; }

    const block_type& operator[] (int i) const
    { return blockVector_[i]; }

    BlockVector& operator= (const field_type& k)
    {
        blockVector_ = k;
        return *this;
    }

    BlockVector& operator*= (const field_type& k)
    {
        blockVector_ *= k;
        return *this;
    }

    BlockVector& operator/= (const field_type& k)
    {
        blockVector_ /= k;
        return *this;
    }

    BlockVector& operator+= (const BlockVector& y)
    {
        blockVector_ += y.blockVector_;
        return *this;
    }

    BlockVector& operator-= (const BlockVector& y)
    {
        blockVector_ -= y.blockVector_;
        return *this;
    }

    BlockVector operator* (const BlockVector& y)
    {
        return blockVector_ * y.blockVector_;
    }

    BlockVector& axpy(const field_type& a, const BlockVector& y)
    {
        blockVector_.axpy(a, y.blockVector_);
        return *this;
    }

    BlockVector& dot(const BlockVector& y)
    {
        blockVector_.dot(y.blockVector_);
        return *this;
    }

    auto one_norm () const
    { return blockVector_.one_norm(); }

    auto two_norm () const
    { return blockVector_.two_norm(); }

    auto two_norm2 () const
    { return blockVector_.two_norm2(); }

    auto infinity_norm () const
    { return blockVector_.infinity_norm(); }

    auto begin() const
    { return blockVector_.begin(); }

    auto end() const
    { return blockVector_.end(); }

    void reserve(size_type capacity)
    { blockVector_.reserve(capacity); }

    size_type capacity() const
    { return blockVector_.capacity(); }

    size_type size() const
    { return blockVector_.size(); }

    void resize(size_type size)
    { blockVector_.resize(size); }

    BlockVectorType& native()
    { return blockVector_; }

    const BlockVectorType& native() const
    { return blockVector_; }

private:
    BlockVectorType blockVector_;
};

// forward declare
template<class BlockType, class PriVarsType, class State>
struct BlockVectorView;

template<class BlockType, class State>
struct ConstBlockVectorView;

template<class BlockType, class State>
struct ConstBlockVectorView
{
    template<class B, class P, class S>
    friend struct BlockVectorView;

    ConstBlockVectorView(const BlockType& priVars, const State& state) : priVars_(&priVars), state_(&state)
    {}

    operator BlockType() const
    { return *priVars_; }

    State state() const
    { return *state_; }

    const auto& operator[] (int i) const
    { return (*priVars_)[i]; }

private:
    const BlockType* priVars_;
    const State* state_;
};

template<class BlockType, class PriVarsType, class State>
struct BlockVectorView
{
    BlockVectorView() : priVars_(nullptr), state_(nullptr) {}

    BlockVectorView(BlockType& priVars, State& state)
    {
        priVars_ = &priVars;
        state_ = &state;
    }

    auto& operator= (const ConstBlockVectorView<BlockType, State>& other)
    {
        if (priVars_)
        {
            *priVars_ = (*other.priVars_);
            *priVars_ = other.state();
        }
        else
        {
            priVarsStored_ = std::make_unique<BlockType>((*other.priVars_));
            stateStored_ = other.state();
        }

        return *this;
    }

    auto& operator= (const PriVarsType& other)
    {
        if (priVars_)
        {
            *priVars_ = other;
            *state_ = other.state();
        }
        else
        {
            priVarsStored_ = std::make_unique<BlockType>(other);
            stateStored_ = other.state();
        }

        return *this;
    }

    void setState(int s)
    { *state_ = s; }

    operator PriVarsType() const
    {
        PriVarsType privars;
        if (priVars_)
        {
            privars = *priVars_;
            privars.setState(*state_);
        }
        else
        {
            privars = *priVarsStored_;
            privars.setState(stateStored_);
        }

        return privars;
    }

    auto& operator[] (int i)
    {
        if (priVars_)
            return (*priVars_)[i];
        else
            return (*priVarsStored_)[i];
    }

    const auto& operator[] (int i) const
    {
        if (priVars_)
            return (*priVars_)[i];
        else
            return (*priVarsStored_)[i];
    }

    static constexpr auto size()
    { return BlockType::size(); }

    State state() const
    {
        if (priVars_)
            return *state_;
        else
            return stateStored_;
    }

private:
    BlockType* priVars_;
    State* state_;
    State stateStored_;
    std::unique_ptr<BlockType> priVarsStored_;
};

template<class BVType, class PriVarsType, class StateVectorType>
class BlockVectorWithState
{
    using BlockVectorType = std::decay_t<BVType>;
    using State = typename std::decay_t<StateVectorType>::value_type;
public:
    using field_type = typename BlockVectorType::field_type;
    using block_type = typename BlockVectorType::block_type;
    using allocator_type = typename BlockVectorType::allocator_type;
    using size_type = typename BlockVectorType::size_type;
    using Iterator = typename BlockVectorType::Iterator;
    using ConstIterator = typename BlockVectorType::ConstIterator;

    BlockVectorWithState() = default;

    //! make vector with _n components
    explicit BlockVectorWithState (size_type n) : blockVector_(n)
    {
        states_.resize(blockVector_.size());
    }

    /** \brief Construct from a std::initializer_list */
    BlockVectorWithState (const std::initializer_list<block_type>& l) : blockVector_(l)
    {
        states_.resize(blockVector_.size());
    }

    template<typename S>
    BlockVectorWithState (size_type n, S capacity) : blockVector_(n, capacity)
    {
        states_.resize(blockVector_.size());
    }

    //! constructor for using this class as a view (member variables should be references)
    //! TODO add some enable_if (also to other ctors) magic to prevent misuse
    BlockVectorWithState(BlockVectorType& otherDofs, std::vector<State>& otherStates) : blockVector_(otherDofs), states_(otherStates)
    {
        static_assert(std::is_lvalue_reference_v<decltype(blockVector_)>);
    }

    auto operator [](size_type i)
    {
        return BlockVectorView<block_type, PriVarsType, State>(blockVector_[i], states_[i]);
    }

    auto operator[] (size_type i) const
    {
        return ConstBlockVectorView<block_type, State>(blockVector_[i], states_[i]);
    }

    BlockVectorWithState& operator= (const field_type& k)
    {
        blockVector_ = k;
        return *this;
    }

    BlockVectorWithState& operator*= (const field_type& k)
    {
        blockVector_ *= k;
        return *this;
    }

    BlockVectorWithState& operator/= (const field_type& k)
    {
        blockVector_ /= k;
        return *this;
    }

    BlockVectorWithState& operator+= (const BlockVectorWithState& y)
    {
        blockVector_ += y.blockVector_;
        return *this;
    }

    BlockVectorWithState& operator-= (const BlockVectorWithState& y)
    {
        blockVector_ -= y.blockVector_;
        return *this;
    }

    BlockVectorWithState operator* (const BlockVectorWithState& y)
    {
        return blockVector_ * y.blockVector_;
    }

    BlockVectorWithState& axpy(const field_type& a, const BlockVectorWithState& y)
    {
        blockVector_.axpy(a, y.blockVector_);
        return *this;
    }

    BlockVectorWithState& dot(const BlockVectorWithState& y)
    {
        blockVector_.dot(y.blockVector_);
        return *this;
    }

    auto one_norm () const
    { return blockVector_.one_norm(); }

    auto two_norm () const
    { return blockVector_.two_norm(); }

    auto two_norm2 () const
    { return blockVector_.two_norm2(); }

    auto infinity_norm () const
    { return blockVector_.infinity_norm(); }

    auto begin() const
    { return blockVector_.begin(); }

    auto end() const
    { return blockVector_.end(); }

    void reserve(size_type capacity)
    {
        blockVector_.reserve(capacity);
        states_.reserve(capacity);
    }

    size_type capacity() const
    { return blockVector_.capacity(); }

    size_type size() const
    { return blockVector_.size(); }

    void resize(size_type size)
    {
        blockVector_.resize(size);
        states_.resize(size);
    }

    //! Returns a reference to the stored DOFs.
    BlockVectorType& native()
    { return blockVector_; }

    //! Returns a reference to the stored DOFs.
    const BlockVectorType& native() const
    { return blockVector_; }

    //! Returns a deep copy of the stored DOFs.
    BlockVectorType nativeDeepCopy() const
    { return blockVector_; }

private:
    BVType blockVector_; // BVType might be a reference
    StateVectorType states_; // StateVectorType might be a reference
};

namespace Detail {

template<class PrimaryVariables>
struct StateVectorHelper
{
    using type = void;
};

template<class PV, class S>
struct StateVectorHelper<Dumux::SwitchablePrimaryVariables<PV, S>>
{
    // using PrimaryVariables = Dumux::SwitchablePrimaryVariables<PV, S>;
    // using State = decltype(std::declval<PrimaryVariables>.state());
    // using type = std::vector<PrimaryVariables, State>;
};

}

template<class MTBVType, class... PriVarsTypes>
class MultiTypeBlockVectorWithState : public MTBVType
{
//     using StateVectors = std::tuple<Detail::StateVectorHelper<PriVarsTypes>...>;

// public:
//     // using size_type = typename MTBVType::size_type;
//     using MTBVType::MTBVType;
//     using MTBVType::operator=;
//     using MTBVType::operator-=;
//     using MTBVType::operator+=;
//     using MTBVType::operator*=;
//     using MTBVType::operator/=;
//     using MTBVType::operator[];


// private:
//     // MTBVType multiTypeBlockVector_;
//     StateVectors stateVectors_;
};

} // end namespace Dumux::Istl

namespace Dumux {

// // forward declare
// template<class PV, class S>
// class SwitchablePrimaryVariables;

// template<class PV>
// class NumEqVectorTraits;

namespace Detail {

template <class T>
using NativeStorageDetector = decltype(std::declval<T>().native());

template<class T>
static constexpr bool hasNativeStorage()
{ return Dune::Std::is_detected<NativeStorageDetector, T>::value; }

template <class T>
using StateDetector = decltype(std::declval<T>().state());

template<class T>
static constexpr bool hasState()
{ return Dune::Std::is_detected<StateDetector, T>::value; }

//! Use a standard Dune::BlockVector for PrimaryVariables without state.
template<class... PrimaryVariables>
struct SolutionVectorHelper
{
    // static_assert(sizeof...(PrimaryVariables) == 1);
    using type = Dune::BlockVector<PrimaryVariables...>;
};

//! Use a wrapped type for PrimaryVariables without state.
template<class PV, class S>
struct SolutionVectorHelper <Dumux::SwitchablePrimaryVariables<PV, S>>
{
    using type = Dumux::Istl::BlockVectorWithState<Dune::BlockVector<PV>, Dumux::SwitchablePrimaryVariables<PV, S>, std::vector<S>>;
};

template<class... PrimaryVariables>
struct MultiTypeSolutionVectorHelper
{
    // TODO make distinction between regular privars and such with state
    using type0 = Dune::MultiTypeBlockVector<Dune::BlockVector<typename NumEqVectorTraits<PrimaryVariables>::type>...>;

    using type = Dumux::Istl::MultiTypeBlockVectorWithState<type0, PrimaryVariables...>;
    // using type = Dune::MultiTypeBlockVector<Dune::BlockVector<PrimaryVariables>...>;
};

template<int i, class... PrimaryVariables>
struct SolutionVectorChooser
{
    using type = typename MultiTypeSolutionVectorHelper<PrimaryVariables...>::type;
};

template<class PrimaryVariables>
struct SolutionVectorChooser<1, PrimaryVariables>
{
    using type = typename SolutionVectorHelper<PrimaryVariables>::type;
};



// Helper struct containing the native function for standard solution vectors.
template<class SolutionVector>
struct NativeChooser
{
    template<class SV>
    static decltype(auto) native(SV&& sol)
    {
        using PrimaryVariable = std::decay_t<decltype(sol[0])>;

        if constexpr (Detail::hasNativeStorage<SolutionVector>())
            return sol.native();
        else if constexpr (!Detail::hasState<PrimaryVariable>())
            return sol;
        else
        {
            using BlockType = typename Dumux::NumEqVectorTraits<PrimaryVariable>::type;
            Dune::BlockVector<BlockType> result(sol.size());
            for (auto i = 0; i < sol.size(); ++i)
                result[i] = sol[i];

            return result;
        }
    }
};

// Helper struct containing the native function for MultiTypeBlockVector solution vectors
template<class... Args>
struct NativeChooser<Dumux::Istl::MultiTypeBlockVectorWithState<Args...>>
{
    template<class SV>
    static decltype(auto) native(SV&& sol)
    {
        return sol;
    }
};


} // end namespace Detail

//! Helper type to determine whether a given type is a Dune::MultiTypeBlockVector
template<class... T>
struct isMultiTypeBlockVector<Dumux::Istl::MultiTypeBlockVectorWithState<T...> > : public std::true_type {};

//! Helper alias to get the correct solution vector based on the chosen PrimaryVariables.
template<class ...PrimaryVariables>
using SolutionVector = typename Detail::SolutionVectorChooser<sizeof...(PrimaryVariables), PrimaryVariables...>::type;

template<class SolutionVector>
struct SolutionVectorTraits
{
    using NativeType = SolutionVector;
};

template<class PV, class S>
struct SolutionVectorTraits<Dune::BlockVector<SwitchablePrimaryVariables<PV, S>>>
{
    using NativeType = Dune::BlockVector<PV>;
};

template<class BlockVectorType, class PrivarsType, class StateVectorType>
struct SolutionVectorTraits<Dumux::Istl::BlockVectorWithState<BlockVectorType, PrivarsType, StateVectorType>>
{
    using NativeType = BlockVectorType;
};

template<class SolutionVector>
decltype(auto) native(SolutionVector&& sol)
{
    return Detail::NativeChooser<std::decay_t<SolutionVector>>::native(sol);
}

} // end namespace Dumux

#endif
