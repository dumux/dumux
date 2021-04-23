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
#ifndef DUMUX_COMMON_BLOCKVECTOR_HH
#define DUMUX_COMMON_BLOCKVECTOR_HH

// // forward declaration
// namespace Dune {
//
// template<class... Args>
// class MultiTypeBlockVector;
//
// } // end namespace Dune

namespace Dumux::Istl {

template<class BVType>
class BlockVector
{
public:
    using field_type = typename BVType::field_type;
    using block_type = typename BVType::block_type;
    using allocator_type = typename BVType::allocator_type;
    using size_type = typename BVType::size_type;
    using Iterator = typename BVType::Iterator;
    using ConstIterator = typename BVType::ConstIterator;

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

    BVType& native()
    { return blockVector_; }

    const BVType& native() const
    { return blockVector_; }

private:
    BVType blockVector_;
};

// forward declare
template<class PV, class State>
struct BlockVectorView;

template<class PV, class State>
struct ConstBlockVectorView;

template<class PV, class State>
struct ConstBlockVectorView
{
    template<typename T, class S>
    friend struct BlockVectorView;

    // friend BlockVectorView<PV>& operator= (const ConstBlockVectorView<PV>& other);

    ConstBlockVectorView(const PV& sol, const State& state) : solution_(&sol), state_(&state)
    {}

   operator PV() const
   { return *solution_; }

   State state() const
   { return *state_; }

   const auto& operator[] (int i) const
   {
       return (*solution_)[i];
   }

private:
   const PV* solution_;
   const State* state_;
};


template<class PV, class State>
struct BlockVectorView
{
    template<typename T, class S>
    friend struct ConstBlockVectorView;

    BlockVectorView() : solution_(nullptr), state_(nullptr) {}

    BlockVectorView(PV& sol, State& state)
    {
        solution_ = &sol;
        state_ = &state;
        // PV::narf();
    }

   auto& operator= (const ConstBlockVectorView<PV, State>& other)
   {
       if (solution_)
       {
           *solution_ = (*other.solution_);
           *state_ = other.state();
       }
       else
       {
           solStored_ = std::make_unique<PV>((*other.solution_));
           stateStored_ = other.state();
       }

       return *this;
   }

   auto& operator= (const SwitchablePrimaryVariables<PV, int>& other)
   {
       if (solution_)
       {
           *solution_ = other;
           *state_ = other.state();
       }
       else
       {
           solStored_ = std::make_unique<PV>(other);
           stateStored_ = other.state();
       }

       return *this;
   }

   void setState(int s)
   { *state_ = s; }

   operator SwitchablePrimaryVariables<PV, int>() const
   {
       SwitchablePrimaryVariables<PV, int> privars;
       if (solution_)
       {
           privars = *solution_;
           privars.setState(*state_);
       }
       else
       {
           privars = *solStored_;
           privars.setState(stateStored_);
       }

       return privars;
   }
   // operator PV() const
   // {
   //     if (solution_)
   //         return *solution_;
   //     else
   //         return *solStored_;
   // }

   auto& operator[] (int i)
   {
       if (solution_)
           return (*solution_)[i];
       else
           return (*solStored_)[i];
   }

   const auto& operator[] (int i) const
   {
       if (solution_)
           return (*solution_)[i];
       else
           return (*solStored_)[i];
   }

   static constexpr auto size()
   { return PV::size(); }



   State state() const
   {
       if (solution_)
           return *state_;
       else
           return stateStored_;
   }

private:
   PV* solution_;
   State* state_;
   State stateStored_;
   std::unique_ptr<PV> solStored_;
};

template<class BVType, class State>
class BlockVectorWithState
{

public:
    using field_type = typename BVType::field_type;
    using block_type = typename BVType::block_type;
    using allocator_type = typename BVType::allocator_type;
    using size_type = typename BVType::size_type;
    using Iterator = typename BVType::Iterator;
    using ConstIterator = typename BVType::ConstIterator;

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

    auto operator [](int i)
    {
        return BlockVectorView<block_type, State>(blockVector_[i], states_[i]);
    }

    const auto operator[] (int i) const
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

    BVType& native()
    { return blockVector_; }

    const BVType& native() const
    { return blockVector_; }

private:
    BVType blockVector_;
    std::vector<State> states_;
};

} // end namespace Dumux::Istl

#endif
