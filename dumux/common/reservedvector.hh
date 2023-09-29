// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup Common
 * \copybrief Dumux::ReservedVector
 */
#ifndef DUMUX_RESERVED_VECTOR_HH
#define DUMUX_RESERVED_VECTOR_HH

#include <bit>
#include <array>
#include <vector>
#include <utility>
#include <algorithm>
#include <iterator>
#include <memory_resource>
#include <initializer_list>

namespace Dumux {

/*!
 * \ingroup Common
 * \brief Reserved vector implementation around std::pmr::vector.
 *        Uses preallocated memory but allows for dynamic insertion
 *        of elements beyond the preallocated memory (in contrast to
 *        Dune::ReservedVector).
 */
template<typename T, std::size_t N>
class ReservedVector {
    using Vector = std::pmr::vector<T>;

 public:
    ReservedVector() = default;

    ReservedVector(const ReservedVector& other)
    : ReservedVector()
    {
        elements_.clear();
        elements_.reserve(other.size());
        std::copy(other.elements_.begin(), other.elements_.end(), std::back_inserter(elements_));
    }

    ReservedVector(std::size_t n, const T& r)
    : ReservedVector()
    {
        elements_.resize(n, r);
    }

    ReservedVector(std::initializer_list<T>&& initList)
    : ReservedVector()
    {
        elements_.resize(initList.size());
        std::copy(initList.begin(), initList.end(), elements_.begin());
    }

    ReservedVector(ReservedVector&& other)
    : ReservedVector()
    {
        elements_.clear();
        elements_.reserve(other.size());
        std::move(other.elements_.begin(), other.elements_.end(), std::back_inserter(elements_));
    }

    ReservedVector& operator=(const ReservedVector& other) {
        elements_ = Vector{typename Vector::allocator_type{&resource_}};
        elements_.reserve(other.size());
        std::copy(other.elements_.begin(), other.elements_.end(), std::back_inserter(elements_));
        return *this;
    };

    void clear() { elements_.clear(); }
    std::size_t size() const { return elements_.size(); }
    void reserve(std::size_t n) { elements_.reserve(n); }
    void resize(std::size_t n) { elements_.resize(n); }
    void resize(std::size_t n, const T& value) { elements_.resize(n, value); }

    template<typename... Args>
    void push_back(Args&&... args) { elements_.push_back(std::forward<Args>(args)...); }
    void push_back(const T& element) { elements_.push_back(element); }
    void push_back(T&& element) { elements_.push_back(std::move(element)); }

    decltype(auto) back() { return elements_.back(); }
    decltype(auto) back() const { return elements_.back(); }

    decltype(auto) begin() { return elements_.begin(); }
    decltype(auto) begin() const { return elements_.begin(); }

    decltype(auto) end() { return elements_.end(); }
    decltype(auto) end() const { return elements_.end(); }

    decltype(auto) operator[](std::size_t i) { return elements_[i]; }
    decltype(auto) operator[](std::size_t i) const { return elements_[i]; }

    decltype(auto) at(std::size_t i) { return elements_.at(i); }
    decltype(auto) at(std::size_t i) const { return elements_.at(i); }

 private:
    std::array<std::byte, N*sizeof(T)> buffer_;
    std::pmr::monotonic_buffer_resource resource_{buffer_.data(), buffer_.size()};
    Vector elements_{typename Vector::allocator_type{&resource_}};
};

}  // namespace Dumux

#endif
