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
 * \ingroup Common
 * \brief A wrapper that can either contain an object of T or be empty.
 *        This might be used as a workaround for non-default constructible classes.
 * \note  Replace this with std::optional when C++17 is available
 */
#ifndef DUMUX_COMMON_OPTIONAL_HH
#define DUMUX_COMMON_OPTIONAL_HH

#include <utility>

#include <dune/common/typeutilities.hh>

namespace Dumux {

/*!
 * \ingroup Common
 * \brief A wrapper that can either contain an object of T or be empty
 * \tparam T Type of wrapped objects
 */
template<class T>
class Optional
{
public:

    Optional() :
        p_(nullptr)
    {}

    template<class TT, Dune::disableCopyMove<Optional, TT> = 0>
    Optional(TT&& t) :
        p_(nullptr)
    {
        emplace(std::forward<TT>(t));
    }

    Optional(Optional&& other)
    {
        if (other)
            p_ = new (buffer_) T(std::move(other.value()));
        else
            p_ = nullptr;
    }

    Optional(const Optional& other)
    {
        if (other)
            p_ = new (buffer_) T(other.value());
        else
            p_ = nullptr;
    }

    ~Optional()
    {
        if (operator bool())
            p_->~T();
    }

    template<class TT, Dune::disableCopyMove<Optional, TT> = 0 >
    Optional& operator=(TT&& t)
    {
        if (operator bool())
            *p_ = std::forward<T>(t);
        else
            p_ = new (buffer_) T(std::forward<T>(t));
        return *this;
    }

    Optional& operator=(const Optional& other)
    {
        if (other)
            *this = other.value();
        else if (operator bool())
        {
            p_->~T();
            p_ = nullptr;
        }
        return *this;
    }

    Optional& operator=(Optional&& other)
    {
        if (other)
            *this = std::move(other.value());
        else if (operator bool())
        {
            p_->~T();
            p_ = nullptr;
        }
        return *this;
    }

    explicit operator bool() const
    {
        return p_;
    }

    const T& value() const
    {
        return *p_;
    }

    T& value()
    {
        return *p_;
    }

    template< class... Args >
    void emplace(Args&&... args)
    {
        if (operator bool())
            p_->~T();
        p_ = new (buffer_) T(std::forward<Args>(args)...);
    }

    void release()
    {
        if (operator bool())
        {
            p_->~T();
            p_ = nullptr;
        }
    }

private:

    alignas(T) char buffer_[sizeof(T)];
    T* p_;
};


} // namespace Dumux

#endif // DUMUX_COMMON_OPTIONAL_HH
