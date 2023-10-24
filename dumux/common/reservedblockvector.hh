// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup Core
 * \brief A arithmetic block vector type based on DUNE's reserved vector
 */
#ifndef DUMUX_RESERVED_BLOCK_VECTOR_HH
#define DUMUX_RESERVED_BLOCK_VECTOR_HH

#include <algorithm>
#include <dune/common/reservedvector.hh>

namespace Dumux {

/*!
 * \ingroup Core
 * \brief A arithmetic block vector type based on DUNE's reserved vector
 */
template<class BlockType, int capacity>
class ReservedBlockVector : public Dune::ReservedVector<BlockType, capacity>
{
    using Base = Dune::ReservedVector<BlockType, capacity>;
public:

    using size_type = typename Base::size_type;
    using value_type = BlockType;

    using Base::Base;

    explicit ReservedBlockVector() : Base() {}
    explicit ReservedBlockVector(size_type size) : Base() { this->resize(size); }

    ReservedBlockVector(const ReservedBlockVector&) = default;
    ReservedBlockVector(ReservedBlockVector&&) = default;

    ReservedBlockVector& operator= (const ReservedBlockVector&) = default;
    ReservedBlockVector& operator= (ReservedBlockVector&&) = default;

    ~ReservedBlockVector() = default;

    //! assignment from scalar
    ReservedBlockVector& operator= (const typename BlockType::field_type& v)
    {
       std::fill(this->begin(), this->end(), v);
       return *this;
    }

    //! vector space addition
    ReservedBlockVector& operator+= (const ReservedBlockVector& other)
    {
      for (size_type i = 0; i < this->size(); ++i)
          (*this)[i] += other[i];
      return *this;
    }

    //! vector space subtraction
    ReservedBlockVector& operator-= (const ReservedBlockVector& other)
    {
      for (size_type i = 0; i < this->size(); ++i)
          (*this)[i] -= other[i];
      return *this;
    }

    //! division by scalar
    ReservedBlockVector& operator/= (const typename BlockType::field_type& v)
    {
      for (size_type i = 0; i < this->size(); ++i)
          (*this)[i] /= v;
      return *this;
    }

    //! multiplication by scalar
    ReservedBlockVector& operator*= (const typename BlockType::field_type& v)
    {
      for (size_type i = 0; i < this->size(); ++i)
          (*this)[i] *= v;
      return *this;
    }
};

} // end namespace Dumux

#endif
