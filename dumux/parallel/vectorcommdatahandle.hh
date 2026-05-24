// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup Parallel
 * \brief Contains a class to exchange entries of a vector
 */
#ifndef DUMUX_VECTOR_COMM_DATA_HANDLE_HH
#define DUMUX_VECTOR_COMM_DATA_HANDLE_HH

#include <algorithm>
#include <bitset>

#include <dune/grid/common/datahandleif.hh>

namespace Dumux {

namespace Detail {

    struct SetEqual
    {
        template<class A, class B>
        static void apply(A& a, const B& b)
        { a = b; }
    };

    struct Sum
    {
        template<class A, class B>
        static void apply(A& a, const B& b)
        { a += b; }
    };

    struct Max
    {
        template<class A, class B>
        static void apply(A& a, const B& b)
        {
            using std::max;
            a = max(a,b);
        }
    };

    struct Min
    {
        template<class A, class B>
        static void apply(A& a, const B& b)
        {
            using std::min;
            a = min(a,b);
        }
    };
} // end namespace Detail

/*!
 * \ingroup Parallel
 * \brief A data handle class to exchange entries of a vector
 */
template<class Mapper, class Vector, int entityCodim,
         class ScatterOperator, class DataT = typename Vector::value_type>
class VectorCommDataHandle
  : public Dune::CommDataHandleIF<VectorCommDataHandle<Mapper, Vector, entityCodim, ScatterOperator, DataT>, DataT>
{
public:
  //! export type of data for message buffer
  using DataType = DataT;

  VectorCommDataHandle(const Mapper& mapper, Vector& vector)
  : mapper_(mapper), vector_(vector)
  {}

  //! returns true if data for this codim should be communicated
  bool contains(int dim, int codim) const
  { return (codim == entityCodim); }

  //! returns true if size per entity of given dim and codim is a constant
  bool fixedSize(int dim, int codim) const
  { return true; }

  /*!
   * \brief how many objects of type DataType have to be sent for a given entity
   * \note Only the sender side needs to know this size.
   */
  template<class Entity>
  std::size_t size(Entity& entity) const
  { return 1; }

  //! pack data from user to message buffer
  template<class MessageBuffer, class Entity>
  void gather(MessageBuffer& buff, const Entity& entity) const
  { buff.write(vector_[mapper_.index(entity)]); }

  /*!
   * \brief unpack data from message buffer to user
   * \note n is the number of objects sent by the sender
   */
  template<class MessageBuffer, class Entity>
  void scatter(MessageBuffer& buff, const Entity& entity, std::size_t n)
  {
      DataType x;
      buff.read(x);
      ScatterOperator::apply(vector_[mapper_.index(entity)], x);
  }

protected:
  const Mapper& mapper_;
  Vector& vector_;
};

/*!
 * \ingroup Parallel
 * \brief A data handle class to exchange entries of a vector for multiple codims in one communication call
 */
template<class Mapper, class Vector, std::size_t numCodims,
         class ScatterOperator, class DataT = typename Vector::value_type>
class MultiCodimVectorCommDataHandle
  : public Dune::CommDataHandleIF<MultiCodimVectorCommDataHandle<Mapper, Vector, numCodims, ScatterOperator, DataT>, DataT>
{
public:
  //! export type of data for message buffer
  using DataType = DataT;

  MultiCodimVectorCommDataHandle(const Mapper& mapper,
                                 Vector& vector,
                                 std::bitset<numCodims> activeCodims)
  : mapper_(mapper), vector_(vector), activeCodims_(std::move(activeCodims))
  {}

  //! returns true if data for this codim should be communicated
  bool contains(int dim, int codim) const
  {
      return codim >= 0
             && codim < static_cast<int>(numCodims)
             && activeCodims_.test(codim);
  }

  //! returns true if size per entity of given dim and codim is a constant
  bool fixedSize(int dim, int codim) const
  { return true; }

  /*!
   * \brief how many objects of type DataType have to be sent for a given entity
   * \note Only the sender side needs to know this size.
   */
  template<class Entity>
  std::size_t size(Entity& entity) const
  { return 1; }

  //! pack data from user to message buffer
  template<class MessageBuffer, class Entity>
  void gather(MessageBuffer& buff, const Entity& entity) const
  { buff.write(vector_[mapper_.index(entity)]); }

  /*!
   * \brief unpack data from message buffer to user
   * \note n is the number of objects sent by the sender
   */
  template<class MessageBuffer, class Entity>
  void scatter(MessageBuffer& buff, const Entity& entity, std::size_t n)
  {
      DataType x;
      buff.read(x);
      ScatterOperator::apply(vector_[mapper_.index(entity)], x);
  }

protected:
  const Mapper& mapper_;
  Vector& vector_;
  std::bitset<numCodims> activeCodims_;
};

template<class Mapper, class Vector, int codim, class DataType = typename Vector::value_type>
using VectorCommDataHandleEqual = VectorCommDataHandle<Mapper, Vector, codim, Detail::SetEqual, DataType>;

template<class Mapper, class Vector, int codim, class DataType = typename Vector::value_type>
using VectorCommDataHandleSum = VectorCommDataHandle<Mapper, Vector, codim, Detail::Sum, DataType>;

template<class Mapper, class Vector, int codim, class DataType = typename Vector::value_type>
using VectorCommDataHandleMin = VectorCommDataHandle<Mapper, Vector, codim, Detail::Min, DataType>;

template<class Mapper, class Vector, int codim, class DataType = typename Vector::value_type>
using VectorCommDataHandleMax = VectorCommDataHandle<Mapper, Vector, codim, Detail::Max, DataType>;

template<class Mapper, class Vector, std::size_t numCodims, class DataType = typename Vector::value_type>
using MultiCodimVectorCommDataHandleEqual = MultiCodimVectorCommDataHandle<Mapper, Vector, numCodims, Detail::SetEqual, DataType>;

template<class Mapper, class Vector, std::size_t numCodims, class DataType = typename Vector::value_type>
using MultiCodimVectorCommDataHandleSum = MultiCodimVectorCommDataHandle<Mapper, Vector, numCodims, Detail::Sum, DataType>;

template<class Mapper, class Vector, std::size_t numCodims, class DataType = typename Vector::value_type>
using MultiCodimVectorCommDataHandleMin = MultiCodimVectorCommDataHandle<Mapper, Vector, numCodims, Detail::Min, DataType>;

template<class Mapper, class Vector, std::size_t numCodims, class DataType = typename Vector::value_type>
using MultiCodimVectorCommDataHandleMax = MultiCodimVectorCommDataHandle<Mapper, Vector, numCodims, Detail::Max, DataType>;

} // end namespace Dumux

#endif
