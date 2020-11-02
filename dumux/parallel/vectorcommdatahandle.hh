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
 * \ingroup Parallel
 * \brief Contains a class to exchange entries of a vector
 */
#ifndef DUMUX_VECTOR_COMM_DATA_HANDLE_HH
#define DUMUX_VECTOR_COMM_DATA_HANDLE_HH

#include <algorithm>

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

template<class Mapper, class Vector, int codim, class DataType = typename Vector::value_type>
using VectorCommDataHandleEqual = VectorCommDataHandle<Mapper, Vector, codim, Detail::SetEqual, DataType>;

template<class Mapper, class Vector, int codim, class DataType = typename Vector::value_type>
using VectorCommDataHandleSum = VectorCommDataHandle<Mapper, Vector, codim, Detail::Sum, DataType>;

template<class Mapper, class Vector, int codim, class DataType = typename Vector::value_type>
using VectorCommDataHandleMin = VectorCommDataHandle<Mapper, Vector, codim, Detail::Min, DataType>;

template<class Mapper, class Vector, int codim, class DataType = typename Vector::value_type>
using VectorCommDataHandleMax = VectorCommDataHandle<Mapper, Vector, codim, Detail::Max, DataType>;

} // end namespace Dumux

#endif
