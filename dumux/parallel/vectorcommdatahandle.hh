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
 * \ingroup Linear
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
 * \ingroup Linear
 * \brief A data handle class to exchange entries of a vector
 */
template<class Mapper, class Vector, int entityCodim, class ScatterOperator>
class VectorCommDataHandle
  : public Dune::CommDataHandleIF<VectorCommDataHandle<Mapper,Vector, entityCodim, ScatterOperator>,
                                  typename Vector::value_type>
{
public:
  //! export type of data for message buffer
  using DataType = typename Vector::value_type;

  VectorCommDataHandle(const Mapper& mapper, Vector& vector)
  : mapper_(mapper), vector_(vector)
  {}

  //! returns true if data for this codim should be communicated
  bool contains(int dim, int codim) const
  { return (codim == entityCodim); }

  //! returns true if size per entity of given dim and codim is a constant
  bool fixedsize(int dim, int codim) const
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

template<class Mapper, class Vector, int codim>
using VectorCommDataHandleEqual = VectorCommDataHandle<Mapper,Vector, codim, Detail::SetEqual>;

template<class Mapper, class Vector, int codim>
using VectorCommDataHandleSum = VectorCommDataHandle<Mapper,Vector, codim, Detail::Sum>;

template<class Mapper, class Vector, int codim>
using VectorCommDataHandleMin = VectorCommDataHandle<Mapper,Vector, codim, Detail::Min>;

template<class Mapper, class Vector, int codim>
using VectorCommDataHandleMax = VectorCommDataHandle<Mapper,Vector, codim, Detail::Max>;

} // end namespace Dumux

#endif // DUMUX_VECTOR_COMM_DATA_HANDLE_HH
