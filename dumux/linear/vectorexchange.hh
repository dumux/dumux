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
#ifndef DUMUX_VECTOR_EXCHANGE_HH
#define DUMUX_VECTOR_EXCHANGE_HH

#include <dune/grid/common/datahandleif.hh>

namespace Dumux {

/*!
 * \ingroup Linear
 * \todo why is this needed? is parallel/vectorexhange.hh not the same?
 * \brief A data handle class to exchange entries of a vector
 */
template<class Mapper, class Vector> // mapper type and vector type
class VectorExchange
  : public Dune::CommDataHandleIF<VectorExchange<Mapper,Vector>,
                                  typename Vector::value_type>
{
public:
  //! export type of data for message buffer
  using DataType = typename Vector::value_type;

  //! returns true if data for this codim should be communicated
  bool contains (int dim, int codim) const
  {
        return (codim == 0);
  }

  //! returns true if size per entity of given dim and codim is a constant
  bool fixedsize (int dim, int codim) const
  {
        return true;
  }

  /*!
   * \brief how many objects of type DataType have to be sent for a given entity
   * \note Only the sender side needs to know this size.
   */
  template<class Entity>
  size_t size (Entity& entity) const
  {
        return 1;
  }

  //! pack data from user to message buffer
  template<class MessageBuffer, class Entity>
  void gather (MessageBuffer& buff, const Entity& entity) const
  {
      buff.write(dataVector_[mapper_.index(entity)]);
  }

  /*!
   * \brief unpack data from message buffer to user
   * \note n is the number of objects sent by the sender
   */
  template<class MessageBuffer, class Entity>
  void scatter (MessageBuffer& buff, const Entity& entity, size_t n)
  {
      DataType x;
      buff.read(x);
      dataVector_[mapper_.index(entity)] = x;
  }

  //! constructor
  VectorExchange (const Mapper& mapper, Vector& dataVector)
        : mapper_(mapper), dataVector_(dataVector)
  {}

private:
  const Mapper& mapper_;
  Vector& dataVector_;
};

} // end namespace Dumux

#endif
