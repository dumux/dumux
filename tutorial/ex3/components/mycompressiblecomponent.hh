// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*****************************************************************************
 *   See the file COPYING for full copying permissions.                      *
 *                                                                           *
 *   This program is free software: you can redistribute it and/or modify    *
 *   it under the terms of the GNU General Public License as published by    *
 *   the Free Software Foundation, either version 2 of the License, or       *
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
 * \ingroup Components
 * \brief A ficitious component to be implemented in exercise 3.
 */
#ifndef DUMUX_MYCOMPRESSIBLECOMPONENT_HH
#define DUMUX_MYCOMPRESSIBLECOMPONENT_HH

#include <dumux/material/idealgas.hh>
#include <dumux/material/components/component.hh>


namespace Dumux
{
/*!
 * \ingroup Components
 * \brief A ficitious component to be implemented in exercise 3.
 *
 * \tparam Scalar The type used for scalar values
 */
template <class Scalar>
class MyCompressibleComponent : public Component<Scalar, MyCompressibleComponent<Scalar> >
{

public:
    /*!
     * \brief A human readable name for MyCompressibleComponent.
     */
    static std::string name()
    { return "MyCompressibleComponent"; }

    /*!
     * TODO: Copy the methods implemented in MyIncompressibleComponent and substitute
     *       the density calculation by the expression given in the exercise description.
     */
};

} // end namespace

#endif
