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
 * \brief A ficitious component to be implemented in tutorial exercise 3.
 */
#ifndef DUMUX_MYINCOMPRESSIBLECOMPONENT_HH
#define DUMUX_MYINCOMPRESSIBLECOMPONENT_HH

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
class MyIncompressibleComponent : public Component<Scalar, MyIncompressibleComponent<Scalar> >
{
public:
    /*!
     * \brief A human readable name for MyIncompressibleComponent.
     */
    static std::string name()
    { return "MyIncompressibleComponent"; }

    /*!
     * TODO: Implement the methods for the component data given in the exercise description.
     */
};

} // end namespace

#endif
