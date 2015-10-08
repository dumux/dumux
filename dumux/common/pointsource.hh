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
 * \brief A point source class,
 *        i.e. sources located at a single point in space
 */

#ifndef DUMUX_POINTSOURCE_HH
#define DUMUX_POINTSOURCE_HH

namespace Dumux
{

namespace Properties
{

// Property forward declarations
NEW_PROP_TAG(Scalar);
NEW_PROP_TAG(GridView);
NEW_PROP_TAG(PrimaryVariables);

} // end namespace Properties

/*!
 * \ingroup Common
 * \brief A point source class
 */
template<class TypeTag>
class PointSource
{
    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
    typedef typename GET_PROP_TYPE(TypeTag, GridView) GridView;
    typedef typename GET_PROP_TYPE(TypeTag, PrimaryVariables) PrimaryVariables;

    static const int dimworld = GridView::dimensionworld;
    typedef typename Dune::FieldVector<Scalar, dimworld> GlobalPosition;

public:
    // Contructor
    PointSource(GlobalPosition pos, PrimaryVariables values)
      : pos_(pos), values_(values) {}

    //! return the source values
    const PrimaryVariables& values() const
    { return values_; }

    //! return the source position
    const GlobalPosition& position() const
    { return pos_; }

    //! Divide values by scalar
    template<class T>
    void divideValues(T&& n)
    { values_ /= std::forward<T>(n); }

    //! Add something to the values
    template<class T>
    void addToValues(T&& n)
    { values_ += std::forward<T>(n); }

private:
    GlobalPosition pos_;
    PrimaryVariables values_;
};

} // end namespace Dumux

#endif
