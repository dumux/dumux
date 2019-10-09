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
 * \ingroup BoundaryTests
 * \brief The spatial parameters class for the test problem using the 1p cc model.
 */

#ifndef DUMUX_1P_TEST_SPATIALPARAMS_HH
#define DUMUX_1P_TEST_SPATIALPARAMS_HH

#include <dumux/material/spatialparams/fv1p.hh>

namespace Dumux
{

/*!
 * \ingroup BoundaryTests
 * \brief The spatial parameters class for the test problem using the
 *        1p cc model.
 */
template<class GridGeometry, class Scalar>
class OnePSpatialParams
: public FVSpatialParamsOneP<GridGeometry, Scalar,
                             OnePSpatialParams<GridGeometry, Scalar>>
{
    using GridView = typename GridGeometry::GridView;
    using ParentType = FVSpatialParamsOneP<GridGeometry, Scalar,
                                           OnePSpatialParams<GridGeometry, Scalar>>;

    using Element = typename GridView::template Codim<0>::Entity;
    using GlobalPosition = typename Element::Geometry::GlobalCoordinate;

public:
    // export permeability type
    using PermeabilityType = Scalar;

    OnePSpatialParams(std::shared_ptr<const GridGeometry> gridGeometry)
        : ParentType(gridGeometry)
    {
        permeability_ = getParam<Scalar>("Darcy.SpatialParams.Permeability");
        alphaBJ_ = getParam<Scalar>("Darcy.SpatialParams.AlphaBeaversJoseph");
    }

    /*!
     * \brief Function for defining the (intrinsic) permeability \f$[m^2]\f$.
     *
     * \param globalPos The global position
     * \return the intrinsic permeability
     */
    PermeabilityType permeabilityAtPos(const GlobalPosition& globalPos) const
    { return permeability_; }

    /*! \brief Defines the porosity in [-].
     *
     * \param globalPos The global position
     */
    Scalar porosityAtPos(const GlobalPosition& globalPos) const
    { return 0.4; }

    /*! \brief Defines the Beavers-Joseph coefficient in [-].
     *
     * \param globalPos The global position
     */
    Scalar beaversJosephCoeffAtPos(const GlobalPosition& globalPos) const
    { return alphaBJ_; }


private:
    Scalar permeability_;
    Scalar alphaBJ_;
};

} // end namespace Dumux

#endif
