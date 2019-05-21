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
 * \ingroup ShallowWaterTests
 * \brief The spatial parameters for the rough channel problem.
 */
#ifndef DUMUX_ROUGH_CHANNEL_SPATIAL_PARAMETERS_HH
#define DUMUX_ROUGH_CHANNEL_SPATIAL_PARAMETERS_HH

#include <dumux/material/spatialparams/fv.hh>

namespace Dumux {

/*!
 * \ingroup ShallowWaterTests
 * \brief The spatial parameters class for the rough channel test.
 *
 */
template<class FVGridGeometry, class Scalar>
class RoughChannelSpatialParams
: public FVSpatialParams<FVGridGeometry, Scalar,
                         RoughChannelSpatialParams<FVGridGeometry, Scalar>>
{
    using ThisType = RoughChannelSpatialParams<FVGridGeometry, Scalar>;
    using ParentType = FVSpatialParams<FVGridGeometry, Scalar, ThisType>;
    using GridView = typename FVGridGeometry::GridView;
    using FVElementGeometry = typename FVGridGeometry::LocalView;
    using SubControlVolume = typename FVElementGeometry::SubControlVolume;
    using Element = typename GridView::template Codim<0>::Entity;
    using GlobalPosition = typename Element::Geometry::GlobalCoordinate;

public:
    RoughChannelSpatialParams(std::shared_ptr<const FVGridGeometry> fvGridGeometry)
    : ParentType(fvGridGeometry)
    {}

    /*! \brief Define the porosity in [-].
   *
   * \param globalPos The global position where we evaluate
   */
    Scalar porosityAtPos(const GlobalPosition& globalPos) const
    { return 1.0; }


    /*! \brief Define the gravitation.
    *
    * \return gravity constant
    */
    Scalar gravity(const GlobalPosition& globalPos) const
    {
        return gravity_;
    }

    /*! \brief Define the friction value.
    *
    * \return friction value
    */
    Scalar frictionValue(const GlobalPosition& globalPos) const
    {
        return frictionValue_;
    }

    /*! \brief Define the friction value.
    *
    * \return friction value
    */
    Scalar frictionValue() const
    {
        return frictionValue_;
    }

    /*! \brief Define the bed surface
    *
    * \param element The current element
    * \param scv The sub-control volume inside the element.
    * \return the bed surface
    */
    Scalar bedSurface(const Element& element,
                      const SubControlVolume& scv) const
    {
        // todo depends on index e.g. eIdx = scv.elementIndex();
        return 10.0 - element.geometry().center()[0] * bedSlope_;
    }

private:
    static constexpr Scalar gravity_ = 9.81;
    static constexpr Scalar bedSlope_ = 0.001;
    static constexpr Scalar frictionValue_ = 0.025;
};

} // end namespace Dumux

#endif
