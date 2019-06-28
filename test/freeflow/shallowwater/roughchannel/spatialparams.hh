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
#include <dumux/common/parameters.hh>
#include <dumux/material/fluidmatrixinteractions/frictionlaws/frictionlaw.hh>
#include <dumux/material/fluidmatrixinteractions/frictionlaws/manning.hh>
#include <dumux/material/fluidmatrixinteractions/frictionlaws/nikuradse.hh>

namespace Dumux {

/*!
 * \ingroup ShallowWaterTests
 * \brief The spatial parameters class for the rough channel test.
 *
 */
template<class FVGridGeometry, class Scalar, class VolumeVariables>
class RoughChannelSpatialParams
: public FVSpatialParams<FVGridGeometry, Scalar,
                         RoughChannelSpatialParams<FVGridGeometry, Scalar, VolumeVariables>>
{
    using ThisType = RoughChannelSpatialParams<FVGridGeometry, Scalar, VolumeVariables>;
    using ParentType = FVSpatialParams<FVGridGeometry, Scalar, ThisType>;
    using GridView = typename FVGridGeometry::GridView;
    using FVElementGeometry = typename FVGridGeometry::LocalView;
    using SubControlVolume = typename FVElementGeometry::SubControlVolume;
    using Element = typename GridView::template Codim<0>::Entity;
    using GlobalPosition = typename Element::Geometry::GlobalCoordinate;

public:
    RoughChannelSpatialParams(std::shared_ptr<const FVGridGeometry> fvGridGeometry)
    : ParentType(fvGridGeometry)
    {
        gravity_ = getParam<Scalar>("Problem.Gravity");
        bedSlope_ = getParam<Scalar>("Problem.BedSlope");
        frictionValue_ = getParam<Scalar>("Problem.FrictionValue");
        frictionLawType_ = getParam<std::string>("Problem.FrictionLaw");
        initFrictionLaw();
    }

    /*!
     * \brief Initialize FrictionLaws
     */
    void initFrictionLaw()
    {
      if (frictionLawType_ == "Manning")
      {
          frictionLaw_ = std::make_unique<FrictionLawManning<VolumeVariables>>(gravity_, frictionValue_);
      }
      if (frictionLawType_ == "Nikuradse")
      {
          frictionLaw_ = std::make_unique<FrictionLawNikuradse<VolumeVariables>>(frictionValue_);
      }
      else
      {
          std::cout<<"The FrictionLaw in params.input is unknown. Valid entries are 'Manning' and 'Nikuradse'!"<<std::endl;
      }

    }

    /*! \brief Define the gravitation.
    *
    * \return gravity constant
    */
    Scalar gravity(const GlobalPosition& globalPos) const
    {
        return gravity_;
    }

    /*! \brief Define the gravitation.
    *
    * \return gravity constant
    */
    Scalar gravity() const
    {
        return gravity_;
    }

    /*! \brief Define the friction value.
    *
    * The unit of the friction value depends on the used friction law.
    * The Mannng friction value is just a coefficient and has no unit $[-]$
    * In the Nikuradse law the friction value is the equivalent sand roughness
    * with the unit $[m]$.
    *
    * \return friction value
    */
    Scalar frictionValue(const GlobalPosition& globalPos) const
    {
        return frictionValue_;
    }

    /*! \brief Define the friction value.
    *
    * The unit of the friction value depends on the used friction law.
    * The Mannng friction value is just a coefficient and has no unit $[-]$
    * In the Nikuradse law the friction value is the equivalent sand roughness
    * with the unit $[m]$.
    *
    * \return friction value
    */
    Scalar frictionValue() const
    {
        return frictionValue_;
    }

    /*! \brief Get the frictionLaw.
    *
    * Get the frictionLaw, which already includes the friction value.
    *
    * \return frictionLaw
    */

    const FrictionLaw<VolumeVariables>& frictionLaw(const Element element) const
    {
        return *frictionLaw_;
    }

    /*! \brief Define the bed surface
    *
    * \param element The current element
    * \param scv The sub-control volume inside the element.
    *
    * \return The bed surface
    */
    Scalar bedSurface(const Element& element,
                      const SubControlVolume& scv) const
    {
        // todo depends on index e.g. eIdx = scv.elementIndex();
        return 10.0 - element.geometry().center()[0] * bedSlope_;
    }

private:
    Scalar gravity_;
    Scalar bedSlope_;
    Scalar frictionValue_;
    std::string frictionLawType_;
    std::unique_ptr<FrictionLaw<VolumeVariables>> frictionLaw_;
};

} // end namespace Dumux

#endif
