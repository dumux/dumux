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

// the header guard
#ifndef DUMUX_ROUGH_CHANNEL_SPATIAL_PARAMETERS_HH
#define DUMUX_ROUGH_CHANNEL_SPATIAL_PARAMETERS_HH

// We include the basic spatial parameters for finite volumes file from which we will inherit
#include <dumux/material/spatialparams/fv.hh>
// The parameters header is needed to retrieve run-time parameters.
#include <dumux/common/parameters.hh>
// We include all friction laws, between we can choose for the calculation of the bottom friction source.
#include <dumux/material/fluidmatrixinteractions/frictionlaws/frictionlaw.hh>
#include <dumux/material/fluidmatrixinteractions/frictionlaws/manning.hh>
#include <dumux/material/fluidmatrixinteractions/frictionlaws/nikuradse.hh>

// We enter the namespace Dumux. All Dumux functions and classes are in a namespace Dumux, to make sure they don`t clash with symbols from other libraries you may want to use in conjunction with Dumux.
namespace Dumux {

//In the RoughChannelSpatialParams class we define all functions needed to describe the spatial distributed parameters.
template<class GridGeometry, class Scalar, class VolumeVariables>
class RoughChannelSpatialParams
: public FVSpatialParams<GridGeometry, Scalar,
                         RoughChannelSpatialParams<GridGeometry, Scalar, VolumeVariables>>
{
    // We introduce using declarations that are derived from the property system which we need in this class
    using ThisType = RoughChannelSpatialParams<GridGeometry, Scalar, VolumeVariables>;
    using ParentType = FVSpatialParams<GridGeometry, Scalar, ThisType>;
    using GridView = typename GridGeometry::GridView;
    using FVElementGeometry = typename GridGeometry::LocalView;
    using SubControlVolume = typename FVElementGeometry::SubControlVolume;
    using Element = typename GridView::template Codim<0>::Entity;
    using GlobalPosition = typename Element::Geometry::GlobalCoordinate;

public:
    // In the constructor be read some values from the `params.input` and initialize the friciton law.
    RoughChannelSpatialParams(std::shared_ptr<const GridGeometry> fvGridGeometry)
    : ParentType(fvGridGeometry)
    {
        gravity_ = getParam<Scalar>("Problem.Gravity");
        bedSlope_ = getParam<Scalar>("Problem.BedSlope");
        frictionLawType_ = getParam<std::string>("Problem.FrictionLaw");
        initFrictionLaw();
    }

    // We initialize the friction law based on the law specified in `params.input`.
    void initFrictionLaw()
    {
      if (frictionLawType_ == "Manning")
      {
          Scalar manningN = getParam<Scalar>("Problem.ManningN");
          frictionLaw_ = std::make_unique<FrictionLawManning<VolumeVariables>>(gravity_, manningN);
      }
      else if (frictionLawType_ == "Nikuradse")
      {
          Scalar ks = getParam<Scalar>("Problem.Ks");
          frictionLaw_ = std::make_unique<FrictionLawNikuradse<VolumeVariables>>(ks);
      }
      else
      {
          std::cout<<"The FrictionLaw in params.input is unknown. Valid entries are `Manning` and `Nikuradse`!"<<std::endl;
      }
    }

    // Use this function, if you want to vary the value for the gravity.
    Scalar gravity(const GlobalPosition& globalPos) const
    {
        return gravity_;
    }

    // Use this function for a constant gravity.
    Scalar gravity() const
    {
        return gravity_;
    }

    // This function returns an object of the friction law class, which is initialized with the appropriate friction values. If you want to use different friciton values or laws, you have to use a vector of unique_ptr for `frictionLaw_` and pick the right friction law instances via the `element` argument.
    const FrictionLaw<VolumeVariables>& frictionLaw(const Element& element,
                                                    const SubControlVolume& scv) const
    {
        return *frictionLaw_;
    }

    // Define the bed surface based on the `bedSlope_`.
    Scalar bedSurface(const Element& element,
                      const SubControlVolume& scv) const
    {
        return 10.0 - element.geometry().center()[0] * bedSlope_;
    }

private:
    Scalar gravity_;
    Scalar bedSlope_;
    std::string frictionLawType_;
    std::unique_ptr<FrictionLaw<VolumeVariables>> frictionLaw_;
};
// end of namespace Dumux.
}

#endif
