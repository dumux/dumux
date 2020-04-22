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

#ifndef DUMUX_ROUGH_CHANNEL_SPATIAL_PARAMETERS_HH
#define DUMUX_ROUGH_CHANNEL_SPATIAL_PARAMETERS_HH

// ## Parameter distributions (`spatialparams.hh`)
//
// This file contains the __spatial parameters class__ which defines the
// the friction law, including it's friction parameter, the acceleration
// due to gravity and the altitude of the channel bed surface. In this example only the bed
// surface has a non constant distribution.
//
// [[content]]
//
// ### Include files
// We include the basic spatial parameters file for finite volumes, from which we will inherit.
#include <dumux/material/spatialparams/fv.hh>
// We include all friction laws.
#include <dumux/material/fluidmatrixinteractions/frictionlaws/frictionlaw.hh>
#include <dumux/material/fluidmatrixinteractions/frictionlaws/manning.hh>
#include <dumux/material/fluidmatrixinteractions/frictionlaws/nikuradse.hh>
#include <dumux/material/fluidmatrixinteractions/frictionlaws/nofriction.hh>

// ### The spatial parameters class
//
// In the `RoughChannelSpatialParams` class, we define all functions needed to describe
// the rough channel for the shallow water problem.
// We inherit from the `FVSpatialParams` class, which is the base class
// for spatial parameters in the context of
// applications using finite volume discretization schemes.
// [[codeblock]]
namespace Dumux {

template<class GridGeometry, class Scalar, class VolumeVariables>
class RoughChannelSpatialParams
: public FVSpatialParams<GridGeometry, Scalar,
                         RoughChannelSpatialParams<GridGeometry, Scalar, VolumeVariables>>
{
    // This convenience aliases will be used throughout this class
    using ThisType = RoughChannelSpatialParams<GridGeometry, Scalar, VolumeVariables>;
    using ParentType = FVSpatialParams<GridGeometry, Scalar, ThisType>;
    using GridView = typename GridGeometry::GridView;
    using FVElementGeometry = typename GridGeometry::LocalView;
    using SubControlVolume = typename FVElementGeometry::SubControlVolume;
    using Element = typename GridView::template Codim<0>::Entity;
    using GlobalPosition = typename Element::Geometry::GlobalCoordinate;
    // [[/codeblock]]

    // In the following, the properties of the the rough channel are set. Namely, these are
    // the friction law, including it's friction parameter, the acceleration
    // due to gravity and the altitude of the channel bed surface.
    // [[codeblock]]
public:
    // In the constructor we read some values from the `params.input` and initialize the friciton law.
    RoughChannelSpatialParams(std::shared_ptr<const GridGeometry> gridGeometry)
    : ParentType(gridGeometry)
    {
        gravity_ = getParam<Scalar>("Problem.Gravity");
        bedSlope_ = getParam<Scalar>("Problem.BedSlope");
        frictionLawType_ = getParam<std::string>("Problem.FrictionLaw");
        initFrictionLaw();
    }

    // This function handles the initialization of the friction law based on the settings
    // specified in `params.input`.
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
      else if (frictionLawType_ == "None")
      {
          frictionLaw_ = std::make_unique<FrictionLawNoFriction<VolumeVariables>>();
      }
      else
      {
          std::cout<<"The FrictionLaw in params.input is unknown. Valid entries are `Manning`,"
                     " `Nikuradse` and `None`!"<<std::endl;
      }
    }

    // This function returns an object of the friction law class, already initialized with a friction value.
    const FrictionLaw<VolumeVariables>& frictionLaw(const Element& element,
                                                    const SubControlVolume& scv) const
    {
        return *frictionLaw_;
    }

    // This function returns the acceleration due to gravity.
    Scalar gravity(const GlobalPosition& globalPos) const
    {
        return gravity_;
    }

    // Define the bed surface based on the bed slope and the bed level at the inflow (10 m).
    Scalar bedSurface(const Element& element,
                      const SubControlVolume& scv) const
    {
        return 10.0 - element.geometry().center()[0] * bedSlope_;
    }

// We declare the private variables of the problem.
private:
    Scalar gravity_;
    Scalar bedSlope_;
    std::string frictionLawType_;
    std::unique_ptr<FrictionLaw<VolumeVariables>> frictionLaw_;
}; // end class definition of RoughChannelSpatialParams
} // end of namespace Dumux.
// [[/codeblock]]
// [[/content]]
#endif
