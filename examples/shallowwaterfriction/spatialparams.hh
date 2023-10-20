// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//

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
// [[details]] includes
// We include the basic spatial parameters file for finite volumes, from which we will inherit.
#include <dumux/freeflow/spatialparams.hh>
// We include all friction laws.
#include <dumux/material/fluidmatrixinteractions/frictionlaws/frictionlaw.hh>
#include <dumux/material/fluidmatrixinteractions/frictionlaws/manning.hh>
#include <dumux/material/fluidmatrixinteractions/frictionlaws/nikuradse.hh>
#include <dumux/material/fluidmatrixinteractions/frictionlaws/nofriction.hh>
// [[/details]]

//
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
: public FreeFlowSpatialParams<GridGeometry, Scalar,
                               RoughChannelSpatialParams<GridGeometry, Scalar, VolumeVariables>>
{
// [[/codeblock]]
    // [[details]] convenience aliases
    using ThisType = RoughChannelSpatialParams<GridGeometry, Scalar, VolumeVariables>;
    using ParentType = FreeFlowSpatialParams<GridGeometry, Scalar, ThisType>;
    using GridView = typename GridGeometry::GridView;
    using FVElementGeometry = typename GridGeometry::LocalView;
    using SubControlVolume = typename FVElementGeometry::SubControlVolume;
    using Element = typename GridView::template Codim<0>::Entity;
    using GlobalPosition = typename Element::Geometry::GlobalCoordinate;
    // [[/details]]

    // In the constructor, the properties of the the rough channel are set. Namely, these are
    // the friction law, including it's friction parameter, the acceleration
    // due to gravity and the altitude of the channel bed surface.
    // [[codeblock]]
public:
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
    // [[/codeblock]]

    // The following functions expose the parameters required by the model.

    // [[codeblock]]
    // This function returns an object of the friction law class, already initialized with a friction value.
    const FrictionLaw<VolumeVariables>& frictionLaw(const Element& element,
                                                    const SubControlVolume& scv) const
    { return *frictionLaw_; }

    // This function returns the acceleration due to gravity.
    Scalar gravity(const GlobalPosition& globalPos) const
    { return gravity_; }

    // Define the bed surface based on the bed slope and the bed level at the inflow (10 m).
    Scalar bedSurface(const Element& element,
                      const SubControlVolume& scv) const
    { return 10.0 - element.geometry().center()[0] * bedSlope_; }
    // [[/codeblock]]

    // [[details]] private variables
private:
    Scalar gravity_;
    Scalar bedSlope_;
    std::string frictionLawType_;
    std::unique_ptr<FrictionLaw<VolumeVariables>> frictionLaw_;
}; // end class definition of RoughChannelSpatialParams
} // end of namespace Dumux.
// [[/details]]
// [[/content]]
#endif
