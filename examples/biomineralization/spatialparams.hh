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

#ifndef DUMUX_MICP_COLUMN_SIMPLE_CHEM_SPATIAL_PARAMS_HH
#define DUMUX_MICP_COLUMN_SIMPLE_CHEM_SPATIAL_PARAMS_HH

// ## Parameter distributions (`spatialparams_1p.hh`)
//
// [[content]]
//
// This file contains the __spatial parameters class__ which defines the
// distributions for the porous medium parameters permeability and porosity
// over the computational grid
//
// ### Include files
// [[details]] includes
// We include the basic spatial parameters for finite volumes file from which we will inherit
#include <dumux/material/spatialparams/fv.hh>
// We include the files for the two-phase laws: the Brooks-Corey pc-Sw and relative permeability laws
#include <dumux/material/fluidmatrixinteractions/2p/brookscorey.hh>
// We include the laws for changing porosity due to precipitation
#include <dumux/material/fluidmatrixinteractions/porosityprecipitation.hh>
// We include the laws for changing permeability based on porosity change according to Kozeny-Carman
#include <dumux/material/fluidmatrixinteractions/permeabilitykozenycarman.hh>
// [[/details]]
//
// ### The spatial parameters class
// In the `ICPSpatialParams` class, we define all functions needed to describe
// the porous medium, e.g. porosity and permeability.
// We inherit from the `FVSpatialParams` class which is the base class for spatial paramters using finite volume discretization schemes.

// [[codeblock]]
namespace Dumux {

// In the ICPSpatialParams class we define all functions needed to describe the spatial distributed parameters.
template<class GridGeometry, class Scalar>
class ICPSpatialParams
: public FVSpatialParams<GridGeometry, Scalar, ICPSpatialParams<GridGeometry, Scalar>>
{
    // We introduce using declarations that are derived from the property system which we need in this class
    using ThisType =  ICPSpatialParams<GridGeometry, Scalar>;
    using ParentType = FVSpatialParams<GridGeometry, Scalar, ThisType>;
    using GridView = typename GridGeometry::GridView;
    using SubControlVolume = typename GridGeometry::SubControlVolume;
    using Element = typename GridView::template Codim<0>::Entity;
    using GlobalPosition = typename Element::Geometry::GlobalCoordinate;
    using PcKrSwCurve = FluidMatrix::BrooksCoreyDefault<Scalar>;

public:
    // type used for the permeability (i.e. tensor or scalar)
    using PermeabilityType = Scalar;
    // [[/codeblock]]
    // #### Using porosity and permeability laws to return the updated values
    // Due due calcium carbonate precipitation the porosity and the permeability change with time. At first, the initial values for the porosity and permeability are defined. Moreover, the functions that return the updated values, based on the chosen laws are defined.
    // [[codeblock]]
    ICPSpatialParams(std::shared_ptr<const GridGeometry> gridGeometry)
    : ParentType(gridGeometry)
    , pcKrSw_("SpatialParams") // initialize from input file
    {
        // We read reference values for porosity and permeability from the input
        referencePorosity_     = getParam<Scalar>("SpatialParams.ReferencePorosity", 0.4);
        referencePermeability_ = getParam<Scalar>("SpatialParams.ReferencePermeability", 2.e-10);
    }

    // We return the reference or initial porosity.
    //This reference porosity is the porosity, for which the permeability is know and set as reference permeability
    Scalar referencePorosity(const Element& element, const SubControlVolume &scv) const
    { return referencePorosity_; }

    // We return the volume fraction of the inert (unreactive) component
    template<class SolidSystem>
    Scalar inertVolumeFractionAtPos(const GlobalPosition& globalPos, int compIdx) const
    { return 1.0-referencePorosity_; }

    // [[codeblock]]
    // We return the updated porosity using the specified porosity law
    template<class ElementSolution>
    Scalar porosity(const Element& element,
                    const SubControlVolume& scv,
                    const ElementSolution& elemSol) const
    { return poroLaw_.evaluatePorosity(element, scv, elemSol, referencePorosity_); }

    // We return the updated permeability using the specified permeability law
    template<class ElementSolution>
    PermeabilityType permeability(const Element& element,
                        const SubControlVolume& scv,
                        const ElementSolution& elemSol) const
    {
        const auto poro = porosity(element, scv, elemSol);
        return permLaw_.evaluatePermeability(referencePermeability_, referencePorosity_, poro);
    }
    // [[/codeblock]]
    // Return the fluid matrix interactions (here only Brooks-Corey curves)
    auto fluidMatrixInteractionAtPos(const GlobalPosition& globalPos) const
    { return makeFluidMatrixInteraction(pcKrSw_); }

    // Define the wetting phase
    template<class FluidSystem>
    int wettingPhaseAtPos(const GlobalPosition& globalPos) const
    { return FluidSystem::H2OIdx; }

    // The remainder of this class contains the data members and defines the porosity law which describes the change of porosity due to calcium carbonate precipitation.
    // Additionally the change of porosity results in a change of permeability. This relation is described in the permeability law, in this case the Kozeny-Carman porosity-permeability relation
    // [[codeblock]]
private:

    const PcKrSwCurve pcKrSw_;
    // Setting porosity precipitation as the porosity law
    PorosityPrecipitation<Scalar, /*numFluidComponents*/9, /*activeComponents*/2> poroLaw_;
    // Setting the Kozeny-Carman porosity-permeability relation as the permeability law
    PermeabilityKozenyCarman<PermeabilityType> permLaw_;
    //The reference porosity is the porosity, for which the permeability is know and set as reference permeability
    Scalar referencePorosity_;
    //The reference permeability is the known (measured) permeability, of the porous medium in the initial condition, before the solid phases change during the simulation
    PermeabilityType referencePermeability_ = 0.0;
};// end class definition of ICPSpatialParams
} //end namespace Dumux
// [[/codeblock]]
// [[/content]]
#endif
