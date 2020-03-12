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
 * \brief The spatial parameters class for the poroelastic domain in the
 *        elastic single-phase facet coupling test.
 */
#ifndef DUMUX_ANALYTIC_CRACK_POROELASTIC_SPATIAL_PARAMS_HH
#define DUMUX_ANALYTIC_CRACK_POROELASTIC_SPATIAL_PARAMS_HH

#include <dumux/geomechanics/lameparams.hh>
#include <dumux/material/spatialparams/fvporoelastic.hh>
#include <dumux/material/fluidmatrixinteractions/porositydeformation.hh>

namespace Dumux {

/*!
 * \brief Definition of the spatial parameters for the poro-elastic
 *        sub-problem in the coupled poro-mechanical el1p problem.
 */
template<class Scalar, class GridGeometry>
class PoroElasticSpatialParams : public FVSpatialParamsPoroElastic< Scalar,
                                                                    GridGeometry,
                                                                    PoroElasticSpatialParams<Scalar, GridGeometry> >
{
    using ThisType = PoroElasticSpatialParams<Scalar, GridGeometry>;
    using ParentType = FVSpatialParamsPoroElastic<Scalar, GridGeometry, ThisType>;

    using GridView = typename GridGeometry::GridView;
    using Element = typename GridView::template Codim<0>::Entity;
    using GlobalPosition = typename Element::Geometry::GlobalCoordinate;

public:
    //! Export the type of the lame parameters
    using LameParams = Dumux::LameParams<Scalar>;

    PoroElasticSpatialParams(std::shared_ptr<const GridGeometry> gridGeometry, const std::string& paramGroup = "")
    : ParentType(gridGeometry)
    , initPorosity_(getParamFromGroup<Scalar>(paramGroup, "SpatialParams.InitialPorosity"))
    , alpha_(getParamFromGroup<Scalar>(paramGroup, "SpatialParams.BiotCoefficient"))
    {
        // Young's modulus [Pa]
       Scalar E = getParamFromGroup<Scalar>(paramGroup, "SpatialParams.EModule");
       // Poisson's ratio [-]
       Scalar nu = getParamFromGroup<Scalar>(paramGroup, "SpatialParams.Nu");
       // Lame parameters [Pa]
       lameParams_.setLambda( (E * nu) / ((1 + nu)*(1 - 2 * nu)) );
       lameParams_.setMu( E / (2 * (1 + nu)) );
    }

    //! Defines the Lame parameters.
    const LameParams& lameParamsAtPos(const GlobalPosition& globalPos) const
    { return lameParams_; }

    //! Returns the porosity of the porous medium.
    Scalar porosityAtPos(const GlobalPosition& elemSol) const
    { return initPorosity_; }

    //! Returns the Biot coefficient of the porous medium.
    Scalar biotCoefficientAtPos(const GlobalPosition& globalPos) const
    { return alpha_; }

private:
    Scalar initPorosity_;
    Scalar alpha_;
    LameParams lameParams_;
};

} // end namespace Dumux

#endif
