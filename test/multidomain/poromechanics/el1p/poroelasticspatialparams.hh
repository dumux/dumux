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
 * \ingroup MultiDomain
 * \ingroup Geomechanics
 * \ingroup PoroElastic
 * \brief Definition of the spatial parameters for the poro-elastic problem
 */
#ifndef DUMUX_POROELASTIC_SPATIAL_PARAMS_HH
#define DUMUX_POROELASTIC_SPATIAL_PARAMS_HH

#include <dumux/geomechanics/lameparams.hh>
#include <dumux/material/spatialparams/fvporoelastic.hh>

namespace Dumux {

/*!
 * \ingroup MultiDomain
 * \ingroup Geomechanics
 * \ingroup PoroElastic
 * \brief Definition of the spatial parameters for the poro-elastic
 *        sub-problem in the coupled poro-mechanical el1p problem.
 */
template<class Scalar, class FVGridGeometry>
class PoroElasticSpatialParams : public FVSpatialParamsPoroElastic< Scalar,
                                                                    FVGridGeometry,
                                                                    PoroElasticSpatialParams<Scalar, FVGridGeometry> >
{
    using ThisType = PoroElasticSpatialParams<Scalar, FVGridGeometry>;
    using ParentType = FVSpatialParamsPoroElastic<Scalar, FVGridGeometry, ThisType>;

    using SubControlVolume = typename FVGridGeometry::SubControlVolume;
    using GridView = typename FVGridGeometry::GridView;
    using Element = typename GridView::template Codim<0>::Entity;
    using GlobalPosition = typename Element::Geometry::GlobalCoordinate;

public:
    //! export the type of the lame parameters
    using LameParams = Dumux::LameParams<Scalar>;

    //! The constructor
    PoroElasticSpatialParams(std::shared_ptr<const FVGridGeometry> fvGridGeometry)
    : ParentType(fvGridGeometry)
    , porosity_(getParam<Scalar>("SpatialParams.Porosity"))
    {
        // Young's modulus [Pa]
       Scalar E = 6.e9;
       // Poisson's ratio [-]
       Scalar nu = 0.2;
       // Lame parameters [Pa]
       lameParams_.setLambda( (E * nu) / ((1 + nu)*(1 - 2 * nu)) );
       lameParams_.setMu( E / (2 * (1 + nu)) );
    }

    //! Define the Lame parameters
    const LameParams& lameParamsAtPos(const GlobalPosition& globalPos) const
    { return lameParams_; }

    //! Return the porosity of the porous medium
    Scalar porosityAtPos(const GlobalPosition& globalPos) const
    { return porosity_; }

    //! Return the biot coefficient of the porous medium
    Scalar biotCoefficientAtPos(const GlobalPosition& globalPos) const
    { return 1.0; }

private:
    Scalar porosity_;
    LameParams lameParams_;
};
} // end namespace Dumux
#endif
