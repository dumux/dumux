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
 * \ingroup GeomechanicsTests
 * \brief Definition of the spatial parameters for the poro-elastic problem.
 */

#ifndef DUMUX_POROELASTIC_SPATIAL_PARAMS_HH
#define DUMUX_POROELASTIC_SPATIAL_PARAMS_HH

#include <dumux/geomechanics/lameparams.hh>
#include <dumux/material/spatialparams/fvporoelastic.hh>
#include <dumux/material/fluidmatrixinteractions/porositydeformation.hh>

namespace Dumux {

/*!
 * \ingroup GeomechanicsTests
 * \brief Definition of the spatial parameters for the poro-elastic problem.
 */
template<class Scalar, class GridGeometry>
class PoroElasticSpatialParams : public FVSpatialParamsPoroElastic< Scalar,
                                                                    GridGeometry,
                                                                    PoroElasticSpatialParams<Scalar, GridGeometry> >
{
    using ThisType = PoroElasticSpatialParams<Scalar, GridGeometry>;
    using ParentType = FVSpatialParamsPoroElastic<Scalar, GridGeometry, ThisType>;

    using SubControlVolume = typename GridGeometry::SubControlVolume;
    using GridView = typename GridGeometry::GridView;
    using Element = typename GridView::template Codim<0>::Entity;
    using GlobalPosition = typename Element::Geometry::GlobalCoordinate;

public:
    //! Export the type of the lame parameters
    using LameParams = Dumux::LameParams<Scalar>;

    PoroElasticSpatialParams(std::shared_ptr<const GridGeometry> gridGeometry)
    : ParentType(gridGeometry)
    {
        lameParams_.setLambda(2);
        lameParams_.setMu(2);
    }

    //! Defines the Lame parameters.
    const LameParams& lameParamsAtPos(const GlobalPosition& globalPos) const
    { return lameParams_; }

    //! Returns the porosity of the porous medium.
    template< class ElemSol >
    Scalar porosity(const Element& element,
                    const SubControlVolume& scv,
                    const ElemSol& elemSol) const
    {
        PorosityDeformation<Scalar> poroLaw;
        return poroLaw.evaluatePorosity(this->gridGeometry(), element, scv, elemSol, /*refPoro*/0.3);
    }

    //! Returns the Biot coefficient of the porous medium.
    Scalar biotCoefficientAtPos(const GlobalPosition& globalPos) const
    { return 1.0; }

private:
    LameParams lameParams_;
};
} // end namespace Dumux
#endif
