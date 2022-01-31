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
#include <dumux/geomechanics/poroelastic/fvspatialparams.hh>
#include <dumux/material/fluidmatrixinteractions/porositydeformation.hh>

namespace Dumux {

/*!
 * \ingroup GeomechanicsTests
 * \brief Definition of the spatial parameters for the poro-elastic problem.
 */
template<class Scalar, class GridGeometry, class FluidSystem, class PrimaryVariables, class Indices>
class PoroElasticSpatialParams : public FVPoroElasticSpatialParams< GridGeometry,
                                                                    Scalar,
                                                                    PoroElasticSpatialParams<Scalar, GridGeometry, FluidSystem, PrimaryVariables, Indices> >
{
    using ThisType = PoroElasticSpatialParams<Scalar, GridGeometry, FluidSystem, PrimaryVariables, Indices>;
    using ParentType = FVPoroElasticSpatialParams<GridGeometry, Scalar, ThisType>;

    using SubControlVolume = typename GridGeometry::SubControlVolume;
    using GridView = typename GridGeometry::GridView;
    using Element = typename GridView::template Codim<0>::Entity;
    using GlobalPosition = typename Element::Geometry::GlobalCoordinate;

    static constexpr Scalar pi = M_PI;
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

    /*!
     * \brief Returns the effective pore pressure
     *
     * \note We use the x-displacement as pressure solution. The shift to
     *       higher values is done to see a mor pronounced effect in stresses.
     *
     * \param globalPos The global position
     */
    Scalar effectivePorePressureAtPos(const GlobalPosition& globalPos) const
    { return exactSolution(globalPos)[0] + 10; }

    /*!
     * \brief Returns the effective fluid density.
     *
     * \param globalPos The global position
     */
    Scalar effectiveFluidDensityAtPos(const GlobalPosition& globalPos) const
    {
        // This test uses the constant component, obtain density only once
        static const Scalar rho = FluidSystem::density(
            effectivePorePressureAtPos(globalPos), this->temperatureAtPos(globalPos)
        );
        return rho;
    }

    //! Returns the Biot coefficient of the porous medium.
    Scalar biotCoefficientAtPos(const GlobalPosition& globalPos) const
    { return 1.0; }

    /*!
     * \brief Evaluates the exact displacement to this problem at a given position.
     */
    PrimaryVariables exactSolution(const GlobalPosition& globalPos) const
    {
        using std::sin;

        const auto x = globalPos[0];
        const auto y = globalPos[1];

        PrimaryVariables exact(0.0);
        exact[Indices::momentum(/*x-dir*/0)] = (x-x*x)*sin(2*pi*y);
        exact[Indices::momentum(/*y-dir*/1)] = sin(2*pi*x)*sin(2*pi*y);
        return exact;
    }

private:
    LameParams lameParams_;
};

} // end namespace Dumux

#endif
