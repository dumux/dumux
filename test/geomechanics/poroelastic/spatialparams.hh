// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
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
     *       higher values is done to see a more pronounced effect in stresses.
     *
     * \param globalPos The global position
     */
    Scalar effectivePorePressureAtPos(const GlobalPosition& globalPos) const
    {
        using std::sin;
        const auto x = globalPos[0];
        const auto y = globalPos[1];
        return (x-x*x)*sin(2.*pi*y) + 10;
    }

    /*!
     * \brief Returns the effective pore pressure gradient
     * \param globalPos The global position
     */
    GlobalPosition effectivePorePressureGradient(const GlobalPosition& globalPos) const
    {
        using std::sin;
        using std::cos;
        const auto x = globalPos[0];
        const auto y = globalPos[1];
        return {{
            (1.-2.*x)*sin(2.*pi*y),
            2.*pi*(x-x*x)*cos(2.*pi*y),
        }};
    }

    /*!
     * \brief Returns the effective fluid density.
     *
     * \param globalPos The global position
     */
    Scalar effectiveFluidDensityAtPos(const GlobalPosition& globalPos) const
    {
        return FluidSystem::density(
            effectivePorePressureAtPos(globalPos),
            this->temperatureAtPos(globalPos)
        );
    }

    //! Returns the Biot coefficient of the porous medium.
    Scalar biotCoefficientAtPos(const GlobalPosition& globalPos) const
    { return 1.0; }

private:
    LameParams lameParams_;
};

} // end namespace Dumux

#endif
