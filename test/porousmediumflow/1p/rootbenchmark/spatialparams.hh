// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup OnePTests
 * \brief Root benchmark
 */

#ifndef DUMUX_ONEP_ROOTBENCHMARK_SPATIALPARAMS_HH
#define DUMUX_ONEP_ROOTBENCHMARK_SPATIALPARAMS_HH

#include <dumux/porousmediumflow/properties.hh>
#include <dumux/porousmediumflow/fvspatialparams1p.hh>

namespace Dumux {

/*!
 * \ingroup OnePTests
 * \brief Root benchmark spatial params
 */
template<class GridGeometry, class Scalar>
class RootBenchmarkSpatialParams
: public FVPorousMediumFlowSpatialParamsOneP<GridGeometry, Scalar,
                                         RootBenchmarkSpatialParams<GridGeometry, Scalar>>
{
    using FVElementGeometry = typename GridGeometry::LocalView;
    using SubControlVolume = typename FVElementGeometry::SubControlVolume;
    using GridView = typename GridGeometry::GridView;
    using Element = typename GridView::template Codim<0>::Entity;

    using ThisType = RootBenchmarkSpatialParams<GridGeometry, Scalar>;
    using ParentType = FVPorousMediumFlowSpatialParamsOneP<GridGeometry, Scalar, ThisType>;

    using GlobalPosition = typename Element::Geometry::GlobalCoordinate;

public:
    // export permeability type
    using PermeabilityType = Scalar;

    RootBenchmarkSpatialParams(std::shared_ptr<const GridGeometry> gridGeometry)
    : ParentType(gridGeometry)
    {
        radius_ = getParam<double>("Problem.Radius");
        const auto rho = getParam<double>("Component.LiquidDensity");
        const auto mu = getParam<double>("Component.LiquidKinematicViscosity")*rho;
        Kx_ = getParam<double>("Problem.AxialConductivity")/(86400.0*1000*9.81)*1e-6; // cm^3/day -> m^4/(Pa*s)
        Kr_ = getParam<double>("Problem.RadialConductivity")/(86400.0*1000*9.81); // 1/day -> m/(Pa*s)
        permeability_ = Kx_*mu/(M_PI*radius_*radius_);
    }

    /*!
     * \brief Function for defining the (intrinsic) permeability \f$[m^2]\f$.
     */
    template<class ElementSolution>
    PermeabilityType permeability(const Element& element,
                                  const SubControlVolume& scv,
                                  const ElementSolution& elemSol) const
    {
        return permeability_;
    }

    /*!
     * \brief Returns the root radius in \f$[m]\f$.
     */
    Scalar radius() const
    { return radius_; }

    /*!
     * \brief The axial conductivity \f$[m^4 Pa^{-1} s^{-1}]\f$.
     */
    Scalar axialConductivity() const
    { return Kx_; }

    /*!
     * \brief The radial conductivity \f$[m Pa^{-1} s^{-1}]\f$.
     */
    Scalar radialConductivity() const
    { return Kr_; }

    /*!
     * \brief The porosity \f$\mathrm{[-]}\f$.
     */
    Scalar porosityAtPos(const GlobalPosition& globalPos) const
    { return 1.0; }

    /*!
     * \brief Returns how much the domain is extruded at a given sub-control volume.
     * Assume circular cross section with given radius
     */
    Scalar extrusionFactorAtPos(const GlobalPosition& globalPos) const
    { return M_PI*radius_*radius_; }

private:
    Scalar radius_, permeability_, Kx_, Kr_;
    static constexpr Scalar eps_ = 1e-8;
};

} // end namespace Dumux

#endif
