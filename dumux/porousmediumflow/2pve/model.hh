// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup TwoPVEModel
 * \brief Two-phase flow in vertical equilibrium, solved on a coarse grid of vertical columns.
 *
 * This model describes immiscible two-phase flow of a wetting phase \f$w\f$ and a less dense
 * nonwetting phase \f$n\f$ in a formation of height \f$H\f$, for instance the injection of CO2 into
 * a saline aquifer. The phases are assumed to be segregated by gravity and in hydrostatic equilibrium
 * in the vertical direction (vertical equilibrium) \cite Nordbotten2012. The mass balance equations
 * are solved on a coarse grid that consists of a single layer of columns spanning the height of the
 * formation, and the vertical distribution of the saturations, pressures and mobilities is
 * reconstructed in each column on a fine grid that resolves the vertical direction.
 * See \cite Buntic2025 for a detailed description of the model.
 *
 * On the coarse level, the mass balance of each phase \f$\alpha \in \{ w, n \}\f$ reads
 \f[
 \frac{\partial (\bar\phi \varrho_\alpha \bar S_\alpha)}{\partial t}
 -
 \nabla \cdot \left\{ \varrho_\alpha \bar\lambda_\alpha \bar K \nabla P_\alpha \right\} - q_\alpha = 0,
 \f]
 * where:
 * * \f$ z \f$ is the height above the bottom of the formation,
 * * \f$ \bar\phi = \frac{1}{H} \int_0^H \phi \, \mathrm{d}z \f$ is the vertical average of the porosity \f$\phi\f$,
 * * \f$ \bar K = \frac{1}{H} \int_0^H k \, \mathrm{d}z \f$ is the vertical average of the scalar permeability \f$k\f$,
 * * \f$ \bar S_\alpha \f$ is the porosity-weighted vertical average of the saturation of phase \f$\alpha\f$,
 * * \f$ \bar\lambda_\alpha = \int_0^H k \lambda_\alpha \, \mathrm{d}z \big/ \int_0^H k \, \mathrm{d}z \f$ is the
 *   permeability-weighted vertical average of the mobility \f$ \lambda_\alpha = k_{r\alpha}/\mu_\alpha \f$ of phase \f$\alpha\f$,
 *   with the relative permeability \f$ k_{r\alpha} \f$ and the dynamic viscosity \f$ \mu_\alpha \f$,
 * * \f$ \varrho_\alpha \f$ is the mass density of phase \f$\alpha\f$,
 * * \f$ P_\alpha \f$ is the pressure of phase \f$\alpha\f$ at the bottom of the column,
 * * \f$ q_\alpha \f$ is a source or sink term.
 *
 * The centers of all coarse-level elements lie at the same height, so gravity only enters the
 * coarse level through the reconstruction. The primary variables are \f$ P_w \f$ and \f$ \bar S_n \f$.
 * Extrapolating the hydrostatic nonwetting-phase pressure to the bottom of the column gives the
 * coarse-level capillary pressure
 \f[
 P_n - P_w = p_e - (\varrho_w - \varrho_n) g z_p,
 \f]
 * where \f$ p_e \f$ is the entry pressure, \f$ g \f$ the norm of the gravitational acceleration and
 * \f$ z_p \f$ the gas plume distance, the height of the lower boundary of the mobile nonwetting phase.
 *
 * Above the gas plume distance, the wetting-phase saturation follows from the Brooks-Corey
 * capillary pressure \cite brooks1964hydrau in hydrostatic equilibrium,
 \f[
 S_w(z) = S_{wr} + (1 - S_{wr} - S_{nr}) \left( 1 + \frac{(\varrho_w - \varrho_n) g (z - z_p)}{p_e} \right)^{-\lambda},
 \f]
 * with the residual saturations \f$ S_{wr} \f$ and \f$ S_{nr} \f$ and the Brooks-Corey parameter \f$ \lambda \f$.
 * Between the minimum gas plume distance of all previous time steps \f$ z_{p,\min} \f$ and \f$ z_p \f$,
 * the nonwetting phase is residually trapped and \f$ S_w = 1 - S_{nr} \f$ \cite Doster2013.
 * Below \f$ z_{p,\min} \f$, \f$ S_w = 1 \f$. The gas plume distance is the height for which the
 * vertical integral of this saturation profile equals \f$ H \bar S_w \f$. The fine-level mobilities
 * are the cell averages of the Brooks-Corey relative permeabilities for this profile divided by the viscosities.
 *
 * The current implementation has the following restrictions:
 * - only two- and three-dimensional, structured, axis-aligned grids are
 *   supported
 * - the fine grid must be uniform in the vertical direction
 * - the coarse and fine grids must have identical domain bounds and matching
 *   subdivisions in all horizontal directions
 * - each coarse column must contain exactly one coarse cell in the vertical
 *   direction and the same number of fine cells
 * - each fine-cell center must be contained unambiguously in exactly one
 *   coarse cell
 * - gravity must be nonzero and aligned with the vertical coordinate axis
 * - the first fluid phase must be the wetting phase, and it must be denser than the nonwetting phase
 * - the reconstruction uses a Brooks-Corey material law
 * - the material law parameters, the densities and the viscosities are constant within a column,
 *   the densities and viscosities are evaluated at the coarse-level wetting-phase pressure
 * - the porosity must be uniform within a column
 * - the permeability must be a scalar
 * - only the \f$ P_w - \bar S_n \f$ primary-variable formulation is supported
 * - only the cell-centered two-point flux approximation is supported
 * - only isothermal flow is supported
 *
 * These are restrictions of the current implementation and are not general
 * restrictions of vertical-equilibrium models.
 *
 * A problem using this model has to provide `fineLevelView()`, returning the
 * Dumux::TwoPVEFineLevelView that connects the coarse level to the fine level,
 * and spatial parameters deriving from Dumux::TwoPVESpatialParams. After each
 * time step, the fine-level solution and the history of the columns have to be
 * updated with Dumux::TwoPVEFineLevelView::updateSol, and cached coarse-level
 * volume variables have to be updated afterwards.
 */

#ifndef DUMUX_TWOPVE_MODEL_HH
#define DUMUX_TWOPVE_MODEL_HH

#include <tuple>

#include <dumux/discretization/method.hh>
#include <dumux/porousmediumflow/2p/model.hh>
#include <dumux/porousmediumflow/2pve/volumevariables.hh>

namespace Dumux::Properties {

namespace TTag {

//! The type tag for the two-phase vertical-equilibrium model, derived from the two-phase model
struct TwoPVE
{
    using InheritsFrom = std::tuple<TwoP>;
};

} // namespace TTag

//! Set the volume variables property, the only property in which the VE model differs from the two-phase model
template<class TypeTag>
struct VolumeVariables<TypeTag, TTag::TwoPVE>
{
private:
    using PrimaryVariables = GetPropType<TypeTag, Properties::PrimaryVariables>;
    using FluidSystem = GetPropType<TypeTag, Properties::FluidSystem>;
    using FluidState = GetPropType<TypeTag, Properties::FluidState>;
    using SolidSystem = GetPropType<TypeTag, Properties::SolidSystem>;
    using SolidState = GetPropType<TypeTag, Properties::SolidState>;
    using ModelTraits = GetPropType<TypeTag, Properties::ModelTraits>;
    using PermeabilityType = typename GetPropType<TypeTag, Properties::SpatialParams>::PermeabilityType;
    using DiscretizationMethod = typename GetPropType<TypeTag, Properties::GridGeometry>::DiscretizationMethod;
    static_assert(DiscretizationMethod{} == DiscretizationMethods::cctpfa, "The two-phase VE model is only implemented for the cell-centered TPFA discretization");
    static constexpr bool enableBoxInterfaceSolver = getPropValue<TypeTag, Properties::EnableBoxInterfaceSolver>();
    using SaturationReconstruction = TwoPScvSaturationReconstruction<DiscretizationMethod, enableBoxInterfaceSolver>;

    using Traits = TwoPVolumeVariablesTraits<
        PrimaryVariables,
        FluidSystem,
        FluidState,
        SolidSystem,
        SolidState,
        PermeabilityType,
        ModelTraits,
        SaturationReconstruction>;

public:
    using type = TwoPVEVolumeVariables<Traits>;
};

} // namespace Dumux::Properties

#endif
