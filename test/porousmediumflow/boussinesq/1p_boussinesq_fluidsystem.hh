// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup FluidSystems
 * \brief Dimensionless two-component single-phase fluid system for Boussinesq
 *        porous media simulations.
 *
 * All flow properties (density, viscosity) are unity — the Boussinesq
 * approximation treats the fluid as incompressible everywhere **except**
 * in the buoyancy term of Darcy's law.  That concentration-dependent
 * buoyancy must be supplied by a matching Darcy-flux implementation
 * (see BoussinesqCVFEDarcyLaw).
 *
 * The Rayleigh number Ra (read from \c DimensionlessNumbers.Ra) enters as the
 * **inverse diffusion coefficient** so that the dimensionless transport equation reads
 *   ∂C/∂t + q·∇C = (1/Ra) ∇²C.
 */
#ifndef DUMUX_FLUIDSYSTEMS_BOUSSINESQ_HH
#define DUMUX_FLUIDSYSTEMS_BOUSSINESQ_HH

#include <string>

#include <dumux/material/fluidsystems/base.hh>
#include <dumux/common/parameters.hh>

namespace Dumux::FluidSystems {

/*!
 * \ingroup FluidSystems
 * \brief Dimensionless Boussinesq fluid system for porous-media dissolution benchmarks.
 *
 * Usage:
 *  - Set \c DimensionlessNumbers.Ra in the input file.
 *  - Combine with BoussinesqCVFEDarcyLaw to get concentration-driven buoyancy
 *    while keeping ρ = 1 in the storage and continuity terms.
 */
template<class Scalar>
class BoussinesqFluid : public Base<Scalar, BoussinesqFluid<Scalar>>
{
    using ThisType = BoussinesqFluid<Scalar>;

public:
    using ParameterCache = NullParameterCache;

    static constexpr int numPhases     = 1;
    static constexpr int numComponents = 2;

    static constexpr int phase0Idx  = 0;
    static constexpr int solventIdx = 0; //!< main component (solvent)
    static constexpr int soluteIdx  = 1; //!< dissolved component (solute)

    //! Read Ra from the input file (key: DimensionlessNumbers.Ra)
    static void init()
    {
        rayleighNumber_ = getParam<Scalar>("DimensionlessNumbers.Ra", 100.0);
    }

    static std::string phaseName(int /*phaseIdx*/ = 0) { return "liquid"; }

    static std::string componentName(int compIdx)
    { return compIdx == solventIdx ? "Solvent" : "Solute"; }

    static std::string name() { return "BoussinesqFluid"; }

    static constexpr bool isMiscible()                         { return false; }
    static constexpr bool isGas(int /*phaseIdx*/ = 0)          { return false; }
    static constexpr bool isIdealMixture(int /*phaseIdx*/ = 0) { return true;  }
    // density is constant (Boussinesq), not compressible w.r.t. pressure
    static constexpr bool isCompressible(int /*phaseIdx*/ = 0) { return false; }
    static constexpr bool viscosityIsConstant(int /*phaseIdx*/ = 0) { return true; }

    //! Both components have unit molar mass → mole fraction == mass fraction
    static Scalar molarMass(int /*compIdx*/) { return 1.0; }

    // ----- density (constant = 1, Boussinesq approximation) ---------------
    using Base<Scalar, ThisType>::density;
    template<class FluidState>
    static Scalar density(const FluidState&, int /*phaseIdx*/)
    { return 1.0; }

    using Base<Scalar, ThisType>::molarDensity;
    template<class FluidState>
    static Scalar molarDensity(const FluidState&, int /*phaseIdx*/)
    { return 1.0; }

    // ----- viscosity (dimensionless = 1) -----------------------------------
    using Base<Scalar, ThisType>::viscosity;
    template<class FluidState>
    static Scalar viscosity(const FluidState&, int /*phaseIdx*/)
    { return 1.0; }

    // ----- fugacity (ideal mixture) ----------------------------------------
    using Base<Scalar, ThisType>::fugacityCoefficient;
    template<class FluidState>
    static Scalar fugacityCoefficient(const FluidState&, int /*phaseIdx*/, int /*compIdx*/)
    { return 1.0; }

    // ----- diffusion: D = 1/Ra so that transport eq. has (1/Ra)∇²C ---------
    using Base<Scalar, ThisType>::diffusionCoefficient;
    template<class FluidState>
    static Scalar diffusionCoefficient(const FluidState&, int /*phaseIdx*/, int /*compIdx*/)
    { return 1.0 / rayleighNumber_; }

    using Base<Scalar, ThisType>::binaryDiffusionCoefficient;
    template<class FluidState>
    static Scalar binaryDiffusionCoefficient(const FluidState&,
                                             int /*phaseIdx*/,
                                             int /*compIIdx*/,
                                             int /*compJIdx*/)
    { return 1.0 / rayleighNumber_; }

    //! Volumetric expansion coefficient β_i for component compIdx.
    //! In the dimensionless formulation β_solute = 1 so that the buoyancy
    //! density (1 + Σ β_i·x_i) equals (1 + C) when C ∈ [0,1].
    static Scalar volumetricExpansionCoeff(int compIdx)
    { return (compIdx == soluteIdx) ? 1.0 : 0.0; }

    //! The current Rayleigh number (for diagnostic output)
    static Scalar rayleighNumber() { return rayleighNumber_; }

    // ----- Boussinesq-specific interface used by BoussinesqVorticityLocalResidual ----
    //
    // The local residual needs fluid properties evaluated at face-averaged
    // concentrations.  The primary interface therefore takes the concentration
    // array; derived fluid systems can override these to implement
    // concentration-dependent viscosity or diffusion.

    //! Reference density ρ₀ (constant under Boussinesq approximation)
    static Scalar referenceDensity() { return 1.0; }

    //! Solutal expansion coefficient for the i-th transported solute (0-based index)
    static Scalar solutalExpansionCoefficient(int /*i*/ = 0) { return 1.0; }

    /*!
     * \brief Dynamic viscosity evaluated at the given concentration vector.
     *
     * The default implementation returns the constant dimensionless value 1.
     * Override in a derived class to implement concentration-dependent viscosity,
     * e.g. μ(C) from an experimental correlation.
     *
     * \param C  array of face-averaged solute concentrations (0-based index)
     */
    template<class ConcentrationArray>
    static Scalar viscosity(const ConcentrationArray& /*C*/) { return 1.0; }

    /*!
     * \brief Molecular diffusion coefficient for solute i at the given concentrations.
     *
     * The default returns 1/Ra (constant).  Override for concentration-dependent D,
     * e.g. D(C) from a Stokes-Einstein or fitted correlation.
     *
     * \param i  0-based solute index
     * \param C  array of face-averaged solute concentrations
     */
    template<class ConcentrationArray>
    static Scalar solutalDiffusionCoefficient(int /*i*/, const ConcentrationArray& /*C*/)
    { return 1.0 / rayleighNumber_; }

private:
    inline static Scalar rayleighNumber_ = 100.0;
};

} // namespace Dumux::FluidSystems

#endif