// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup NavierStokesNCTests
 * \brief The properties of the channel flow test for the multi-component staggered grid (Navier-)Stokes model.
 */
#ifndef DUMUX_CHANNEL_MAXWELL_STEPHAN_TEST_PROPERTIES_HH
#define DUMUX_CHANNEL_MAXWELL_STEPHAN_TEST_PROPERTIES_HH

#include <dune/grid/yaspgrid.hh>

#include <dumux/freeflow/navierstokes/mass/1pnc/model.hh>
#include <dumux/freeflow/navierstokes/mass/problem.hh>

#include <dumux/freeflow/navierstokes/momentum/problem.hh>
#include <dumux/freeflow/navierstokes/momentum/model.hh>

#include <dumux/flux/maxwellstefanslaw.hh>
#include <dumux/material/fluidsystems/base.hh>

#include <dumux/discretization/fcstaggered.hh>
#include <dumux/discretization/cctpfa.hh>

#include <dumux/multidomain/traits.hh>
#include <dumux/multidomain/freeflow/couplingmanager.hh>

#include "problem.hh"

namespace Dumux::Properties {

// Create new type tags
namespace TTag {
struct MaxwellStefanNCTest {};
struct MaxwellStefanTestMomentum { using InheritsFrom = std::tuple<MaxwellStefanNCTest, NavierStokesMomentum, FaceCenteredStaggeredModel>; };
struct MaxwellStefanTestMass { using InheritsFrom = std::tuple<MaxwellStefanNCTest, NavierStokesMassOnePNC, CCTpfaModel>; };
} // end namespace TTag

template<class TypeTag>
struct ReplaceCompEqIdx<TypeTag, TTag::MaxwellStefanNCTest> { static constexpr int value = 0; };

// Set the grid type
template<class TypeTag>
struct Grid<TypeTag, TTag::MaxwellStefanNCTest>
{ using type = Dune::YaspGrid<2>; };

// Set the problem property
template<class TypeTag>
struct Problem<TypeTag, TTag::MaxwellStefanTestMomentum>
{ using type = MaxwellStefanNCTestProblem<TypeTag, Dumux::NavierStokesMomentumProblem<TypeTag>>; };

template<class TypeTag>
struct Problem<TypeTag, TTag::MaxwellStefanTestMass>
{ using type = MaxwellStefanNCTestProblem<TypeTag, Dumux::NavierStokesMassProblem<TypeTag>>; };

template<class TypeTag>
struct EnableGridGeometryCache<TypeTag, TTag::MaxwellStefanNCTest> { static constexpr bool value = true; };
template<class TypeTag>
struct EnableGridFluxVariablesCache<TypeTag, TTag::MaxwellStefanNCTest> { static constexpr bool value = true; };
template<class TypeTag>
struct EnableGridVolumeVariablesCache<TypeTag, TTag::MaxwellStefanNCTest> { static constexpr bool value = true; };

template<class TypeTag>
struct UseMoles<TypeTag, TTag::MaxwellStefanNCTest> { static constexpr bool value = true; };

//! Here we set FicksLaw or MaxwellStefansLaw
template<class TypeTag>
struct MolecularDiffusionType<TypeTag, TTag::MaxwellStefanNCTest> { using type = MaxwellStefansLaw<TypeTag>; };

/*!
 * \ingroup NavierStokesNCTests
 * \brief  A simple fluid system with three components for testing the  multi-component diffusion with the Maxwell-Stefan formulation.
 */
template<class TypeTag>
class MaxwellStefanFluidSystem
: public FluidSystems::Base<GetPropType<TypeTag, Properties::Scalar>, MaxwellStefanFluidSystem<TypeTag>>

{
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using ThisType = MaxwellStefanFluidSystem<TypeTag>;
    using Base = FluidSystems::Base<Scalar, ThisType>;

public:
    //! The number of phases
    static constexpr int numPhases = 1;
    static constexpr int numComponents = 3;

    static constexpr int H2Idx = 0;//first major component
    static constexpr int N2Idx = 1;//second major component
    static constexpr int CO2Idx = 2;//secondary component

    //! Human readable component name (index compIdx) (for vtk output)
    static std::string componentName(int compIdx)
    { return "MaxwellStefan_" + std::to_string(compIdx); }

    //! Human readable phase name (index phaseIdx) (for velocity vtk output)
    static std::string phaseName(int phaseIdx = 0)
    { return "Gas"; }

    //! Molar mass in kg/mol of the component with index compIdx
    static Scalar molarMass(unsigned int compIdx)
    { return 0.02896; }

    //! Returns whether the fluids are miscible
    static constexpr bool isMiscible()
    { return false; }

    //! Returns whether the fluids are compressible
    static constexpr bool isCompressible(int phaseIdx = 0)
    { return false; }

    using Base::binaryDiffusionCoefficient;
   /*!
     * \brief Given a phase's composition, temperature and pressure,
     *        returns the binary diffusion coefficient \f$\mathrm{[m^2/s]}\f$ for components
     *        \f$i\f$ and \f$j\f$ in this phase.
     *
     * \param fluidState An arbitrary fluid state
     * \param phaseIdx The index of the fluid phase to consider
     * \param compIIdx The index of the first component to consider
     * \param compJIdx The index of the second component to consider
     */
    template <class FluidState>
    static Scalar binaryDiffusionCoefficient(const FluidState &fluidState,
                                             int phaseIdx,
                                             int compIIdx,
                                             int compJIdx)
    {
        if (compIIdx > compJIdx)
        {
            using std::swap;
            swap(compIIdx, compJIdx);
        }

        if (compIIdx == H2Idx && compJIdx == N2Idx)
                return 83.3e-6;
        if (compIIdx == H2Idx && compJIdx == CO2Idx)
                return 68.0e-6;
        if (compIIdx == N2Idx && compJIdx == CO2Idx)
                return 16.8e-6;
        DUNE_THROW(Dune::InvalidStateException,
                       "Binary diffusion coefficient of components "
                       << compIIdx << " and " << compJIdx << " is undefined!\n");
    }
    using Base::density;
   /*!
     * \brief Given a phase's composition, temperature, pressure, and
     *        the partial pressures of all components, returns its
     *        density \f$\mathrm{[kg/m^3]}\f$.
     *
     * \param phaseIdx index of the phase
     * \param fluidState the fluid state
     *
     */
    template <class FluidState>
    static Scalar density(const FluidState &fluidState,
                          const int phaseIdx)
    {
      return 1;
    }

    using Base::viscosity;
   /*!
     * \brief Calculates the dynamic viscosity of a fluid phase \f$\mathrm{[Pa*s]}\f$
     *
     * \param fluidState An arbitrary fluid state
     * \param phaseIdx The index of the fluid phase to consider
     */
    template <class FluidState>
    static Scalar viscosity(const FluidState &fluidState,
                            int phaseIdx)
    {
        return 1e-6;
    }

    using Base::molarDensity;
    /*!
     * \brief The molar density \f$\rho_{mol,\alpha}\f$
     *   of a fluid phase \f$\alpha\f$ in \f$\mathrm{[mol/m^3]}\f$
     *
     * The molar density for the simple relation is defined by the
     * mass density \f$\rho_\alpha\f$ and the molar mass of the main component \f$M_\kappa\f$:
     *
     * \f[\rho_{mol,\alpha} = \frac{\rho_\alpha}{M_\kappa} \;.\f]
     */
    template <class FluidState>
    static Scalar molarDensity(const FluidState &fluidState, int phaseIdx)
    {
        return density(fluidState, phaseIdx)/molarMass(0);
    }
};

template<class TypeTag>
struct FluidSystem<TypeTag, TTag::MaxwellStefanNCTest> { using type = MaxwellStefanFluidSystem<TypeTag>; };

template<class TypeTag>
struct CouplingManager<TypeTag, TTag::MaxwellStefanNCTest>
{
    using Traits = MultiDomainTraits<TTag::MaxwellStefanTestMomentum, TTag::MaxwellStefanTestMass>;
    using type = FreeFlowCouplingManager<Traits>;
};


} // end namespace Dumux::Properties

#endif
