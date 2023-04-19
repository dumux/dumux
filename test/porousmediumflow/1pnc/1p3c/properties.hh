// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/**
 * \file
 * \ingroup OnePNCTests
 * \brief Definition of a problem for a 1p3c problem:
 *        Component transport of N2, CO2 and H2 using the Maxwell-Stefan diffusion law.
 */

#ifndef DUMUX_1P3C_TEST_PROBLEM_PROPERTIES_HH
#define DUMUX_1P3C_TEST_PROBLEM_PROPERTIES_HH

#if HAVE_DUNE_UGGRID
#include <dune/grid/uggrid.hh>
#endif
#include <dune/grid/yaspgrid.hh>

#include <dumux/discretization/cctpfa.hh>
#include <dumux/discretization/box.hh>
#include <dumux/discretization/evalsolution.hh>
#include <dumux/discretization/evalgradients.hh>
#include <dumux/porousmediumflow/1pnc/model.hh>

#include <dumux/material/idealgas.hh>
#include <dumux/material/fluidsystems/base.hh>

#include "problem.hh"
#include "../1p2c/spatialparams.hh"


#include <dumux/flux/maxwellstefanslaw.hh>

namespace Dumux::Properties {

// Create new type tags
namespace TTag {
struct MaxwellStefanOnePThreeCTest { using InheritsFrom = std::tuple<OnePNC>; };
struct MaxwellStefanOnePThreeCTestBox { using InheritsFrom = std::tuple<MaxwellStefanOnePThreeCTest, BoxModel>; };
struct MaxwellStefanOnePThreeCTestCCTpfa { using InheritsFrom = std::tuple<MaxwellStefanOnePThreeCTest, CCTpfaModel>; };
} // end namespace TTag

// Set the grid type
#if HAVE_DUNE_UGGRID
template<class TypeTag>
struct Grid<TypeTag, TTag::MaxwellStefanOnePThreeCTest> { using type = Dune::UGGrid<2>; };
#else
template<class TypeTag>
struct Grid<TypeTag, TTag::MaxwellStefanOnePThreeCTest> { using type = Dune::YaspGrid<2>; };
#endif

// Set the problem property
template<class TypeTag>
struct Problem<TypeTag, TTag::MaxwellStefanOnePThreeCTest> { using type = MaxwellStefanOnePThreeCTestProblem<TypeTag>; };

/*!
 * \ingroup OnePNCTests
 * \brief  A simple fluid system with three components for testing the  multi-component diffusion with the Maxwell-Stefan formulation.
 */
template<class TypeTag>
class H2N2CO2FluidSystem
: public FluidSystems::Base<GetPropType<TypeTag, Properties::Scalar>, H2N2CO2FluidSystem<TypeTag>>

{
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using ThisType = H2N2CO2FluidSystem<TypeTag>;
    using Base = FluidSystems::Base<Scalar, ThisType>;
    using IdealGas = Dumux::IdealGas<Scalar>;

public:
    //! The number of phases
    static constexpr int numPhases = 1;
    static constexpr int numComponents = 3;

    static constexpr int H2Idx = 0; //first major component
    static constexpr int N2Idx = 1; //second major component
    static constexpr int CO2Idx = 2; //secondary component

    //! Human readable component name (index compIdx) (for vtk output)
    static std::string componentName(int compIdx)
    { return "MaxwellStefan_" + std::to_string(compIdx); }

    //! Human readable phase name (index phaseIdx) (for velocity vtk output)
    static std::string phaseName(int phaseIdx = 0)
    { return "Gas"; }

    //! Molar mass in kg/mol of the component with index compIdx
    static Scalar molarMass(unsigned int compIdx)
    {
        switch (compIdx)
        {
        case H2Idx: return 0.002;
        case N2Idx: return 0.028;
        case CO2Idx:return 0.044;
        }
        DUNE_THROW(Dune::InvalidStateException, "Invalid component index " << compIdx);;
    }

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
     */
    template <class FluidState>
    static Scalar density(const FluidState &fluidState,
                          const int phaseIdx)
    {
        Scalar T = fluidState.temperature(phaseIdx);
        Scalar p = fluidState.pressure(phaseIdx);
        return IdealGas::molarDensity(T, p) * fluidState.averageMolarMass(0);
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
        Scalar T = fluidState.temperature(phaseIdx);
        Scalar p = fluidState.pressure(phaseIdx);
        return IdealGas::molarDensity(T,p);
    }
};

// Set fluid configuration
template<class TypeTag>
struct FluidSystem<TypeTag, TTag::MaxwellStefanOnePThreeCTest>
{using type = H2N2CO2FluidSystem<TypeTag>; };

// Set the spatial parameters
template<class TypeTag>
struct SpatialParams<TypeTag, TTag::MaxwellStefanOnePThreeCTest>
{
    using GridGeometry = GetPropType<TypeTag, Properties::GridGeometry>;
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using type = OnePNCTestSpatialParams<GridGeometry, Scalar>;
};

// Define whether mole(true) or mass (false) fractions are used
template<class TypeTag>
struct UseMoles<TypeTag, TTag::MaxwellStefanOnePThreeCTest> { static constexpr bool value = true; };

//! Here we set FicksLaw or MaxwellStefansLaw
template<class TypeTag>
struct MolecularDiffusionType<TypeTag, TTag::MaxwellStefanOnePThreeCTest> { using type = MaxwellStefansLaw<TypeTag>; };

} // end namespace Properties

#endif
