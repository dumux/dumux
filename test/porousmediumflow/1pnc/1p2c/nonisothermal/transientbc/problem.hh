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
/**
 * \file
 * \ingroup OnePNCTests
 * \brief Test for the OnePNCModel in combination with the NI model with transient boundary conditions.
 *
 * The simulation domain is a tube with an elevated temperature on the left hand side.
 */
#ifndef DUMUX_1P2CNI_TRANSIENT_BC_TEST_PROBLEM_HH
#define DUMUX_1P2CNI_TRANSIENT_BC_TEST_PROBLEM_HH

#if HAVE_UG
#include <dune/grid/uggrid.hh>
#endif
#include <dune/grid/yaspgrid.hh>

#include <dumux/common/boundarytypes.hh>

#include <dumux/discretization/elementsolution.hh>
#include <dumux/discretization/cctpfa.hh>
#include <dumux/discretization/ccmpfa.hh>
#include <dumux/discretization/box.hh>
#include <dumux/porousmediumflow/1pnc/model.hh>
#include <dumux/porousmediumflow/problem.hh>

#include <dumux/material/fluidsystems/1padapter.hh>
#include <dumux/material/fluidsystems/h2on2.hh>
#include <dumux/material/components/h2o.hh>
#include "../../spatialparams.hh"

namespace Dumux {

template <class TypeTag>
class OnePTwoCNITransientBCProblem;

namespace Properties {
// Create new type tags
namespace TTag {
struct OnePTwoCNITransientBC { using InheritsFrom = std::tuple<OnePNCNI>; };
struct OnePTwoCNITransientBCCCTpfa { using InheritsFrom = std::tuple<OnePTwoCNITransientBC, CCTpfaModel>; };
struct OnePTwoCNITransientBCCCMpfa { using InheritsFrom = std::tuple<OnePTwoCNITransientBC, CCMpfaModel>; };
struct OnePTwoCNITransientBCBox { using InheritsFrom = std::tuple<OnePTwoCNITransientBC, BoxModel>; };
} // end namespace TTag

// Set the grid type
#if HAVE_UG
template<class TypeTag>
struct Grid<TypeTag, TTag::OnePTwoCNITransientBC> { using type = Dune::UGGrid<2>; };
#else
template<class TypeTag>
struct Grid<TypeTag, TTag::OnePTwoCNITransientBC> { using type = Dune::YaspGrid<2>; };
#endif

// Set the problem property
template<class TypeTag>
struct Problem<TypeTag, TTag::OnePTwoCNITransientBC> { using type = OnePTwoCNITransientBCProblem<TypeTag>; };

// Set fluid configuration
template<class TypeTag>
struct FluidSystem<TypeTag, TTag::OnePTwoCNITransientBC>
{
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using H2ON2 = FluidSystems::H2ON2<Scalar, FluidSystems::H2ON2DefaultPolicy</*simplified=*/true>>;
    using type = FluidSystems::OnePAdapter<H2ON2, H2ON2::liquidPhaseIdx>;
};

// Set the spatial parameters
template<class TypeTag>
struct SpatialParams<TypeTag, TTag::OnePTwoCNITransientBC>
{
    using GridGeometry = GetPropType<TypeTag, Properties::GridGeometry>;
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using type = OnePNCTestSpatialParams<GridGeometry, Scalar>;
};

// Define whether mole(true) or mass (false) fractions are used
template<class TypeTag>
struct UseMoles<TypeTag, TTag::OnePTwoCNITransientBC> { static constexpr bool value = true; };

#ifndef ENABLECACHING
#define ENABLECACHING false
#endif

// Switch on/off caching
template<class TypeTag>
struct EnableGridGeometryCache<TypeTag, TTag::OnePTwoCNITransientBC> { static constexpr bool value = ENABLECACHING; };
template<class TypeTag>
struct EnableGridFluxVariablesCache<TypeTag, TTag::OnePTwoCNITransientBC> { static constexpr bool value = ENABLECACHING; };
template<class TypeTag>
struct EnableGridVolumeVariablesCache<TypeTag, TTag::OnePTwoCNITransientBC> { static constexpr bool value = ENABLECACHING; };
}

/*!
 * \ingroup OnePNCTests
 * \brief Test for the OnePTwoCModel in combination with the NI model for a convection problem.
 *
 * The simulation domain is a tube where water with an elevated temperature is injected
 * at a constant rate on the left hand side.
 *
 * Initially, the domain is fully saturated with water at a constant temperature.
 * On the left hand side water is injected at a constant rate and on the right hand side
 * a Dirichlet boundary with constant pressure, saturation and temperature is applied.
 *
 * The results are compared to an analytical solution where a retarded front velocity is calculated as follows:
  \f[
     v_{Front}=\frac{q S_{water}}{\phi S_{total}}
 \f]
 *
 * The result of the analytical solution is written into the vtu files.
 *
 * This problem uses the \ref OnePModel and \ref NIModel model.
 *
 * To run the simulation execute the following line in shell: <br>
 * <tt>./test_box1p2cniconvection -ParameterFile ./test_box1p2cniconvection.input</tt> or <br>
 * <tt>./test_cc1p2cniconvection -ParameterFile ./test_cc1p2cniconvection.input</tt>
 */

template <class TypeTag>
class OnePTwoCNITransientBCProblem : public PorousMediumFlowProblem<TypeTag>
{
    using ParentType = PorousMediumFlowProblem<TypeTag>;

    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using Indices = typename GetPropType<TypeTag, Properties::ModelTraits>::Indices;
    using GridView = typename GetPropType<TypeTag, Properties::GridGeometry>::GridView;
    using FluidSystem = GetPropType<TypeTag, Properties::FluidSystem>;
    using BoundaryTypes = Dumux::BoundaryTypes<GetPropType<TypeTag, Properties::ModelTraits>::numEq()>;
    using PrimaryVariables = GetPropType<TypeTag, Properties::PrimaryVariables>;
    using GridGeometry = GetPropType<TypeTag, Properties::GridGeometry>;
    using FVElementGeometry = typename GetPropType<TypeTag, Properties::GridGeometry>::LocalView;
    using NumEqVector = GetPropType<TypeTag, Properties::NumEqVector>;
    using ElementVolumeVariables = typename GetPropType<TypeTag, Properties::GridVolumeVariables>::LocalView;
    using SubControlVolumeFace = typename FVElementGeometry::SubControlVolumeFace;
    using Element = typename GridView::template Codim<0>::Entity;
    using SolutionVector = GetPropType<TypeTag, Properties::SolutionVector>;
    using VolumeVariables = GetPropType<TypeTag, Properties::VolumeVariables>;
    using IapwsH2O = Components::H2O<Scalar>;

    // copy some indices for convenience
    enum
    {
        // indices of the primary variables
        pressureIdx = Indices::pressureIdx,
        temperatureIdx = Indices::temperatureIdx,

        // component indices
        H2OIdx = FluidSystem::compIdx(FluidSystem::MultiPhaseFluidSystem::H2OIdx),
        N2Idx = FluidSystem::compIdx(FluidSystem::MultiPhaseFluidSystem::N2Idx),

        // indices of the equations
        contiH2OEqIdx = Indices::conti0EqIdx + H2OIdx,
        contiN2EqIdx = Indices::conti0EqIdx + N2Idx,
        energyEqIdx = Indices::energyEqIdx
    };

    //! Property that defines whether mole or mass fractions are used
    static constexpr bool useMoles = getPropValue<TypeTag, Properties::UseMoles>();
    static const int dimWorld = GridView::dimensionworld;
    using GlobalPosition = typename SubControlVolumeFace::GlobalPosition;

public:
    OnePTwoCNITransientBCProblem(std::shared_ptr<const GridGeometry> gridGeometry)
    : ParentType(gridGeometry)
    {
        //initialize fluid system
        FluidSystem::init();

        // stating in the console whether mole or mass fractions are used
        if(useMoles)
            std::cout<<"problem uses mole fractions"<<std::endl;
        else
            std::cout<<"problem uses mass fractions"<<std::endl;
    }

    /*!
     * \name Problem parameters
     */
    // \{

    /*!
     * \brief Returns the temperature within the domain [K].
     *
     * This problem assumes a temperature of 20 degrees Celsius.
     */
    Scalar temperature() const
    { return 273.15 + 20; } // in [K]

    // \}

    /*!
     * \name Boundary conditions
     */
    // \{

    /*!
     * \brief Specifies which kind of boundary condition should be
     *        used for which equation on a given boundary segment.
     *
     * \param globalPos The position for which the bc type should be evaluated
     */
    BoundaryTypes boundaryTypesAtPos(const GlobalPosition &globalPos) const
    {
        BoundaryTypes values;

        if(globalPos[0] < this->gridGeometry().bBoxMin()[0] + eps_ || globalPos[0] > this->gridGeometry().bBoxMax()[0] - eps_)
            values.setAllDirichlet();
        else
            values.setAllNeumann();

        return values;
    }

    /*!
     * \brief Evaluates the boundary conditions for a Dirichlet boundary segment.
     *
     * \param globalPos The position for which the bc type should be evaluated
     */
    PrimaryVariables dirichletAtPos(const GlobalPosition &globalPos) const
    {
        PrimaryVariables values = initial_(globalPos);

        // make the BCs on the left border time-dependent
        if (globalPos[0] < this->gridGeometry().bBoxMin()[0] + eps_)
        {
            values[pressureIdx] += time_ * 1.0;
            values[N2Idx] += time_ * 1e-8;
            values[temperatureIdx] += time_ * 1e-3;
        }

        return values;
    }

    /*!
     * \brief Evaluates the boundary conditions for a Neumann boundary segment.
     */
    NumEqVector neumannAtPos(const GlobalPosition& globalPos) const
    {
        return NumEqVector(0.0);
    }

    // \}

    /*!
     * \name Volume terms
     */
    // \{

    /*!
     * \brief Evaluates the source term for all phases within a given
     *        sub-control volume.
     *
     * For this method, the \a priVars parameter stores the rate mass
     * of a component is generated or annihilated per volume
     * unit. Positive values mean that mass is created, negative ones
     * mean that it vanishes.
     *
     * The units must be according to either using mole or mass fractions (mole/(m^3*s) or kg/(m^3*s)).
     */
    NumEqVector sourceAtPos(const GlobalPosition &globalPos) const
    { return NumEqVector(0.0); }

    /*!
     * \brief Evaluates the initial value for a control volume.
     *
     * \param globalPos The position for which the initial condition should be evaluated
     *
     * For this method, the \a values parameter stores primary
     * variables.
     */
    PrimaryVariables initialAtPos(const GlobalPosition &globalPos) const
    { return initial_(globalPos); }

    /*!
    * \brief Set the simulation time.
    *
    * \param t The current time.
    */
    void setTime(Scalar t)
    { time_ = t; }

    // \}
private:

    // the internal method for the initial condition
    PrimaryVariables initial_(const GlobalPosition& globalPos) const
    {
        PrimaryVariables priVars;

        if (globalPos[0] < this->gridGeometry().bBoxMin()[0] + eps_)
        {
            priVars[pressureIdx] = 1.1e5; // initial condition for the pressure
            priVars[N2Idx] = 2e-10;  // initial condition for the N2 molefraction
            priVars[temperatureIdx] = 300.00;
        }
        else
        {
            priVars[pressureIdx] = 1.0e5;
            priVars[N2Idx] = 0.0;
            priVars[temperatureIdx] = 285.00;
        }
        return priVars;
    }

        static constexpr Scalar eps_ = 1e-6;
        Scalar time_;
    };

} // end namespace Dumux

#endif
