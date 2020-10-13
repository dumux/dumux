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
 * \ingroup NavierStokesNCTests
 * \brief Channel flow test for the multi-component staggered grid (Navier-)Stokes model.
 */

#ifndef DUMUX_CHANNEL_NC_TEST_PROBLEM_HH
#define DUMUX_CHANNEL_NC_TEST_PROBLEM_HH

#ifndef ENABLECACHING
#define ENABLECACHING 1
#endif

#include <dune/grid/yaspgrid.hh>

#include <dumux/freeflow/navierstokes/problem.hh>

#include <dumux/material/components/simpleh2o.hh>
#include <dumux/material/fluidsystems/1padapter.hh>
#include <dumux/material/fluidsystems/h2oair.hh>

#include <dumux/freeflow/navierstokes/momentum/fluxhelper.hh>
#include <dumux/freeflow/navierstokes/scalarfluxhelper.hh>
#include <dumux/freeflow/navierstokes/momentum/model.hh>
#include <dumux/freeflow/navierstokes/mass/1pnc/model.hh>
#include <dumux/freeflow/navierstokes/problem.hh>
#include <dumux/discretization/fcstaggered.hh>
#include <dumux/discretization/cctpfa.hh>

namespace Dumux {

template <class TypeTag>
class ChannelNCTestProblem;

namespace Properties {

// Create new type tags
namespace TTag {
struct ChannelNCTest {};
struct ChannelNCTestMomentum { using InheritsFrom = std::tuple<ChannelNCTest, NavierStokesMomentum, FaceCenteredStaggeredModel>; };
#if !NONISOTHERMAL
struct ChannelNCTestMass { using InheritsFrom = std::tuple<ChannelNCTest, NavierStokesMassOnePNC, CCTpfaModel>; };
#else
struct ChannelNCTestMass { using InheritsFrom = std::tuple<ChannelNCTest, NavierStokesMassOnePNCNI, CCTpfaModel>; };
#endif
} // end namespace TTag

// Select the fluid system
template<class TypeTag>
struct FluidSystem<TypeTag, TTag::ChannelNCTest>
{
    using H2OAir = FluidSystems::H2OAir<GetPropType<TypeTag, Properties::Scalar>>;
    static constexpr int phaseIdx = H2OAir::liquidPhaseIdx;
    using type = FluidSystems::OnePAdapter<H2OAir, phaseIdx>;
};

template<class TypeTag>
struct ReplaceCompEqIdx<TypeTag, TTag::ChannelNCTest> { static constexpr int value = 0; };

// Set the grid type
template<class TypeTag>
struct Grid<TypeTag, TTag::ChannelNCTest> { using type = Dune::YaspGrid<2>; };

// Set the problem property
template<class TypeTag>
struct Problem<TypeTag, TTag::ChannelNCTest> { using type = Dumux::ChannelNCTestProblem<TypeTag> ; };

template<class TypeTag>
struct EnableGridGeometryCache<TypeTag, TTag::ChannelNCTest> { static constexpr bool value = ENABLECACHING; };
template<class TypeTag>
struct EnableGridFluxVariablesCache<TypeTag, TTag::ChannelNCTest> { static constexpr bool value = ENABLECACHING; };
template<class TypeTag>
struct EnableGridVolumeVariablesCache<TypeTag, TTag::ChannelNCTest> { static constexpr bool value = ENABLECACHING; };
template<class TypeTag>
struct EnableGridFaceVariablesCache<TypeTag, TTag::ChannelNCTest> { static constexpr bool value = ENABLECACHING; };

// Use mole fraction formulation
#if USE_MASS
template<class TypeTag>
struct UseMoles<TypeTag, TTag::ChannelNCTest> { static constexpr bool value = false; };
#else
template<class TypeTag>
struct UseMoles<TypeTag, TTag::ChannelNCTest> { static constexpr bool value = true; };
#endif
} // end namespace Properties

/*!
 * \ingroup NavierStokesNCTests
 * \brief  Test problem for the one-phase (Navier-)Stokes model.
 *
 * Flow from left to right in a channel is considered. A parabolic velocity
 * profile is set at the left boundary, while the pressure is set to
 * a fixed value on the right boundary. The top and bottom boundaries
 * represent solid walls with no-slip/no-flow conditions.
 */
template <class TypeTag>
class ChannelNCTestProblem : public NavierStokesProblem<TypeTag>
{
    using ParentType = NavierStokesProblem<TypeTag>;
    using BoundaryTypes = typename ParentType::BoundaryTypes;
    using GridGeometry = GetPropType<TypeTag, Properties::GridGeometry>;
    using FVElementGeometry = typename GridGeometry::LocalView;
    using SubControlVolume = typename FVElementGeometry::SubControlVolume;
    using SubControlVolumeFace = typename FVElementGeometry::SubControlVolumeFace;
    using Indices = typename GetPropType<TypeTag, Properties::ModelTraits>::Indices;
    using NumEqVector = typename ParentType::NumEqVector;
    using ModelTraits = GetPropType<TypeTag, Properties::ModelTraits>;
    using PrimaryVariables = typename ParentType::PrimaryVariables;
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using SolutionVector = GetPropType<TypeTag, Properties::SolutionVector>;
    using FluidSystem = GetPropType<TypeTag, Properties::FluidSystem>;

    static constexpr auto dimWorld = GridGeometry::GridView::dimensionworld;
    using Element = typename GridGeometry::GridView::template Codim<0>::Entity;
    using GlobalPosition = typename Element::Geometry::GlobalCoordinate;
    using VelocityVector = Dune::FieldVector<Scalar, dimWorld>;

    using CouplingManager = GetPropType<TypeTag, Properties::CouplingManager>;

    using TimeLoopPtr = std::shared_ptr<CheckPointTimeLoop<Scalar>>;

    static constexpr auto compIdx = 1;

public:
    ChannelNCTestProblem(std::shared_ptr<const GridGeometry> gridGeometry, std::shared_ptr<CouplingManager> couplingManager)
    : ParentType(gridGeometry, couplingManager), eps_(1e-6)
    {
        isChannelTest_ = getParam<bool>("Problem.IsChannelTest");

        if (isChannelTest_)
            inletVelocity_ = getParam<Scalar>("Problem.InletVelocity");
    }

   /*!
     * \name Problem parameters
     */
    // \{

   /*!
     * \brief Returns the temperature within the domain in [K].
     *
     * This problem assumes a temperature of 10 degrees Celsius.
     */
    Scalar temperature() const
    { return 273.15 + 10; } // 10C

    /*!
     * \brief Returns a reference pressure at a given sub control volume face.
     *        This pressure is substracted from the actual pressure for the momentum balance
     *        which potentially helps to improve numerical accuracy by avoiding issues related do floating point arithmetic.
     */
    Scalar referencePressure(const Element& element,
                             const FVElementGeometry& fvGeometry,
                             const SubControlVolumeFace& scvf) const
    { return 1.1e5; }

    // \}
   /*!
     * \name Boundary conditions
     */
    // \{

   /*!
     * \brief Specifies which kind of boundary condition should be
     *        used for which equation on a given boundary control volume.
     *
     * \param globalPos The position of the center of the finite volume
     */
    BoundaryTypes boundaryTypesAtPos(const GlobalPosition& globalPos) const
    {
        if (isChannelTest_)
            return boundaryTypesForChannelTest_(globalPos);
        else
            return boundaryTypesForDensityDrivenFlowTest_(globalPos);
    }

   /*!
     * \brief Evaluates the boundary conditions for a Dirichlet control volume.
     *
     * \param globalPos The center of the finite volume which ought to be set.
     */
    PrimaryVariables dirichlet(const Element& element, const SubControlVolumeFace& scvf) const
    {
        if (isChannelTest_)
            return dirichletForChannelTest_(element, scvf);
        else
            return dirichletForDensityDrivenFlowTest_(element, scvf);
    }

    /*!
     * \brief Evaluates the boundary conditions for a Neumann control volume.
     *
     * \param element The element for which the Neumann boundary condition is set
     * \param fvGeometry The fvGeometry
     * \param elemVolVars The element volume variables
     * \param elemFaceVars The element face variables
     * \param scvf The boundary sub control volume face
     */
    template<class ElementVolumeVariables, class ElementFluxVariablesCache>
    NumEqVector neumann(const Element& element,
                        const FVElementGeometry& fvGeometry,
                        const ElementVolumeVariables& elemVolVars,
                        const ElementFluxVariablesCache& elemFluxVarsCache,
                        const SubControlVolumeFace& scvf) const
    {
        if (isChannelTest_)
            return neumannForChannelTest_(element, fvGeometry, elemVolVars, elemFluxVarsCache, scvf);
        else
            return neumannForDensityDrivenFlowTest_(element, fvGeometry, elemVolVars, elemFluxVarsCache, scvf);
    }

    // \}

   /*!
     * \name Volume terms
     */
    // \{

   /*!
     * \brief Evaluates the initial value for a control volume.
     *
     * \param globalPos The global position
     */
    PrimaryVariables initialAtPos(const GlobalPosition& globalPos) const
    {
        if (isChannelTest_)
            return initialForChannelTest_(globalPos);
        else
            return initialForDensityDrivenFlowTest_(globalPos);
    }

    // \}

    void setTimeLoop(TimeLoopPtr timeLoop)
    {
        timeLoop_ = timeLoop;
    }

    Scalar time() const
    {
        return timeLoop_->time();
    }

    //! Enable internal Dirichlet constraints
    static constexpr bool enableInternalDirichletConstraints()
    { return !ParentType::isMomentumProblem(); }

    /*!
     * \brief Tag a degree of freedom to carry internal Dirichlet constraints.
     *        If true is returned for a dof, the equation for this dof is replaced
     *        by the constraint that its primary variable values must match the
     *        user-defined values obtained from the function internalDirichlet(),
     *        which must be defined in the problem.
     *
     * \param element The finite element
     * \param scv The sub-control volume
     */
    std::bitset<PrimaryVariables::dimension> hasInternalDirichletConstraint(const Element& element, const SubControlVolume& scv) const
    {
        std::bitset<PrimaryVariables::dimension> values;

        if (isChannelTest_)
            return values;

        // the pure Neumann problem is only defined up to a constant
        // we create a well-posed problem by fixing the pressure at one dof
        if (scv.dofIndex() == 0)
            values.set(0);

        return values;
    }

    /*!
     * \brief Define the values of internal Dirichlet constraints for a degree of freedom.
     * \param element The finite element
     * \param scv The sub-control volume
     */
    PrimaryVariables internalDirichlet(const Element& element, const SubControlVolume& scv) const
    { return PrimaryVariables(1.1e5); }

private:

    bool isChannelTest_;

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //////////// Problem setup for channel test ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    PrimaryVariables initialForChannelTest_(const GlobalPosition &globalPos) const
    {
        PrimaryVariables values;

        if constexpr (ParentType::isMomentumProblem())
        {
            // parabolic velocity profile
            values[Indices::velocityXIdx] =  inletVelocity_*(globalPos[1] - this->gridGeometry().bBoxMin()[1])*(this->gridGeometry().bBoxMax()[1] - globalPos[1])
            / (0.25*(this->gridGeometry().bBoxMax()[1] - this->gridGeometry().bBoxMin()[1])*(this->gridGeometry().bBoxMax()[1] - this->gridGeometry().bBoxMin()[1]));

            values[Indices::velocityYIdx] = 0.0;
        }
        else
        {
            values[Indices::pressureIdx] = 1.1e+5;
            values[Indices::conti0EqIdx+compIdx] = 0.0;
#if NONISOTHERMAL
            values[Indices::temperatureIdx] = 283.15;
#endif
        }

        return values;
    }

    BoundaryTypes boundaryTypesForChannelTest_(const GlobalPosition &globalPos) const
    {
        BoundaryTypes values;

        if constexpr (ParentType::isMomentumProblem())
        {
            values.setDirichlet(Indices::velocityXIdx);
            values.setDirichlet(Indices::velocityYIdx);

            if (isOutlet_(globalPos))
                values.setAllNeumann();
        }
        else
        {
            values.setAllNeumann();

            if (isInlet_(globalPos))
                values.setAllDirichlet();
        }

        return values;
    }

    PrimaryVariables dirichletForChannelTest_(const Element& element, const SubControlVolumeFace& scvf) const
    {
        const auto& globalPos = scvf.ipGlobal();
        PrimaryVariables values = initialAtPos(globalPos);

        if constexpr (!ParentType::isMomentumProblem())
        {
            // give the system some time so that the pressure can equilibrate, then start the injection of the tracer
            if (isInlet_(globalPos))
            {
                values[Indices::pressureIdx] = this->couplingManager().cellPressure(element, scvf);

                if (time() >= 10.0 || inletVelocity_  < eps_)
                {
                    Scalar moleFracTransportedComp = 1e-3;
#if USE_MASS
                    Scalar averageMolarMassPhase = moleFracTransportedComp * FluidSystem::molarMass(compIdx)
                                                   + (1. - moleFracTransportedComp)  * FluidSystem::molarMass(1-compIdx);
                    values[Indices::conti0EqIdx+compIdx] = moleFracTransportedComp * FluidSystem::molarMass(compIdx)
                                               / averageMolarMassPhase;
#else
                    values[Indices::conti0EqIdx+compIdx] = moleFracTransportedComp;
#endif
#if NONISOTHERMAL
                values[Indices::temperatureIdx] = 293.15;
#endif
                }
            }
        }

        return values;
    }

    template<class ElementVolumeVariables, class ElementFluxVariablesCache>
    NumEqVector neumannForChannelTest_(const Element& element,
                                       const FVElementGeometry& fvGeometry,
                                       const ElementVolumeVariables& elemVolVars,
                                       const ElementFluxVariablesCache& elemFluxVarsCache,
                                       const SubControlVolumeFace& scvf) const
    {
        NumEqVector values(0.0);

        if constexpr (ParentType::isMomentumProblem())
        {
            using FluxHelper = NavierStokesMomentumBoundaryFluxHelper;
            values = FluxHelper::fixedPressureMomentumFlux(*this, element, fvGeometry, scvf, elemVolVars, elemFluxVarsCache, 1.1e+5, true /*zeroNormalVelocityGradient*/);
        }
        else
        {
            using FluxHelper = NavierStokesScalarBoundaryFluxHelper<AdvectiveFlux<ModelTraits>>;
            if (isOutlet_(scvf.ipGlobal()))
                values = FluxHelper::scalarOutflowFlux(*this, element, fvGeometry, scvf, elemVolVars);
        }

        return values;
    }

    bool isInlet_(const GlobalPosition& globalPos) const
    {
        return globalPos[0] < eps_;
    }

    bool isOutlet_(const GlobalPosition& globalPos) const
    {
        return globalPos[0] > this->gridGeometry().bBoxMax()[0] - eps_;
    }

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //////////// Problem setup for density driven flow test ///////////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    PrimaryVariables initialForDensityDrivenFlowTest_(const GlobalPosition &globalPos) const
    {
        PrimaryVariables values(0.0);

        if constexpr (!ParentType::isMomentumProblem())
            values[Indices::pressureIdx] = 1.1e5;
        return values;
    }

    BoundaryTypes boundaryTypesForDensityDrivenFlowTest_(const GlobalPosition& globalPos) const
    {
        BoundaryTypes values;

        if constexpr (ParentType::isMomentumProblem())
            values.setAllDirichlet();
        else
        {
            values.setAllNeumann();
            if (globalPos[1] > this->gridGeometry().bBoxMax()[1] - eps_)
            {
                static const bool useWholeLength = getParam<bool>("Problem.UseWholeLength");
                if (useWholeLength || (globalPos[0] > 0.4 && globalPos[0] < 0.6))
                    values.setAllDirichlet();
            }

        }

        return values;
    }

    PrimaryVariables dirichletForDensityDrivenFlowTest_(const Element& element, const SubControlVolumeFace& scvf) const
    {
        PrimaryVariables values(0.0);

        if constexpr (!ParentType::isMomentumProblem())
        {
            values[Indices::pressureIdx] = 1.1e+5;
            values[Indices::conti0EqIdx+compIdx] = 1e-3;
        }

        return values;
    }

    template<class ElementVolumeVariables, class ElementFluxVariablesCache>
    NumEqVector neumannForDensityDrivenFlowTest_(const Element& element,
                                                 const FVElementGeometry& fvGeometry,
                                                 const ElementVolumeVariables& elemVolVars,
                                                 const ElementFluxVariablesCache& elemFluxVarsCache,
                                                 const SubControlVolumeFace& scvf) const
    {
        // no flow everywhere
        return NumEqVector(0.0);
    }

    const Scalar eps_;
    Scalar inletVelocity_;
    TimeLoopPtr timeLoop_;
};
} // end namespace Dumux

#endif
