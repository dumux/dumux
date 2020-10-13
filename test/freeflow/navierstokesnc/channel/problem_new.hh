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

#include <dumux/freeflow/navierstokes/fluxhelper.hh>
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
// #if !NONISOTHERMAL
struct ChannelNCTestMass { using InheritsFrom = std::tuple<ChannelNCTest, NavierStokesMassOnePNC, CCTpfaModel>; };
// #else
// struct ChannelNCTest { using InheritsFrom = std::tuple<NavierStokesNCNI, StaggeredFreeFlowModel>; };
// #endif
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
        inletVelocity_ = getParam<Scalar>("Problem.InletVelocity");
        FluidSystem::init();

        if constexpr (!ParentType::isMomentumProblem())
            deltaP_.resize(this->gridGeometry().numDofs());
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
    BoundaryTypes boundaryTypesAtPos(const GlobalPosition &globalPos) const
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

   /*!
     * \brief Evaluates the boundary conditions for a Dirichlet control volume.
     *
     * \param globalPos The center of the finite volume which ought to be set.
     */
    PrimaryVariables dirichlet(const Element& element, const SubControlVolumeFace& scvf) const
    {
        const auto& globalPos = scvf.ipGlobal();
        PrimaryVariables values = initialAtPos(globalPos);

        if constexpr (!ParentType::isMomentumProblem())
        {
            // give the system some time so that the pressure can equilibrate, then start the injection of the tracer
            if (isInlet_(globalPos))
            {
                values[Indices::pressureIdx] = values[Indices::pressureIdx] = this->couplingManager().cellPressure(element, scvf);

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
        NumEqVector values(0.0);

        if constexpr (ParentType::isMomentumProblem())
        {
            values = NavierStokesBoundaryFluxHelper<ModelTraits>::fixedPressureMomentumFlux(*this, element, fvGeometry, scvf, elemVolVars, elemFluxVarsCache, 1.1e+5, true /*zeroNormalVelocityGradient*/);
        }
        else
        {
            if (isOutlet_(scvf.ipGlobal()))
                values = NavierStokesBoundaryFluxHelper<ModelTraits>::scalarOutflowFlux(*this, element, fvGeometry, scvf, elemVolVars);
        }

        return values;
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
    PrimaryVariables initialAtPos(const GlobalPosition &globalPos) const
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

   /*!
     * \brief Adds additional VTK output data to the VTKWriter.
     *
     * Function is called by the output module on every write.
     *
     * \param gridVariables The grid variables
     * \param sol The solution vector
     */
    template<class GridVariables, class SolutionVector>
    void calculateDeltaP(const GridVariables& gridVariables, const SolutionVector& sol)
    {
        for (const auto& element : elements(this->gridGeometry().gridView()))
        {
            auto fvGeometry = localView(this->gridGeometry());
            fvGeometry.bindElement(element);
            for (auto&& scv : scvs(fvGeometry))
            {
                auto ccDofIdx = scv.dofIndex();

                auto elemVolVars = localView(gridVariables.curGridVolVars());
                elemVolVars.bindElement(element, fvGeometry, sol);

                deltaP_[ccDofIdx] = elemVolVars[scv].pressure() - 1.1e5;
            }
        }
    }

    auto& getDeltaP() const
    { return deltaP_; }

    // \}

    void setTimeLoop(TimeLoopPtr timeLoop)
    {
        timeLoop_ = timeLoop;
    }

    Scalar time() const
    {
        return timeLoop_->time();
    }

private:

    bool isInlet_(const GlobalPosition& globalPos) const
    {
        return globalPos[0] < eps_;
    }

    bool isOutlet_(const GlobalPosition& globalPos) const
    {
        return globalPos[0] > this->gridGeometry().bBoxMax()[0] - eps_;
    }

    const Scalar eps_;
    Scalar inletVelocity_;
    TimeLoopPtr timeLoop_;
    std::vector<Scalar> deltaP_;
};
} // end namespace Dumux

#endif
