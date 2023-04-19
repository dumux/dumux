// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup NavierStokesTests
 * \brief Channel flow test for the staggered grid (Navier-)Stokes model.
 */

#ifndef DUMUX_CHANNEL_TEST_PROBLEM_HH
#define DUMUX_CHANNEL_TEST_PROBLEM_HH

#include <dumux/common/parameters.hh>
#include <dumux/common/properties.hh>
#include <dumux/common/timeloop.hh>

#include <dumux/freeflow/navierstokes/boundarytypes.hh>

#include <dumux/freeflow/navierstokes/momentum/fluxhelper.hh>
#include <dumux/freeflow/navierstokes/scalarfluxhelper.hh>
#include <dumux/freeflow/navierstokes/mass/1p/advectiveflux.hh>

namespace Dumux {

/*!
 * \ingroup NavierStokesTests
 * \brief  Test problem for the one-phase (Navier-) Stokes problem in a channel.
 *
 * Flow from left to right in a two-dimensional channel is considered. At the inlet (left),
 * fixed values for velocity are set, while at the outlet (right), a fixed pressure
 * boundary condition is used. The channel is confined by solid walls at the top and bottom
 * of the domain which corresponds to no-slip/no-flow conditions.
 * For the non-isothermal test, water of increased temperature is injected at the inlet
 * while the walls are fully isolating.
 */
template <class TypeTag, class BaseProblem>
class ChannelTestProblem : public BaseProblem
{
    using ParentType = BaseProblem;

    using BoundaryTypes = typename ParentType::BoundaryTypes;
    using GridGeometry = GetPropType<TypeTag, Properties::GridGeometry>;
    using FVElementGeometry = typename GridGeometry::LocalView;
    using SubControlVolumeFace = typename GridGeometry::SubControlVolumeFace;
    using InitialValues = typename ParentType::InitialValues;
    using Sources = typename ParentType::Sources;
    using DirichletValues = typename ParentType::DirichletValues;
    using BoundaryFluxes = typename ParentType::BoundaryFluxes;
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using ModelTraits = GetPropType<TypeTag, Properties::ModelTraits>;
    using SolutionVector = GetPropType<TypeTag, Properties::SolutionVector>;

    static constexpr auto dimWorld = GridGeometry::GridView::dimensionworld;
    using Element = typename FVElementGeometry::Element;
    using GlobalPosition = typename Element::Geometry::GlobalCoordinate;

    using CouplingManager = GetPropType<TypeTag, Properties::CouplingManager>;

    // the types of outlet boundary conditions
    enum class OutletCondition
    {
        outflow, unconstrainedOutflow, neumannXdirichletY, neumannXneumannY
    };

public:
    using Indices = typename GetPropType<TypeTag, Properties::ModelTraits>::Indices;

    ChannelTestProblem(std::shared_ptr<const GridGeometry> gridGeometry, std::shared_ptr<CouplingManager> couplingManager)
    : ParentType(gridGeometry, couplingManager)
    {
        inletVelocity_ = getParam<Scalar>("Problem.InletVelocity");
        dynamicViscosity_ = getParam<Scalar>("Component.LiquidKinematicViscosity", 1.0)
                            * getParam<Scalar>("Component.LiquidDensity", 1.0);

        const auto outletBC = getParam<std::string>("Problem.OutletCondition", "Outflow");
        if (outletBC == "Outflow")
            outletCondition_ = OutletCondition::outflow;
        else if (outletBC == "UnconstrainedOutflow")
            outletCondition_ = OutletCondition::unconstrainedOutflow;
        else if (outletBC == "NeumannX_DirichletY")
            outletCondition_ = OutletCondition::neumannXdirichletY;
        else if (outletBC == "NeumannX_NeumannY")
            outletCondition_ = OutletCondition::neumannXneumannY;
        else
            DUNE_THROW(Dune::InvalidStateException, outletBC + " is not a valid outlet boundary condition");

        useVelocityProfile_ = getParam<bool>("Problem.UseVelocityProfile", false);
        outletPressure_ = getParam<Scalar>("Problem.OutletPressure", 1.1e5);
    }

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
        BoundaryTypes values;

        if constexpr (ParentType::isMomentumProblem())
        {
            values.setDirichlet(Indices::velocityXIdx);
            values.setDirichlet(Indices::velocityYIdx);

            if (isOutlet_(globalPos))
            {
                if (outletCondition_ == OutletCondition::neumannXdirichletY)
                {
                    values.setNeumann(Indices::momentumXBalanceIdx);
                    values.setDirichlet(Indices::velocityYIdx);
                }
                else
                    values.setAllNeumann();
            }
        }
        else
        {
            values.setAllNeumann();

            if (isInlet_(globalPos))
            {
                values.setDirichlet(Indices::pressureIdx);
#if NONISOTHERMAL
                values.setDirichlet(Indices::temperatureIdx);
#endif
            }
        }

        return values;
    }


    /*!
     * \brief Evaluates the boundary conditions for a Dirichlet control volume.
     *
     * \param globalPos The center of the finite volume which ought to be set.
     */
    DirichletValues dirichlet(const Element& element, const SubControlVolumeFace& scvf) const
    {
        const auto& globalPos = scvf.ipGlobal();
        DirichletValues values = initialAtPos(globalPos);

        if constexpr (ParentType::isMomentumProblem())
        {
            values[Indices::velocityXIdx] = parabolicProfile(globalPos[1], inletVelocity_);
            if (!useVelocityProfile_ && isInlet_(globalPos))
            values[Indices::velocityXIdx] = inletVelocity_;
        }
        else
        {
            values[Indices::pressureIdx] = this->couplingManager().cellPressure(element, scvf);

#if NONISOTHERMAL
        // give the system some time so that the pressure can equilibrate, then start the injection of the hot liquid
        if (time() >= 200.0)
            values[Indices::temperatureIdx] = 293.15;
#endif
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
    BoundaryFluxes neumann(const Element& element,
                           const FVElementGeometry& fvGeometry,
                           const ElementVolumeVariables& elemVolVars,
                           const ElementFluxVariablesCache& elemFluxVarsCache,
                           const SubControlVolumeFace& scvf) const
    {
        BoundaryFluxes values(0.0);

        if constexpr (ParentType::isMomentumProblem())
        {
            using FluxHelper = NavierStokesMomentumBoundaryFluxHelper;

            if (outletCondition_ == OutletCondition::unconstrainedOutflow)
                values = FluxHelper::fixedPressureMomentumFlux(*this, fvGeometry, scvf, elemVolVars, elemFluxVarsCache, outletPressure_, false /*zeroNormalVelocityGradient*/);
            else if (outletCondition_ == OutletCondition::outflow)
                values = FluxHelper::fixedPressureMomentumFlux(*this, fvGeometry, scvf, elemVolVars, elemFluxVarsCache, outletPressure_, true /*zeroNormalVelocityGradient*/);
            else
            {
                assert(outletCondition_ == OutletCondition::neumannXneumannY || outletCondition_ == OutletCondition::neumannXdirichletY);
                values = FluxHelper::fixedPressureMomentumFlux(*this, fvGeometry, scvf, elemVolVars, elemFluxVarsCache, outletPressure_, false /*zeroNormalVelocityGradient*/);

                Dune::FieldMatrix<Scalar, dimWorld, dimWorld> shearRate(0.0); // gradV + gradV^T
                shearRate[0][1] = dudy(scvf.ipGlobal()[1], inletVelocity_); // we assume du/dx = dv/dy = dv/dx = 0
                shearRate[1][0] = shearRate[0][1];

                const auto normal = scvf.unitOuterNormal();
                BoundaryFluxes normalGradient(0.0);
                shearRate.mv(normal, normalGradient);
                values -= this->effectiveViscosity(element, fvGeometry, scvf)*normalGradient;
            }
        }
        else
        {
            using FluxHelper = NavierStokesScalarBoundaryFluxHelper<AdvectiveFlux<ModelTraits>>;

            if (isOutlet_(scvf.ipGlobal()))
                values = FluxHelper::scalarOutflowFlux(*this, element, fvGeometry, scvf, elemVolVars);
        }

        return values;
    }

    /*!
     * \brief A parabolic velocity profile.
     *
     * \param y The position where the velocity is evaluated.
     * \param vMax The profile's maximum velocity.
     */
    Scalar parabolicProfile(const Scalar y, const Scalar vMax) const
    {
        const Scalar yMin = this->gridGeometry().bBoxMin()[1];
        const Scalar yMax = this->gridGeometry().bBoxMax()[1];
        return  vMax * (y - yMin)*(yMax - y) / (0.25*(yMax - yMin)*(yMax - yMin));
    }

    /*!
     * \brief The partial derivative of the horizontal velocity (following a parabolic profile for
     *         Stokes flow) w.r.t. to the y-coordinate (du/dy).
     *
     * \param y The position where the derivative is evaluated.
     * \param vMax The profile's maximum velocity.
     */
    Scalar dudy(const Scalar y, const Scalar vMax) const
    {
        const Scalar yMin = this->gridGeometry().bBoxMin()[1];
        const Scalar yMax = this->gridGeometry().bBoxMax()[1];
        return vMax * (4.0*yMin + 4*yMax - 8.0*y) / ((yMin-yMax)*(yMin-yMax));
    }

    /*!
     * \brief Return the analytical solution of the problem at a given position
     *
     * \param globalPos The global position
     */
    DirichletValues analyticalSolution(const GlobalPosition& globalPos, Scalar time = 0.0) const
    {
        DirichletValues values;

        if constexpr (ParentType::isMomentumProblem())
        {
            values[Indices::velocityXIdx] = parabolicProfile(globalPos[1], inletVelocity_);
            values[Indices::velocityYIdx] = 0.0;
        }
        else
        {
            const Scalar yMin = this->gridGeometry().bBoxMin()[1];
            const Scalar yMax = this->gridGeometry().bBoxMax()[1];
            const Scalar velocityQuadraticCoefficient = - inletVelocity_ / (0.25*(yMax - yMin)*(yMax - yMin));
            const Scalar pressureLinearCoefficient = 2.0 * velocityQuadraticCoefficient * dynamicViscosity_;
            const Scalar pressureConstant = -pressureLinearCoefficient * this->gridGeometry().bBoxMax()[0]  + outletPressure_;
            values[Indices::pressureIdx] =  pressureLinearCoefficient * globalPos[0] + pressureConstant;
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
    InitialValues initialAtPos(const GlobalPosition& globalPos) const
    {
        InitialValues values;

        if constexpr (ParentType::isMomentumProblem())
        {
            values[Indices::velocityYIdx] = 0.0;
            if (useVelocityProfile_)
                values[Indices::velocityXIdx] = parabolicProfile(globalPos[1], inletVelocity_);
            else
                values[Indices::velocityXIdx] = inletVelocity_;
        }
        else
        {
            values[Indices::pressureIdx] = outletPressure_;
#if NONISOTHERMAL
            values[Indices::temperatureIdx] = 283.15;
#endif
        }



        return values;
    }

    // \}

    /*!
     * \brief Returns a reference pressure at a given sub control volume face.
     *        This pressure is subtracted from the actual pressure for the momentum balance
     *        which potentially helps to improve numerical accuracy by avoiding issues related do floating point arithmetic.
     */
    Scalar referencePressure(const Element& element,
                             const FVElementGeometry& fvGeometry,
                             const SubControlVolumeFace& scvf) const
    { return outletPressure_; }

    void setTime(Scalar time)
    {
        time_ = time;
    }

    Scalar time() const
    {
        return time_;
    }

    bool hasAnalyticalSolution() const
    { return outletCondition_ == OutletCondition::neumannXneumannY; }

private:
    bool isInlet_(const GlobalPosition& globalPos) const
    {
        return globalPos[0] < eps_;
    }

    bool isOutlet_(const GlobalPosition& globalPos) const
    {
        return globalPos[0] > this->gridGeometry().bBoxMax()[0] - eps_;
    }

    static constexpr Scalar eps_=1e-6;
    Scalar inletVelocity_;
    Scalar dynamicViscosity_;
    Scalar outletPressure_;
    OutletCondition outletCondition_;
    bool useVelocityProfile_;
    Scalar time_;
};
} // end namespace Dumux

#endif
