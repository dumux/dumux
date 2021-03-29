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
 * \ingroup NavierStokesTests
 * \brief Channel flow test for the staggered grid (Navier-)Stokes model.
 */

#ifndef DUMUX_CHANNEL_TEST_PROBLEM_HH
#define DUMUX_CHANNEL_TEST_PROBLEM_HH

#include <dune/grid/yaspgrid.hh>



#include <dumux/material/fluidsystems/1pliquid.hh>
#include <dumux/material/components/constant.hh>
#include <dumux/material/components/simpleh2o.hh>

#include <dumux/freeflow/navierstokes/fluxhelper.hh>
#include <dumux/freeflow/navierstokes/momentum/model.hh>
#include <dumux/freeflow/navierstokes/mass/1p/model.hh>
#include <dumux/freeflow/navierstokes/problem.hh>
#include <dumux/discretization/fcstaggered.hh>
#include <dumux/discretization/cctpfa.hh>

namespace Dumux {
template <class TypeTag>
class ChannelTestProblem;

namespace Properties {
// Create new type tags
namespace TTag {
struct ChannelTest {};
struct ChannelTestMomentum { using InheritsFrom = std::tuple<ChannelTest, NavierStokesMomentum, FaceCenteredStaggeredModel>; };
#if !NONISOTHERMAL
struct ChannelTestMass { using InheritsFrom = std::tuple<ChannelTest, NavierStokesMassOneP, CCTpfaModel>; };
#else
struct ChannelTestMass { using InheritsFrom = std::tuple<ChannelTest, NavierStokesMassOnePNI, CCTpfaModel>; };
#endif
} // end namespace TTag

// Set the problem property
template<class TypeTag>
struct Problem<TypeTag, TTag::ChannelTest>
{
    using type = Dumux::ChannelTestProblem<TypeTag> ;
};

// the fluid system
template<class TypeTag>
struct FluidSystem<TypeTag, TTag::ChannelTest>
{
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
#if NONISOTHERMAL
    using type = FluidSystems::OnePLiquid<Scalar, Components::SimpleH2O<Scalar> >;
#else
    using type = FluidSystems::OnePLiquid<Scalar, Components::Constant<1, Scalar> >;
#endif
};

// Set the grid type
template<class TypeTag>
struct Grid<TypeTag, TTag::ChannelTest> { using type = Dune::YaspGrid<2>; };


template<class TypeTag>
struct EnableGridGeometryCache<TypeTag, TTag::ChannelTest> { static constexpr bool value = true; };

template<class TypeTag>
struct EnableGridFluxVariablesCache<TypeTag, TTag::ChannelTest> { static constexpr bool value = true; };
template<class TypeTag>
struct EnableGridVolumeVariablesCache<TypeTag, TTag::ChannelTest> { static constexpr bool value = true; };
} // end namespace Properties


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
template <class TypeTag>
class ChannelTestProblem : public NavierStokesProblem<TypeTag>
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

    static constexpr auto dimWorld = GridGeometry::GridView::dimensionworld;
    using Element = typename GridGeometry::GridView::template Codim<0>::Entity;
    using GlobalPosition = typename Element::Geometry::GlobalCoordinate;
    using VelocityVector = Dune::FieldVector<Scalar, dimWorld>;

    using CouplingManager = GetPropType<TypeTag, Properties::CouplingManager>;
    using TimeLoopPtr = std::shared_ptr<CheckPointTimeLoop<Scalar>>;

    // the types of outlet boundary conditions
    enum class OutletCondition
    {
        outflow, doNothing, neumannXdirichletY, neumannXneumannY
    };

public:
    ChannelTestProblem(std::shared_ptr<const GridGeometry> gridGeometry, std::shared_ptr<CouplingManager> couplingManager)
    : ParentType(gridGeometry, couplingManager)
    {
        inletVelocity_ = getParam<Scalar>("Problem.InletVelocity");
        const auto tmp = getParam<std::string>("Problem.OutletCondition", "Outflow");
        if (tmp == "Outflow")
            outletCondition_ = OutletCondition::outflow;
        else if (tmp == "DoNothing")
            outletCondition_ = OutletCondition::doNothing;
        else if (tmp == "NeumannX_DirichletY")
            outletCondition_ = OutletCondition::neumannXdirichletY;
        else if (tmp == "NeumannX_NeumannY")
            outletCondition_ = OutletCondition::neumannXneumannY;
        else
            DUNE_THROW(Dune::InvalidStateException, tmp + " is not a valid outlet boundary condition");

        useVelocityProfile_ = getParam<bool>("Problem.UseVelocityProfile", false);
        outletPressure_ = getParam<Scalar>("Problem.OutletPressure", 1.1e5);
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
                values.setDirichlet(Indices::phi1Idx);
                values.setDirichlet(Indices::phi2Idx);
                values.setDirichlet(Indices::phi3Idx);
                values.setDirichlet(Indices::u1Idx);
                values.setDirichlet(Indices::u2Idx);
                values.setDirichlet(Indices::u3Idx);
            }
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

         const static auto bcType = getParam<std::string>("Phasefield.Scenario");
         const static auto delayTime = getParam<Scalar>("Problem.InletDelay");
         const static auto windupTime = getParam<Scalar>("Problem.InletWindup");
         if constexpr (ParentType::isMomentumProblem())
         {
             const static Scalar rad = (bcType == "ThinChannel") ?
                 getParam<Scalar>("Phasefield.DirichletRadius") : 1.0;
             values[Indices::velocityXIdx] = parabolicProfile(globalPos[1], inletVelocity_, rad);
             if (!useVelocityProfile_ && isInlet_(globalPos))
                values[Indices::velocityXIdx] = inletVelocity_;
             if (time() < delayTime)
             {
                 values[Indices::velocityXIdx] = 0.0;
             }
             else if (time() < delayTime + windupTime)
             {
                 values[Indices::velocityXIdx] *= (time() - delayTime)/windupTime;
             }
         }
         else
         {
             values[Indices::pressureIdx] = this->couplingManager().cellPressure(element, scvf);

#if NONISOTHERMAL
            // give the system some time so that the pressure can equilibrate, then start the injection of the hot liquid
            if (time() >= 200.0)
                values[Indices::temperatureIdx] = 293.15;
#endif
            const static Scalar xi = getParam<Scalar>("Phasefield.xi");
            const static Scalar rad = getParam<Scalar>("Phasefield.DirichletRadius");
            const static Scalar S = getParam<Scalar>("Phasefield.DirichletScaling");
            const static Scalar a_in = getParam<Scalar>("Phasefield.InletAConcentration");
            const static Scalar b_in = getParam<Scalar>("Phasefield.InletBConcentration");
            const static Scalar c_in = getParam<Scalar>("Phasefield.InletCConcentration");
            if (bcType == "ThinChannel")
            {
                const Scalar s = (globalPos[1]-1.0)*(globalPos[1]-1.0)-rad*rad;
                values[Indices::phi1Idx] = 1.0/(1.0 + std::exp(S*s/xi));
                values[Indices::phi2Idx] = 1.0 - values[Indices::phi1Idx];
                values[Indices::phi3Idx] = 0.0;
            }
            else
            {
                values[Indices::phi1Idx] = 1.0;
                values[Indices::phi2Idx] = 0.0;
                values[Indices::phi3Idx] = 0.0;
            }
            //values[Indices::u1Idx] = a_in;
            //values[Indices::u2Idx] = b_in;
            //values[Indices::u3Idx] = c_in;
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
            if (outletCondition_ == OutletCondition::doNothing)
                values = NavierStokesBoundaryFluxHelper<ModelTraits>::fixedPressureMomentumFlux(*this, element, fvGeometry, scvf, elemVolVars, elemFluxVarsCache, outletPressure_, false /*zeroNormalVelocityGradient*/);
            else if (outletCondition_ == OutletCondition::outflow)
                values = NavierStokesBoundaryFluxHelper<ModelTraits>::fixedPressureMomentumFlux(*this, element, fvGeometry, scvf, elemVolVars, elemFluxVarsCache, outletPressure_, true /*zeroNormalVelocityGradient*/);
            else
            {
                assert(outletCondition_ == OutletCondition::neumannXneumannY);
                values = NavierStokesBoundaryFluxHelper<ModelTraits>::fixedPressureMomentumFlux(*this, element, fvGeometry, scvf, elemVolVars, elemFluxVarsCache, outletPressure_, false /*zeroNormalVelocityGradient*/);
                values[Indices::momentumYBalanceIdx] = -dudy(scvf.ipGlobal()[1], inletVelocity_) * this->effectiveViscosity(element, fvGeometry, scvf) * scvf.directionSign();
            }
        }
        else
        {
            if (isOutlet_(scvf.ipGlobal()))
                values = NavierStokesBoundaryFluxHelper<ModelTraits>::scalarOutflowFlux(*this, element, fvGeometry, scvf, elemVolVars);
            if (isInlet_(scvf.ipGlobal()))
            {
                const static Scalar a_in = getParam<Scalar>("Phasefield.InletAConcentration");
                const static Scalar b_in = getParam<Scalar>("Phasefield.InletBConcentration");
                const static Scalar c_in = getParam<Scalar>("Phasefield.InletCConcentration");
                const auto volumeFlux = this->faceVelocity(element,fvGeometry, scvf) *
                    scvf.unitOuterNormal();
                const auto insideVolVars = elemVolVars[scvf.insideScvIdx()];
                values[Indices::u1Idx] = volumeFlux * a_in;
                values[Indices::u2Idx] = volumeFlux * b_in;
                values[Indices::u3Idx] = volumeFlux * c_in;
            }
        }

        return values;
    }

    /*!
     * \brief A parabolic velocity profile.
     *
     * \param y The position where the velocity is evaluated.
     * \param vMax The profile's maxmium velocity.
     */
    Scalar parabolicProfile(const Scalar y, const Scalar vMax, const Scalar radius = 1.0) const
    {
        const Scalar bMin = this->gridGeometry().bBoxMin()[1];
        const Scalar bMax = this->gridGeometry().bBoxMax()[1];
        const Scalar yMin = bMin/2.0 * (1.0 + radius) + bMax/2.0 * (1.0 - radius);
        const Scalar yMax = bMax/2.0 * (1.0 + radius) + bMin/2.0 * (1.0 - radius);
        if (y < yMin || y > yMax)
            return 0.0;
        return  vMax * (y - yMin)*(yMax - y) / (0.25*(yMax - yMin)*(yMax - yMin));
    }

    /*!
     * \brief The partial dervivative of the horizontal velocity (following a parabolic profile for
     *         Stokes flow) w.r.t. to the y-coordinate (du/dy).
     *
     * \param y The position where the derivative is evaluated.
     * \param vMax The profile's maxmium velocity.
     */
    Scalar dudy(const Scalar y, const Scalar vMax) const
    {
        const Scalar yMin = this->gridGeometry().bBoxMin()[1];
        const Scalar yMax = this->gridGeometry().bBoxMax()[1];
        return vMax * (4.0*yMin + 4*yMax - 8.0*y) / ((yMin-yMax)*(yMin-yMax));
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
        PrimaryVariables values;
        const static auto icType = getParam<std::string>("Phasefield.Scenario");

        if constexpr (ParentType::isMomentumProblem())
        {
            values[Indices::velocityYIdx] = 0.0;
            if (icType == "Grain")
                return values;
            if (useVelocityProfile_)
            {
                const static Scalar rad = (icType == "ThinChannel") ?
                    getParam<Scalar>("Phasefield.DirichletRadius") : 1.0;
                values[Indices::velocityXIdx] = parabolicProfile(globalPos[1], inletVelocity_, rad);
            }
            else
                values[Indices::velocityXIdx] = inletVelocity_;
        }
        else
        {
            values[Indices::pressureIdx] = outletPressure_;
#if NONISOTHERMAL
            values[Indices::temperatureIdx] = 283.15;
#endif
            const static Scalar xi = getParam<Scalar>("Phasefield.xi");
            const static Scalar rad = getParam<Scalar>("Phasefield.StartingRadius");
            const static Scalar S = getParam<Scalar>("Phasefield.StartingScaling");
            const static Scalar grainX = getParam<Scalar>("Phasefield.GrainX");
            const static Scalar grainY = getParam<Scalar>("Phasefield.GrainY");
            if (icType == "ThinChannel")
            {
                const Scalar s = (globalPos[1]-1.0)*(globalPos[1]-1.0) -rad*rad;
                values[Indices::phi1Idx] = 1.0/(1.0 + std::exp(S*s/xi));
                values[Indices::phi2Idx] = 1.0 - values[Indices::phi1Idx];
                values[Indices::phi3Idx] = 0.0;
            }
            else if (icType == "Grain")
            {
                const Scalar s = (globalPos[1]-grainY)*(globalPos[1]-grainY) +
                    (globalPos[0]-grainX)*(globalPos[0]-grainX) -rad*rad;
                const Scalar s2 = (globalPos[0] - grainX);
                values[Indices::phi1Idx] = 1.0/(1.0 + std::exp(S*-s/xi));
                values[Indices::phi2Idx] = (1.0 - values[Indices::phi1Idx])/(1.0 + std::exp(S*-s2/xi));
                values[Indices::phi3Idx] = 1-values[Indices::phi1Idx]-values[Indices::phi2Idx];
            }
            else
            {
                values[Indices::phi1Idx] = 1.0;
                values[Indices::phi2Idx] = 0.0;
                values[Indices::phi3Idx] = 0.0;
            }
            const auto a0 = getParam<Scalar>("Phasefield.InitialAConcentration");
            const auto b0 = getParam<Scalar>("Phasefield.InitialBConcentration");
            const auto c0 = getParam<Scalar>("Phasefield.InitialCConcentration");
            values[Indices::u1Idx] = a0;
            values[Indices::u2Idx] = b0;
            values[Indices::u3Idx] = c0;
        }



        return values;
    }

    /*!
     * \brief Returns a reference pressure at a given sub control volume face.
     *        This pressure is substracted from the actual pressure for the momentum balance
     *        which potentially helps to improve numerical accuracy by avoiding issues related do floating point arithmetic.
     */
    Scalar referencePressure(const Element& element,
                             const FVElementGeometry& fvGeometry,
                             const SubControlVolumeFace& scvf) const
    { return outletPressure_; }

    // \}
    void setTimeLoop(TimeLoopPtr timeLoop)
    {
        timeLoop_ = timeLoop;
    }

    Scalar time() const
    {
        return timeLoop_->time();
    }

    Scalar timestepsize() const
    {
        return timeLoop_->timeStepSize();
    }

    /*!
     * \brief Return the analytical solution of the problem at a given position
     *
     * \param globalPos The global position
     */
    PrimaryVariables analyticalSolution(const GlobalPosition& globalPos) const
    {
        PrimaryVariables values;

        if constexpr (ParentType::isMomentumProblem())
        {
            values[Indices::velocityXIdx] = parabolicProfile(globalPos[1], inletVelocity_);
            values[Indices::velocityYIdx] = 0.0;
        }
        else
        {
            const Scalar yMin = this->gridGeometry().bBoxMin()[1];
            const Scalar yMax = this->gridGeometry().bBoxMax()[1];
            static const Scalar kinematicViscosity = getParam<Scalar>("Component.LiquidKinematicViscosity");
            const Scalar velocityQuadraticCoefficient = - inletVelocity_ / (0.25*(yMax - yMin)*(yMax - yMin));
            const Scalar pressureLinearCoefficient = 2.0 * velocityQuadraticCoefficient * kinematicViscosity;
            const Scalar pressureConstant = -pressureLinearCoefficient * this->gridGeometry().bBoxMax()[0]  + outletPressure_;
            values[Indices::pressureIdx] =  pressureLinearCoefficient * globalPos[0] + pressureConstant;
        }

        return values;
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

    static constexpr Scalar eps_ = 1e-6;
    Scalar inletVelocity_;
    Scalar outletPressure_;
    OutletCondition outletCondition_;
    bool useVelocityProfile_;
    TimeLoopPtr timeLoop_;
};
} // end namespace Dumux

#endif
