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

#include <dumux/discretization/staggered/freeflow/properties.hh>

#include <dumux/freeflow/navierstokes/boundarytypes.hh>
#include <dumux/freeflow/navierstokes/model.hh>
#include <dumux/freeflow/navierstokes/problem.hh>

#include <dumux/material/fluidsystems/1pliquid.hh>
#include <dumux/material/components/constant.hh>
#include <dumux/material/components/simpleh2o.hh>

#include "../../l2error.hh"

namespace Dumux {
template <class TypeTag>
class ChannelTestProblem;

namespace Properties {
// Create new type tags
namespace TTag {
#if !NONISOTHERMAL
struct ChannelTest { using InheritsFrom = std::tuple<NavierStokes, StaggeredFreeFlowModel>; };
#else
struct ChannelTest { using InheritsFrom = std::tuple<NavierStokesNI, StaggeredFreeFlowModel>; };
#endif
} // end namespace TTag

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

// Set the problem property
template<class TypeTag>
struct Problem<TypeTag, TTag::ChannelTest> { using type = Dumux::ChannelTestProblem<TypeTag> ; };

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

    using BoundaryTypes = Dumux::NavierStokesBoundaryTypes<GetPropType<TypeTag, Properties::ModelTraits>::numEq()>;
    using GridGeometry = GetPropType<TypeTag, Properties::GridGeometry>;
    using FVElementGeometry = typename GridGeometry::LocalView;
    using SubControlVolumeFace = typename GridGeometry::SubControlVolumeFace;
    using Indices = typename GetPropType<TypeTag, Properties::ModelTraits>::Indices;
    using NumEqVector = GetPropType<TypeTag, Properties::NumEqVector>;
    using PrimaryVariables = GetPropType<TypeTag, Properties::PrimaryVariables>;
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using ModelTraits = GetPropType<TypeTag, Properties::ModelTraits>;
    using SolutionVector = GetPropType<TypeTag, Properties::SolutionVector>;

    using Element = typename GridGeometry::GridView::template Codim<0>::Entity;
    using GlobalPosition = typename Element::Geometry::GlobalCoordinate;
    static constexpr auto dimWorld = GetPropType<TypeTag, Properties::GridGeometry>::GridView::dimensionworld;
    using VelocityVector = Dune::FieldVector<Scalar, dimWorld>;

    using TimeLoopPtr = std::shared_ptr<CheckPointTimeLoop<Scalar>>;

    // the types of outlet boundary conditions
    enum class OutletCondition
    {
        outflow, doNothing, neumannXdirichletY, neumannXneumannY
    };

public:
    ChannelTestProblem(std::shared_ptr<const GridGeometry> gridGeometry)
    : ParentType(gridGeometry)
    {
        inletVelocity_ = getParam<Scalar>("Problem.InletVelocity");
        kinematicViscosity_ = getParam<Scalar>("Component.LiquidKinematicViscosity", 1.0);

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

        hasAnalyticalSolution_ = (outletCondition_ == OutletCondition::neumannXneumannY);
        if (hasAnalyticalSolution_)
            createAnalyticalSolution_();
    }

   /*!
     * \name Problem parameters
     */
    // \{

    void printL2Error(const SolutionVector& curSol) const
    {
        using L2Error = NavierStokesTestL2Error<Scalar, ModelTraits, PrimaryVariables>;
        const auto l2error = L2Error::calculateL2Error(*this, curSol);
        const int numCellCenterDofs = this->gridGeometry().numCellCenterDofs();
        const int numFaceDofs = this->gridGeometry().numFaceDofs();
        std::cout << std::setprecision(8) << "** L2 error (abs/rel) for "
                << std::setw(6) << numCellCenterDofs << " cc dofs and " << numFaceDofs << " face dofs (total: " << numCellCenterDofs + numFaceDofs << "): "
                << std::scientific
                << "L2(p) = " << l2error.first[Indices::pressureIdx] << " / " << l2error.second[Indices::pressureIdx]
                << " , L2(vx) = " << l2error.first[Indices::velocityXIdx] << " / " << l2error.second[Indices::velocityXIdx]
                << " , L2(vy) = " << l2error.first[Indices::velocityYIdx] << " / " << l2error.second[Indices::velocityYIdx]
                << std::endl;
    }

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
    BoundaryTypes boundaryTypesAtPos(const GlobalPosition &globalPos) const
    {
        BoundaryTypes values;

        if (isInlet_(globalPos))
        {
            values.setDirichlet(Indices::velocityXIdx);
            values.setDirichlet(Indices::velocityYIdx);
#if NONISOTHERMAL
            values.setDirichlet(Indices::temperatureIdx);
#endif
        }
        else if (isOutlet_(globalPos))
        {
            if (outletCondition_ == OutletCondition::outflow)
                values.setDirichlet(Indices::pressureIdx);
            else if (outletCondition_ == OutletCondition::doNothing ||
                     outletCondition_ == OutletCondition::neumannXneumannY)
            {
                values.setNeumann(Indices::momentumXBalanceIdx);
                values.setNeumann(Indices::momentumYBalanceIdx);
            }
            else // (outletCondition_ == OutletCondition::neumannXdirichletY)
            {
                values.setDirichlet(Indices::velocityYIdx);
                values.setNeumann(Indices::momentumXBalanceIdx);
            }
#if NONISOTHERMAL
            values.setOutflow(Indices::energyEqIdx);
#endif
        }
        else
        {
            values.setDirichlet(Indices::velocityXIdx);
            values.setDirichlet(Indices::velocityYIdx);
#if NONISOTHERMAL
            values.setNeumann(Indices::energyEqIdx);
#endif
        }

        return values;
    }

    /*!
      * \brief Evaluates the boundary conditions for a Dirichlet control volume.
      *
      * \param globalPos The center of the finite volume which ought to be set.
      */
    PrimaryVariables dirichletAtPos(const GlobalPosition &globalPos) const
    {
        PrimaryVariables values = initialAtPos(globalPos);

        if(isInlet_(globalPos))
        {
#if NONISOTHERMAL
            // give the system some time so that the pressure can equilibrate, then start the injection of the hot liquid
            if(time() >= 200.0)
                 values[Indices::temperatureIdx] = 293.15;
#endif
        }
        else
            values[Indices::velocityXIdx] = 0.0;

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
    template<class ElementVolumeVariables, class ElementFaceVariables>
    NumEqVector neumann(const Element& element,
                        const FVElementGeometry& fvGeometry,
                        const ElementVolumeVariables& elemVolVars,
                        const ElementFaceVariables& elemFaceVars,
                        const SubControlVolumeFace& scvf) const
    {
        NumEqVector values(0.0);

        values[scvf.directionIndex()] = initialAtPos(scvf.center())[Indices::pressureIdx];

        // make sure to normalize the pressure if the property is set true
        if (getPropValue<TypeTag, Properties::NormalizePressure>())
            values[scvf.directionIndex()] -= initialAtPos(scvf.center())[Indices::pressureIdx];

        if (outletCondition_ != OutletCondition::doNothing)
            values[1] = -dudy(scvf.center()[1], inletVelocity_) * elemVolVars[scvf.insideScvIdx()].viscosity();

        return values;
    }

    /*!
     * \brief A parabolic velocity profile.
     *
     * \param y The position where the velocity is evaluated.
     * \param vMax The profile's maxmium velocity.
     */
    Scalar parabolicProfile(const Scalar y, const Scalar vMax) const
    {
        const Scalar yMin = this->gridGeometry().bBoxMin()[1];
        const Scalar yMax = this->gridGeometry().bBoxMax()[1];
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

    /*!
     * \brief Return the analytical solution of the problem at a given position
     *
     * \param globalPos The global position
     */
    PrimaryVariables analyticalSolution(const GlobalPosition& globalPos) const
    {
        Scalar x = globalPos[0];
        Scalar y = globalPos[1];

        const Scalar yMin = this->gridGeometry().bBoxMin()[1];
        const Scalar yMax = this->gridGeometry().bBoxMax()[1];

        Scalar velocityQuadraticCoefficient = - inletVelocity_ / (0.25*(yMax - yMin)*(yMax - yMin));
        Scalar pressureLinearCoefficient = 2.0 * velocityQuadraticCoefficient * kinematicViscosity_;
        Scalar pressureConstant = -pressureLinearCoefficient * this->gridGeometry().bBoxMax()[0]  + outletPressure_;

        PrimaryVariables values;
        values[Indices::pressureIdx] =  pressureLinearCoefficient * x + pressureConstant;
        values[Indices::velocityXIdx] = parabolicProfile(y, inletVelocity_);
        values[Indices::velocityYIdx] = 0;

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
    PrimaryVariables initialAtPos(const GlobalPosition& globalPos) const
    {
        PrimaryVariables values;

        values[Indices::pressureIdx] = outletPressure_;
        values[Indices::velocityYIdx] = 0.0;

        if (useVelocityProfile_)
            values[Indices::velocityXIdx] = parabolicProfile(globalPos[1], inletVelocity_);
        else
            values[Indices::velocityXIdx] = inletVelocity_;

#if NONISOTHERMAL
        values[Indices::temperatureIdx] = 283.15;
#endif

        return values;
    }

    // \}
    void setTimeLoop(TimeLoopPtr timeLoop)
    {
        timeLoop_ = timeLoop;
        if(inletVelocity_ > eps_)
            timeLoop_->setCheckPoint({200.0, 210.0});
    }

    Scalar time() const
    {
        return timeLoop_->time();
    }

/*!
     * \brief Returns the analytical solution for the pressure
     */
    auto& getAnalyticalPressureSolution() const
    {
        return analyticalPressure_;
    }

   /*!
     * \brief Returns the analytical solution for the velocity
     */
    auto& getAnalyticalVelocitySolution() const
    {
        return analyticalVelocity_;
    }

   /*!
     * \brief Returns the analytical solution for the velocity at the faces
     */
    auto& getAnalyticalVelocitySolutionOnFace() const
    {
        return analyticalVelocityOnFace_;
    }

    bool hasAnalyticalSolution()
    {
        return hasAnalyticalSolution_;
    }

private:

   /*!
     * \brief Adds additional VTK output data to the VTKWriter. Function is called by the output module on every write.
     */
    void createAnalyticalSolution_()
    {
        analyticalPressure_.resize(this->gridGeometry().numCellCenterDofs());
        analyticalVelocity_.resize(this->gridGeometry().numCellCenterDofs());
        analyticalVelocityOnFace_.resize(this->gridGeometry().numFaceDofs());

        for (const auto& element : elements(this->gridGeometry().gridView()))
        {
            auto fvGeometry = localView(this->gridGeometry());
            fvGeometry.bindElement(element);
            for (auto&& scv : scvs(fvGeometry))
            {
                auto ccDofIdx = scv.dofIndex();
                auto ccDofPosition = scv.dofPosition();
                auto analyticalSolutionAtCc = analyticalSolution(ccDofPosition);

                // velocities on faces
                for (auto&& scvf : scvfs(fvGeometry))
                {
                    const auto faceDofIdx = scvf.dofIndex();
                    const auto faceDofPosition = scvf.center();
                    const auto dirIdx = scvf.directionIndex();
                    const auto analyticalSolutionAtFace = analyticalSolution(faceDofPosition);
                    analyticalVelocityOnFace_[faceDofIdx][dirIdx] = analyticalSolutionAtFace[Indices::velocity(dirIdx)];
                }

                analyticalPressure_[ccDofIdx] = analyticalSolutionAtCc[Indices::pressureIdx];

                for(int dirIdx = 0; dirIdx < ModelTraits::dim(); ++dirIdx)
                    analyticalVelocity_[ccDofIdx][dirIdx] = analyticalSolutionAtCc[Indices::velocity(dirIdx)];
            }
        }
     }

    bool isInlet_(const GlobalPosition& globalPos) const
    {
        return globalPos[0] < eps_;
    }

    bool isOutlet_(const GlobalPosition& globalPos) const
    {
        return globalPos[0] > this->gridGeometry().bBoxMax()[0] - eps_;
    }

    static constexpr Scalar eps_=1e-6;
    bool printL2Error_;
    Scalar inletVelocity_;
    Scalar kinematicViscosity_;
    Scalar outletPressure_;
    OutletCondition outletCondition_;
    bool useVelocityProfile_;
    bool hasAnalyticalSolution_;
    TimeLoopPtr timeLoop_;

    std::vector<Scalar> analyticalPressure_;
    std::vector<VelocityVector> analyticalVelocity_;
    std::vector<VelocityVector> analyticalVelocityOnFace_;
};
} // end namespace Dumux

#endif
