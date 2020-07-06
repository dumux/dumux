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

#include <dumux/common/boundarytypes.hh>
#include <dumux/freeflow/navierstokes/momentum/model.hh>
#include <dumux/freeflow/navierstokes/momentum/problem.hh>
#include <dumux/freeflow/navierstokes/massandenergy/model.hh>
#include <dumux/freeflow/navierstokes/massandenergy/problem.hh>
#include <dumux/discretization/fcstaggered.hh>
#include <dumux/discretization/cctpfa.hh>

namespace Dumux {
template <class TypeTag>
class ChannelTestProblem;

namespace Properties {
// // Create new type tags
// namespace TTag {
// #if !NONISOTHERMAL
// struct ChannelTest { using InheritsFrom = std::tuple<NavierStokes, StaggeredFreeFlowModel>; };
// #else
// struct ChannelTest { using InheritsFrom = std::tuple<NavierStokesNI, StaggeredFreeFlowModel>; };
// #endif
// } // end namespace TTag
// Create new type tags
namespace TTag {
struct ChannelTest {};
struct ChannelTestMomentum { using InheritsFrom = std::tuple<ChannelTest, NavierStokesMomentum, FaceCenteredStaggeredModel>; };
struct ChannelTestMass { using InheritsFrom = std::tuple<ChannelTest, NavierStokesMassAndEnergy, CCTpfaModel>; };
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

namespace Impl {

template<class TypeTag>
constexpr bool isMomentumProblem()
{
    return GetPropType<TypeTag, Properties::GridGeometry>::discMethod == DiscretizationMethod::fcstaggered;
};

template<class TypeTag>
using BaseProblem = std::conditional_t<isMomentumProblem<TypeTag>(),
                                       NavierStokesMomentumProblem<TypeTag>,
                                       NavierStokesMassAndEnergyProblem<TypeTag>>;

template<class TypeTag>
using BoundaryTypes = std::conditional_t<isMomentumProblem<TypeTag>(),
                                         Dumux::BoundaryTypes<GetPropType<TypeTag, Properties::ModelTraits>::dim()>,
                                         Dumux::BoundaryTypes<GetPropType<TypeTag, Properties::ModelTraits>::numEq()>>;
template<class TypeTag>
using NumEqVector = std::conditional_t<isMomentumProblem<TypeTag>(),
                                       Dune::FieldVector<GetPropType<TypeTag, Properties::Scalar>, GetPropType<TypeTag, Properties::ModelTraits>::dim()>,
                                       GetPropType<TypeTag, Properties::NumEqVector>>;

template<class TypeTag>
using PrimaryVariables = NumEqVector<TypeTag>;

}

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
class ChannelTestProblem : public Impl::BaseProblem<TypeTag>
{
    using ParentType = Impl::BaseProblem<TypeTag>;
    using BoundaryTypes = Impl::BoundaryTypes<TypeTag>;
    using GridGeometry = GetPropType<TypeTag, Properties::GridGeometry>;
    using FVElementGeometry = typename GridGeometry::LocalView;
    using SubControlVolume = typename FVElementGeometry::SubControlVolume;
    using SubControlVolumeFace = typename FVElementGeometry::SubControlVolumeFace;
    using Indices = typename GetPropType<TypeTag, Properties::ModelTraits>::Indices;
    using NumEqVector = Impl::NumEqVector<TypeTag>;
    using ModelTraits = GetPropType<TypeTag, Properties::ModelTraits>;
    using PrimaryVariables = Impl::PrimaryVariables<TypeTag>;
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using SolutionVector = GetPropType<TypeTag, Properties::SolutionVector>;
    using Element = typename GridGeometry::GridView::template Codim<0>::Entity;
    using GlobalPosition = typename Element::Geometry::GlobalCoordinate;

    using TimeLoopPtr = std::shared_ptr<CheckPointTimeLoop<Scalar>>;
    using CouplingManager = GetPropType<TypeTag, Properties::CouplingManager>;

    // the types of outlet boundary conditions
    enum class OutletCondition
    {
        outflow, doNothing, neumannXdirichletY, neumannXneumannY
    };

public:
    ChannelTestProblem(std::shared_ptr<const GridGeometry> gridGeometry, std::shared_ptr<CouplingManager> couplingManager)
    : ParentType(gridGeometry, couplingManager)
    , couplingManager_(couplingManager)
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
    BoundaryTypes boundaryTypesAtPos(const GlobalPosition &globalPos) const
    {
        BoundaryTypes values;



        // if(isInlet_(globalPos))
        // {
        //     if constexpr (Impl::isMomentumProblem<TypeTag>())
        //     {
        //         values.setDirichlet(Indices::velocityXIdx);
        //         values.setDirichlet(Indices::velocityYIdx);
        //     }
        //     else
        //     {
        //         values.setNeumann(0); // TODO idx
        //     }

        // }
        // else if(isOutlet_(globalPos))
        // {

        //     if constexpr (Impl::isMomentumProblem<TypeTag>())
        //     {
        //         // values.setNeumann(Indices::momentumXBalanceIdx);
        //         // values.setNeumann(Indices::momentumYBalanceIdx);
        //         values.setDirichlet(Indices::velocityXIdx);
        //         values.setDirichlet(Indices::velocityYIdx);
        //     }
        //     else
        //     {
        //         values.setNeumann(0); // TODO idx
        //         // values.setDirichlet(0); // TODO idx
        //     }

        // }
        // else
        // {
            if constexpr (Impl::isMomentumProblem<TypeTag>())
            {
                values.setDirichlet(Indices::velocityXIdx);
                values.setDirichlet(Indices::velocityYIdx);

                // if (isOutlet_(globalPos))
                //     values.setAllNeumann();
            }
            else
                values.setDirichlet(0);
                // values.setNeumann(0);

        // }

        return values;
    }

        //! Enable internal Dirichlet constraints
    static constexpr bool enableInternalDirichletConstraints()
    { return !Impl::isMomentumProblem<TypeTag>(); }

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

        if (scv.dofIndex() == 0)
            values.set(0);
        // the pure Neumann problem is only defined up to a constant
        // we create a well-posed problem by fixing the pressure at one dof
        return values;
    }

    /*!
     * \brief Define the values of internal Dirichlet constraints for a degree of freedom.
     * \param element The finite element
     * \param scv The sub-control volume
     */
    PrimaryVariables internalDirichlet(const Element& element, const SubControlVolume& scv) const
    { return PrimaryVariables(0); }


   /*!
     * \brief Evaluates the boundary conditions for a Dirichlet control volume.
     *
     * \param globalPos The center of the finite volume which ought to be set.
     */
    PrimaryVariables dirichletAtPos(const GlobalPosition &globalPos) const
    {
        PrimaryVariables values = initialAtPos(globalPos);

//         if(isInlet_(globalPos))
//         {
// #if NONISOTHERMAL
//             // give the system some time so that the pressure can equilibrate, then start the injection of the hot liquid
//             if(time() >= 200.0)
//                 values[Indices::temperatureIdx] = 293.15;
// #endif
//         }
//         else
//         {
//             if constexpr (Impl::isMomentumProblem<TypeTag>())
//             {
//                 values[Indices::velocityXIdx] = 0.0;
//             }
//         }
        if constexpr (Impl::isMomentumProblem<TypeTag>())
            values[Indices::velocityXIdx] = parabolicProfile(globalPos[1],1);

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

        if constexpr (Impl::isMomentumProblem<TypeTag>())
        {}
        else
        {
            if (isInlet_(scvf.ipGlobal()) || isOutlet_(scvf.ipGlobal()))
            {
                const auto insideDensity = elemVolVars[scvf.insideScvIdx()].density();
                return elemFluxVarsCache[scvf].normalVelocity() * insideDensity * directionIndex(scvf.unitOuterNormal());
            }
        }



        // values[scvf.directionIndex()] = initialAtPos(scvf.center())[Indices::pressureIdx];

        // // make sure to normalize the pressure if the property is set true
        // if (getPropValue<TypeTag, Properties::NormalizePressure>())
        //     values[scvf.directionIndex()] -= initialAtPos(scvf.center())[Indices::pressureIdx];

        // if (outletCondition_ != OutletCondition::doNothing)
        //     values[1] = -dudy(scvf.center()[1], inletVelocity_) * elemVolVars[scvf.insideScvIdx()].viscosity();

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

        if constexpr (Impl::isMomentumProblem<TypeTag>())
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
        }

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
    Scalar outletPressure_;
    OutletCondition outletCondition_;
    bool useVelocityProfile_;
    TimeLoopPtr timeLoop_;
    std::shared_ptr<CouplingManager> couplingManager_;
};
} // end namespace Dumux

#endif
