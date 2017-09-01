// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*****************************************************************************
 *   See the file COPYING for full copying permissions.                      *
 *                                                                           *
 *   This program is free software: you can redistribute it and/or modify    *
 *   it under the terms of the GNU General Public License as published by    *
 *   the Free Software Foundation, either version 2 of the License, or       *
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
 * \brief Isothermal two-component stokes subproblem with air flowing
 *        from the left to the right and coupling at the bottom.
 */
#ifndef DUMUX_STOKES2C_SUBPROBLEM_HH
#define DUMUX_STOKES2C_SUBPROBLEM_HH

#include <dumux/implicit/staggered/properties.hh>
#include <dumux/freeflow/staggerednc/model.hh>
#include <dumux/implicit/problem.hh>

// coupling-specific includes
#include <dumux/multidomain/subproblemproperties.hh>

// TODO necessary?
//#include <dumux/material/components/simpleh2o.hh>
//#include <dumux/material/fluidsystems/liquidphase.hh>
//#include <dumux/material/components/constant.hh>

//#include <dumux/multidomain/2cstokes2p2c/stokesnccouplinglocalresidual.hh>
//#include <dumux/multidomain/subdomainpropertydefaults.hh>

namespace Dumux
{

template <class TypeTag>
class Stokes2cSubProblem;

namespace Properties
{
NEW_TYPE_TAG(Stokes2cSubProblem, INHERITS_FROM(StaggeredModel, NavierStokesNC));

// Set the problem property
SET_TYPE_PROP(Stokes2cSubProblem, Problem, Stokes2cSubProblem<TypeTag>);

// Set the grid type
#if ENABLE_3D
SET_TYPE_PROP(Stokes2cSubProblem, Grid, Dune::YaspGrid<3, Dune::TensorProductCoordinates<typename GET_PROP_TYPE(TypeTag, Scalar), 3> >);
#else
SET_TYPE_PROP(Stokes2cSubProblem, Grid, Dune::YaspGrid<2, Dune::TensorProductCoordinates<typename GET_PROP_TYPE(TypeTag, Scalar), 2> >);
#endif

// Used the fluid system from the coupled problem
SET_TYPE_PROP(Stokes2cSubProblem, FluidSystem,
//              typename GET_PROP_TYPE(typename GET_PROP_TYPE(TypeTag, MultiDomainTypeTag), FluidSystem));
              typename GET_PROP_TYPE(typename GET_PROP_TYPE(TypeTag, GlobalProblemTypeTag), FluidSystem));

// Use formulation based on mass fractions
SET_BOOL_PROP(Stokes2cSubProblem, UseMoles, false);

// Disable gravity
SET_BOOL_PROP(Stokes2cSubProblem, ProblemEnableGravity, false);

// Switch inertia term off
SET_BOOL_PROP(Stokes2cSubProblem, EnableInertiaTerms, false);

// TODO not needed in old multidomain, used in stokestestproblem 1p
//SET_PROP(StokesTestProblem, Fluid)
//{
//private:
//    using Scalar = typename GET_PROP_TYPE(TypeTag, Scalar);
//public:
//    using type = FluidSystems::LiquidPhase<Scalar, Dumux::Constant<TypeTag, Scalar> >;
//};
//
SET_BOOL_PROP(Stokes2cSubProblem, EnableGlobalFVGeometryCache, true);
//
SET_BOOL_PROP(Stokes2cSubProblem, EnableGlobalFluxVariablesCache, true);
SET_BOOL_PROP(Stokes2cSubProblem, EnableGlobalVolumeVariablesCache, true);
//
//
//// Set the grid parameter group
SET_STRING_PROP(Stokes2cSubProblem, GridParameterGroup, "StokesGrid");

NEW_PROP_TAG(GlobalProblemTypeTag);
NEW_PROP_TAG(CouplingManager);
}

/*!
 * \ingroup ImplicitTestProblems
 * \ingroup TwoPTwoCStokesTwoCModel
 * \brief Isothermal two-component stokes subproblem with air flowing
 *        from the left to the right and coupling at the bottom.
 *
 * The Stokes subdomain is sized 0.25m times 0.25m. The boundary conditions
 * for the momentum balances are all set to Dirichlet, except on the right
 * boundary, where outflow conditions are set. The mass balance receives
 * outflow BCs, which are replaced in the localresidual by the sum
 * of the two momentum balances. On the right boundary Dirichlet BCs are
 * set for the pressure.
 *
 * This sub problem uses the \ref NavierStokesNCModel. It is part of the
 * 2cstokes2p2c model and is combined with the 2p2csubproblem for
 * the Darcy domain.
 */
template <class TypeTag>
class Stokes2cSubProblem : public NavierStokesProblem<TypeTag>
{
    using ParentType = NavierStokesProblem<TypeTag>;
    using GridView = typename GET_PROP_TYPE(TypeTag, GridView);
    using TimeManager = typename GET_PROP_TYPE(TypeTag, TimeManager);

    using Indices = typename GET_PROP_TYPE(TypeTag, Indices); // TODO from original 2cstokes2p2c
    enum {
        // Number of equations and grid dimension
        numEq = GET_PROP_VALUE(TypeTag, NumEq),
        dim = GridView::dimension,
        dimWorld = GridView::dimensionworld
    };
    enum { // equation indices
        massBalanceIdx = Indices::massBalanceIdx,

        momentumXIdx = Indices::momentumXBalanceIdx, //!< Index of the x-component of the momentum balance
        momentumYIdx = Indices::momentumYBalanceIdx, //!< Index of the y-component of the momentum balance
        momentumZIdx = Indices::momentumZBalanceIdx, //!< Index of the z-component of the momentum balance

        transportEqIdx = Indices::transportEqIdx //!< Index of the transport equation (massfraction)
    };
    enum { // primary variable indices
            pressureIdx = Indices::pressureIdx,
            velocityXIdx = Indices::velocityXIdx,
            velocityYIdx = Indices::velocityYIdx,
            velocityZIdx = Indices::velocityZIdx,
            massOrMoleFracIdx = Indices::massOrMoleFracIdx
    };
    enum { phaseIdx = Indices::phaseIdx };
    enum { numComponents = Indices::numComponents };
    enum {
        transportCompIdx = Indices::transportCompIdx, //!< water component index
        phaseCompIdx = Indices::phaseCompIdx          //!< air component index
    };

    using PrimaryVariables = typename GET_PROP_TYPE(TypeTag, PrimaryVariables);
    using BoundaryTypes = typename GET_PROP_TYPE(TypeTag, BoundaryTypes);

    using Element = typename GridView::template Codim<0>::Entity;
    using Vertex = typename GridView::template Codim<dim>::Entity;
    using Intersection = typename GridView::Intersection;
    using CoordScalar = typename GridView::ctype;

    using SubControlVolumeFace =  typename GET_PROP_TYPE(TypeTag, SubControlVolumeFace);
    using Fluid = typename GET_PROP_TYPE(TypeTag, Fluid);
    using Scalar = typename GET_PROP_TYPE(TypeTag, Scalar);
    using GlobalPosition = Dune::FieldVector<CoordScalar, dimWorld>;

    using GlobalTypeTag = typename GET_PROP_TYPE(TypeTag, GlobalProblemTypeTag);
    using CouplingManager = typename GET_PROP_TYPE(GlobalTypeTag, CouplingManager);

    using BoundaryValues = typename GET_PROP_TYPE(TypeTag, BoundaryValues);
    using InitialValues = typename GET_PROP_TYPE(TypeTag, BoundaryValues);

    using FVElementGeometry = typename GET_PROP_TYPE(TypeTag, FVElementGeometry);
    using FluidSystem = typename GET_PROP_TYPE(TypeTag, FluidSystem);
    using FluidState = typename GET_PROP_TYPE(TypeTag, FluidState);


public:
    /*!
     * \brief The sub-problem for the Stokes subdomain
     *
     * \param timeManager The TimeManager which is used by the simulation
     * \param gridView The simulation's idea about physical space
     */
    Stokes2cSubProblem(TimeManager &timeManager, const GridView &gridView)
        : ParentType(timeManager, gridView)
    {
        std::vector<Scalar> positions0 = GET_RUNTIME_PARAM_FROM_GROUP(TypeTag, std::vector<Scalar>, Grid, Positions0);
        std::vector<Scalar> positions1 = GET_RUNTIME_PARAM_FROM_GROUP(TypeTag, std::vector<Scalar>, Grid, Positions1);

        bBoxMin_[0] = positions0.front();
        bBoxMax_[0] = positions0.back();
        bBoxMin_[1] = GET_RUNTIME_PARAM_FROM_GROUP(TypeTag, Scalar, Grid, InterfacePosY);
        bBoxMax_[1] = positions1.back();
        runUpDistanceX_ = GET_RUNTIME_PARAM_FROM_GROUP(TypeTag, Scalar, Grid, RunUpDistanceX); // first part of the interface without coupling

        refVelocity_ = GET_RUNTIME_PARAM_FROM_GROUP(TypeTag, Scalar, FreeFlow, RefVelocity);
        refPressure_ = GET_RUNTIME_PARAM_FROM_GROUP(TypeTag, Scalar, FreeFlow, RefPressure);
        refMassfrac_ = GET_RUNTIME_PARAM_FROM_GROUP(TypeTag, Scalar, FreeFlow, RefMassfrac);
        refTemperature_ = GET_RUNTIME_PARAM_FROM_GROUP(TypeTag, Scalar, FreeFlow, RefTemperature);

        sinusVAmplitude_ = GET_RUNTIME_PARAM_FROM_GROUP(TypeTag, Scalar, FreeFlow, SinusVelAmplitude);
        sinusVPeriod_ = GET_RUNTIME_PARAM_FROM_GROUP(TypeTag, Scalar, FreeFlow, SinusVelPeriod);
        sinusPAmplitude_ = GET_RUNTIME_PARAM_FROM_GROUP(TypeTag, Scalar, FreeFlow, SinusPressureAmplitude);
        sinusPPeriod_ = GET_RUNTIME_PARAM_FROM_GROUP(TypeTag, Scalar, FreeFlow, SinusPressurePeriod);
        sinusXAmplitude_ = GET_RUNTIME_PARAM_FROM_GROUP(TypeTag, Scalar, FreeFlow, SinusConcentrationAmplitude);
        sinusXPeriod_ = GET_RUNTIME_PARAM_FROM_GROUP(TypeTag, Scalar, FreeFlow, SinusConcentrationPeriod);
        sinusTAmplitude_ = GET_RUNTIME_PARAM_FROM_GROUP(TypeTag, Scalar, FreeFlow, SinusTemperatureAmplitude);
        sinusTPeriod_ = GET_RUNTIME_PARAM_FROM_GROUP(TypeTag, Scalar, FreeFlow, SinusTemperaturePeriod);

        initializationTime_ = GET_RUNTIME_PARAM_FROM_GROUP(TypeTag, Scalar, TimeManager, InitTime);
    }

    // functions have to be overwritten, otherwise they remain uninitialized
    //! \copydoc ImplicitProblem::bBoxMin()
    const GlobalPosition &bBoxMin() const
    { return bBoxMin_; }

    //! \copydoc ImplicitProblem::bBoxMax()
    const GlobalPosition &bBoxMax() const
    { return bBoxMax_; }

    /*!
     * \name Problem parameters
     */
    // \{

    /*!
     * \brief Returns the problem name
     *
     * This is used as a prefix for files generated by the simulation.
     */
    const std::string &name() const
    { return GET_RUNTIME_PARAM_FROM_GROUP(TypeTag, std::string, Output, NameFF); }

    /*!
     * \brief Returns the temperature within the domain.
     *
     * This problem assumes a constant temperature, which can
     * be set in the parameter file.
     */
    Scalar temperature() const
    {
        return refTemperature_;
    };

    // \}

    /*!
     * \name Boundary conditions
     */
    // \{

    //! \copydoc ImplicitProblem::boundaryTypesAtPos()
    BoundaryTypes boundaryTypesAtPos(const GlobalPosition &globalPos) const
    {
        BoundaryTypes values;
        const Scalar time = this->timeManager().time();

        values.setAllDirichlet();

        if (onUpperBoundary_(globalPos))
            values.setNeumann(transportEqIdx);

        // Left inflow boundaries should be Neumann, otherwise the
        // evaporative fluxes are much more grid dependent
        if (onLeftBoundary_(globalPos))
        {
            values.setNeumann(transportEqIdx);

            if (onUpperBoundary_(globalPos)) // corner point
                values.setAllDirichlet();
        }

        if (onRightBoundary_(globalPos))
        {
            values.setAllOutflow();

            if (onUpperBoundary_(globalPos)) // corner point
                values.setAllDirichlet();
        }

        if (onLowerBoundary_(globalPos))
        {
            values.setAllDirichlet();
            values.setNeumann(transportEqIdx);

            if (globalPos[0] > runUpDistanceX_-eps_ && time > initializationTime_)
            {
                values.setAllCouplingDirichlet();
                values.setCouplingNeumann(momentumXIdx);
                values.setCouplingNeumann(momentumYIdx);
            }
        }

        // the mass balance has to be of type outflow
        // it does not get a coupling condition, since pn is a condition for stokes
        values.setOutflow(massBalanceIdx);

        // set pressure at one point, do NOT specify this
        // if the Darcy domain has a Dirichlet condition for pressure
        if (onRightBoundary_(globalPos))
        {
            if (time > initializationTime_)
                values.setDirichlet(pressureIdx);
            else
                if (!onLowerBoundary_(globalPos) && !onUpperBoundary_(globalPos))
                    values.setDirichlet(pressureIdx);
        }

        return values;
    }

    /*!
     * \brief Return Dirichlet boundary values at a given position
     *
     * \param globalPos The global position
     */
    BoundaryValues dirichletAtPos(const GlobalPosition &globalPos) const
    {
        BoundaryValues values;
        FluidState fluidState;
        updateFluidStateForBC_(fluidState);

        const Scalar density =
                FluidSystem::density(fluidState, phaseIdx);

        values[velocityXIdx] = xVelocity_(globalPos);
        values[velocityYIdx] = 0.0;
        values[pressureIdx] = refPressure()
                + density*this->gravity()[1]*(globalPos[1] - bBoxMin_[1]);
        values[massOrMoleFracIdx] = refMassfrac();

        return values;
    }

    /*!
     * \brief Evaluate the boundary conditions for a Neumann
     *        boundary segment.
     *
     * \param values The Neumann values for the conservation equations in units of
     *                 \f$ [ \textnormal{unit of conserved quantity} / (m^{\textrm{dim}-1} \cdot s )] \f$
     * \param globalPos The global position
     */
    PrimaryVariables neumannAtPos(const GlobalPosition &globalPos) const
    {
        PrimaryVariables values(0.0);

        FluidState fluidState;
        updateFluidStateForBC_(fluidState);

        const Scalar density =
                FluidSystem::density(fluidState, phaseIdx);
        const Scalar xVelocity = xVelocity_(globalPos);

        if (onLeftBoundary_(globalPos)
                && globalPos[1] > bBoxMin_[1] - eps_ && globalPos[1] < bBoxMax_[1] + eps_)
        {
            // rho*v*X at inflow
            values[transportEqIdx] = -xVelocity * density * refMassfrac();
        }

        return values;
    }

    // \}

    /*!
     * \name Volume terms
     */
    // \{

    /*!
     * \brief Returns the source term
     *
     * \param globalPos The global position
     */
    PrimaryVariables sourceAtPos(const GlobalPosition &globalPos) const
    {
        PrimaryVariables values;
        // The source term of the mass balance has to be chosen as
        // div (q_momentum) in the problem file
        values = Scalar(0);

        return values;
    }

    /*!
     * \brief Evaluate the initial value for a control volume.
     *
     * \param globalPos The global position
     */
    InitialValues initialAtPos(const GlobalPosition &globalPos) const
    {
        InitialValues values(0.0);
        initial_(values, globalPos);

        return values;
    }
    // \}

    //! \brief Returns the reference velocity.
    const Scalar refVelocity() const
    { return refVelocity_ + variation_(sinusVAmplitude_, sinusVPeriod_); }

    //! \brief Returns the reference pressure.
    const Scalar refPressure() const
    { return refPressure_ + variation_(sinusPAmplitude_, sinusPPeriod_); }

    //! \brief Returns the reference mass fraction.
    const Scalar refMassfrac() const
    { return refMassfrac_ + variation_(sinusXAmplitude_, sinusXPeriod_); }

    //! \brief Returns the reference temperature.
    const Scalar refTemperature() const
    { return refTemperature_+ variation_(sinusTAmplitude_, sinusTPeriod_); }


    /*!
     * \brief Set the coupling manager
     * \param couplingManager The coupling manager
     *
     */
    void setCouplingManager(std::shared_ptr<CouplingManager> couplingManager)
    {
        couplingManager_ = couplingManager;
    }

    /*!
     * \brief Get the coupling manager
     *
     */
    CouplingManager& couplingManager() const
    { return *couplingManager_; }

    bool onCouplingInterface(const GlobalPosition &globalPos) const
    {return onLowerBoundary_(globalPos); }

private:
    /*!
     * \brief Internal method for the initial condition
     *        (reused for the dirichlet conditions!)
     */
    void initial_(InitialValues &values,
                  const GlobalPosition &globalPos) const
    {
        FluidState fluidState;
        updateFluidStateForBC_(fluidState);

        const Scalar density =
                FluidSystem::density(fluidState, phaseIdx);

        values[velocityXIdx] = xVelocity_(globalPos);
        values[velocityYIdx] = 0.;

        values[pressureIdx] = refPressure()
                + density*this->gravity()[1]*(globalPos[1] - bBoxMin_[1]);
        values[massOrMoleFracIdx] = refMassfrac();
    }

    //! \brief set the profile of the inflow velocity (horizontal direction)
    const Scalar xVelocity_(const GlobalPosition &globalPos) const
    {
        const Scalar vmax = refVelocity();
        return  4*vmax*(globalPos[1] - bBoxMin_[1])*(bBoxMax_[1] - globalPos[1])
                / (height_()*height_()) + 0.00134;
    }

    //! \brief updates the fluid state to obtain required quantities for IC/BC
    void updateFluidStateForBC_(FluidState& fluidState) const
    {
        fluidState.setTemperature(refTemperature());
        fluidState.setPressure(phaseIdx, refPressure());
        // setMassFraction() has only to be called 1-numComponents times
        fluidState.setMassFraction(phaseIdx, transportCompIdx, refMassfrac());
    }

    // can be used for the variation of a boundary condition
    const Scalar variation_(const Scalar amplitude, const Scalar period) const
    { return sin(2*M_PI*this->timeManager().time()/period) * amplitude; }

    bool onLeftBoundary_(const GlobalPosition &globalPos) const
    { return globalPos[0] < bBoxMin_[0] + eps_; }

    bool onRightBoundary_(const GlobalPosition &globalPos) const
    { return globalPos[0] > bBoxMax_[0] - eps_; }

    bool onLowerBoundary_(const GlobalPosition &globalPos) const
    { return globalPos[1] < bBoxMin_[1] + eps_; }

    bool onUpperBoundary_(const GlobalPosition &globalPos) const
    { return globalPos[1] > bBoxMax_[1] - eps_; }

    // the height of the free-flow domain
    const Scalar height_() const
    { return bBoxMax_[1] - bBoxMin_[1]; }

    static constexpr Scalar eps_ = 1e-8;
    GlobalPosition bBoxMin_;
    GlobalPosition bBoxMax_;

    Scalar refVelocity_;
    Scalar refPressure_;
    Scalar refMassfrac_;
    Scalar refTemperature_;

    Scalar sinusVAmplitude_;
    Scalar sinusVPeriod_;
    Scalar sinusPAmplitude_;
    Scalar sinusPPeriod_;
    Scalar sinusXAmplitude_;
    Scalar sinusXPeriod_;
    Scalar sinusTAmplitude_;
    Scalar sinusTPeriod_;

    Scalar runUpDistanceX_;
    Scalar initializationTime_;

    std::shared_ptr<CouplingManager> couplingManager_;

};
} //end namespace Dumux

#endif // DUMUX_STOKES2C_SUBPROBLEM_HH
