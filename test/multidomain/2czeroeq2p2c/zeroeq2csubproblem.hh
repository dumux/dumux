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
 * \brief Isothermal two-component ZeroEq subproblem with air flowing
 *        from the left to the right and coupling at the bottom.
 */
#ifndef DUMUX_ZEROEQTWOCSUBPROBLEM_HH
#define DUMUX_ZEROEQTWOCSUBPROBLEM_HH

#include <dumux/freeflow/zeroeqnc/model.hh>
#include <dumux/multidomain/subdomainpropertydefaults.hh>
#include <dumux/multidomain/2cstokes2p2c/stokesnccouplinglocalresidual.hh>

namespace Dumux
{

template <class TypeTag>
class ZeroEq2cSubProblem;

namespace Properties
{
NEW_TYPE_TAG(ZeroEq2cSubProblem,
             INHERITS_FROM(BoxZeroEqnc, SubDomain));

// Set the problem property
SET_TYPE_PROP(ZeroEq2cSubProblem, Problem, ZeroEq2cSubProblem<TypeTag>);

// Use the StokencCouplingLocalResidual for the computation of the local residual in the ZeroEq domain
SET_TYPE_PROP(ZeroEq2cSubProblem, LocalResidual,
              StokesncCouplingLocalResidual<TypeTag>);

// Used the fluid system from the coupled problem
SET_TYPE_PROP(ZeroEq2cSubProblem,
              FluidSystem,
              typename GET_PROP_TYPE(typename GET_PROP_TYPE(TypeTag, MultiDomainTypeTag), FluidSystem));

// Disable use of mole formulation
SET_BOOL_PROP(ZeroEq2cSubProblem, UseMoles, false);

// Disable gravity
SET_BOOL_PROP(ZeroEq2cSubProblem, ProblemEnableGravity, false);

// Enable Navier-Stokes
SET_BOOL_PROP(ZeroEq2cSubProblem, EnableNavierStokes, true);

// Set the properties for variable inflow BC
NEW_PROP_TAG(FreeFlowSinusVelocityAmplitude);
NEW_PROP_TAG(FreeFlowSinusVelocityPeriod);
SET_SCALAR_PROP(ZeroEq2cSubProblem, FreeFlowSinusVelocityAmplitude, 0.0);
SET_SCALAR_PROP(ZeroEq2cSubProblem, FreeFlowSinusVelocityPeriod, 3600.0);
NEW_PROP_TAG(FreeFlowSinusPressureAmplitude);
NEW_PROP_TAG(FreeFlowSinusPressurePeriod);
SET_SCALAR_PROP(ZeroEq2cSubProblem, FreeFlowSinusPressureAmplitude, 0.0);
SET_SCALAR_PROP(ZeroEq2cSubProblem, FreeFlowSinusPressurePeriod, 3600.0);
NEW_PROP_TAG(FreeFlowSinusConcentrationAmplitude);
NEW_PROP_TAG(FreeFlowSinusConcentrationPeriod);
SET_SCALAR_PROP(ZeroEq2cSubProblem, FreeFlowSinusConcentrationAmplitude, 0.0);
SET_SCALAR_PROP(ZeroEq2cSubProblem, FreeFlowSinusConcentrationPeriod, 3600.0);
NEW_PROP_TAG(FreeFlowSinusTemperatureAmplitude);
NEW_PROP_TAG(FreeFlowSinusTemperaturePeriod);
SET_SCALAR_PROP(ZeroEq2cSubProblem, FreeFlowSinusTemperatureAmplitude, 0.0);
SET_SCALAR_PROP(ZeroEq2cSubProblem, FreeFlowSinusTemperaturePeriod, 3600.0);
}

/*!
 * \ingroup ImplicitTestProblems
 * \ingroup TwoPTwoCZeroEqTwoCModel
 * \brief Isothermal two-component ZeroEq subproblem with air flowing
 *        from the left to the right and coupling at the bottom.
 *
 * The free-flow subdomain is sized 0.5m times 0.5m. Dry air is flowing from left (Dirichlet)
 * to right (outflow), at the right half of the bottom the coupling conditions
 * are applied to all balance equations. They handle the exchange to the porous-medium
 * subdomain.
 *
 * This subproblem uses the \ref ZeroEqncModel. It is part of a multidomain model and
 * combined with the 2p2csubproblem for the porous-medium domain.
 */
template <class TypeTag>
class ZeroEq2cSubProblem : public ZeroEqProblem<TypeTag>
{
    typedef ZeroEq2cSubProblem<TypeTag> ThisType;
    typedef ZeroEqProblem<TypeTag> ParentType;

    typedef typename GET_PROP_TYPE(TypeTag, GridView) GridView;
    typedef typename GET_PROP_TYPE(TypeTag, TimeManager) TimeManager;

    typedef typename GET_PROP_TYPE(TypeTag, Indices) Indices;

    enum {
        dim = GridView::dimension
    };
    enum { // equation indices
        massBalanceIdx = Indices::massBalanceIdx,
        momentumXIdx = Indices::momentumXIdx, // Index of the x-component of the momentum balance
        momentumYIdx = Indices::momentumYIdx, // Index of the y-component of the momentum balance
        momentumZIdx = Indices::momentumZIdx, // Index of the z-component of the momentum balance
        transportEqIdx = Indices::transportEqIdx // Index of the transport equation (massfraction)
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
        transportCompIdx = Indices::transportCompIdx, // water component index
        phaseCompIdx = Indices::phaseCompIdx          // air component index
    };

    typedef typename GET_PROP_TYPE(TypeTag, PrimaryVariables) PrimaryVariables;
    typedef typename GET_PROP_TYPE(TypeTag, BoundaryTypes) BoundaryTypes;

    typedef typename GridView::template Codim<0>::Entity Element;
    typedef typename GridView::template Codim<dim>::Entity Vertex;
    typedef typename GridView::ctype CoordScalar;
    typedef typename GridView::Intersection Intersection;

    typedef typename GET_PROP_TYPE(TypeTag, FVElementGeometry) FVElementGeometry;
    typedef typename GET_PROP_TYPE(TypeTag, FluidSystem) FluidSystem;
    typedef typename GET_PROP_TYPE(TypeTag, FluidState) FluidState;
    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
    typedef Dune::FieldVector<CoordScalar, dim> GlobalPosition;


public:
    /*!
     * \brief The sub-problem for the ZeroEq subdomain
     *
     * \param timeManager The TimeManager which is used by the simulation
     * \param gridView The simulation's idea about physical space
     */
    ZeroEq2cSubProblem(TimeManager &timeManager, const GridView gridView)
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
    }

    // functions have to be overwritten, otherwise they remain uninitialized
    //! \copydoc ImplicitProblem::bBoxMin()
    const GlobalPosition &bBoxMin() const
    { return bBoxMin_; }

    //! \copydoc ImplicitProblem::bBoxMax()
    const GlobalPosition &bBoxMax() const
    { return bBoxMax_; }

    /*!
     * \brief Returns the problem name
     *
     * This is used as a prefix for files generated by the simulation.
     */
    const std::string &name() const
    { return GET_RUNTIME_PARAM_FROM_GROUP(TypeTag, std::string, Output, NameFF); }

    /*!
     * \brief Returns the temperature \f$ K \f$
     *
     * \param globalPos The global position
     */
    Scalar temperatureAtPos(const GlobalPosition &globalPos) const
    {
        return refTemperature_;
    }

    /*!
     * \name Boundary conditions
     */
    // \{

    //! \copydoc ImplicitProblem::boundaryTypesAtPos()
    void boundaryTypesAtPos(BoundaryTypes &values,
                            const GlobalPosition &globalPos) const
    {
        values.setAllDirichlet();

        if (onUpperBoundary_(globalPos))
        {
            values.setNeumann(transportEqIdx);
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

            if (globalPos[0] > runUpDistanceX_-eps_)
            {
                values.setAllCouplingDirichlet();
                values.setCouplingNeumann(momentumXIdx);
                values.setCouplingNeumann(momentumYIdx);
            }
        }
        if (onLeftBoundary_(globalPos))
        {
            // Left inflow boundaries should be Neumann, otherwise the
            // evaporative fluxes are much more grid dependent
            values.setNeumann(transportEqIdx);

            if (onUpperBoundary_(globalPos) || onLowerBoundary_(globalPos)) // corner point
                values.setAllDirichlet();
        }

        // the mass balance has to be of type outflow
        // it does not get a coupling condition, since pn is a condition for ZeroEq
        values.setOutflow(massBalanceIdx);

        if (onRightBoundary_(globalPos))
            values.setDirichlet(pressureIdx);
    }

    //! \copydoc ImplicitProblem::dirichletAtPos()
    void dirichletAtPos(PrimaryVariables &values,
                        const GlobalPosition &globalPos) const
    {
        initial_(values, globalPos);
    }

    //! \copydoc ImplicitProblem::neumannAtPos()
    void neumannAtPos(PrimaryVariables &values,
                      const GlobalPosition &globalPos) const
    {
        values = 0.0;

        FluidState fluidState;
        updateFluidStateForBC_(fluidState);
        const Scalar density = FluidSystem::density(fluidState, phaseIdx);

        // rho*v*X at inflow
        if (onLeftBoundary_(globalPos)
                && globalPos[1] > bBoxMin_[1] - eps_ && globalPos[1] < bBoxMax_[1] + eps_)
        {
            values[transportEqIdx] = -xVelocity_(globalPos) * density * refMassfrac();
        }
    }

    // \}

    /*!
     * \name Volume terms
     */
    // \{

    //! \copydoc ImplicitProblem::sourceAtPos()
    void sourceAtPos(PrimaryVariables &values,
                     const GlobalPosition &globalPos) const
    {
        // The source term of the mass balance has to be chosen as
        // div (q_momentum) in the problem file
        values = Scalar(0);
    }

    //! \copydoc ImplicitProblem::initialAtPos()
    void initialAtPos(PrimaryVariables &values,
                      const GlobalPosition &globalPos) const
    {
        initial_(values, globalPos);
    }

    // \}

    //! \brief Returns the velocity at the inflow.
    const Scalar refVelocity() const
    {
        return refVelocity_ + variation_(GET_PARAM_FROM_GROUP(TypeTag, Scalar, FreeFlow, SinusVelocityAmplitude),
                                         GET_PARAM_FROM_GROUP(TypeTag, Scalar, FreeFlow, SinusVelocityPeriod));
    }

    //! \brief Returns the pressure at the inflow.
    const Scalar refPressure() const
    {
        return refPressure_ + variation_(GET_PARAM_FROM_GROUP(TypeTag, Scalar, FreeFlow, SinusPressureAmplitude),
                                         GET_PARAM_FROM_GROUP(TypeTag, Scalar, FreeFlow, SinusPressurePeriod));
    }

    //! \brief Returns the mass fraction at the inflow.
    const Scalar refMassfrac() const
    {
        return refMassfrac_ + variation_(GET_PARAM_FROM_GROUP(TypeTag, Scalar, FreeFlow, SinusConcentrationAmplitude),
                                         GET_PARAM_FROM_GROUP(TypeTag, Scalar, FreeFlow, SinusConcentrationPeriod));
    }

    //! \brief Returns the temperature at the inflow.
    const Scalar refTemperature() const
    {
        return refTemperature_ + variation_(GET_PARAM_FROM_GROUP(TypeTag, Scalar, FreeFlow, SinusTemperatureAmplitude),
                                            GET_PARAM_FROM_GROUP(TypeTag, Scalar, FreeFlow, SinusTemperaturePeriod));
    }

private:
    // Internal method for the initial condition (reused for the dirichlet conditions!)
    void initial_(PrimaryVariables &values,
                  const GlobalPosition &globalPos) const
    {
        FluidState fluidState;
        updateFluidStateForBC_(fluidState);
        const Scalar density = FluidSystem::density(fluidState, phaseIdx);

        values[velocityXIdx] = xVelocity_(globalPos);
        values[velocityYIdx] = 0.0;
        values[pressureIdx] = refPressure()
                              + density * this->gravity()[1] * (globalPos[1] - bBoxMin_[1]);
        values[massOrMoleFracIdx] = refMassfrac();
    }

    // Set the profile of the inflow velocity (horizontal direction)
    const Scalar xVelocity_(const GlobalPosition &globalPos) const
    {
        if (onUpperBoundary_(globalPos) || onLowerBoundary_(globalPos))
            return 0.0;
        return refVelocity();
    }

    // Updates the fluid state to obtain required quantities for IC/BC
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

    static constexpr Scalar eps_ = 1e-8;
    GlobalPosition bBoxMin_;
    GlobalPosition bBoxMax_;
    Scalar runUpDistanceX_;

    Scalar refVelocity_;
    Scalar refPressure_;
    Scalar refMassfrac_;
    Scalar refTemperature_;
};
} //end namespace Dumux

#endif // DUMUX_ZEROEQ2C_SUBPROBLEM_HH
