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
 * \brief Non-isothermal two-component ZeroEq subproblem with air flowing
 *        from the left to the right and coupling at the bottom.
 */
#ifndef DUMUX_ZEROEQTWOCNI_SUBPROBLEM_HH
#define DUMUX_ZEROEQTWOCNI_SUBPROBLEM_HH

#include <dumux/freeflow/zeroeqncni/zeroeqncnimodel.hh>
#include <dumux/multidomain/common/subdomainpropertydefaults.hh>
#include <dumux/multidomain/2cnistokes2p2cni/stokesncnicouplinglocalresidual.hh>

namespace Dumux
{

template <class TypeTag>
class ZeroEq2cniSubProblem;

namespace Properties
{
NEW_TYPE_TAG(ZeroEq2cniSubProblem,
             INHERITS_FROM(BoxZeroEqncni, SubDomain));

// Set the problem property
SET_TYPE_PROP(ZeroEq2cniSubProblem, Problem, Dumux::ZeroEq2cniSubProblem<TypeTag>);

// Use the StokesncniCouplingLocalResidual for the computation of the local residual in the ZeroEq domain
SET_TYPE_PROP(ZeroEq2cniSubProblem, LocalResidual, StokesncniCouplingLocalResidual<TypeTag>);

// Set the property for the material parameters by extracting it from the material law.
SET_TYPE_PROP(ZeroEq2cniSubProblem,
              MaterialLawParams,
              typename GET_PROP_TYPE(TypeTag, MaterialLaw)::Params);

// Use the fluid system from the coupled problem
SET_TYPE_PROP(ZeroEq2cniSubProblem,
              FluidSystem,
              typename GET_PROP_TYPE(typename GET_PROP_TYPE(TypeTag, MultiDomainTypeTag), FluidSystem));

// Disable use of mole formulation
SET_BOOL_PROP(ZeroEq2cniSubProblem, UseMoles, false);

// Disable gravity
SET_BOOL_PROP(ZeroEq2cniSubProblem, ProblemEnableGravity, false);

// Enable Navier-Stokes
SET_BOOL_PROP(ZeroEq2cniSubProblem, EnableNavierStokes, true);

// Set the properties for variable inflow BC
NEW_PROP_TAG(FreeFlowSinusVelocityAmplitude);
NEW_PROP_TAG(FreeFlowSinusVelocityPeriod);
SET_SCALAR_PROP(ZeroEq2cniSubProblem, FreeFlowSinusVelocityAmplitude, 0.0);
SET_SCALAR_PROP(ZeroEq2cniSubProblem, FreeFlowSinusVelocityPeriod, 3600.0);
NEW_PROP_TAG(FreeFlowSinusPressureAmplitude);
NEW_PROP_TAG(FreeFlowSinusPressurePeriod);
SET_SCALAR_PROP(ZeroEq2cniSubProblem, FreeFlowSinusPressureAmplitude, 0.0);
SET_SCALAR_PROP(ZeroEq2cniSubProblem, FreeFlowSinusPressurePeriod, 3600.0);
NEW_PROP_TAG(FreeFlowSinusConcentrationAmplitude);
NEW_PROP_TAG(FreeFlowSinusConcentrationPeriod);
SET_SCALAR_PROP(ZeroEq2cniSubProblem, FreeFlowSinusConcentrationAmplitude, 0.0);
SET_SCALAR_PROP(ZeroEq2cniSubProblem, FreeFlowSinusConcentrationPeriod, 3600.0);
NEW_PROP_TAG(FreeFlowSinusTemperatureAmplitude);
NEW_PROP_TAG(FreeFlowSinusTemperaturePeriod);
SET_SCALAR_PROP(ZeroEq2cniSubProblem, FreeFlowSinusTemperatureAmplitude, 0.0);
SET_SCALAR_PROP(ZeroEq2cniSubProblem, FreeFlowSinusTemperaturePeriod, 3600.0);
}

/*!
 * \ingroup ImplicitTestProblems
 * \ingroup TwoPTwoCNIZeroEqTwoCNIModel
 * \brief Non-isothermal two-component ZeroEq subproblem with air flowing
 *        from the left to the right and coupling at the bottom.
 *
 * The free-flow subdomain is sized 0.75m times 0.5m. Dry and hot air is flowing from left (Dirichlet)
 * to right (outflow), at the middle third of the bottom the coupling conditions
 * are applied to all balance equations. They handle the exchange to the porous-medium
 * subdomain.
 *
 * This subproblem uses the \ref ZeroEqncniModel. It is part of a multidomain model and
 * combined with the 2p2cnisubproblem for the porous-medium domain.
 */
template <class TypeTag>
class ZeroEq2cniSubProblem : public ZeroEqProblem<TypeTag>
{
    typedef ZeroEq2cniSubProblem<TypeTag> ThisType;
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
        transportEqIdx = Indices::transportEqIdx, // Index of the transport equation (massfraction)
        energyEqIdx =    Indices::energyEqIdx // Index of the energy equation (temperature)
    };
    enum { // primary variable indices
        pressureIdx = Indices::pressureIdx,
        velocityXIdx = Indices::velocityXIdx,
        velocityYIdx = Indices::velocityYIdx,
        velocityZIdx = Indices::velocityZIdx,
        massOrMoleFracIdx = Indices::massOrMoleFracIdx,
        temperatureIdx = Indices::temperatureIdx
    };
    enum {
        transportCompIdx = Indices::transportCompIdx, // water component index
        phaseCompIdx = Indices::phaseCompIdx // air component index
    };

    typedef typename GET_PROP_TYPE(TypeTag, PrimaryVariables) PrimaryVariables;
    typedef typename GET_PROP_TYPE(TypeTag, BoundaryTypes) BoundaryTypes;

    typedef typename GridView::template Codim<0>::Entity Element;
    typedef typename GridView::template Codim<dim>::Entity Vertex;
    typedef typename GridView::ctype CoordScalar;

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
    ZeroEq2cniSubProblem(TimeManager &timeManager, const GridView &gridView)
        : ParentType(timeManager, gridView)
    {
        bBoxMin_[0] = GET_RUNTIME_PARAM_FROM_GROUP(TypeTag, Scalar, Grid, LowerLeftX);
        bBoxMax_[0] = GET_RUNTIME_PARAM_FROM_GROUP(TypeTag, Scalar, Grid, UpperRightX);
        bBoxMin_[1] = GET_RUNTIME_PARAM_FROM_GROUP(TypeTag, Scalar, Grid, InterfacePosY);
        bBoxMax_[1] = GET_RUNTIME_PARAM_FROM_GROUP(TypeTag, Scalar, Grid, UpperRightY);
        runUpDistanceX1_ = GET_RUNTIME_PARAM_FROM_GROUP(TypeTag, Scalar, Grid, RunUpDistanceX1); // first part of the interface without coupling
        runUpDistanceX2_ = GET_RUNTIME_PARAM_FROM_GROUP(TypeTag, Scalar, Grid, RunUpDistanceX2); // second part of the interface without coupling

        refVelocity_ = GET_RUNTIME_PARAM_FROM_GROUP(TypeTag, Scalar, FreeFlow, RefVelocity);
        refPressure_ = GET_RUNTIME_PARAM_FROM_GROUP(TypeTag, Scalar, FreeFlow, RefPressure);
        refMassfrac_ = GET_RUNTIME_PARAM_FROM_GROUP(TypeTag, Scalar, FreeFlow, RefMassfrac);
        refTemperature_ = GET_RUNTIME_PARAM_FROM_GROUP(TypeTag, Scalar, FreeFlow, RefTemperature);
    }

    // functions have to be overwritten, otherwise they remain uninitialized
    //! \copydoc Dumux::ImplicitProblem::bBoxMin()
    const GlobalPosition &bBoxMin() const
    { return bBoxMin_; }

    //! \copydoc Dumux::ImplicitProblem::bBoxMax()
    const GlobalPosition &bBoxMax() const
    { return bBoxMax_; }

    /*!
     * \name Problem parameters
     */
    // \{

    //! \copydoc Dumux::ImplicitProblem::name()
    const std::string &name() const
    { return GET_RUNTIME_PARAM_FROM_GROUP(TypeTag, std::string, Output, NameFF); }

    // \}

    /*!
     * \name Boundary conditions
     */
    // \{

    //! \copydoc Dumux::ImplicitProblem::boundaryTypesAtPos()
    void boundaryTypesAtPos(BoundaryTypes &values,
                            const GlobalPosition &globalPos) const
    {
        values.setAllDirichlet();


        if (onUpperBoundary_(globalPos))
        {
            values.setNeumann(transportEqIdx);
            values.setDirichlet(temperatureIdx);
        }

        if (onLowerBoundary_(globalPos))
        {
            values.setNeumann(transportEqIdx);
            values.setDirichlet(temperatureIdx);

            if (globalPos[0] > runUpDistanceX1_ - eps_
                && globalPos[0] < runUpDistanceX2_ + eps_)
            {
                values.setAllCouplingDirichlet();
            }
        }

        if (onRightBoundary_(globalPos))
        {
            values.setAllOutflow();

            if (onUpperBoundary_(globalPos) || onLowerBoundary_(globalPos)) // corner points
                values.setAllDirichlet();
        }

        if (onLeftBoundary_(globalPos))
        {
            values.setAllDirichlet();
        }

        // the mass balance has to be of type outflow
        // it does not get a coupling condition, since pn is a condition for stokes
        values.setOutflow(massBalanceIdx);

        if (onRightBoundary_(globalPos))
        {
            values.setAllOutflow();
            values.setDirichlet(pressureIdx);
        }
    }

    //! \copydoc Dumux::ImplicitProblem::dirichletAtPos()
    void dirichletAtPos(PrimaryVariables &values, const GlobalPosition &globalPos) const
    {
        values = 0.0;

        values[velocityXIdx] = xVelocity_(globalPos);
        values[velocityYIdx] = 0.0;
        values[pressureIdx] = refPressure()
                              + 1.189 * this->gravity()[1] * (globalPos[1] - bBoxMin_[1]);
        values[massOrMoleFracIdx] = refMassfrac();
        values[temperatureIdx] = refTemperature();
    }

    //! \copydoc Dumux::ImplicitProblem::neumannAtPos()
    void neumannAtPos(PrimaryVariables &values, const GlobalPosition &globalPos) const
    {
        values = 0.;
    }

    //! \copydoc Dumux::ImplicitProblem::sourceAtPos()
    void sourceAtPos(PrimaryVariables &values,
                     const GlobalPosition &globalPos) const
    {
        values = 0.0;
    }

    //! \copydoc Dumux::ImplicitProblem::initialAtPos()
    void initialAtPos(PrimaryVariables &values, const GlobalPosition &globalPos) const
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
    // Internal method for the initial and Dirichlet conditions
    void initial_(PrimaryVariables &values,
                  const GlobalPosition &globalPos) const
    {
        values[velocityXIdx] = xVelocity_(globalPos);
        values[velocityYIdx] = 0.;

        values[pressureIdx] = refPressure() + 1.189 * this->gravity()[1] * (globalPos[1] - bBoxMin_[1]);
        values[massOrMoleFracIdx] = refMassfrac();
        values[temperatureIdx] = refTemperature();
    }

    // returns the inflow velocity profile
    const Scalar xVelocity_(const GlobalPosition &globalPos) const
    {
        if (onUpperBoundary_(globalPos) || onLowerBoundary_(globalPos))
            return 0.0;
        return refVelocity();
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

    Scalar refVelocity_;
    Scalar refPressure_;
    Scalar refMassfrac_;
    Scalar refTemperature_;

    Scalar runUpDistanceX1_;
    Scalar runUpDistanceX2_;
};
} //end namespace

#endif // DUMUX_ZEROEQTWOCNI_SUBPROBLEM_HH
