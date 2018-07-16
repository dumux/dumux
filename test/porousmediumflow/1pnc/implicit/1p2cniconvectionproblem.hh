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
 * \ingroup OnePNCTests
 * \brief Test for the OnePNCModel in combination with the NI model for a conduction problem:
 * The simulation domain is a tube where with an elevated temperature on the left hand side.
 */
#ifndef DUMUX_1P2CNI_CONVECTION_TEST_PROBLEM_HH
#define DUMUX_1P2CNI_CONVECTION_TEST_PROBLEM_HH

#if HAVE_UG
#include <dune/grid/uggrid.hh>
#endif
#include <dune/grid/yaspgrid.hh>

#include <dumux/discretization/elementsolution.hh>
#include <dumux/discretization/cellcentered/tpfa/properties.hh>
#include <dumux/discretization/cellcentered/mpfa/properties.hh>
#include <dumux/discretization/box/properties.hh>
#include <dumux/porousmediumflow/1pnc/model.hh>
#include <dumux/porousmediumflow/problem.hh>

#include <dumux/material/fluidsystems/h2on2.hh>
#include <dumux/material/components/h2o.hh>
#include "1pnctestspatialparams.hh"

namespace Dumux
{

template <class TypeTag>
class OnePTwoCNIConvectionProblem;

namespace Properties
{
NEW_TYPE_TAG(OnePTwoCNIConvectionTypeTag, INHERITS_FROM(OnePNCNI));
NEW_TYPE_TAG(OnePTwoCNIConvectionCCTpfaTypeTag, INHERITS_FROM(CCTpfaModel, OnePTwoCNIConvectionTypeTag));
NEW_TYPE_TAG(OnePTwoCNIConvectionCCMpfaTypeTag, INHERITS_FROM(CCMpfaModel, OnePTwoCNIConvectionTypeTag));
NEW_TYPE_TAG(OnePTwoCNIConvectionBoxTypeTag, INHERITS_FROM(BoxModel, OnePTwoCNIConvectionTypeTag));

// Set the grid type
#if HAVE_UG
SET_TYPE_PROP(OnePTwoCNIConvectionTypeTag, Grid, Dune::UGGrid<2>);
#else
SET_TYPE_PROP(OnePTwoCNIConvectionTypeTag, Grid, Dune::YaspGrid<2>);
#endif

// Set the problem property
SET_TYPE_PROP(OnePTwoCNIConvectionTypeTag, Problem, OnePTwoCNIConvectionProblem<TypeTag>);

// Set fluid configuration
SET_TYPE_PROP(OnePTwoCNIConvectionTypeTag,
              FluidSystem,
              FluidSystems::H2ON2<typename GET_PROP_TYPE(TypeTag, Scalar)>);

// Set the spatial parameters
SET_TYPE_PROP(OnePTwoCNIConvectionTypeTag, SpatialParams, OnePNCTestSpatialParams<TypeTag>);

// Define whether mole(true) or mass (false) fractions are used
SET_BOOL_PROP(OnePTwoCNIConvectionTypeTag, UseMoles, true);
}


/*!
 * \ingroup OnePNCTests
 * \brief Test for the OnePTwoCModel in combination with the NI model for a convection problem:
 * The simulation domain is a tube where water with an elevated temperature is injected
 * at a constant rate on the left hand side.
 *
 * Initially the domain is fully saturated with water at a constant temperature.
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
class OnePTwoCNIConvectionProblem : public PorousMediumFlowProblem<TypeTag>
{
    using ParentType = PorousMediumFlowProblem<TypeTag>;

    using Scalar = typename GET_PROP_TYPE(TypeTag, Scalar);
    using Indices = typename GET_PROP_TYPE(TypeTag, ModelTraits)::Indices;
    using GridView = typename GET_PROP_TYPE(TypeTag, GridView);
    using FluidSystem = typename GET_PROP_TYPE(TypeTag, FluidSystem);
    using BoundaryTypes = typename GET_PROP_TYPE(TypeTag, BoundaryTypes);
    using PrimaryVariables = typename GET_PROP_TYPE(TypeTag, PrimaryVariables);
    using FVGridGeometry = typename GET_PROP_TYPE(TypeTag, FVGridGeometry);
    using FVElementGeometry = typename GET_PROP_TYPE(TypeTag, FVGridGeometry)::LocalView;
    using NumEqVector = typename GET_PROP_TYPE(TypeTag, NumEqVector);
    using ElementVolumeVariables = typename GET_PROP_TYPE(TypeTag, GridVolumeVariables)::LocalView;
    using SubControlVolumeFace = typename FVElementGeometry::SubControlVolumeFace;
    using Element = typename GridView::template Codim<0>::Entity;
    using SolutionVector = typename GET_PROP_TYPE(TypeTag, SolutionVector);
    using VolumeVariables = typename GET_PROP_TYPE(TypeTag, VolumeVariables);
    using IapwsH2O = Components::H2O<Scalar>;

    // copy some indices for convenience
    enum
    {
        // indices of the primary variables
        pressureIdx = Indices::pressureIdx,
        temperatureIdx = Indices::temperatureIdx,

        // equation indices
        contiH2OIdx = Indices::conti0EqIdx + FluidSystem::H2OIdx,
        contiN2Idx = Indices::conti0EqIdx + FluidSystem::N2Idx,
        energyEqIdx = Indices::energyEqIdx
    };

    //! property that defines whether mole or mass fractions are used
    static constexpr bool useMoles = GET_PROP_VALUE(TypeTag, UseMoles);
    static const int dimWorld = GridView::dimensionworld;
    using GlobalPosition = typename SubControlVolumeFace::GlobalPosition;

public:
    OnePTwoCNIConvectionProblem(std::shared_ptr<const FVGridGeometry> fvGridGeometry)
    : ParentType(fvGridGeometry)
    {
        //initialize fluid system
        FluidSystem::init();

        // stating in the console whether mole or mass fractions are used
        if(useMoles)
            std::cout<<"problem uses mole fractions"<<std::endl;
        else
            std::cout<<"problem uses mass fractions"<<std::endl;

        temperatureExact_.resize(fvGridGeometry->numDofs(), 290.0);

        darcyVelocity_ = getParam<Scalar>("Problem.DarcyVelocity");

        temperatureHigh_ = 291.;
        temperatureLow_ = 290.;
        pressureHigh_ = 2e5;
        pressureLow_ = 1e5;
    }

    //! get the analytical temperature
    const std::vector<Scalar>& getExactTemperature()
    {
        return temperatureExact_;
    }

    //! udpate the analytical temperature
    void updateExactTemperature(const SolutionVector& curSol, Scalar time)
    {
        const auto someElement = *(elements(this->fvGridGeometry().gridView()).begin());

        auto someElemSol = elementSolution(someElement, curSol, this->fvGridGeometry());
        const auto someInitSol = initialAtPos(someElement.geometry().center());

        auto fvGeometry = localView(this->fvGridGeometry());
        fvGeometry.bindElement(someElement);
        const auto someScv = *(scvs(fvGeometry).begin());

        VolumeVariables volVars;
        volVars.update(someElemSol, *this, someElement, someScv);

        const auto porosity = this->spatialParams().porosity(someElement, someScv, someElemSol);
        const auto densityW = volVars.density();
        const auto heatCapacityW = IapwsH2O::liquidHeatCapacity(someInitSol[temperatureIdx], someInitSol[pressureIdx]);
        const auto storageW =  densityW*heatCapacityW*porosity;
        const auto densityS = volVars.solidDensity();
        const auto heatCapacityS = volVars.solidHeatCapacity();
        const auto storageTotal = storageW + densityS*heatCapacityS*(1 - porosity);
        std::cout << "storage: " << storageTotal << '\n';

        using std::max;
        time = max(time, 1e-10);
        const Scalar retardedFrontVelocity = darcyVelocity_*storageW/storageTotal/porosity;
        std::cout << "retarded velocity: " << retardedFrontVelocity << '\n';

        for (const auto& element : elements(this->fvGridGeometry().gridView()))
        {
            auto fvGeometry = localView(this->fvGridGeometry());
            fvGeometry.bindElement(element);
            for (auto&& scv : scvs(fvGeometry))
            {
                auto dofIdxGlobal = scv.dofIndex();
                auto dofPosition = scv.dofPosition();
                temperatureExact_[dofIdxGlobal] = (dofPosition[0] < retardedFrontVelocity*time) ? temperatureHigh_ : temperatureLow_;
            }
        }
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

        if(globalPos[0] > this->fvGridGeometry().bBoxMax()[0] - eps_)
        {
            values.setAllDirichlet();
        }
        else
        {
            values.setAllNeumann();
        }

        return values;
    }

    /*!
     * \brief Evaluate the boundary conditions for a dirichlet
     *        boundary segment.
     *
     * \param globalPos The position for which the bc type should be evaluated
     */
    PrimaryVariables dirichletAtPos(const GlobalPosition &globalPos) const
    {
        PrimaryVariables values = initial_(globalPos);
        return values;
    }

    /*!
     * \brief Evaluate the boundary conditions for a neumann
     *        boundary segment.
     *
     * This is the method for the case where the Neumann condition is
     * potentially solution dependent and requires some quantities that
     * are specific to the fully-implicit method.
     *
     * \param values The neumann values for the conservation equations in units of
     *                 \f$ [ \textnormal{unit of conserved quantity} / (m^2 \cdot s )] \f$
     * \param element The finite element
     * \param fvGeometry The finite-volume geometry
     * \param elemVolVars All volume variables for the element
     * \param scvf The sub control volume face
     *
     * For this method, the \a values parameter stores the flux
     * in normal direction of each phase. Negative values mean influx.
     * E.g. for the mass balance that would the mass flux in \f$ [ kg / (m^2 \cdot s)] \f$.
     */
    NumEqVector neumann(const Element& element,
                           const FVElementGeometry& fvGeometry,
                           const ElementVolumeVariables& elemVolVars,
                           const SubControlVolumeFace& scvf) const
    {
        NumEqVector flux(0.0);
        const auto& globalPos = scvf.ipGlobal();
        const auto& scv = fvGeometry.scv(scvf.insideScvIdx());

        if(globalPos[0] < eps_)
        {
             flux[contiH2OIdx] = -darcyVelocity_*elemVolVars[scv].molarDensity();
             flux[contiN2Idx] = -darcyVelocity_*elemVolVars[scv].molarDensity()*elemVolVars[scv].moleFraction(0, FluidSystem::N2Idx);
             flux[energyEqIdx] = -darcyVelocity_
                                 *elemVolVars[scv].density()
                                 *IapwsH2O::liquidEnthalpy(temperatureHigh_, elemVolVars[scv].pressure());
        }

        return flux;
    }

    // \}

    /*!
     * \name Volume terms
     */
    // \{

    /*!
     * \brief Evaluate the source term for all phases within a given
     *        sub-control-volume.
     *
     * For this method, the \a priVars parameter stores the rate mass
     * of a component is generated or annihilate per volume
     * unit. Positive values mean that mass is created, negative ones
     * mean that it vanishes.
     *
     * The units must be according to either using mole or mass fractions. (mole/(m^3*s) or kg/(m^3*s))
     */
    NumEqVector sourceAtPos(const GlobalPosition &globalPos) const
    { return NumEqVector(0.0); }

    /*!
     * \brief Evaluate the initial value for a control volume.
     *
     * \param globalPos The position for which the initial condition should be evaluated
     *
     * For this method, the \a values parameter stores primary
     * variables.
     */
    PrimaryVariables initialAtPos(const GlobalPosition &globalPos) const
    { return initial_(globalPos); }

    // \}
private:

    // the internal method for the initial condition
    PrimaryVariables initial_(const GlobalPosition &globalPos) const
    {
        PrimaryVariables priVars;
        priVars[pressureIdx] = pressureLow_; // initial condition for the pressure
        priVars[FluidSystem::N2Idx] = 1e-10;  // initial condition for the N2 molefraction
        priVars[temperatureIdx] = temperatureLow_;
        return priVars;
    }
        static constexpr Scalar eps_ = 1e-6;
        Scalar temperatureHigh_;
        Scalar temperatureLow_;
        Scalar pressureHigh_;
        Scalar pressureLow_;
        Scalar darcyVelocity_;

        std::vector<Scalar> temperatureExact_;
    };

} //end namespace Dumux

#endif
