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
 * \ingroup RichardsTests
 * \brief Test for the RichardsModel in combination with the NI model for a convection problem:
 * The simulation domain is a tube where water with an elevated temperature is injected
 * at a constant rate on the left hand side.
 */
#ifndef DUMUX_RICHARDS_CONVECTION_PROBLEM_HH
#define DUMUX_RICHARDS_CONVECTION_PROBLEM_HH

#include <math.h>

#include <dumux/discretization/cellcentered/tpfa/properties.hh>
#include <dumux/discretization/box/properties.hh>

#include <dumux/porousmediumflow/problem.hh>
#include <dumux/porousmediumflow/richards/model.hh>
#include <dumux/material/fluidmatrixinteractions/2p/thermalconductivitysomerton.hh>
#include <dumux/material/fluidsystems/h2on2.hh>
#include "richardsnispatialparams.hh"

namespace Dumux
{
/**
 * \ingroup RichardsTests
 * \brief Test for the RichardsModel in combination with the NI model for a convection problem:
 * The simulation domain is a tube where water with an elevated temperature is injected
 * at a constant rate on the left hand side.
 */
template <class TypeTag>
class RichardsNIConvectionProblem;

namespace Properties
{
NEW_TYPE_TAG(RichardsNIConvectionProblem, INHERITS_FROM(RichardsNI, RichardsNISpatialParams));
NEW_TYPE_TAG(RichardsNIConvectionBoxProblem, INHERITS_FROM(BoxModel, RichardsNIConvectionProblem));
NEW_TYPE_TAG(RichardsNIConvectionCCProblem, INHERITS_FROM(CCTpfaModel, RichardsNIConvectionProblem));

// Set the grid type
SET_TYPE_PROP(RichardsNIConvectionProblem, Grid, Dune::YaspGrid<2>);

// Set the problem property
SET_TYPE_PROP(RichardsNIConvectionProblem, Problem,
              RichardsNIConvectionProblem<TypeTag>);

// Set the fluid system
SET_TYPE_PROP(RichardsNIConvectionProblem, FluidSystem, FluidSystems::H2ON2<typename GET_PROP_TYPE(TypeTag, Scalar), false>);

// Set the spatial parameters
SET_TYPE_PROP(RichardsNIConvectionProblem,
              SpatialParams,
              RichardsNISpatialParams<TypeTag>);
}


/*!
 * \ingroup RichardsModel
 * \ingroup ImplicitTestProblems
 *
 * \brief Test for the RichardsModel in combination with the NI model for a convection problem:
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
 * This problem uses the \ref RichardsModel and \ref NIModel model.
 *
 * To run the simulation execute the following line in shell: <br>
 * <tt>./test_boxrichardsniconvection -ParameterFile ./test_boxrichardsniconvection.input</tt> or <br>
 * <tt>./test_ccrichardsniconvection -ParameterFile ./test_ccrichardsniconvection.input</tt>
 */
template <class TypeTag>
class RichardsNIConvectionProblem : public PorousMediumFlowProblem<TypeTag>
{
    using ParentType = PorousMediumFlowProblem<TypeTag>;

    using GridView = typename GET_PROP_TYPE(TypeTag, GridView);
    using Scalar = typename GET_PROP_TYPE(TypeTag, Scalar);
    using FVElementGeometry = typename GET_PROP_TYPE(TypeTag, FVGridGeometry)::LocalView;
    using FVGridGeometry = typename GET_PROP_TYPE(TypeTag, FVGridGeometry);
    using PrimaryVariables = typename GET_PROP_TYPE(TypeTag, PrimaryVariables);
    using FluidSystem = typename GET_PROP_TYPE(TypeTag, FluidSystem);
    using BoundaryTypes = typename GET_PROP_TYPE(TypeTag, BoundaryTypes);
    using ThermalConductivityModel = typename GET_PROP_TYPE(TypeTag, ThermalConductivityModel);
    using VolumeVariables = typename GET_PROP_TYPE(TypeTag, VolumeVariables);
    using ElementVolumeVariables = typename GET_PROP_TYPE(TypeTag, ElementVolumeVariables);
    using ElementSolutionVector = typename GET_PROP_TYPE(TypeTag, ElementSolutionVector);
    using SolutionVector = typename GET_PROP_TYPE(TypeTag, SolutionVector);
    using SubControlVolumeFace = typename FVElementGeometry::SubControlVolumeFace;
    using IapwsH2O = H2O<Scalar>;

    // copy some indices for convenience
    using Indices = typename GET_PROP_TYPE(TypeTag, Indices);
    enum { dimWorld = GridView::dimensionworld };

    enum {
        pressureIdx = Indices::pressureIdx,
        wPhaseOnly = Indices::wPhaseOnly,
        wPhaseIdx = Indices::wPhaseIdx,
        temperatureIdx = Indices::temperatureIdx
    };

    using Element = typename GridView::template Codim<0>::Entity;

    using GlobalPosition = Dune::FieldVector<Scalar, dimWorld>;

public:
    RichardsNIConvectionProblem(std::shared_ptr<const FVGridGeometry> fvGridGeometry)
    : ParentType(fvGridGeometry)
    {
        //initialize fluid system
        FluidSystem::init();

        name_ =  getParam<std::string>("Problem.Name");
        darcyVelocity_ =  getParam<Scalar>("Problem.DarcyVelocity");
        temperatureHigh_ = 291.;
        temperatureLow_ = 290.;
        pressureHigh_ = 2e5;
        pressureLow_ = 1e5;
        temperatureExact_.resize(fvGridGeometry->numDofs());
    }


   /*!
     * \brief Append all quantities of interest which can be derived
     *        from the solution of the current time step to the VTK
     *        writer.
     */
    //! get the analytical temperature
    const std::vector<Scalar>& getExactTemperature()
    {
        return temperatureExact_;
    }

  //! udpate the analytical temperature
    void updateExactTemperature(const SolutionVector& curSol, Scalar time)
    {
        const auto someElement = *(elements(this->fvGridGeometry().gridView()).begin());

        ElementSolutionVector someElemSol(someElement, curSol, this->fvGridGeometry());
        const auto someInitSol = initialAtPos(someElement.geometry().center());

        auto fvGeometry = localView(this->fvGridGeometry());
        fvGeometry.bindElement(someElement);
        const auto someScv = *(scvs(fvGeometry).begin());

        VolumeVariables volVars;
        volVars.update(someElemSol, *this, someElement, someScv);

        const auto porosity = this->spatialParams().porosity(someElement, someScv, someElemSol);
        const auto densityW = volVars.density(wPhaseIdx);
        const auto heatCapacityW = IapwsH2O::liquidHeatCapacity(someInitSol[temperatureIdx], someInitSol[pressureIdx]);
        const auto densityS = this->spatialParams().solidDensity(someElement, someScv, someElemSol);
        const auto heatCapacityS = this->spatialParams().solidHeatCapacity(someElement, someScv, someElemSol);
        const auto storage = densityW*heatCapacityW*porosity + densityS*heatCapacityS*(1 - porosity);
        const auto effectiveThermalConductivity = ThermalConductivityModel::effectiveThermalConductivity(volVars, this->spatialParams(),
                                                                                                         someElement, fvGeometry, someScv);
        using std::max;
        time = max(time, 1e-10);
        for (const auto& element : elements(this->fvGridGeometry().gridView()))
        {
            auto fvGeometry = localView(this->fvGridGeometry());
            fvGeometry.bindElement(element);

            for (auto&& scv : scvs(fvGeometry))
            {
               auto globalIdx = scv.dofIndex();
               const auto& globalPos = scv.dofPosition();
               using std::erf;
               using std::sqrt;
               temperatureExact_[globalIdx] = temperatureHigh_ + (someInitSol[temperatureIdx] - temperatureHigh_)
                                              *erf(0.5*sqrt(globalPos[0]*globalPos[0]*storage/time/effectiveThermalConductivity));

            }
        }
    }
    /*!
     * \name Problem parameters
     */
    // \{

    /*!
     * \brief The problem name.
     *
     * This is used as a prefix for files generated by the simulation.
     */
    const std::string& name() const
    {
        return name_;
    }

    // \}

    /*!
     * \name Boundary conditions
     */
    // \{

    /*!
     * \brief Specifies which kind of boundary condition should be
     *        used for which equation on a given boundary segment.
     *
     * \param globalPos The position for which the boundary type is set
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
     * \param values The dirichlet values for the primary variables
     * \param globalPos The position for which the bc type should be evaluated
     *
     * For this method, the \a values parameter stores primary variables.
     */
    PrimaryVariables dirichletAtPos(const GlobalPosition &globalPos) const
    {
        return initial_(globalPos);
    }

    /*!
     * \brief Evaluate the boundary conditions for a neumann
     *        boundary segment.
     *
     * \param element The finite element
     * \param fvGeometry The finite-volume geometry in the box scheme
     * \param elemVolVars The element volume variables
     * \param scvf The subcontrolvolume face
     *  Negative values mean influx.
     */
    PrimaryVariables neumann(const Element &element,
                             const FVElementGeometry& fvGeometry,
                             const ElementVolumeVariables& elemVolVars,
                             const SubControlVolumeFace& scvf) const
    {
        PrimaryVariables values(0.0);
        const auto globalPos = scvf.ipGlobal();

        if(globalPos[0] < eps_)
        {
            values[pressureIdx] = -darcyVelocity_*elemVolVars[scvf.insideScvIdx()].density(wPhaseIdx);
            values[temperatureIdx] = -darcyVelocity_*elemVolVars[scvf.insideScvIdx()].density(wPhaseIdx)
                                     *IapwsH2O::liquidEnthalpy(temperatureHigh_, elemVolVars[scvf.insideScvIdx()].pressure(wPhaseIdx));
        }
        return values;
    }

    // \}

    /*!
     * \name Volume terms
     */
    // \{


    /*!
     * \brief Returns the reference pressure [Pa] of the non-wetting
     *        fluid phase within a finite volume
     *
     * This problem assumes a constant reference pressure of 1 bar.
     *
     * \param element The DUNE Codim<0> entity which intersects with
     *                the finite volume in question
     * \param fvGeometry The finite volume geometry of the element
     * \param scvIdx The sub control volume index inside the finite
     *               volume geometry
     */
    Scalar nonWettingReferencePressure() const
    { return 1e5; };

    /*!
     * \brief Evaluate the initial value for a control volume.
     *
     * \param values The initial values for the primary variables
     * \param globalPos The position for which the initial condition should be evaluated
     *
     * For this method, the \a values parameter stores primary
     * variables.
     */
   PrimaryVariables initialAtPos(const GlobalPosition &globalPos) const
    {
        return initial_(globalPos);
    }

    // \}

private:
    PrimaryVariables initial_(const GlobalPosition &globalPos) const
    {
        PrimaryVariables priVars(0.0);
        priVars.setState(wPhaseOnly);
        priVars[pressureIdx] = pressureLow_; // initial condition for the pressure
        priVars[temperatureIdx] = temperatureLow_;
        return priVars;
    }

    Scalar temperatureHigh_;
    Scalar temperatureLow_;
    Scalar pressureHigh_;
    Scalar pressureLow_;
    Scalar darcyVelocity_;
    static constexpr Scalar eps_ = 1e-6;
    std::string name_;
    std::vector<Scalar> temperatureExact_;
};

} //end namespace
#endif // DUMUX_RICHARDSNI_CONVECTION_PROBLEM_HH
