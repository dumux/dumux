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
/*!
 * \file
 * \ingroup TwoPNCTests
 * \brief Definition of a problem for water management in PEM fuel cells.
 */
#ifndef DUMUX_FUELCELL_PROBLEM_HH
#define DUMUX_FUELCELL_PROBLEM_HH

#include <dune/grid/yaspgrid.hh>

#include <dumux/discretization/elementsolution.hh>
#include <dumux/discretization/cellcentered/tpfa/properties.hh>
#include <dumux/discretization/box/properties.hh>
#include <dumux/porousmediumflow/2pnc/model.hh>
#include <dumux/porousmediumflow/problem.hh>
#include <dumux/material/fluidsystems/h2on2o2.hh>
#ifdef NONISOTHERMAL
#include <dumux/material/chemistry/electrochemistry/electrochemistryni.hh>
#else
#include <dumux/material/chemistry/electrochemistry/electrochemistry.hh>
#endif
#include "fuelcellspatialparams.hh"

namespace Dumux {

/*!
 * \ingroup TwoPNCTests
 * \brief Definition of a problem for water management in PEM fuel cells.
 */
template <class TypeTag>
class FuelCellProblem;

namespace Properties {
#ifdef NONISOTHERMAL
NEW_TYPE_TAG(FuelCell, INHERITS_FROM(TwoPNCNI));
NEW_TYPE_TAG(FuelCellNIBox, INHERITS_FROM(BoxModel, FuelCell));
#else
NEW_TYPE_TAG(FuelCell, INHERITS_FROM(TwoPNC));
NEW_TYPE_TAG(FuelCellBox, INHERITS_FROM(BoxModel, FuelCell));
NEW_TYPE_TAG(FuelCellCCTpfa, INHERITS_FROM(CCTpfaModel, FuelCell));
#endif

// Set the grid type
SET_TYPE_PROP(FuelCell, Grid, Dune::YaspGrid<2>);

// Set the problem property
SET_TYPE_PROP(FuelCell, Problem, FuelCellProblem<TypeTag>);

// Set the spatial parameters
SET_PROP(FuelCell, SpatialParams)
{
    using FVGridGeometry = typename GET_PROP_TYPE(TypeTag, FVGridGeometry);
    using Scalar = typename GET_PROP_TYPE(TypeTag, Scalar);
    using type = FuelCellSpatialParams<FVGridGeometry, Scalar>;
};

// Set the primary variable combination for the 2pnc model
SET_PROP(FuelCell, Formulation)
{ static constexpr auto value = TwoPFormulation::p1s0; };

// Set fluid configuration
SET_PROP(FuelCell, FluidSystem)
{
private:
    using Scalar = typename GET_PROP_TYPE(TypeTag, Scalar);
public:
    using type = FluidSystems::H2ON2O2<Scalar>;
};
} // end namespace Properties

/*!
 * \ingroup TwoPNCModel
 * \ingroup ImplicitTestProblems
 * \brief Problem or water management in PEM fuel cells.
 *
 * To run the simulation execute the following line in shell:
 * <tt>./test_box2pnc</tt>
 */
template <class TypeTag>
class FuelCellProblem : public PorousMediumFlowProblem<TypeTag>
{
    using ParentType = PorousMediumFlowProblem<TypeTag>;

    using Scalar = typename GET_PROP_TYPE(TypeTag, Scalar);
    using Indices = typename GET_PROP_TYPE(TypeTag, ModelTraits)::Indices;
    using FluidSystem = typename GET_PROP_TYPE(TypeTag, FluidSystem);
    using BoundaryTypes = typename GET_PROP_TYPE(TypeTag, BoundaryTypes);
    using PrimaryVariables = typename GET_PROP_TYPE(TypeTag, PrimaryVariables);
    using NumEqVector = typename GET_PROP_TYPE(TypeTag, NumEqVector);
    using ElementVolumeVariables = typename GET_PROP_TYPE(TypeTag, GridVolumeVariables)::LocalView;
    using FVElementGeometry = typename GET_PROP_TYPE(TypeTag, FVGridGeometry)::LocalView;
    using SubControlVolume = typename FVElementGeometry::SubControlVolume;
    using GridView = typename GET_PROP_TYPE(TypeTag, GridView);
    using Element = typename GridView::template Codim<0>::Entity;
    using FVGridGeometry = typename GET_PROP_TYPE(TypeTag, FVGridGeometry);
    using SolutionVector = typename GET_PROP_TYPE(TypeTag, SolutionVector);
    using VolumeVariables = typename GET_PROP_TYPE(TypeTag, VolumeVariables);
    // Select the electrochemistry method
#ifdef NONISOTHERMAL
    using ElectroChemistry = typename Dumux::ElectroChemistryNI<Scalar, Indices, FluidSystem, FVGridGeometry, ElectroChemistryModel::Ochs>;
#else
    using ElectroChemistry = typename Dumux::ElectroChemistry<Scalar, Indices, FluidSystem, FVGridGeometry, ElectroChemistryModel::Ochs>;
#endif
    static constexpr int dim = GridView::dimension;
    static constexpr int dimWorld = GridView::dimensionworld;
    static constexpr bool isBox = FVGridGeometry::discMethod == DiscretizationMethod::box;
    using GlobalPosition = typename SubControlVolume::GlobalPosition;

    enum { dofCodim = isBox ? dim : 0 };
public:
    /*!
     * \brief The constructor
     *
     * \param timeManager The time manager
     * \param gridView The grid view
     */
    FuelCellProblem(std::shared_ptr<const FVGridGeometry> fvGridGeometry)
    : ParentType(fvGridGeometry)
    {
        nTemperature_       = getParam<int>("Problem.NTemperature");
        nPressure_          = getParam<int>("Problem.NPressure");
        pressureLow_        = getParam<Scalar>("Problem.PressureLow");
        pressureHigh_       = getParam<Scalar>("Problem.PressureHigh");
        temperatureLow_     = getParam<Scalar>("Problem.TemperatureLow");
        temperatureHigh_    = getParam<Scalar>("Problem.TemperatureHigh");
        temperature_        = getParam<Scalar>("Problem.InitialTemperature");

        name_               = getParam<std::string>("Problem.Name");

        pO2Inlet_           = getParam<Scalar>("ElectroChemistry.pO2Inlet");

        FluidSystem::init(/*Tmin=*/temperatureLow_,
                          /*Tmax=*/temperatureHigh_,
                          /*nT=*/nTemperature_,
                          /*pmin=*/pressureLow_,
                          /*pmax=*/pressureHigh_,
                          /*np=*/nPressure_);
    }

    /*!
     * \name Problem parameters
     */

    /*!
     * \brief The problem name.
     *
     * This is used as a prefix for files generated by the simulation.
     */
    const std::string& name() const
    { return name_; }

    /*!
     * \brief Returns the temperature within the domain.
     *
     * This problem assumes a temperature of 10 degrees Celsius.
     */
    Scalar temperature() const
    { return temperature_; }

    //! \copydoc Dumux::FVProblem::source()
    NumEqVector source(const Element &element,
                   const FVElementGeometry& fvGeometry,
                   const ElementVolumeVariables& elemVolVars,
                   const SubControlVolume &scv) const
    {
        NumEqVector values(0.0);
        const auto& globalPos = scv.dofPosition();

        //reaction sources from electro chemistry
        if(inReactionLayer_(globalPos))
        {
            const auto& volVars = elemVolVars[scv];
            auto currentDensity = ElectroChemistry::calculateCurrentDensity(volVars);
            ElectroChemistry::reactionSource(values, currentDensity);
        }

        return values;
    }


    /*!
     * \name Boundary conditions
     */
    // \{

    /*!
     * \brief Specifies which kind of boundary condition should be
     *        used for which equation on a given boundary segment
     *
     * \param globalPos The global position
     */
    BoundaryTypes boundaryTypesAtPos(const GlobalPosition &globalPos) const
    {
        BoundaryTypes bcTypes;
        bcTypes.setAllNeumann();

        if (onUpperBoundary_(globalPos)){
            bcTypes.setAllDirichlet();
        }
        return bcTypes;
    }

    /*!
     * \brief Evaluates the boundary conditions for a Dirichlet
     *        boundary segment
     *
     * \param globalPos The global position
     */
    PrimaryVariables dirichletAtPos(const GlobalPosition &globalPos) const
    {
        auto priVars = initial_(globalPos);

        if(onUpperBoundary_(globalPos))
        {
            Scalar pn = 1.0e5;
            priVars[Indices::pressureIdx] = pn;
            priVars[Indices::switchIdx] = 0.3;//Sw for bothPhases
            priVars[Indices::switchIdx+1] = pO2Inlet_/4.315e9; //moleFraction xlO2 for bothPhases
#ifdef NONISOTHERMAL
            priVars[Indices::temperatureIdx] = 293.15;
#endif
        }

        return priVars;
    }

    /*!
     * \name Volume terms
     */


    /*!
     * \brief Evaluates the initial values for a control volume
     *
     * \param globalPos The global position
     */
    PrimaryVariables initialAtPos(const GlobalPosition &globalPos) const
    { return initial_(globalPos); }


    /*!
     * \brief Adds additional VTK output data to the VTKWriter. Function is called by the output module on every write.
     */
    template<class VTKWriter>
    void addVtkFields(VTKWriter& vtk)
    {
        const auto& gridView = this->fvGridGeometry().gridView();
        currentDensity_.resize(gridView.size(dofCodim));
        reactionSourceH2O_.resize(gridView.size(dofCodim));
        reactionSourceO2_.resize(gridView.size(dofCodim));
        Kxx_.resize(gridView.size(dofCodim));
        Kyy_.resize(gridView.size(dofCodim));

        vtk.addField(currentDensity_, "currentDensity [A/cm^2]");
        vtk.addField(reactionSourceH2O_, "reactionSourceH2O [mol/(sm^2)]");
        vtk.addField(reactionSourceO2_, "reactionSourceO2 [mol/(sm^2)]");
        vtk.addField(Kxx_, "Kxx");
        vtk.addField(Kyy_, "Kyy");
    }

    void updateVtkFields(const SolutionVector& curSol)
    {
        for (const auto& element : elements(this->fvGridGeometry().gridView()))
        {
            auto elemSol = elementSolution(element, curSol, this->fvGridGeometry());

            auto fvGeometry = localView(this->fvGridGeometry());
            fvGeometry.bindElement(element);

            for (auto&& scv : scvs(fvGeometry))
            {
                VolumeVariables volVars;
                volVars.update(elemSol, *this, element, scv);
                const auto& globalPos = scv.dofPosition();
                const auto dofIdxGlobal = scv.dofIndex();

                if(inReactionLayer_(globalPos))
                {
                    //reactionSource Output
                    PrimaryVariables source;
                    auto i = ElectroChemistry::calculateCurrentDensity(volVars);
                    ElectroChemistry::reactionSource(source, i);

                    reactionSourceH2O_[dofIdxGlobal] = source[Indices::conti0EqIdx + FluidSystem::H2OIdx];
                    reactionSourceO2_[dofIdxGlobal] = source[Indices::conti0EqIdx + FluidSystem::O2Idx];

                    //Current Output in A/cm^2
                    currentDensity_[dofIdxGlobal] = i/10000;
                }
                else
                {
                    reactionSourceH2O_[dofIdxGlobal] = 0.0;
                    reactionSourceO2_[dofIdxGlobal] = 0.0;
                    currentDensity_[dofIdxGlobal] = 0.0;
                }
                Kxx_[dofIdxGlobal] = volVars.permeability()[0][0];
                Kyy_[dofIdxGlobal] = volVars.permeability()[1][1];
            }
        }
    }

private:

    PrimaryVariables initial_(const GlobalPosition &globalPos) const
    {
        PrimaryVariables priVars(0.0);
        priVars.setState(Indices::bothPhases);

        Scalar pn = 1.0e5;
        priVars[Indices::pressureIdx] = pn;
        priVars[Indices::switchIdx] = 0.3;//Sw for bothPhases
        priVars[Indices::switchIdx+1] = pO2Inlet_/4.315e9; //moleFraction xlO2 for bothPhases
#ifdef NONISOTHERMAL
        priVars[Indices::temperatureIdx] = 293.15;
#endif

        return priVars;
    }

    bool onUpperBoundary_(const GlobalPosition &globalPos) const
    { return globalPos[1] > this->fvGridGeometry().bBoxMax()[1] - eps_; }

    bool inReactionLayer_(const GlobalPosition& globalPos) const
    { return globalPos[1] < 0.1*(this->fvGridGeometry().bBoxMax()[1] - this->fvGridGeometry().bBoxMin()[1]) + eps_; }

    Scalar temperature_;
    static constexpr Scalar eps_ = 1e-6;
    int nTemperature_;
    int nPressure_;
    std::string name_ ;
    Scalar pressureLow_, pressureHigh_;
    Scalar temperatureLow_, temperatureHigh_;
    Scalar pO2Inlet_;
    std::vector<double> currentDensity_;
    std::vector<double> reactionSourceH2O_;
    std::vector<double> reactionSourceO2_;
    std::vector<double> Kxx_;
    std::vector<double> Kyy_;
};

} //end namespace Dumux

#endif
