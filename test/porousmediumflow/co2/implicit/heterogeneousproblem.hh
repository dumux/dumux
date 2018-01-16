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
 * \ingroup CO2Tests
 * \brief Definition of a problem, where CO2 is injected in a reservoir.
 */
#ifndef DUMUX_HETEROGENEOUS_PROBLEM_HH
#define DUMUX_HETEROGENEOUS_PROBLEM_HH

#include <dumux/discretization/cellcentered/tpfa/properties.hh>
#include <dumux/discretization/box/properties.hh>

#include <dumux/porousmediumflow/problem.hh>
#include <dumux/porousmediumflow/co2/model.hh>

#include <dumux/material/fluidsystems/brineco2.hh>
#include <dumux/discretization/box/scvftoscvboundarytypes.hh>

#include "heterogeneousspatialparameters.hh"
#include "heterogeneousco2tables.hh"

// per default use isothermal model
#ifndef ISOTHERMAL
#define ISOTHERMAL 1
#endif

namespace Dumux
{
/*!
 * \ingroup CO2Tests
 * \brief Definition of a problem, where CO2 is injected in a reservoir.
 */
template <class TypeTag>
class HeterogeneousProblem;

namespace Properties
{
NEW_TYPE_TAG(HeterogeneousTypeTag, INHERITS_FROM(TwoPTwoCCO2, HeterogeneousSpatialParams));
NEW_TYPE_TAG(HeterogeneousBoxTypeTag, INHERITS_FROM(BoxModel, HeterogeneousTypeTag));
NEW_TYPE_TAG(HeterogeneousCCTpfaTypeTag, INHERITS_FROM(CCTpfaModel, HeterogeneousTypeTag));

// Set the grid type
SET_TYPE_PROP(HeterogeneousTypeTag, Grid, Dune::ALUGrid<2, 2, Dune::cube, Dune::nonconforming>);

// Set the problem property
SET_TYPE_PROP(HeterogeneousTypeTag, Problem, HeterogeneousProblem<TypeTag>);

// Set fluid configuration
SET_TYPE_PROP(HeterogeneousTypeTag, FluidSystem, FluidSystems::BrineCO2<typename GET_PROP_TYPE(TypeTag, Scalar),
                                                                        HeterogeneousCO2Tables::CO2Tables>);

// Use Moles
SET_BOOL_PROP(HeterogeneousTypeTag, UseMoles, false);

#if !ISOTHERMAL
NEW_TYPE_TAG(HeterogeneousNITypeTag, INHERITS_FROM(TwoPTwoCCO2NI, HeterogeneousSpatialParams));
NEW_TYPE_TAG(HeterogeneousNIBoxTypeTag, INHERITS_FROM(BoxModel, HeterogeneousNITypeTag));
NEW_TYPE_TAG(HeterogeneousNICCTpfaTypeTag, INHERITS_FROM(CCTpfaModel, HeterogeneousNITypeTag));

// Set the grid type
SET_TYPE_PROP(HeterogeneousNITypeTag, Grid, Dune::ALUGrid<2, 2, Dune::cube, Dune::nonconforming>);

// Set the problem property
SET_TYPE_PROP(HeterogeneousNITypeTag, Problem, HeterogeneousProblem<TypeTag>);

// Set fluid configuration
SET_TYPE_PROP(HeterogeneousNITypeTag, FluidSystem, FluidSystems::BrineCO2<typename GET_PROP_TYPE(TypeTag, Scalar),
                                                                        HeterogeneousCO2Tables::CO2Tables>);

// Use Moles
SET_BOOL_PROP(HeterogeneousNITypeTag, UseMoles, false);
#endif
}


/*!
 * \ingroup CO2Model
 * \ingroup ImplicitTestProblems
 * \brief Definition of a problem, where CO2 is injected in a reservoir.
 *
 * The domain is sized 200m times 100m and consists of four layers, a
 * permeable reservoir layer at the bottom, a barrier rock layer with reduced permeability, another reservoir layer
 * and at the top a barrier rock layer with a very low permeablility.
 *
 * CO2 is injected at the permeable bottom layer
 * from the left side. The domain is initially filled with brine.
 *
 * The grid is unstructered and permeability and porosity for the elements are read in from the grid file. The grid file
 * also contains so-called boundary ids which can be used assigned during the grid creation in order to differentiate
 * between different parts of the boundary.
 * These boundary ids can be imported into the problem where the boundary conditions can then be assigned accordingly.
 *
 * The model is able to use either mole or mass fractions. The property useMoles can be set to either true or false in the
 * problem file. Make sure that the according units are used in the problem setup. The default setting for useMoles is false.
 *
 * To run the simulation execute the following line in shell (works with the box and cell centered spatial discretization method):
 * <tt>./test_ccco2 </tt> or <tt>./test_boxco2 </tt>
 */
template <class TypeTag >
class HeterogeneousProblem : public PorousMediumFlowProblem<TypeTag>
{
    using ParentType = PorousMediumFlowProblem<TypeTag>;
    using GridView = typename GET_PROP_TYPE(TypeTag, GridView);
    using Scalar = typename GET_PROP_TYPE(TypeTag, Scalar);
    using FluidSystem = typename GET_PROP_TYPE(TypeTag, FluidSystem);
    using VolumeVariables = typename GET_PROP_TYPE(TypeTag, VolumeVariables);
    using ElementVolumeVariables = typename GET_PROP_TYPE(TypeTag, ElementVolumeVariables);
    using Indices = typename GET_PROP_TYPE(TypeTag, Indices);

    enum { dimWorld = GridView::dimensionworld };

    // copy some indices for convenience
    enum {
        pressureIdx = Indices::pressureIdx,
        switchIdx = Indices::switchIdx,

        wPhaseIdx = Indices::wPhaseIdx,
        nPhaseIdx = Indices::nPhaseIdx,

        wPhaseOnly = Indices::wPhaseOnly,

        nCompIdx = FluidSystem::nCompIdx,
        BrineIdx = FluidSystem::BrineIdx,
        CO2Idx = FluidSystem::CO2Idx
    };
    enum {
        conti0EqIdx = Indices::conti0EqIdx,
        contiCO2EqIdx = conti0EqIdx + CO2Idx
    };

#if !ISOTHERMAL
    enum {
        temperatureIdx = Indices::temperatureIdx,
        energyEqIdx = Indices::energyEqIdx,
    };
#endif

    using PrimaryVariables = typename GET_PROP_TYPE(TypeTag, PrimaryVariables);
    using NeumannFluxes = typename GET_PROP_TYPE(TypeTag, NumEqVector);
    using BoundaryTypes = typename GET_PROP_TYPE(TypeTag, BoundaryTypes);
    using Element = typename GridView::template Codim<0>::Entity;
    using FVGridGeometry = typename GET_PROP_TYPE(TypeTag, FVGridGeometry);
    using FVElementGeometry = typename GET_PROP_TYPE(TypeTag, FVElementGeometry);
    using SubControlVolumeFace = typename FVElementGeometry::SubControlVolumeFace;
    using SubControlVolume = typename FVElementGeometry::SubControlVolume;
    using GlobalPosition = Dune::FieldVector<Scalar, dimWorld>;
    using CO2 = Dumux::CO2<Scalar, HeterogeneousCO2Tables::CO2Tables>;

    //! property that defines whether mole or mass fractions are used
    static constexpr bool useMoles = GET_PROP_VALUE(TypeTag, UseMoles);

    // the discretization method we are using
    static constexpr auto discMethod = GET_PROP_VALUE(TypeTag, DiscretizationMethod);

public:
    /*!
     * \brief The constructor
     *
     * \param timeManager The time manager
     * \param gridView The grid view
     */
    HeterogeneousProblem(std::shared_ptr<const FVGridGeometry> fvGridGeometry)
    : ParentType(fvGridGeometry)
    , injectionTop_(1)
    , injectionBottom_(2)
    , dirichletBoundary_(3)
    , noFlowBoundary_(4)
    {
        nTemperature_       = getParam<int>("FluidSystem.NTemperature");
        nPressure_          = getParam<int>("FluidSystem.NPressure");
        pressureLow_        = getParam<Scalar>("FluidSystem.PressureLow");
        pressureHigh_       = getParam<Scalar>("FluidSystem.PressureHigh");
        temperatureLow_     = getParam<Scalar>("FluidSystem.TemperatureLow");
        temperatureHigh_    = getParam<Scalar>("FluidSystem.TemperatureHigh");
        depthBOR_           = getParam<Scalar>("Problem.DepthBOR");
        name_               = getParam<std::string>("Problem.Name");
        injectionRate_      = getParam<Scalar>("Problem.InjectionRate");

#if !ISOTHERMAL
        injectionPressure_ = getParam<Scalar>("Problem.InjectionPressure");
        injectionTemperature_ = getParam<Scalar>("Problem.InjectionTemperature");
#endif

        // set the spatial parameters by reading the DGF grid file
        this->spatialParams().getParamsFromGrid();

        // initialize the tables of the fluid system
        FluidSystem::init(/*Tmin=*/temperatureLow_,
                          /*Tmax=*/temperatureHigh_,
                          /*nT=*/nTemperature_,
                          /*pmin=*/pressureLow_,
                          /*pmax=*/pressureHigh_,
                          /*np=*/nPressure_);

        //stating in the console whether mole or mass fractions are used
        if(useMoles)
            std::cout<<"problem uses mole fractions"<<std::endl;
        else
            std::cout<<"problem uses mass fractions"<<std::endl;

        // precompute the boundary types for the box method from the cell-centered boundary types
        scvfToScvBoundaryTypes_.computeBoundaryTypes(*this);
    }

    /*!
     * \brief Append all quantities of interest which can be derived
     *        from the solution of the current time step to the VTK
     *        writer.
     */
    template<class VTKWriter>
    void addFieldsToWriter(VTKWriter& vtk)
    {
        const auto numElements = this->fvGridGeometry().gridView().size(0);
        const auto numDofs = this->fvGridGeometry().numDofs();

        vtkKxx_.resize(numElements);
        vtkPorosity_.resize(numElements);
        vtkBoxVolume_.resize(numDofs, 0.0);

        vtk.addField(vtkKxx_, "Kxx");
        vtk.addField(vtkPorosity_, "cellwisePorosity");
        vtk.addField(vtkBoxVolume_, "boxVolume");

#if !ISOTHERMAL
        vtk.addVolumeVariable([](const VolumeVariables& v){ return v.enthalpy(wPhaseIdx); }, "enthalpyW");
        vtk.addVolumeVariable([](const VolumeVariables& v){ return v.enthalpy(nPhaseIdx); }, "enthalpyN");
#else
        vtkTemperature_.resize(numDofs, 0.0);
        vtk.addField(vtkTemperature_, "temperature");
#endif

        const auto& gridView = this->fvGridGeometry().gridView();
        for (const auto& element : elements(gridView))
        {
            const auto eIdx = this->fvGridGeometry().elementMapper().index(element);
            auto fvGeometry = localView(this->fvGridGeometry());
            fvGeometry.bindElement(element);

            for (const auto& scv : scvs(fvGeometry))
            {
                const auto dofIdxGlobal = scv.dofIndex();
                vtkBoxVolume_[dofIdxGlobal] += scv.volume();
#if ISOTHERMAL
                vtkTemperature_[dofIdxGlobal] = initialTemperatureField_(scv.dofPosition());
#endif
            }

            vtkKxx_[eIdx] = this->spatialParams().permeability(eIdx);
            vtkPorosity_[eIdx] = this->spatialParams().porosity(eIdx);
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
    { return name_; }

    /*!
     * \brief Returns the temperature within the domain.
     *
     * \param globalPos The position
     *
     * This problem assumes a geothermal gradient with
     * a surface temperature of 10 degrees Celsius.
     */
    Scalar temperatureAtPos(const GlobalPosition &globalPos) const
    { return initialTemperatureField_(globalPos); }

    // \}

    /*!
     * \name Boundary conditions
     */
    // \{

    /*!
     * \brief Specifies which kind of boundary condition should be
     *        used for which equation on a given boundary segment.
     *
     * \param element The finite element
     * \param scv The sub control volume
     */
    BoundaryTypes boundaryTypes(const Element &element,
                                const SubControlVolume &scv) const
    { return scvfToScvBoundaryTypes_.boundaryTypes(scv); }

    /*!
     * \brief Specifies which kind of boundary condition should be
     *        used for which equation on a given boundary segment.
     *
     * \param element The finite element
     * \param scvf The sub control volume face
     */
    BoundaryTypes boundaryTypes(const Element &element,
                                const SubControlVolumeFace &scvf) const
    {
        BoundaryTypes bcTypes;
        const auto boundaryId = scvf.boundaryFlag();

        if (boundaryId < 1 || boundaryId > 4)
            DUNE_THROW(Dune::InvalidStateException, "Invalid boundary ID: " << boundaryId);

        if (boundaryId == dirichletBoundary_)
            bcTypes.setAllDirichlet();
        else
            bcTypes.setAllNeumann();

        return bcTypes;
    }

    /*!
     * \brief Evaluates the boundary conditions for a Dirichlet
     *        boundary segment
     *
     * \param returns the Dirichlet values for the conservation equations in
     *               \f$ [ \textnormal{unit of primary variable} ] \f$
     * \param globalPos The global position
     */
    PrimaryVariables dirichletAtPos(const GlobalPosition &globalPos) const
    { return initial_(globalPos); }

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
    NeumannFluxes neumann(const Element& element,
                          const FVElementGeometry& fvGeometry,
                          const ElementVolumeVariables& elemVolvars,
                          const SubControlVolumeFace& scvf) const
    {
        const auto boundaryId = scvf.boundaryFlag();

        NeumannFluxes fluxes(0.0);
         // kg/(m^2*s) or mole/(m^2*s) depending on useMoles
        if (boundaryId == injectionBottom_)
        {
            fluxes[contiCO2EqIdx] = useMoles ? -injectionRate_/FluidSystem::molarMass(nCompIdx) : -injectionRate_;
#if !ISOTHERMAL
            // energy fluxes are always mass specific
            fluxes[energyEqIdx] = -injectionRate_/*kg/(m^2 s)*/*CO2::gasEnthalpy(
                                    injectionTemperature_, injectionPressure_)/*J/kg*/; // W/(m^2)
#endif
        }

        return fluxes;
    }

    // \}

    /*!
     * \name Volume terms
     */
    // \{

    /*!
     * \brief Evaluates the initial values at a position
     *
     * \returns the initial values for the conservation equations in
     *           \f$ [ \textnormal{unit of primary variables} ] \f$
     * \param globalPos The global position
     */
    PrimaryVariables initialAtPos(const GlobalPosition &globalPos) const
    {
        return initial_(globalPos);
    }

    // \}

private:
    /*!
     * \brief Evaluates the initial values for a control volume
     *
     * The internal method for the initial condition
     *
     * \param values Stores the initial values for the conservation equations in
     *               \f$ [ \textnormal{unit of primary variables} ] \f$
     * \param globalPos The global position
     */
    PrimaryVariables initial_(const GlobalPosition &globalPos) const
    {
        PrimaryVariables values(0.0);
        values.setState(wPhaseOnly);

        const Scalar temp = initialTemperatureField_(globalPos);
        const Scalar densityW = FluidSystem::Brine::liquidDensity(temp, 1e7);

        const Scalar moleFracLiquidCO2 = 0.00;
        const Scalar moleFracLiquidBrine = 1.0 - moleFracLiquidCO2;

        const Scalar meanM = FluidSystem::molarMass(BrineIdx)*moleFracLiquidBrine
                             + FluidSystem::molarMass(CO2Idx)*moleFracLiquidCO2;

        if(useMoles) // mole-fraction formulation
            values[switchIdx] = moleFracLiquidCO2;
        else // mass-fraction formulation
            values[switchIdx] = moleFracLiquidCO2*FluidSystem::molarMass(CO2Idx)/meanM;

        values[pressureIdx] = 1.0e5 - densityW*this->gravity()[dimWorld-1]*(depthBOR_ - globalPos[dimWorld-1]);

#if !ISOTHERMAL
        values[temperatureIdx] = temp;
#endif
        return values;
    }

    Scalar initialTemperatureField_(const GlobalPosition globalPos) const
    {
        return 283.0 + (depthBOR_ - globalPos[dimWorld-1])*0.03;
    }

    Scalar depthBOR_;
    Scalar injectionRate_;

#if !ISOTHERMAL
    Scalar injectionPressure_, injectionTemperature_;
#endif

    int nTemperature_;
    int nPressure_;

    std::string name_ ;

    Scalar pressureLow_, pressureHigh_;
    Scalar temperatureLow_, temperatureHigh_;

    int injectionTop_;
    int injectionBottom_;
    int dirichletBoundary_;
    int noFlowBoundary_;

    // vtk output
    std::vector<Scalar> vtkKxx_, vtkPorosity_, vtkBoxVolume_, vtkTemperature_;
    ScvfToScvBoundaryTypes<BoundaryTypes, discMethod> scvfToScvBoundaryTypes_;
};

} //end namespace Dumux

#endif
