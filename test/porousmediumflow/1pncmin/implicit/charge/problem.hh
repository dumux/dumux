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
 * \ingroup OnePNCMinTests
 * \brief Definition of a problem for thermochemical heat storage using \f$ \textnormal{CaO},   \textnormal{Ca} \left( \textnormal{OH} \right)_2\f$.
 */
#ifndef DUMUX_CHARGE_PROBLEM_HH
#define DUMUX_CHARGE_PROBLEM_HH

#include <dune/grid/yaspgrid.hh>
#include <dumux/porousmediumflow/1pncmin/model.hh>
#include <dumux/discretization/elementsolution.hh>
#include <dumux/discretization/box/properties.hh>
#include <dumux/discretization/cellcentered/tpfa/properties.hh>
#include <dumux/discretization/cellcentered/mpfa/properties.hh>
#include <dumux/porousmediumflow/problem.hh>
#include <dumux/material/fluidsystems/1padapter.hh>
#include <dumux/material/fluidsystems/h2on2.hh>
#include <dumux/material/fluidmatrixinteractions/1p/thermalconductivityaverage.hh>
#include <dumux/material/components/cao2h2.hh>
#include <dumux/material/solidsystems/compositionalsolidphase.hh>
#include <dumux/material/fluidsystems/1padapter.hh>
// box solution dependent Neumann, outflow
#include <dumux/discretization/evalgradients.hh>
#include <dumux/discretization/evalsolution.hh>
#include <dumux/material/fluidsystems/h2on2.hh>
#include "../discharge/spatialparams.hh"
#include "../discharge/reaction.hh"
#include "../discharge/modifiedcao.hh"

namespace Dumux {

template <class TypeTag>
class ChargeProblem;

namespace Properties {
// Create new type tags
namespace TTag {
struct Charge { using InheritsFrom = std::tuple<OnePNCMinNI>; };
struct ChargeBox { using InheritsFrom = std::tuple<Charge, BoxModel>; };
} // end namespace TTag

// Set the grid type
template<class TypeTag>
struct Grid<TypeTag, TTag::Charge> { using type = Dune::YaspGrid<2>; };
// Set the problem property
template<class TypeTag>
struct Problem<TypeTag, TTag::Charge> { using type = ChargeProblem<TypeTag>; };

// The fluid system
template<class TypeTag>
struct FluidSystem<TypeTag, TTag::Charge>
{
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using H2ON2 = FluidSystems::H2ON2<Scalar>;
    static constexpr auto phaseIdx = H2ON2::gasPhaseIdx; // simulate the gas phase
    using type = FluidSystems::OnePAdapter<H2ON2, phaseIdx>;
};

template<class TypeTag>
struct SolidSystem<TypeTag, TTag::Charge>
{
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using ComponentOne = Components::ModifiedCaO<Scalar>;
    using ComponentTwo = Components::CaO2H2<Scalar>;
    using type = SolidSystems::CompositionalSolidPhase<Scalar, ComponentOne, ComponentTwo>;
};

// // Enable velocity output
// template<class TypeTag>
// struct VtkAddVelocity<TypeTag, TTag::Charge> { static constexpr bool value = false; };

// Set the spatial parameters
template<class TypeTag>
struct SpatialParams<TypeTag, TTag::Charge>
{
    using FVGridGeometry = GetPropType<TypeTag, Properties::FVGridGeometry>;
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using type = ThermoChemSpatialParams<FVGridGeometry, Scalar>;
};

// Define whether mole(true) or mass (false) fractions are used
template<class TypeTag>
struct UseMoles<TypeTag, TTag::Charge> { static constexpr bool value = true; };
}
/*!
 * \ingroup OnePNCMinTests
 *
 * \brief Test for the 1pncmin model in combination with the NI model for a quasi batch
 * reaction of Calciumoxyde to Calciumhydroxide.
 *
 * The boundary conditions of the batch test are such, that there are no gradients for temperature, pressure and gas water concentration within the reactor.
 *
 * To run the simulation execute the following line in shell:
 * <tt>./test_1pncminni_box -ParameterFile </tt>
 * The test only runs for the box discretization.
 */
template <class TypeTag>
class ChargeProblem : public PorousMediumFlowProblem<TypeTag>
{
    using ParentType = PorousMediumFlowProblem<TypeTag>;
    using GridView = GetPropType<TypeTag, Properties::GridView>;
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using FluidSystem = GetPropType<TypeTag, Properties::FluidSystem>;
    using SolidSystem = GetPropType<TypeTag, Properties::SolidSystem>;
    using VolumeVariables = GetPropType<TypeTag, Properties::VolumeVariables>;
    using Indices = typename GetPropType<TypeTag, Properties::ModelTraits>::Indices;
    using ElementVolumeVariables = typename GetPropType<TypeTag, Properties::GridVolumeVariables>::LocalView;
    using PrimaryVariables = GetPropType<TypeTag, Properties::PrimaryVariables>;
    using BoundaryTypes = GetPropType<TypeTag, Properties::BoundaryTypes>;
    using Element = typename GridView::template Codim<0>::Entity;
    using FVElementGeometry = typename GetPropType<TypeTag, Properties::FVGridGeometry>::LocalView;
    using SubControlVolume = typename FVElementGeometry::SubControlVolume;
    using SubControlVolumeFace = typename FVElementGeometry::SubControlVolumeFace;
    using FVGridGeometry = GetPropType<TypeTag, Properties::FVGridGeometry>;
    using NumEqVector = GetPropType<TypeTag, Properties::NumEqVector>;
    using SolutionVector = GetPropType<TypeTag, Properties::SolutionVector>;
    using ReactionRate = ThermoChemReaction;
    using FluidState = GetPropType<TypeTag, Properties::FluidState>;

    enum { dim = GridView::dimension };
    enum { dimWorld = GridView::dimensionworld };

    enum
    {
        // Indices of the primary variables
        pressureIdx = Indices::pressureIdx, //gas-phase pressure

        //firstMoleFracIdx = FluidSystem::compIdx(FluidSystem::MultiPhaseFluidSystem::H2OIdx), // mole fraction water
        H2OIdx = FluidSystem::compIdx(FluidSystem::MultiPhaseFluidSystem::H2OIdx), // mole fraction water

        CaOIdx = FluidSystem::numComponents,
        CaO2H2Idx = FluidSystem::numComponents+1,

        // Equation Indices
        conti0EqIdx = Indices::conti0EqIdx,

        // Phase Indices
        phaseIdx = FluidSystem::phase0Idx,
        cPhaseIdx = SolidSystem::comp0Idx,
        hPhaseIdx = SolidSystem::comp0Idx+1,

        temperatureIdx = Indices::temperatureIdx,
        energyEqIdx = Indices::energyEqIdx
    };

    using GlobalPosition = typename SubControlVolumeFace::GlobalPosition;

public:
    /*!
     * \brief The constructor
     *
     * \param fvGridGeometry The finite volume grid geometry
     */
    ChargeProblem(std::shared_ptr<const FVGridGeometry> fvGridGeometry)
        : ParentType(fvGridGeometry)
    {
        name_      = getParam<std::string>("Problem.Name");
        FluidSystem::init(/*tempMin=*/473.15,
                          /*tempMax=*/623.0,
                          /*numTemptempSteps=*/25,
                          /*startPressure=*/0,
                          /*endPressure=*/9e6,
                          /*pressureSteps=*/200);


        unsigned int codim = GET_PROP_TYPE(TypeTag, FVGridGeometry)::discMethod == DiscretizationMethod::box ? dim : 0;
        permeability_.resize(fvGridGeometry->gridView().size(codim));
        porosity_.resize(fvGridGeometry->gridView().size(codim));
        reactionRate_.resize(fvGridGeometry->gridView().size(codim));
    }

    /*!
     * \brief The problem name.
     *
     * This is used as a prefix for files generated by the simulation.
     */
    const std::string name() const
    { return name_; }

    /*!
     * \brief Sets the currently used time step size.
     *
     * This is necessary to limit the source terms to the maximum possible rate.
     */
    void setTimeStepSize( Scalar timeStepSize )
     {
        timeStepSize_ = timeStepSize;
     }

    /*!
     * \name Boundary conditions
     *
     * \brief Specifies which kind of boundary condition should be
     *        used for which equation on a given boundary segment
     *
     * \param globalPos The global position
     */
    BoundaryTypes boundaryTypesAtPos( const GlobalPosition &globalPos) const
    {
        BoundaryTypes values;

            values.setNeumann(pressureIdx);
            values.setNeumann(H2OIdx);
            values.setNeumann(temperatureIdx);
            values.setNeumann(CaO2H2Idx);
            values.setNeumann(CaOIdx);
        return values;
    }
    /*!
     * \brief Evaluates the boundary conditions for a Neumann
     *        boundary segment in dependency on the current solution.
     *
     * \param element The element
     * \param fvGeometry The finite volume geometry
     * \param elemVolVars The element volume variables
     * \param scvf The subcontrolvolume face
     *
     * \f$ [ \textnormal{unit of conserved quantity} / (m^(dim-1) \cdot s )] \f$
     * Negative values indicate an inflow.
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
           Scalar InFlowAir = 4.64;
           Scalar InFlowH2O = 0.01;
           Scalar tIn = 773.15;

           FluidState fluidstateBorder;

           fluidstateBorder.setTemperature(tIn);
           fluidstateBorder.setPressure(phaseIdx, elemVolVars[scv].pressure(phaseIdx));

           Scalar T = elemVolVars[scv].temperature();

           Scalar deltaH = 0.0; // if temperature at the right border > 573.15 K, cool down
                                //temperature of injected fluid via the enthalpyflux deltaH

           deltaH = -(T-tIn)* 1e5;

           Scalar hInAir = InFlowAir*FluidSystem::molarMass(H2OIdx-1)
                           *FluidSystem::componentEnthalpy(fluidstateBorder, phaseIdx, H2OIdx-1);

           Scalar hInH2O = InFlowH2O*FluidSystem::molarMass(H2OIdx)
                           *FluidSystem::componentEnthalpy(fluidstateBorder, phaseIdx, H2OIdx);
           flux[pressureIdx] = - InFlowAir; //[mol/s] gas inflow of the air component
           flux[H2OIdx] = - InFlowH2O;//[mol/s] gas inflow of the water component
           flux[temperatureIdx] = - (hInAir + hInH2O -deltaH); //[J/s] enthalpy inflow
        }
        // outflow BC
        if(globalPos[0] > this->fvGridGeometry().bBoxMax()[0] - eps_ )
        {
            // set a fixed pressure on the right side of the domain
            const Scalar dirichletPressure = 1.0e5;
            const Scalar K = elemVolVars[scv].permeability();
            const Scalar molarDensity = elemVolVars[scv].molarDensity(phaseIdx);
            const Scalar density = elemVolVars[scv].density(phaseIdx);

            // construct the element solution
            const auto elemSol = [&]()
            {
                auto sol = elementSolution(element, elemVolVars, fvGeometry);

                for(auto&& scvf : scvfs(fvGeometry))
                    if(scvf.center()[0] > this->fvGridGeometry().bBoxMax()[0] - eps_)
                        sol[fvGeometry.scv(scvf.insideScvIdx()).localDofIndex()][pressureIdx] = dirichletPressure;

                return sol;
            }();
            // evaluate the gradient
            const auto gradient = [&]()->GlobalPosition
            {
                const auto grads = evalGradients(element, element.geometry(), fvGeometry.fvGridGeometry(), elemSol, globalPos);
                return grads[pressureIdx];

            }();

            // calculate the flux
            Scalar tpfaFlux = gradient * scvf.unitOuterNormal();
            tpfaFlux *= -1.0  * K;

            if(tpfaFlux < 0) tpfaFlux = 0.0; //make sure that there is no influx from the right

            Scalar tpfaFluxMole = tpfaFlux * molarDensity * elemVolVars[scv].mobility(phaseIdx);
            Scalar tpfaFluxMass = tpfaFlux * density * elemVolVars[scv].mobility(phaseIdx);

            // emulate an outflow condition for the component transport on the right side
            flux[pressureIdx] = tpfaFluxMole * (elemVolVars[scv].moleFraction(phaseIdx, (H2OIdx-1)));
            flux[H2OIdx] = tpfaFluxMole * elemVolVars[scv].moleFraction(phaseIdx, H2OIdx);
            flux[temperatureIdx] = tpfaFluxMass * (FluidSystem::enthalpy(elemVolVars[scv].fluidState(), phaseIdx));
        }
       return flux;
    }
    /*!
     *
     *
     * \brief Evaluates the initial values for a control volume in
     *               \f$ [ \textnormal{unit of primary variables} ] \f$
     *
     * \param globalPos The global position
     */
    PrimaryVariables initialAtPos(const GlobalPosition& globalPos) const
    {
        PrimaryVariables priVars(0.0);

        Scalar pInit;
        Scalar tInit;
        Scalar h2oInit;
        Scalar CaOInit;
        Scalar CaO2H2Init;

        pInit = getParam<Scalar>("Problem.PressureInitial");
        tInit = getParam<Scalar>("Problem.TemperatureInitial");
        h2oInit = getParam<Scalar>("Problem.VaporInitial");
        CaOInit = getParam<Scalar>("Problem.CaOInitial");
        CaO2H2Init = getParam<Scalar>("Problem.CaO2H2Initial");

        priVars[pressureIdx] = pInit;
        priVars[H2OIdx]   = h2oInit;
        priVars[temperatureIdx] = tInit;
        priVars[CaOIdx] = CaOInit;
        priVars[CaO2H2Idx]   = CaO2H2Init;

        return priVars;
    }

    /*!
     * \brief Evaluate the source term for all phases within a given
     *        sub-control-volume in units of \f$ [ \textnormal{unit of conserved quantity} / (m^3 \cdot s )] \f$.
     *
     * This is the method for the case where the source term is
     * potentially solution dependent and requires some quantities that
     * are specific to the fully-implicit method.
     *
     * \param element The finite element
     * \param fvGeometry The finite-volume geometry
     * \param elemVolVars All volume variables for the element
     * \param scv The subcontrolvolume
     *
     * For this method, the \a values parameter stores the conserved quantity rate
     * generated or annihilate per volume unit. Positive values mean
     * that the conserved quantity is created, negative ones mean that it vanishes.
     * E.g. for the mass balance that would be a mass rate in \f$ [ kg / (m^3 \cdot s)] \f$.
     */
     PrimaryVariables source(const Element &element,
                            const FVElementGeometry& fvGeometry,
                            const ElementVolumeVariables& elemVolVars,
                            const SubControlVolume &scv) const
    {

        PrimaryVariables source(0.0);
        const auto& volVars = elemVolVars[scv];

        Scalar qMass = rrate_.thermoChemReactionSimple(volVars);

        const auto elemSol = elementSolution(element, elemVolVars, fvGeometry);
        Scalar qMole = qMass/FluidSystem::molarMass(H2OIdx)*(1-volVars.porosity());

        // make sure not more solid reacts than present
        // In this test, we only consider charge. Therefore, we use the hPhaseIdx for CaO2H2.
        if (-qMole*timeStepSize_ + volVars.solidVolumeFraction(hPhaseIdx)* volVars.solidComponentMolarDensity(hPhaseIdx) < 0 + eps_)
        {
            qMole = -volVars.solidVolumeFraction(hPhaseIdx)* volVars.solidComponentMolarDensity(hPhaseIdx)/timeStepSize_;
        }

        source[conti0EqIdx+CaO2H2Idx] = qMole;
        source[conti0EqIdx+CaOIdx] = - qMole;
        source[conti0EqIdx+H2OIdx] = - qMole;

        Scalar deltaH = 108.3e3; // J/mol
        source[energyEqIdx] = qMole * (deltaH - 4*(volVars.pressure(phaseIdx)/volVars.molarDensity(phaseIdx))) ;

        return source;
    }
   /*!
     * \brief Return the permeability
     */
    const std::vector<Scalar>& getPerm()
    {
        return permeability_;
    }

   /*!
     * \brief Return the porosity
     */
    const std::vector<Scalar>& getPoro()
    {
        return porosity_;
    }

     /*!
     * \brief Return the reaction rate
     */
    const std::vector<Scalar>& getRRate()
    {
        return reactionRate_;
    }

    /*!
     * \brief Adds additional VTK output data to the VTKWriter. Function is called by the output module on every write.
     */
    void updateVtkOutput(const SolutionVector& curSol)
    {
        for (const auto& element : elements(this->fvGridGeometry().gridView()))
        {
            const auto elemSol = elementSolution(element, curSol, this->fvGridGeometry());

            auto fvGeometry = localView(this->fvGridGeometry());
            fvGeometry.bindElement(element);

            for (auto&& scv : scvs(fvGeometry))
            {
                VolumeVariables volVars;
                volVars.update(elemSol, *this, element, scv);
                const auto dofIdxGlobal = scv.dofIndex();
                permeability_[dofIdxGlobal] = this->spatialParams().permeability(element, scv, elemSol);
                porosity_[dofIdxGlobal] = volVars.porosity();
                reactionRate_[dofIdxGlobal] = rrate_.thermoChemReaction(volVars);
            }
        }
    }


private:
    std::string name_;

    static constexpr Scalar eps_ = 1e-6;
    std::vector<double> permeability_;
    std::vector<double> porosity_;
    std::vector<double> reactionRate_;
    std::ofstream outputFile_;
    ReactionRate rrate_;
    Scalar timeStepSize_;
};

} //end namespace

#endif
