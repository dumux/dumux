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
/**
 * \file
 * \ingroup OnePNCTests
 * \brief Definition of a problem, for the 1pnc problem:
 * Component transport of nitrogen dissolved in the water phase.
 */
#ifndef DUMUX_1P2C_TEST_PROBLEM_HH
#define DUMUX_1P2C_TEST_PROBLEM_HH

#if HAVE_UG
#include <dune/grid/uggrid.hh>
#endif
#include <dune/grid/yaspgrid.hh>

#include <dumux/discretization/cctpfa.hh>
#include <dumux/discretization/ccmpfa.hh>
#include <dumux/discretization/box.hh>
#include <dumux/discretization/evalsolution.hh>
#include <dumux/discretization/evalgradients.hh>
#include <dumux/porousmediumflow/1pnc/model.hh>
#include <dumux/porousmediumflow/problem.hh>

#include <dumux/material/fluidsystems/h2on2.hh>
#include <dumux/material/fluidsystems/1padapter.hh>

#include "../spatialparams.hh"

namespace Dumux {

template <class TypeTag>
class OnePTwoCTestProblem;

namespace Properties {

// Create new type tags
namespace TTag {
struct OnePTwoCTest { using InheritsFrom = std::tuple<OnePNC>; };
struct OnePTwoCTestBox { using InheritsFrom = std::tuple<OnePTwoCTest, BoxModel>; };
struct OnePTwoCTestCCTpfa { using InheritsFrom = std::tuple<OnePTwoCTest, CCTpfaModel>; };
struct OnePTwoCTestCCMpfa { using InheritsFrom = std::tuple<OnePTwoCTest, CCMpfaModel>; };
} // end namespace TTag

// Set the grid type
#if HAVE_UG
template<class TypeTag>
struct Grid<TypeTag, TTag::OnePTwoCTest> { using type = Dune::UGGrid<2>; };
#else
template<class TypeTag>
struct Grid<TypeTag, TTag::OnePTwoCTest> { using type = Dune::YaspGrid<2>; };
#endif

// Set the problem property
template<class TypeTag>
struct Problem<TypeTag, TTag::OnePTwoCTest> { using type = OnePTwoCTestProblem<TypeTag>; };

// Set fluid configuration
template<class TypeTag>
struct FluidSystem<TypeTag, TTag::OnePTwoCTest>
{
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using H2ON2 = FluidSystems::H2ON2<Scalar, FluidSystems::H2ON2DefaultPolicy</*simplified=*/true>>;
    using type = FluidSystems::OnePAdapter<H2ON2, H2ON2::liquidPhaseIdx>;
};

// Set the spatial parameters
template<class TypeTag>
struct SpatialParams<TypeTag, TTag::OnePTwoCTest>
{
    using FVGridGeometry = GetPropType<TypeTag, Properties::FVGridGeometry>;
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using type = OnePNCTestSpatialParams<FVGridGeometry, Scalar>;
};

// Define whether mole(true) or mass (false) fractions are used
template<class TypeTag>
struct UseMoles<TypeTag, TTag::OnePTwoCTest> { static constexpr bool value = true; };
}


/*!
 * \ingroup OnePNCTests
 * \brief Definition of a problem, for the 1pnc problem:
 * Nitrogen is dissolved in the water phase and
 * is transported with the water flow from the left side to the right.
 *
 * The model domain is 1 m times 1 m with a discretization length of 0.05 m
 * and homogeneous soil properties (\f$ \mathrm{K=10e-10, \Phi=0.4, \tau=0.28}\f$).
 * Initially the domain is filled with pure water.
 *
 * At the left side, a Dirichlet condition defines a nitrogen mole fraction
 * of 0.3 mol/mol.
 * The water phase flows from the left side to the right due to the applied pressure
 * gradient of 1e5 Pa/m. The nitrogen is transported with the water flow
 * and leaves the domain at the right boundary
 * where again Dirichlet boundary conditions are applied. Here, the nitrogen mole
 * fraction is set to 0.0 mol/mol.
 *
 * This problem uses the \ref OnePNCModel model.
 *
 * To run the simulation execute the following line in shell:
 * <tt>./test_box1pnc -parameterFile ./test_box1pnc.input</tt> or
 * <tt>./test_cc1pnc -parameterFile ./test_cc1pnc.input</tt>
 */
template <class TypeTag>
class OnePTwoCTestProblem : public PorousMediumFlowProblem<TypeTag>
{
    using ParentType = PorousMediumFlowProblem<TypeTag>;

    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using Indices = typename GetPropType<TypeTag, Properties::ModelTraits>::Indices;
    using GridView = GetPropType<TypeTag, Properties::GridView>;
    using FluidSystem = GetPropType<TypeTag, Properties::FluidSystem>;
    using BoundaryTypes = GetPropType<TypeTag, Properties::BoundaryTypes>;
    using PrimaryVariables = GetPropType<TypeTag, Properties::PrimaryVariables>;
    using FVGridGeometry = GetPropType<TypeTag, Properties::FVGridGeometry>;
    using FVElementGeometry = typename GetPropType<TypeTag, Properties::FVGridGeometry>::LocalView;
    using NumEqVector = GetPropType<TypeTag, Properties::NumEqVector>;
    using ElementVolumeVariables = typename GetPropType<TypeTag, Properties::GridVolumeVariables>::LocalView;
    using SubControlVolumeFace = typename FVElementGeometry::SubControlVolumeFace;
    using Element = typename GridView::template Codim<0>::Entity;

    // copy some indices for convenience
    enum
    {
        // indices of the primary variables
        pressureIdx = Indices::pressureIdx,

        // component indices
        H2OIdx = FluidSystem::compIdx(FluidSystem::MultiPhaseFluidSystem::H2OIdx),
        N2Idx = FluidSystem::compIdx(FluidSystem::MultiPhaseFluidSystem::N2Idx),

        // indices of the equations
        contiH2OEqIdx = Indices::conti0EqIdx + H2OIdx,
        contiN2EqIdx = Indices::conti0EqIdx + N2Idx
    };

    //! property that defines whether mole or mass fractions are used
    static constexpr bool useMoles = getPropValue<TypeTag, Properties::UseMoles>();
    static const bool isBox = GetPropType<TypeTag, Properties::FVGridGeometry>::discMethod == DiscretizationMethod::box;

    static const int dimWorld = GridView::dimensionworld;
    using GlobalPosition = typename SubControlVolumeFace::GlobalPosition;

public:
    OnePTwoCTestProblem(std::shared_ptr<const FVGridGeometry> fvGridGeometry)
    : ParentType(fvGridGeometry), useNitscheTypeBc_(getParam<bool>("Problem.UseNitscheTypeBc", false))
    {
        //initialize fluid system
        FluidSystem::init();

        // stating in the console whether mole or mass fractions are used
        if(useMoles)
            std::cout<<"problem uses mole fractions"<<std::endl;
        else
            std::cout<<"problem uses mass fractions"<<std::endl;
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

        if(globalPos[0] < eps_)
            values.setAllDirichlet();
        else
            values.setAllNeumann();

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

        // condition for the N2 molefraction at left boundary
        if (globalPos[0] < eps_ )
        {
            values[N2Idx] = 2.0e-5;
        }

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
        // set a fixed pressure on the right side of the domain
        const Scalar dirichletPressure = 1e5;

        NumEqVector flux(0.0);
        const auto& ipGlobal = scvf.ipGlobal();
        const auto& volVars = elemVolVars[scvf.insideScvIdx()];

        // no-flow everywhere except at the right boundary
        if(ipGlobal[0] < this->fvGridGeometry().bBoxMax()[0] - eps_)
            return flux;

        // if specified in the input file, use a Nitsche type boundary condition for the box model,
        // otherwise compute the acutal fluxes explicitly
        if(isBox && useNitscheTypeBc_)
        {
            flux[contiH2OEqIdx] = (volVars.pressure() - dirichletPressure) * 1e7;
            flux[contiN2EqIdx] = flux[contiH2OEqIdx]  * (useMoles ? volVars.moleFraction(0, N2Idx) :
                                                                    volVars.massFraction(0, N2Idx));
            return flux;
        }

        // construct the element solution
        const auto elemSol = [&]()
        {
            auto sol = elementSolution(element, elemVolVars, fvGeometry);

            if(isBox)
                for(auto&& scvf : scvfs(fvGeometry))
                    if(scvf.center()[0] > this->fvGridGeometry().bBoxMax()[0] - eps_)
                        sol[fvGeometry.scv(scvf.insideScvIdx()).localDofIndex()][pressureIdx] = dirichletPressure;

            return sol;
        }();

        // evaluate the gradient
        const auto gradient = [&]()->GlobalPosition
        {
            if(isBox)
            {
                const auto grads = evalGradients(element, element.geometry(), fvGeometry.fvGridGeometry(), elemSol, ipGlobal);
                return grads[pressureIdx];
            }

            else
            {
                const auto& scvCenter = fvGeometry.scv(scvf.insideScvIdx()).center();
                const Scalar scvCenterPresureSol = elemSol[0][pressureIdx];
                auto grad = ipGlobal - scvCenter;
                grad /= grad.two_norm2();
                grad *= (dirichletPressure - scvCenterPresureSol);
                return grad;
            }
        }();

        const Scalar K = volVars.permeability();
        const Scalar density = useMoles ? volVars.molarDensity() : volVars.density();

        // calculate the flux
        Scalar tpfaFlux = gradient * scvf.unitOuterNormal();
        tpfaFlux *= -1.0  * K;
        tpfaFlux *=  density * volVars.mobility();
        flux[contiH2OEqIdx] = tpfaFlux;

        // emulate an outflow condition for the component transport on the right side
        flux[contiN2EqIdx] = tpfaFlux  * (useMoles ? volVars.moleFraction(0, N2Idx) : volVars.massFraction(0, N2Idx));

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
        priVars[pressureIdx] = 2e5 - 1e5*globalPos[0]; // initial condition for the pressure
        priVars[N2Idx] = 0.0;  // initial condition for the N2 molefraction
        return priVars;
    }
        static constexpr Scalar eps_ = 1e-6;
        bool useNitscheTypeBc_;
    };

} //end namespace Dumux

#endif
