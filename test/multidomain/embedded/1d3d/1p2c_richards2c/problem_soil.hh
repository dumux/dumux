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
 * \ingroup EmbeddedTests
 * \brief Definition of a problem, for the 1p2c problem:
 * Component transport of oxygen in interstitial fluid.
 */

#ifndef DUMUX_TISSUE_PROBLEM_HH
#define DUMUX_TISSUE_PROBLEM_HH

#include <dune/geometry/quadraturerules.hh>
#include <dune/grid/yaspgrid.hh>
#include <dune/grid/uggrid.hh>
#include <dune/localfunctions/lagrange/pqkfactory.hh>

#include <dumux/common/boundarytypes.hh>
#include <dumux/common/math.hh>
#include <dumux/common/parameters.hh>
#include <dumux/common/properties.hh>
#include <dumux/discretization/cctpfa.hh>

#include <dumux/porousmediumflow/richardsnc/model.hh>
#include <dumux/porousmediumflow/problem.hh>
#include <dumux/material/components/simpleh2o.hh>
#include <dumux/material/components/constant.hh>
#include <dumux/material/fluidsystems/liquidphase2c.hh>

#include "spatialparams_soil.hh"

namespace Dumux {

template <class TypeTag>
class SoilProblem;

namespace Properties {

// Create new type tags
namespace TTag {
struct Soil { using InheritsFrom = std::tuple<RichardsNC, CCTpfaModel>; };
} // end namespace TTag

// Set the grid type
#if HAVE_UG
template<class TypeTag>
struct Grid<TypeTag, TTag::Soil> { using type = Dune::UGGrid<3>; };
#else
template<class TypeTag>
struct Grid<TypeTag, TTag::Soil> { using type = Dune::YaspGrid<3, Dune::EquidistantOffsetCoordinates<double, 3>>; };
#endif

template<class TypeTag>
struct EnableGridGeometryCache<TypeTag, TTag::Soil> { static constexpr bool value = true; };
template<class TypeTag>
struct EnableGridVolumeVariablesCache<TypeTag, TTag::Soil> { static constexpr bool value = true; };
template<class TypeTag>
struct EnableGridFluxVariablesCache<TypeTag, TTag::Soil> { static constexpr bool value = true; };
template<class TypeTag>
struct SolutionDependentAdvection<TypeTag, TTag::Soil> { static constexpr bool value = false; };
template<class TypeTag>
struct SolutionDependentMolecularDiffusion<TypeTag, TTag::Soil> { static constexpr bool value = false; };
template<class TypeTag>
struct SolutionDependentHeatConduction<TypeTag, TTag::Soil> { static constexpr bool value = false; };

// Set the problem property
template<class TypeTag>
struct Problem<TypeTag, TTag::Soil> { using type = SoilProblem<TypeTag>; };

// Set the spatial parameters
template<class TypeTag>
struct SpatialParams<TypeTag, TTag::Soil>
{
    using type = SoilSpatialParams<GetPropType<TypeTag, Properties::GridGeometry>,
                                   GetPropType<TypeTag, Properties::Scalar>>;
};

// Set the fluid system
template<class TypeTag>
struct FluidSystem<TypeTag, TTag::Soil>
{
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using type = FluidSystems::LiquidPhaseTwoC<Scalar, Components::SimpleH2O<Scalar>,
                                                       Components::Constant<1, Scalar>>;
};

template<class TypeTag>
struct UseMoles<TypeTag, TTag::Soil> { static constexpr bool value = true; };

} // end namespace Properties


/*!
 * \ingroup EmbeddedTests
 * \brief Definition of a problem, for the 1p2c problem:
 * Component transport of oxygen in interstitial fluid.
 */
template <class TypeTag>
class SoilProblem : public PorousMediumFlowProblem<TypeTag>
{
    using ParentType = PorousMediumFlowProblem<TypeTag>;
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using GridGeometry = GetPropType<TypeTag, Properties::GridGeometry>;
    using FVElementGeometry = typename GridGeometry::LocalView;
    using SubControlVolume = typename GridGeometry::SubControlVolume;
    using GridView = typename GridGeometry::GridView;
    using GlobalPosition = typename GridGeometry::GlobalCoordinate;
    using PrimaryVariables = GetPropType<TypeTag, Properties::PrimaryVariables>;
    using NumEqVector = GetPropType<TypeTag, Properties::NumEqVector>;
    using SolutionVector = GetPropType<TypeTag, Properties::SolutionVector>;
    using GridVariables = GetPropType<TypeTag, Properties::GridVariables>;
    using BoundaryTypes = Dumux::BoundaryTypes<GetPropType<TypeTag, Properties::ModelTraits>::numEq()>;
    using PointSource = GetPropType<TypeTag, Properties::PointSource>;
    using Element = typename GridView::template Codim<0>::Entity;

    using CouplingManager = GetPropType<TypeTag, Properties::CouplingManager>;

public:
    using FluidSystem = GetPropType<TypeTag, Properties::FluidSystem>;
    enum Indices {
        // world dimension
        dim = GridView::dimension,
        dimWorld = GridView::dimensionworld,

        pressureIdx = 0,
        transportCompIdx = 1,

        conti0EqIdx = 0,
        transportEqIdx = 1,

        liquidPhaseIdx = FluidSystem::liquidPhaseIdx
    };

    SoilProblem(std::shared_ptr<const GridGeometry> gridGeometry,
                std::shared_ptr<CouplingManager> couplingManager)
    : ParentType(gridGeometry, "Soil")
    , couplingManager_(couplingManager)
    {
        //read parameters from input file
        name_  =  getParam<std::string>("Vtk.OutputName") + "_" + getParamFromGroup<std::string>(this->paramGroup(), "Problem.Name");
        contaminantMoleFraction_ = getParam<Scalar>("Problem.ContaminantMoleFraction");

        // for initial conditions
        const Scalar sw = getParam<Scalar>("Problem.InitTopSaturation", 0.3); // start with 30% saturation on top
        pcTop_ = this->spatialParams().fluidMatrixInteractionAtPos(gridGeometry->bBoxMax()).pc(sw);
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
     * \brief Returns the temperature within the domain [K].
     *
     * This problem assumes a temperature of 10 degrees Celsius.
     */
    Scalar temperature() const
    { return 273.15 + 10; } // in [K]

    /*
      * \brief Returns the reference pressure [Pa] of the nonwetting
     *        fluid phase within a finite volume.
     *
     * This problem assumes a constant reference pressure of 1 bar.
     */
    Scalar nonwettingReferencePressure() const
    { return 1.0e5; }


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
        values.setAllNeumann();
        return values;
    }

    /*!
     * \brief Evaluates the boundary conditions for a Dirichlet boundary segment.
     *
     * \param globalPos The position for which the bc type should be evaluated
     *
     * For this method, the \a values parameter stores primary variables.
     */
    PrimaryVariables dirichletAtPos(const GlobalPosition &globalPos) const
    { return initialAtPos(globalPos); }

    // \}

    /*!
     * \name Volume terms
     */
    // \{

    /*!
     * \brief Applies a vector of point sources which are possibly solution dependent.
     *
     * \param pointSources A vector of Dumux::PointSource s that contain
              source values for all phases and space positions.
     *
     * For this method, the \a values method of the point source
     * has to return the absolute mass rate in kg/s. Positive values mean
     * that mass is created, negative ones mean that it vanishes.
     */
    void addPointSources(std::vector<PointSource>& pointSources) const
    { pointSources = this->couplingManager().bulkPointSources(); }

    /*!
     * \brief Evaluates the point sources (added by addPointSources)
     *        for all phases within a given sub control volume.
     *
     * This is the method for the case where the point source is
     * solution dependent and requires some quantities that
     * are specific to the fully-implicit method.
     *
     * \param source A single point source
     * \param element The finite element
     * \param fvGeometry The finite-volume geometry
     * \param elemVolVars All volume variables for the element
     * \param scv The sub control volume within the element
     *
     * For this method, the \a values() method of the point sources returns
     * the absolute rate mass generated or annihilated in kg/s. Positive values mean
     * that mass is created, negative ones mean that it vanishes.
     */
    template<class ElementVolumeVariables>
    void pointSource(PointSource& source,
                     const Element &element,
                     const FVElementGeometry& fvGeometry,
                     const ElementVolumeVariables& elemVolVars,
                     const SubControlVolume &scv) const
    {
        NumEqVector sourceValues;

        // compute source at every integration point
        const auto priVars3D = this->couplingManager().bulkPriVars(source.id());
        const auto priVars1D = this->couplingManager().lowDimPriVars(source.id());
        const Scalar pressure3D = priVars3D[pressureIdx];
        const Scalar pressure1D = priVars1D[pressureIdx];

        const auto& spatialParams = this->couplingManager().problem(Dune::index_constant<1>{}).spatialParams();
        const auto lowDimElementIdx = this->couplingManager().pointSourceData(source.id()).lowDimElementIdx();
        const Scalar Kr = spatialParams.Kr(lowDimElementIdx);
        const Scalar rootRadius = spatialParams.radius(lowDimElementIdx);

        // sink defined as radial flow Jr * density [m^2 s-1]* [kg m-3]
        const auto molarDensityH20 = 1000 / 0.018;
        sourceValues[conti0EqIdx] = 2 * M_PI * rootRadius * Kr * (pressure1D - pressure3D) * molarDensityH20;

        const Scalar x3D = priVars3D[transportCompIdx];
        const Scalar x1D = priVars1D[transportCompIdx];

        //! Advective transport over root wall
        // compute correct upwind concentration
        if (sourceValues[conti0EqIdx] > 0)
            sourceValues[transportEqIdx] = sourceValues[conti0EqIdx]*x1D;
        else
            sourceValues[transportEqIdx] = sourceValues[conti0EqIdx]*x3D;

        //! Diffusive transport over root wall
        const auto molarDensityD20 = 1000 / 0.020;
        sourceValues[transportEqIdx] += 2 * M_PI * rootRadius * 1.0e-8 * (x1D - x3D) * molarDensityD20;

        sourceValues *= source.quadratureWeight()*source.integrationElement();
        source = sourceValues;
    }

    /*!
     * \brief Evaluates the initial value for a control volume.
     *
     * \param globalPos The position for which the initial condition should be evaluated
     *
     * For this method, the \a values parameter stores primary
     * variables.
     */
    PrimaryVariables initialAtPos(const GlobalPosition &globalPos) const
    {
        const auto& gg = this->gridGeometry();
        static const Scalar extend = 0.15*(gg.bBoxMax()[0]-gg.bBoxMin()[0]);
        const auto xTracer = [&]()
        {
            auto contaminationPos = gg.bBoxMax()-gg.bBoxMin();
            contaminationPos[0] *= 0.26;
            contaminationPos[1] *= 0.56;
            contaminationPos[2] *= 0.26;
            contaminationPos += gg.bBoxMin();

            if ((globalPos - contaminationPos).infinity_norm() <  extend + eps_)
                return contaminantMoleFraction_;
            else
                return 0.0;
        }();

        PrimaryVariables priVars(0.0);
        //! Hydrostatic pressure profile
        priVars[pressureIdx] = (nonwettingReferencePressure() - pcTop_)
                                -9.81*1000*(globalPos[dimWorld-1] - gg.bBoxMax()[dimWorld-1]);
        priVars[transportCompIdx] = xTracer;
        return priVars;
    }

    //! Get the coupling manager
    const CouplingManager& couplingManager() const
    { return *couplingManager_; }

private:
    Scalar pcTop_, contaminantMoleFraction_;

    static constexpr Scalar eps_ = 1.5e-7;
    std::string name_;

    std::shared_ptr<CouplingManager> couplingManager_;
};

} // end namespace Dumux

#endif
