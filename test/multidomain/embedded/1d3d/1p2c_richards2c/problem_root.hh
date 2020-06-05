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
/*!
 * \file
 * \ingroup EmbeddedTests
 * \brief A test problem for the one-phase root model:
 * Sap is flowing through a 1d network root xylem.
 */

#ifndef DUMUX_ROOT_PROBLEM_HH
#define DUMUX_ROOT_PROBLEM_HH

#include <dune/foamgrid/foamgrid.hh>

#include <dumux/common/boundarytypes.hh>
#include <dumux/common/parameters.hh>
#include <dumux/common/properties.hh>
#include <dumux/discretization/cctpfa.hh>

#include <dumux/porousmediumflow/1pnc/model.hh>
#include <dumux/porousmediumflow/problem.hh>
#include <dumux/material/components/simpleh2o.hh>
#include <dumux/material/components/constant.hh>
#include <dumux/material/fluidsystems/liquidphase2c.hh>

#include "spatialparams_root.hh"

namespace Dumux {
// forward declaration
template <class TypeTag> class RootProblem;

namespace Properties {

// Create new type tags
namespace TTag {
struct Root { using InheritsFrom = std::tuple<OnePNC, CCTpfaModel>; };
} // end namespace TTag

// Set the grid type
template<class TypeTag>
struct Grid<TypeTag, TTag::Root> { using type = Dune::FoamGrid<1, 3>; };

template<class TypeTag>
struct EnableGridGeometryCache<TypeTag, TTag::Root> { static constexpr bool value = true; };
template<class TypeTag>
struct EnableGridVolumeVariablesCache<TypeTag, TTag::Root> { static constexpr bool value = true; };
template<class TypeTag>
struct EnableGridFluxVariablesCache<TypeTag, TTag::Root> { static constexpr bool value = true; };
template<class TypeTag>
struct SolutionDependentAdvection<TypeTag, TTag::Root> { static constexpr bool value = false; };
template<class TypeTag>
struct SolutionDependentMolecularDiffusion<TypeTag, TTag::Root> { static constexpr bool value = false; };
template<class TypeTag>
struct SolutionDependentHeatConduction<TypeTag, TTag::Root> { static constexpr bool value = false; };

// Set the problem property
template<class TypeTag>
struct Problem<TypeTag, TTag::Root> { using type = RootProblem<TypeTag>; };

// Set the fluid system
template<class TypeTag>
struct FluidSystem<TypeTag, TTag::Root>
{
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using type = FluidSystems::LiquidPhaseTwoC<Scalar, Components::SimpleH2O<Scalar>,
                                                       Components::Constant<1, Scalar>>;
};

// Set the spatial parameters
template<class TypeTag>
struct SpatialParams<TypeTag, TTag::Root>
{
    using type = RootSpatialParams<GetPropType<TypeTag, Properties::GridGeometry>,
                                   GetPropType<TypeTag, Properties::Scalar>>;
};

template<class TypeTag>
struct UseMoles<TypeTag, TTag::Root> { static constexpr bool value = true; };

} // end namespace Properties

/*!
 * \ingroup EmbeddedTests
 * \brief Exact solution 1D-3D
 */
template <class TypeTag>
class RootProblem : public PorousMediumFlowProblem<TypeTag>
{
    using ParentType = PorousMediumFlowProblem<TypeTag>;
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using PointSource = GetPropType<TypeTag, Properties::PointSource>;
    using PrimaryVariables = GetPropType<TypeTag, Properties::PrimaryVariables>;
    using NeumannFluxes = GetPropType<TypeTag, Properties::NumEqVector>;
    using SourceValues = GetPropType<TypeTag, Properties::NumEqVector>;
    using BoundaryTypes = Dumux::BoundaryTypes<GetPropType<TypeTag, Properties::ModelTraits>::numEq()>;
    using GridGeometry = GetPropType<TypeTag, Properties::GridGeometry>;
    using GridView = typename GridGeometry::GridView;
    using FVElementGeometry = typename GridGeometry::LocalView;
    using SubControlVolume = typename GridGeometry::SubControlVolume;
    using SubControlVolumeFace = typename GridGeometry::SubControlVolumeFace;
    using GlobalPosition = typename GridGeometry::GlobalCoordinate;
    using Element = typename GridView::template Codim<0>::Entity;
    using SolutionVector = GetPropType<TypeTag, Properties::SolutionVector>;
    using GridVariables = GetPropType<TypeTag, Properties::GridVariables>;

    using CouplingManager = GetPropType<TypeTag, Properties::CouplingManager>;

public:
    using FluidSystem = GetPropType<TypeTag, Properties::FluidSystem>;
    enum Indices {
        // Grid and world dimension
        dim = GridView::dimension,
        dimworld = GridView::dimensionworld,

        pressureIdx = 0,
        transportCompIdx = 1,

        conti0EqIdx = 0,
        transportEqIdx = 1,

        liquidPhaseIdx = 0
    };

    template<class SpatialParams>
    RootProblem(std::shared_ptr<const GridGeometry> gridGeometry,
                std::shared_ptr<SpatialParams> spatialParams,
                std::shared_ptr<CouplingManager> couplingManager)
    : ParentType(gridGeometry, spatialParams, "Root")
    , couplingManager_(couplingManager)
    {
        // read parameters from input file
        name_  =  getParam<std::string>("Vtk.OutputName") + "_" + getParamFromGroup<std::string>(this->paramGroup(), "Problem.Name");
        transpirationRate_ = getParam<Scalar>("BoundaryConditions.TranspirationRate");
        initPressure_ = getParam<Scalar>("BoundaryConditions.InitialRootPressure");
    }

    /*!
     * \brief Returns how much the domain is extruded at a given sub-control volume.
     *
     * The extrusion factor here extrudes the 1d line to a circular tube with
     * cross-section area pi*r^2.
     */
    template<class ElementSolution>
    Scalar extrusionFactor(const Element &element,
                           const SubControlVolume &scv,
                           const ElementSolution& elemSol) const
    {
        const auto eIdx = this->gridGeometry().elementMapper().index(element);
        const auto radius = this->spatialParams().radius(eIdx);
        return M_PI*radius*radius;
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
     * \brief Returns the temperature within the domain in [K].
     */
    Scalar temperature() const
    { return 273.15 + 10.0; }

    // \}
    /*!
     * \name Boundary conditions
     */
    // \{

   /*!
     * \brief Specifies which kind of boundary condition should be
     *        used for which equation on a given boundary segment.
     *
     * \param globalPos The global position
     */
    BoundaryTypes boundaryTypesAtPos(const GlobalPosition &globalPos) const
    {
        BoundaryTypes values;
        values.setAllNeumann();
        return values;
    }

    /*!
     * \brief Evaluates the boundary conditions for a Dirichlet control volume.
     *
     * \param globalPos The global position
     *
     * For this method, the \a values parameter stores primary variables.
     */
    PrimaryVariables dirichletAtPos(const GlobalPosition& globalPos) const
    { return initialAtPos(globalPos); }


    /*!
     * \brief Evaluates the boundary conditions for a Neumann boundary segment.
     *
     * For this method, the \a priVars parameter stores the mass flux
     * in normal direction of each component. Negative values mean
     * influx.
     */
    template<class ElementVolumeVariables, class ElementFluxVarsCache>
    NeumannFluxes neumann(const Element& element,
                          const FVElementGeometry& fvGeometry,
                          const ElementVolumeVariables& elemVolvars,
                          const ElementFluxVarsCache& elemFluxVarsCache,
                          const SubControlVolumeFace& scvf) const
    {
        NeumannFluxes values(0.0);
        if (scvf.center()[2] + eps_ > this->gridGeometry().bBoxMax()[2])
        {
            const auto& volVars = elemVolvars[scvf.insideScvIdx()];
            const Scalar value = transpirationRate_ * volVars.molarDensity(liquidPhaseIdx)/volVars.density(liquidPhaseIdx);

            values[conti0EqIdx] = value / volVars.extrusionFactor() / scvf.area();
            // use upwind mole fraction to get outflow condition for the tracer
            values[transportEqIdx] = values[conti0EqIdx] * volVars.moleFraction(liquidPhaseIdx, transportCompIdx);
        }
        return values;

    }

    // \}

    /*!
     * \name Volume terms
     */
    // \{

     /*!
     * \brief Applies a vector of point sources which are possibly solution dependent.
     *
     * \param pointSources A vector of PointSource s that contain
              source values for all phases and space positions.
     *
     * For this method, the \a values method of the point source
     * has to return the absolute mass rate in kg/s. Positive values mean
     * that mass is created, negative ones mean that it vanishes.
     */
    void addPointSources(std::vector<PointSource>& pointSources) const
    { pointSources = this->couplingManager().lowDimPointSources(); }

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
        SourceValues sourceValues;

        // compute source at every integration point
        const auto priVars3D = this->couplingManager().bulkPriVars(source.id());
        const auto priVars1D = this->couplingManager().lowDimPriVars(source.id());
        const Scalar pressure3D = priVars3D[pressureIdx];
        const Scalar pressure1D = priVars1D[pressureIdx];

        const auto lowDimElementIdx = this->couplingManager().pointSourceData(source.id()).lowDimElementIdx();
        const Scalar Kr = this->spatialParams().Kr(lowDimElementIdx);
        const Scalar rootRadius = this->spatialParams().radius(lowDimElementIdx);

        // sink defined as radial flow Jr * density [m^2 s-1]* [kg m-3]
        const auto molarDensityH20 = 1000 / 0.018;
        sourceValues[conti0EqIdx] = 2 * M_PI * rootRadius * Kr * (pressure3D - pressure1D) * molarDensityH20;

        const Scalar x3D = priVars3D[transportCompIdx];
        const Scalar x1D = priVars1D[transportCompIdx];

        //! Advective transport over root wall
        // compute correct upwind concentration
        if (sourceValues[conti0EqIdx] > 0)
            sourceValues[transportEqIdx] = sourceValues[conti0EqIdx]*x3D;
        else
            sourceValues[transportEqIdx] = sourceValues[conti0EqIdx]*x1D;

        //! Diffusive transport over root wall
        const auto molarDensityD20 = 1000 / 0.020;
        sourceValues[transportEqIdx] += 2 * M_PI * rootRadius * 1.0e-8 * (x3D - x1D) * molarDensityD20;

        sourceValues *= source.quadratureWeight()*source.integrationElement();
        source = sourceValues;
    }

    /*!
     * \brief Evaluates the initial value for a control volume.
     *
     * For this method, the \a priVars parameter stores primary
     * variables.
     */
    PrimaryVariables initialAtPos(const GlobalPosition &globalPos) const
    { return PrimaryVariables({initPressure_, 0.0}); }

    // \}

    /*!
     * \brief Adds additional VTK output data to the VTKWriter.
     *
     * Function is called by the output module on every write.
     */
    template<class VtkOutputModule>
    void addVtkOutputFields(VtkOutputModule& vtk) const
    {
        vtk.addField(this->spatialParams().getRadii(), "radius");
    }

    //! Get the coupling manager
    const CouplingManager& couplingManager() const
    { return *couplingManager_; }

private:
    Scalar transpirationRate_, initPressure_;

    static constexpr Scalar eps_ = 1.5e-7;
    std::string name_;

    std::shared_ptr<CouplingManager> couplingManager_;
};

} // end namespace Dumux

#endif
