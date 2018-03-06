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
 *
 * \brief A test problem for the one-phase root model:
 * Sap is flowing through a 1d network root xylem.
 */
#ifndef DUMUX_ROOT_PROBLEM_HH
#define DUMUX_ROOT_PROBLEM_HH

#include <dumux/common/parameters.hh>
#include <dumux/common/properties.hh>
#include <dumux/discretization/cellcentered/tpfa/properties.hh>

#include <dumux/porousmediumflow/1pnc/model.hh>
#include <dumux/porousmediumflow/problem.hh>
#include <dumux/material/components/simpleh2o.hh>
#include <dumux/material/components/constant.hh>
#include <dumux/material/fluidsystems/liquidphase2c.hh>

#include <dumux/material/components/simpleh2o.hh>
#include <dumux/material/fluidsystems/liquidphase.hh>

#include "rootspatialparams.hh"

namespace Dumux {
// forward declaration
template <class TypeTag> class RootProblem;

namespace Properties {

NEW_TYPE_TAG(RootTypeTag, INHERITS_FROM(CCTpfaModel, OnePNC));

// Set the grid type
SET_TYPE_PROP(RootTypeTag, Grid, Dune::FoamGrid<1, 3>);

SET_BOOL_PROP(RootTypeTag, EnableFVGridGeometryCache, true);
SET_BOOL_PROP(RootTypeTag, EnableGridVolumeVariablesCache, true);
SET_BOOL_PROP(RootTypeTag, EnableGridFluxVariablesCache, true);
SET_BOOL_PROP(RootTypeTag, SolutionDependentAdvection, false);
SET_BOOL_PROP(RootTypeTag, SolutionDependentMolecularDiffusion, false);
SET_BOOL_PROP(RootTypeTag, SolutionDependentHeatConduction, false);

// Set the problem property
SET_TYPE_PROP(RootTypeTag, Problem, RootProblem<TypeTag>);

// Set the fluid system
SET_PROP(RootTypeTag, FluidSystem)
{
    using Scalar = typename GET_PROP_TYPE(TypeTag, Scalar);
    using type = FluidSystems::LiquidPhaseTwoC<Scalar, Components::SimpleH2O<Scalar>,
                                                       Components::Constant<1, Scalar>>;
};

// Set the spatial parameters
SET_TYPE_PROP(RootTypeTag, SpatialParams, RootSpatialParams<TypeTag>);

SET_BOOL_PROP(RootTypeTag, UseMoles, true);

} // end namespace Properties

/*!
 * \ingroup OnePTests
 * \brief Exact solution 1D-3D
 */
template <class TypeTag>
class RootProblem : public PorousMediumFlowProblem<TypeTag>
{
    using ParentType = PorousMediumFlowProblem<TypeTag>;
    using GridView = typename GET_PROP_TYPE(TypeTag, GridView);
    using Scalar = typename GET_PROP_TYPE(TypeTag, Scalar);
    using ElementVolumeVariables = typename GET_PROP_TYPE(TypeTag, ElementVolumeVariables);
    using VolumeVariables = typename GET_PROP_TYPE(TypeTag, VolumeVariables);
    using PointSource = typename GET_PROP_TYPE(TypeTag, PointSource);
    using PrimaryVariables = typename GET_PROP_TYPE(TypeTag, PrimaryVariables);
    using NeumannFluxes = typename GET_PROP_TYPE(TypeTag, NumEqVector);
    using SourceValues = typename GET_PROP_TYPE(TypeTag, NumEqVector);
    using BoundaryTypes = typename GET_PROP_TYPE(TypeTag, BoundaryTypes);
    using FVGridGeometry = typename GET_PROP_TYPE(TypeTag, FVGridGeometry);
    using FVElementGeometry = typename FVGridGeometry::LocalView;
    using SubControlVolume = typename FVGridGeometry::SubControlVolume;
    using SubControlVolumeFace = typename FVGridGeometry::SubControlVolumeFace;
    using ElementSolutionVector = typename GET_PROP_TYPE(TypeTag, ElementSolutionVector);
    using SolutionVector = typename GET_PROP_TYPE(TypeTag, SolutionVector);
    using GridVariables = typename GET_PROP_TYPE(TypeTag, GridVariables);
    using Element = typename GridView::template Codim<0>::Entity;
    using GlobalPosition = Dune::FieldVector<Scalar, GridView::dimensionworld>;

    using CouplingManager = typename GET_PROP_TYPE(TypeTag, CouplingManager);

public:
    using FluidSystem = typename GET_PROP_TYPE(TypeTag, FluidSystem);
    enum Indices {
        // Grid and world dimension
        dim = GridView::dimension,
        dimworld = GridView::dimensionworld,

        pressureIdx = 0,
        transportCompIdx = 1,

        conti0EqIdx = 0,
        transportEqIdx = 1,

        wPhaseIdx = 0
    };

    RootProblem(std::shared_ptr<const FVGridGeometry> fvGridGeometry,
                std::shared_ptr<CouplingManager> couplingManager)
    : ParentType(fvGridGeometry)
    , couplingManager_(couplingManager)
    {
        //read parameters from input file
        name_ = getParam<std::string>("Problem.Name") + "_1d";
        transpirationRate_ = getParam<Scalar>("BoundaryConditions.TranspirationRate");
        initPressure_ = getParam<Scalar>("BoundaryConditions.InitialRootPressure");
    }

    /*!
     * \brief Return how much the domain is extruded at a given sub-control volume.
     *
     * The extrusion factor here makes extrudes the 1d line to a circular tube with
     * cross-section area pi*r^2.
     */
    Scalar extrusionFactor(const Element &element,
                           const SubControlVolume &scv,
                           const ElementSolutionVector& elemSol) const
    {
        const auto eIdx = this->fvGridGeometry().elementMapper().index(element);
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
     * \brief Return the temperature within the domain in [K].
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
     * \brief Evaluate the boundary conditions for a dirichlet
     *        control volume.
     *
     * \param globalPos The global position
     *
     * For this method, the \a values parameter stores primary variables.
     */
    PrimaryVariables dirichletAtPos(const GlobalPosition& globalPos) const
    { return initialAtPos(globalPos); }


    /*!
     * \brief Evaluate the boundary conditions for a neumann
     *        boundary segment.
     *
     * For this method, the \a priVars parameter stores the mass flux
     * in normal direction of each component. Negative values mean
     * influx.
     */
    NeumannFluxes neumann(const Element& element,
                          const FVElementGeometry& fvGeometry,
                          const ElementVolumeVariables& elemVolvars,
                          const SubControlVolumeFace& scvf) const
    {
        NeumannFluxes values(0.0);
        if (scvf.center()[2] + eps_ > this->fvGridGeometry().bBoxMax()[2])
        {
            const auto& volVars = elemVolvars[scvf.insideScvIdx()];
            const Scalar value = transpirationRate_ * volVars.molarDensity(wPhaseIdx)/volVars.density(wPhaseIdx);

            values[conti0EqIdx] = value / volVars.extrusionFactor() / scvf.area();
            // use upwind mole fraction to get outflow condition for the tracer
            values[transportEqIdx] = values[conti0EqIdx] * volVars.moleFraction(wPhaseIdx, transportCompIdx);
        }
        return values;

    }

    // \}

    /*!
     * \name Volume terms
     */
    // \{

     /*!
     * \brief Applies a vector of point sources. The point sources
     *        are possibly solution dependent.
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
     * \brief Evaluate the point sources (added by addPointSources)
     *        for all phases within a given sub-control-volume.
     *
     * This is the method for the case where the point source is
     * solution dependent and requires some quantities that
     * are specific to the fully-implicit method.
     *
     * \param pointSource A single point source
     * \param element The finite element
     * \param fvGeometry The finite-volume geometry
     * \param elemVolVars All volume variables for the element
     * \param scv The sub-control volume within the element
     *
     * For this method, the \a values() method of the point sources returns
     * the absolute rate mass generated or annihilate in kg/s. Positive values mean
     * that mass is created, negative ones mean that it vanishes.
     */
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

        //! advective transport over root wall
        // compute correct upwind concentration
        if (sourceValues[conti0EqIdx] > 0)
            sourceValues[transportEqIdx] = sourceValues[conti0EqIdx]*x3D;
        else
            sourceValues[transportEqIdx] = sourceValues[conti0EqIdx]*x1D;

        //! diffusive transport over root wall
        const auto molarDensityD20 = 1000 / 0.020;
        sourceValues[transportEqIdx] += 2 * M_PI * rootRadius * 1.0e-8 * (x3D - x1D) * molarDensityD20;

        sourceValues *= source.quadratureWeight()*source.integrationElement();
        source = sourceValues;
    }

    /*!
     * \brief Evaluate the initial value for a control volume.
     *
     * For this method, the \a priVars parameter stores primary
     * variables.
     */
    PrimaryVariables initialAtPos(const GlobalPosition &globalPos) const
    { return PrimaryVariables({initPressure_, 0.0}); }

    // \}

    /*!
     * \brief Adds additional VTK output data to the VTKWriter. Function is called by the output module on every write.
     */
    template<class VtkOutputModule>
    void addVtkOutputFields(VtkOutputModule& vtk) const
    {
        vtk.addField(this->spatialParams().getRadii(), "radius");
    }

    //! Set the coupling manager
    void setCouplingManager(std::shared_ptr<CouplingManager> cm)
    { couplingManager_ = cm; }

    //! Get the coupling manager
    const CouplingManager& couplingManager() const
    { return *couplingManager_; }

private:
    Scalar transpirationRate_, initPressure_;

    static constexpr Scalar eps_ = 1.5e-7;
    std::string name_;

    std::shared_ptr<CouplingManager> couplingManager_;
};

} //end namespace Dumux

#endif
