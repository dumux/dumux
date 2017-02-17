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
 * \brief A sub problem for the rootsystem
 */
#ifndef DUMUX_ROOTSYSTEM_TEST_PROBLEM_HH
#define DUMUX_ROOTSYSTEM_TEST_PROBLEM_HH

#include <cmath>

#include <dumux/implicit/cellcentered/tpfa/properties.hh>
#include <dumux/porousmediumflow/implicit/problem.hh>
#include <dumux/porousmediumflow/1p2c/implicit/model.hh>
#include <dumux/material/components/simpleh2o.hh>
#include <dumux/material/components/constant.hh>
#include <dumux/material/fluidsystems/liquidphase2c.hh>

//! get the properties needed for subproblems
#include <dumux/mixeddimension/subproblemproperties.hh>

#include "rootsystemtestspatialparams.hh"

namespace Dumux
{
template <class TypeTag>
class RootsystemTestProblem;

namespace Properties
{
NEW_TYPE_TAG(RootsystemTestProblem, INHERITS_FROM(OnePTwoC));
NEW_TYPE_TAG(RootsystemTestBoxProblem, INHERITS_FROM(BoxModel, RootsystemTestProblem));
NEW_TYPE_TAG(RootsystemTestCCProblem, INHERITS_FROM(CCTpfaModel, RootsystemTestProblem));

SET_PROP(RootsystemTestProblem, FluidSystem)
{
    using Scalar = typename GET_PROP_TYPE(TypeTag, Scalar);
    using type = FluidSystems::LiquidPhaseTwoC<TypeTag, SimpleH2O<Scalar>, Constant<TypeTag,Scalar>>;
};

// Set the grid type
SET_TYPE_PROP(RootsystemTestProblem, Grid, Dune::FoamGrid<1, 3>);

SET_BOOL_PROP(RootsystemTestProblem, EnableGlobalFVGeometryCache, true);
SET_BOOL_PROP(RootsystemTestProblem, EnableGlobalVolumeVariablesCache, true);
SET_BOOL_PROP(RootsystemTestProblem, EnableGlobalFluxVariablesCache, true);
SET_BOOL_PROP(RootsystemTestProblem, SolutionDependentAdvection, false);
SET_BOOL_PROP(RootsystemTestProblem, SolutionDependentMolecularDiffusion, false);
SET_BOOL_PROP(RootsystemTestProblem, SolutionDependentHeatConduction, false);

// Set the problem property
SET_TYPE_PROP(RootsystemTestProblem, Problem, RootsystemTestProblem<TypeTag>);

// Set the spatial parameters
SET_TYPE_PROP(RootsystemTestProblem, SpatialParams, RootsystemTestSpatialParams<TypeTag>);

// Enable gravity
SET_BOOL_PROP(RootsystemTestProblem, ProblemEnableGravity, true);

// Enable velocity output
SET_BOOL_PROP(RootsystemTestProblem, VtkAddVelocity, true);

// Use mole fractions
SET_BOOL_PROP(RootsystemTestProblem, UseMoles, true);
}

/*!
 * \ingroup OneDRootSystem
 * \ingroup ImplicitTestProblems
 * \brief TODO
 */
template <class TypeTag>
class RootsystemTestProblem : public ImplicitPorousMediaProblem<TypeTag>
{
    using ParentType = ImplicitPorousMediaProblem<TypeTag>;
    using GridView = typename GET_PROP_TYPE(TypeTag, GridView);
    using Scalar = typename GET_PROP_TYPE(TypeTag, Scalar);
    using ElementVolumeVariables = typename GET_PROP_TYPE(TypeTag, ElementVolumeVariables);
    using PointSource = typename GET_PROP_TYPE(TypeTag, PointSource);
    using Indices = typename GET_PROP_TYPE(TypeTag, Indices);
    // copy some indices for convenience
    enum {
        // Grid and world dimension
        dim = GridView::dimension,
        dimworld = GridView::dimensionworld
    };
    enum {
        // indices of the primary variables
        conti0EqIdx = Indices::conti0EqIdx,
        transportEqIdx = Indices::transportEqIdx,
        pressureIdx = Indices::pressureIdx,
        massOrMoleFracIdx = Indices::massOrMoleFracIdx,
        transportCompIdx = Indices::transportCompIdx
    };

    static const int wPhaseIdx = GET_PROP_VALUE(TypeTag, PhaseIdx);

    using PrimaryVariables = typename GET_PROP_TYPE(TypeTag, PrimaryVariables);
    using BoundaryTypes = typename GET_PROP_TYPE(TypeTag, BoundaryTypes);
    using TimeManager = typename GET_PROP_TYPE(TypeTag, TimeManager);
    using FVElementGeometry = typename GET_PROP_TYPE(TypeTag, FVElementGeometry);
    using SubControlVolume = typename GET_PROP_TYPE(TypeTag, SubControlVolume);
    using SubControlVolumeFace = typename GET_PROP_TYPE(TypeTag, SubControlVolumeFace);
    using ElementSolutionVector = typename GET_PROP_TYPE(TypeTag, ElementSolutionVector);
    using Element = typename GridView::template Codim<0>::Entity;
    using GlobalPosition = Dune::FieldVector<Scalar, dimworld>;

    using GlobalProblemTypeTag = typename GET_PROP_TYPE(TypeTag, GlobalProblemTypeTag);
    using CouplingManager = typename GET_PROP_TYPE(GlobalProblemTypeTag, CouplingManager);

public:
    RootsystemTestProblem(TimeManager &timeManager, const GridView &gridView)
    : ParentType(timeManager, gridView)
    {
        name_ = GET_RUNTIME_PARAM_FROM_GROUP(TypeTag, std::string, Problem, Name) + "-root";
    }
    /*!
     * \name Problem parameters
     */
    // \{

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
        const auto eIdx = this->gridView().indexSet().index(element);
        const auto radius = this->spatialParams().radius(eIdx);
        return M_PI*radius*radius;
    }

    /*!
     * \brief The problem name.
     *
     * This is used as a prefix for files generated by the simulation.
     */
    const std::string& name() const
    {
        return name_;
    }

    /*!
     * \brief Return the temperature within the domain.
     *
     * This problem assumes a temperature of 10 degrees Celsius.
     */
    Scalar temperature() const
    { return 273.15 + 10; } // 10C


    // \}
    /*!
     * \name Boundary conditions
     */
    // \{

    BoundaryTypes boundaryTypesAtPos (const GlobalPosition &globalPos ) const
    {
        BoundaryTypes bcTypes;
        bcTypes.setAllNeumann();
        return bcTypes;
    }


    /*!
     * \brief Specifies which kind of boundary condition should be
     *        used for which equation on a given boundary control volume.
     *
     * \param values The boundary types for the conservation equations
     * \param globalPos The position of the center of the finite volume
     */

    /*!
     * \brief Evaluate the boundary conditions for a dirichlet
     *        control volume.
     *
     * \param values The dirichlet values for the primary variables
     * \param globalPos The center of the finite volume which ought to be set.
     *
     * For this method, the \a values parameter stores primary variables.
     */
    PrimaryVariables dirichletAtPos(const GlobalPosition& globalPos) const
    {
        PrimaryVariables values(0.0);
        if (globalPos[2] + eps_ >  this->bBoxMax()[2] )
        {
            values[pressureIdx] = GET_RUNTIME_PARAM(TypeTag, Scalar, BoundaryConditions.CriticalCollarPressure);
            values[massOrMoleFracIdx] = 0.0;
        }

        return values;
    }

    /*!
     * \brief Evaluate the boundary conditions for a neumann
     *        boundary segment.
     *
     * For this method, the \a priVars parameter stores the mass flux
     * in normal direction of each component. Negative values mean
     * influx.
     */
    PrimaryVariables neumann(const Element& element,
                             const FVElementGeometry& fvGeometry,
                             const ElementVolumeVariables& elemVolvars,
                             const SubControlVolumeFace& scvf) const
    {
        PrimaryVariables values(0.0);
        if (scvf.center()[2] + eps_ > this->bBoxMax()[2])
        {
            const auto& volVars = elemVolvars[scvf.insideScvIdx()];
            // conversion of transpiration rate from kg/s to mol/s
            static const Scalar tr = GET_RUNTIME_PARAM(TypeTag, Scalar, BoundaryConditions.TranspirationRate);
            const Scalar value = tr * volVars.molarDensity(wPhaseIdx)/volVars.density(wPhaseIdx);
            values[conti0EqIdx] = value / volVars.extrusionFactor() / scvf.area();
            // use upwind mole fraction to get outflow condition for the tracer
            values[transportEqIdx] = values[conti0EqIdx]
                                     * volVars.moleFraction(wPhaseIdx, transportCompIdx);
        }
        return values;

    }

    // \}

    /*!
     * \name Volume terms
     */
    // \{

    /*!
     * \brief Evaluate the initial value for a control volume.
     *
     * For this method, the \a priVars parameter stores primary
     * variables.
     */
    PrimaryVariables initialAtPos(const GlobalPosition& globalPos) const
    {
        PrimaryVariables values(0.0);
        values[pressureIdx] =  GET_RUNTIME_PARAM(TypeTag,
                                                 Scalar,
                                                 BoundaryConditions.InitialRootPressure);
        return values;
    }

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
        // compute source at every integration point
        const auto& bulkVolVars = this->couplingManager().bulkVolVars(source.id());
        const auto& lowDimVolVars = elemVolVars[scv];

        const unsigned int lowDimElementIdx = this->couplingManager().pointSourceData(source.id()).lowDimElementIdx();
        const Scalar Kr = this->spatialParams().Kr(lowDimElementIdx);
        const Scalar rootRadius = this->spatialParams().radius(lowDimElementIdx);

        // sink defined as radial flow Jr * density [m^2 s-1]* [kg m-3]
        PrimaryVariables sourceValues(0.0);
        sourceValues[conti0EqIdx] = 2* M_PI *rootRadius * Kr *(bulkVolVars.pressure(wPhaseIdx) - lowDimVolVars.pressure(wPhaseIdx))
                                    *bulkVolVars.molarDensity(wPhaseIdx);

        //! advective transport over root wall
        // compute correct upwind concentration
        if (sourceValues[conti0EqIdx] < 0)
            sourceValues[transportEqIdx] = sourceValues[conti0EqIdx]
                                            * lowDimVolVars.moleFraction(wPhaseIdx, transportCompIdx);
        else
            sourceValues[transportEqIdx] = sourceValues[conti0EqIdx]
                                            * bulkVolVars.moleFraction(wPhaseIdx, transportCompIdx);

        //! diffusive transport over root wall
        sourceValues[transportEqIdx] += 2* M_PI *rootRadius * 1.0e-8
                                        *(bulkVolVars.moleFraction(wPhaseIdx, transportCompIdx) - lowDimVolVars.moleFraction(wPhaseIdx, transportCompIdx))
                                        *0.5*(bulkVolVars.molarDensity(wPhaseIdx) + lowDimVolVars.molarDensity(wPhaseIdx));

        sourceValues *= source.quadratureWeight()*source.integrationElement();
        source = sourceValues;
    }

    //! Called after every time step
    //! Output the total global exchange term
    void postTimeStep()
    {
        ParentType::postTimeStep();

        PrimaryVariables source(0.0);

        if (!(this->timeManager().time() < 0.0))
        {
            for (const auto& element : elements(this->gridView()))
            {
                auto fvGeometry = localView(this->model().globalFvGeometry());
                fvGeometry.bindElement(element);

                auto elemVolVars = localView(this->model().curGlobalVolVars());
                elemVolVars.bindElement(element, fvGeometry, this->model().curSol());

                for (auto&& scv : scvs(fvGeometry))
                {
                    const auto& volVars = elemVolVars[scv];
                    auto pointSources = this->scvPointSources(element, fvGeometry, elemVolVars, scv);
                    // conversion to kg/s
                    pointSources *= scv.volume()*volVars.extrusionFactor()
                                    * volVars.density(wPhaseIdx) / volVars.molarDensity(wPhaseIdx);
                    source += pointSources;
                }
            }
        }

        std::cout << "Global integrated source (root): " << source[conti0EqIdx] << " (kg/s) / "
                  <<                           source[conti0EqIdx]*3600*24*1000 << " (g/day)" << '\n';
    }

    bool shouldWriteRestartFile() const
    {
        return false;
    }

    //! Set the coupling manager
    void setCouplingManager(std::shared_ptr<CouplingManager> cm)
    { couplingManager_ = cm; }

    //! Get the coupling manager
    const CouplingManager& couplingManager() const
    { return *couplingManager_; }

private:
    std::string name_;
    const Scalar eps_ = 1e-9;
    std::shared_ptr<CouplingManager> couplingManager_;
};

} //end namespace Dumux

#endif
