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
 * \brief A sub problem for the richards problem
 */
#ifndef DUMUX_RICHARDS_TEST_PROBLEM_HH
#define DUMUX_RICHARDS_TEST_PROBLEM_HH

#include <cmath>

#include <dumux/implicit/cellcentered/tpfa/properties.hh>
#include <dumux/porousmediumflow/implicit/problem.hh>
#include <dumux/porousmediumflow/richardsnc/implicit/model.hh>
#include <dumux/material/components/simpleh2o.hh>
#include <dumux/material/components/constant.hh>
#include <dumux/material/fluidsystems/liquidphase2c.hh>

//! get the properties needed for subproblems
#include <dumux/mixeddimension/subproblemproperties.hh>

#include "richardstestspatialparams.hh"

namespace Dumux
{
template <class TypeTag>
class RichardsTestProblem;

namespace Properties
{
NEW_TYPE_TAG(RichardsTestProblem, INHERITS_FROM(RichardsNC, RichardsTestSpatialParams));
NEW_TYPE_TAG(RichardsTestBoxProblem, INHERITS_FROM(BoxModel, RichardsTestProblem));
NEW_TYPE_TAG(RichardsTestCCProblem, INHERITS_FROM(CCTpfaModel, RichardsTestProblem));

// Set the grid type
//SET_TYPE_PROP(RichardsTestProblem, Grid, Dune::YaspGrid<3, Dune::EquidistantOffsetCoordinates<double, 3> >);
SET_TYPE_PROP(RichardsTestProblem, Grid, Dune::UGGrid<3>);

SET_BOOL_PROP(RichardsTestProblem, EnableGlobalFVGeometryCache, true);
SET_BOOL_PROP(RichardsTestProblem, EnableGlobalVolumeVariablesCache, true);
SET_BOOL_PROP(RichardsTestProblem, EnableGlobalFluxVariablesCache, true);
SET_BOOL_PROP(RichardsTestProblem, SolutionDependentAdvection, false);
SET_BOOL_PROP(RichardsTestProblem, SolutionDependentMolecularDiffusion, false);
SET_BOOL_PROP(RichardsTestProblem, SolutionDependentHeatConduction, false);

// Set the problem property
SET_TYPE_PROP(RichardsTestProblem, Problem, RichardsTestProblem<TypeTag>);

// Set the fluid system
SET_PROP(RichardsTestProblem, FluidSystem)
{
    using Scalar = typename GET_PROP_TYPE(TypeTag, Scalar);
    using type = FluidSystems::LiquidPhaseTwoC<TypeTag, SimpleH2O<Scalar>, Constant<TypeTag,Scalar>>;
};

// Set the spatial parameters
SET_TYPE_PROP(RichardsTestProblem, SpatialParams, RichardsTestSpatialParams<TypeTag>);

// Enable gravity
SET_BOOL_PROP(RichardsTestProblem, ProblemEnableGravity, true);

// Enable velocity output
SET_BOOL_PROP(RichardsTestProblem, VtkAddVelocity, true);

// Set the grid parameter group
SET_STRING_PROP(RichardsTestProblem, GridParameterGroup, "SoilGrid");

// Use mole fractions
SET_BOOL_PROP(RichardsTestProblem, UseMoles, true);
}

/*!
 * \ingroup OnePBoxModel
 * \ingroup ImplicitTestProblems
 * \brief  Test problem for the one-phase model:
 * water is flowing from bottom to top through and around a low permeable lens.
 *
 * The domain is box shaped. All sides are closed (Neumann 0 boundary)
 * except the top and bottom boundaries (Dirichlet), where water is
 * flowing from bottom to top.
 *
 * In the middle of the domain, a lens with low permeability (\f$K=10e-12\f$)
 * compared to the surrounding material (\f$ K=10e-10\f$) is defined.
 *
 * To run the simulation execute the following line in shell:
 * <tt>./test_box1p -parameterFile test_box1p.input</tt> or
 * <tt>./test_cc1p -parameterFile test_cc1p.input</tt>
 *
 * The same parameter file can be also used for 3d simulation but you need to change line
 * <tt>typedef Dune::SGrid<2,2> type;</tt> to
 * <tt>typedef Dune::SGrid<3,3> type;</tt> in the problem file
 * and use <tt>1p_3d.dgf</tt> in the parameter file.
 */
template <class TypeTag>
class RichardsTestProblem : public ImplicitPorousMediaProblem<TypeTag>
{
    using ParentType = ImplicitPorousMediaProblem<TypeTag>;
    using GridView = typename GET_PROP_TYPE(TypeTag, GridView);
    using Scalar = typename GET_PROP_TYPE(TypeTag, Scalar);
    using FVElementGeometry = typename GET_PROP_TYPE(TypeTag, FVElementGeometry);
    using SubControlVolume = typename GET_PROP_TYPE(TypeTag, SubControlVolume);
    using PrimaryVariables = typename GET_PROP_TYPE(TypeTag, PrimaryVariables);
    using ElementVolumeVariables = typename GET_PROP_TYPE(TypeTag, ElementVolumeVariables);
    using BoundaryTypes = typename GET_PROP_TYPE(TypeTag, BoundaryTypes);
    using TimeManager = typename GET_PROP_TYPE(TypeTag, TimeManager);
    using PointSource = typename GET_PROP_TYPE(TypeTag, PointSource);
    using MaterialLaw = typename GET_PROP_TYPE(TypeTag, MaterialLaw);
    using Indices = typename GET_PROP_TYPE(TypeTag, Indices);

    enum {
        // Grid and world dimension
        dim = GridView::dimension,
        dimWorld = GridView::dimensionworld
    };
    enum { conti0EqIdx = Indices::conti0EqIdx,
           transportEqIdx = Indices::conti0EqIdx + 1 };
    enum { pressureIdx = Indices::pressureIdx };
    enum { wPhaseIdx = Indices::wPhaseIdx };
    enum { transportCompIdx = Indices::compMainIdx + 1 };

    using Element = typename GridView::template Codim<0>::Entity;
    using GlobalPosition = Dune::FieldVector<Scalar, dimWorld>;

    using GlobalProblemTypeTag = typename GET_PROP_TYPE(TypeTag, GlobalProblemTypeTag);
    using CouplingManager = typename GET_PROP_TYPE(GlobalProblemTypeTag, CouplingManager);

public:
    RichardsTestProblem(TimeManager &timeManager, const GridView &gridView)
    : ParentType(timeManager, gridView)
    {
        name_ = GET_RUNTIME_PARAM_FROM_GROUP(TypeTag, std::string, Problem, Name) + "-soil";
        contaminantMoleFraction_ = GET_RUNTIME_PARAM_FROM_GROUP(TypeTag, Scalar, Problem, ContaminantMoleFraction);

        // for initial conditions
        const Scalar sw = 0.3; // start with 30% saturation on top
        pcTop_ = MaterialLaw::pc(this->spatialParams().materialLawParamsAtPos(this->bBoxMax()), sw);
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

    /*!
     * \brief Return the temperature within the domain.
     *
     * This problem assumes a temperature of 10 degrees Celsius.
     */
    Scalar temperature() const
    { return 273.15 + 10; } // 10C
    /*
      * \brief Returns the reference pressure [Pa] of the non-wetting
     *        fluid phase within a finite volume
     *
     * This problem assumes a constant reference pressure of 1 bar.
     */
    Scalar nonWettingReferencePressure() const
    { return 1.0e5; }

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
    { pointSources = this->couplingManager().bulkPointSources(); }

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
        const auto& lowDimVolVars = this->couplingManager().lowDimVolVars(source.id());

        const auto& spatialParams = this->couplingManager().lowDimProblem().spatialParams();
        const unsigned int lowDimElementIdx = this->couplingManager().pointSourceData(source.id()).lowDimElementIdx();
        const Scalar Kr = spatialParams.Kr(lowDimElementIdx);
        const Scalar rootRadius = spatialParams.radius(lowDimElementIdx);

        // sink defined as radial flow Jr * density [m^2 s-1]* [kg m-3]
        PrimaryVariables sourceValues(0.0);
        sourceValues[conti0EqIdx] = 2* M_PI *rootRadius * Kr *(lowDimVolVars.pressure(wPhaseIdx) - bulkVolVars.pressure(wPhaseIdx))
                                    *bulkVolVars.molarDensity(wPhaseIdx);

        //! advective transport over root wall
        // compute correct upwind concentration
        if (sourceValues[conti0EqIdx] > 0)
            sourceValues[transportEqIdx] = sourceValues[conti0EqIdx]
                                            * lowDimVolVars.moleFraction(wPhaseIdx, transportCompIdx);
        else
            sourceValues[transportEqIdx] = sourceValues[conti0EqIdx]
                                            * bulkVolVars.moleFraction(wPhaseIdx, transportCompIdx);

        //! diffusive transport over root wall
        sourceValues[transportEqIdx] += 2* M_PI *rootRadius * 1.0e-8
                                        *(lowDimVolVars.moleFraction(wPhaseIdx, transportCompIdx) - bulkVolVars.moleFraction(wPhaseIdx, transportCompIdx))
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

        std::cout << "Global integrated source (soil): " << source[conti0EqIdx] << " (kg/s) / "
                  <<                           source[conti0EqIdx]*3600*24*1000 << " (g/day)" << '\n';
    }

    // \}
    /*!
     * \name Boundary conditions
     */
    // \{

    /*!
     * \brief Specifies which kind of boundary condition should be
     *        used for which equation on a given boundary control volume.
     *
     * \param globalPos The position of the center of the finite volume
     */
    BoundaryTypes boundaryTypesAtPos(const GlobalPosition &globalPos) const
    {
        BoundaryTypes bcTypes;
        bcTypes.setAllNeumann();
        return bcTypes;
    }

    /*!
     * \brief Evaluate the boundary conditions for a dirichlet
     *        control volume.
     *
     * \param globalPos The center of the finite volume which ought to be set.
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
     * \brief Evaluate the initial values for a control volume.
     *
     * For this method, the \a values parameter stores primary
     * variables.
     *
     * \param globalPos The position for which the boundary type is set
     */
    PrimaryVariables initialAtPos(const GlobalPosition &globalPos) const
    {
        const auto xTracer = [&,this]()
        {
            auto contaminationPos = this->bBoxMax()-this->bBoxMin();
            contaminationPos[0] *= 0.25;
            contaminationPos[1] *= 0.55;
            contaminationPos[2] *= 0.25;
            contaminationPos += this->bBoxMin();

            static const Scalar extend = 0.15*(this->bBoxMax()[0]-this->bBoxMin()[0]);
            if ((globalPos - contaminationPos).infinity_norm() <  extend + eps_)
                return contaminantMoleFraction_;
            else
                return 0.0;
        }();

        PrimaryVariables values(0.0);
        //! hydrostatic pressure profile
        values[pressureIdx] = (nonWettingReferencePressure() - pcTop_)
                              -9.81*1000*(globalPos[dimWorld-1] - this->bBoxMax()[dimWorld-1]);
        values[transportCompIdx] = xTracer;
        return values;

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
    std::shared_ptr<CouplingManager> couplingManager_;
    Scalar contaminantMoleFraction_;
    Scalar pcTop_;
    static constexpr Scalar eps_ = 1e-7;
};

} //end namespace Dumux

#endif
