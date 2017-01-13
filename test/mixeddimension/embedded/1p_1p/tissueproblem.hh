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
 * \brief Definition of a problem, for the 1p2c problem:
 * Component transport of oxygen in interstitial fluid.
 */
#ifndef DUMUX_TISSUE_PROBLEM_HH
#define DUMUX_TISSUE_PROBLEM_HH

#include <dune/geometry/quadraturerules.hh>
#include <dune/geometry/referenceelements.hh>
#include <dune/localfunctions/lagrange/pqkfactory.hh>

#include <dumux/common/math.hh>
#include <dumux/implicit/cellcentered/tpfa/properties.hh>
#include <dumux/porousmediumflow/1p/implicit/model.hh>
#include <dumux/porousmediumflow/implicit/problem.hh>
#include <dumux/material/components/constant.hh>
#include <dumux/material/fluidsystems/liquidphase.hh>

#include <dumux/mixeddimension/subproblemproperties.hh>

#include "tissuespatialparams.hh"

namespace Dumux
{

template <class TypeTag>
class TissueProblem;

namespace Properties
{
NEW_TYPE_TAG(TissueProblem, INHERITS_FROM(OneP));
NEW_TYPE_TAG(TissueCCProblem, INHERITS_FROM(CCTpfaModel, TissueProblem));

// Set the grid type
SET_TYPE_PROP(TissueProblem, Grid, Dune::YaspGrid<3, Dune::EquidistantOffsetCoordinates<typename GET_PROP_TYPE(TypeTag, Scalar), 3> >);

SET_BOOL_PROP(TissueProblem, EnableGlobalFVGeometryCache, true);
SET_BOOL_PROP(TissueProblem, EnableGlobalVolumeVariablesCache, true);
SET_BOOL_PROP(TissueProblem, EnableGlobalFluxVariablesCache, true);
SET_BOOL_PROP(TissueProblem, SolutionDependentAdvection, false);
SET_BOOL_PROP(TissueProblem, SolutionDependentMolecularDiffusion, false);
SET_BOOL_PROP(TissueProblem, SolutionDependentHeatConduction, false);

// Set the grid parameter group
SET_STRING_PROP(TissueProblem, GridParameterGroup, "TissueGrid");

// Set the problem property
SET_TYPE_PROP(TissueProblem, Problem, TissueProblem<TypeTag>);

SET_PROP(TissueProblem, Fluid)
{
private:
    using Scalar = typename GET_PROP_TYPE(TypeTag, Scalar);
public:
    using type = FluidSystems::LiquidPhase<Scalar, Constant<TypeTag, Scalar>>;
};

// Set the spatial parameters
SET_TYPE_PROP(TissueProblem, SpatialParams, Dumux::TissueSpatialParams<TypeTag>);


// Enable velocity output
SET_BOOL_PROP(TissueProblem, VtkAddVelocity, false);

// Disable gravity
SET_BOOL_PROP(TissueProblem, ProblemEnableGravity, false);
}


/*!
 * \ingroup OnePTwoCModel
 * \ingroup ImplicitTestProblems
 */
template <class TypeTag>
class TissueProblem : public ImplicitPorousMediaProblem<TypeTag>
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
    using Indices = typename GET_PROP_TYPE(TypeTag, Indices);

    // copy some indices for convenience
    enum {
        // world dimension
        dim = GridView::dimension,
        dimWorld = GridView::dimensionworld
    };
    enum { pressureIdx = Indices::pressureIdx };
    enum { conti0EqIdx = Indices::conti0EqIdx };

    enum { isBox = GET_PROP_VALUE(TypeTag, ImplicitIsBox) };
    enum { dofCodim = isBox ? dim : 0 };

    using Element = typename GridView::template Codim<0>::Entity;
    using GlobalPosition = Dune::FieldVector<Scalar, dimWorld>;

    using GlobalProblemTypeTag = typename GET_PROP_TYPE(TypeTag, GlobalProblemTypeTag);
    using CouplingManager = typename GET_PROP_TYPE(GlobalProblemTypeTag, CouplingManager);

public:
    TissueProblem(TimeManager &timeManager, const GridView &gridView)
        : ParentType(timeManager, gridView)
    {
        //read parameters from input file
        name_ = GET_RUNTIME_PARAM_FROM_GROUP(TypeTag, std::string, Problem, Name) + "_3d";
    }

    /*!
     * \brief Called by the TimeManager in order to
     *        initialize the problem.
     */
    void init()
    {
        ParentType::init();
        exactPressure_.resize(this->model().numDofs());
        for (const auto& element : elements(this->gridView()))
        {
            auto fvGeometry = localView(this->model().globalFvGeometry());
            fvGeometry.bindElement(element);

            for (auto&& scv : scvs(fvGeometry))
                exactPressure_[scv.dofIndex()] = exactSolution(scv.dofPosition());
        }
    }

    /*!
     * \brief Returns true if a restart file should be written to
     *        disk.
     */
    bool shouldWriteRestartFile() const
    { return false; }

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
     * This problem assumes a temperature of 37 degrees Celsius.
     */
    Scalar temperature() const
    { return 273.15 + 37; } // in [K]

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
        values.setAllDirichlet();
        return values;
    }

    /*!
     * \brief Evaluate the boundary conditions for a dirichlet
     *        boundary segment.
     *
     * \param globalPos The position for which the bc type should be evaluated
     *
     * For this method, the \a values parameter stores primary variables.
     */
    PrimaryVariables dirichletAtPos(const GlobalPosition &globalPos) const
    {
        PrimaryVariables values;
        values = exactSolution(globalPos);
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
        //const auto& bulkVolVars = elemVolVars[scv];
        const auto& bulkVolVars = this->couplingManager().bulkVolVars(source.id());
        const Scalar pressure1D = this->couplingManager().lowDimPriVars(source.id())[pressureIdx];

        // calculate the source
        const Scalar radius = this->couplingManager().radius(source.id());
        const Scalar beta = 2*M_PI/(2*M_PI + std::log(radius));
        const Scalar sourceValue = beta*(pressure1D - bulkVolVars.pressure())*bulkVolVars.density();
        source = sourceValue*source.quadratureWeight()*source.integrationElement();
    }

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
    { return PrimaryVariables(0.0); }

    //! The exact solution
    Scalar exactSolution(const GlobalPosition &globalPos) const
    {
        Dune::FieldVector<double, 2> xy({globalPos[0], globalPos[1]});
        return -1.0*(1+globalPos[2])/(2*M_PI)*std::log(xy.two_norm());

        // use this instead if you are using the coupling manager with circle sources
        // auto R = GET_RUNTIME_PARAM_FROM_GROUP(TypeTag, Scalar, SpatialParams, Radius);
        // if (xy.two_norm() > R)
        //     return -1.0*(1+globalPos[2])/(2*M_PI)*std::log(xy.two_norm());
        // else
        //    return -1.0*(1+globalPos[2])/(2*M_PI)*std::log(R);
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
                    auto pointSources = this->scvPointSources(element, fvGeometry, elemVolVars, scv);
                    pointSources *= scv.volume()*elemVolVars[scv].extrusionFactor();
                    source += pointSources;
                }
            }
        }

        std::cout << "Global integrated source (3D): " << source << '\n';
    }

    /*!
     * \brief Adds additional VTK output data to the VTKWriter. Function is called by the output module on every write.
     */
    void addVtkOutputFields(VtkOutputModule<TypeTag>& outputModule) const
    {
        auto& p = outputModule.createScalarField("exact pressure", dofCodim);
        p = exactPressure_;
    }

    //! Set the coupling manager
    void setCouplingManager(std::shared_ptr<CouplingManager> cm)
    { couplingManager_ = cm; }

    //! Get the coupling manager
    const CouplingManager& couplingManager() const
    { return *couplingManager_; }

    const std::vector<unsigned int>& getAdditionalDofDependencies(unsigned int dofGlobalIdx) const
    { return couplingManager().getAdditionalDofDependencies(dofGlobalIdx); }

private:
    std::vector<Scalar> exactPressure_;

    static constexpr Scalar eps_ = 1.5e-7;
    std::string name_;

    std::shared_ptr<CouplingManager> couplingManager_;
};

} //end namespace Dumux

#endif
