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
 * \ingroup OnePTests
 * \brief Definition of a problem, for the 1p2c problem:
 * Component transport of oxygen in interstitial fluid.
 */
#ifndef DUMUX_TISSUE_PROBLEM_HH
#define DUMUX_TISSUE_PROBLEM_HH

#include <dune/geometry/quadraturerules.hh>
#include <dune/geometry/referenceelements.hh>
#include <dune/localfunctions/lagrange/pqkfactory.hh>

#include <dumux/common/math.hh>
#include <dumux/common/parameters.hh>
#include <dumux/common/properties.hh>
#include <dumux/discretization/cellcentered/tpfa/properties.hh>

#include <dumux/porousmediumflow/1p/model.hh>
#include <dumux/porousmediumflow/problem.hh>
#include <dumux/porousmediumflow/1p/incompressiblelocalresidual.hh>

#include <dumux/material/components/constant.hh>
#include <dumux/material/fluidsystems/1pliquid.hh>

#include "tissuespatialparams.hh"
#include "darcyslaw.hh"

namespace Dumux {

template <class TypeTag>
class TissueProblem;

namespace Properties {

NEW_TYPE_TAG(TissueTypeTag, INHERITS_FROM(OneP));
NEW_TYPE_TAG(TissueCCTypeTag, INHERITS_FROM(CCTpfaModel, TissueTypeTag));

// Set the grid type
SET_TYPE_PROP(TissueTypeTag, Grid, Dune::YaspGrid<3, Dune::EquidistantOffsetCoordinates<typename GET_PROP_TYPE(TypeTag, Scalar), 3> >);

SET_BOOL_PROP(TissueTypeTag, EnableFVGridGeometryCache, true);
SET_BOOL_PROP(TissueTypeTag, EnableGridVolumeVariablesCache, true);
SET_BOOL_PROP(TissueTypeTag, EnableGridFluxVariablesCache, true);
SET_BOOL_PROP(TissueTypeTag, SolutionDependentAdvection, false);
SET_BOOL_PROP(TissueTypeTag, SolutionDependentMolecularDiffusion, false);
SET_BOOL_PROP(TissueTypeTag, SolutionDependentHeatConduction, false);

// Set the problem property
SET_TYPE_PROP(TissueTypeTag, Problem, TissueProblem<TypeTag>);

SET_TYPE_PROP(TissueTypeTag, AdvectionType, MyCCTpfaDarcysLaw<typename GET_PROP_TYPE(TypeTag, Scalar), typename GET_PROP_TYPE(TypeTag, FVGridGeometry)>);

// Set the problem property
//SET_TYPE_PROP(TissueTypeTag, LocalResidual, OnePIncompressibleLocalResidual<TypeTag>);

// the fluid system
SET_PROP(TissueTypeTag, FluidSystem)
{
    using Scalar = typename GET_PROP_TYPE(TypeTag, Scalar);
    using type = FluidSystems::OnePLiquid<Scalar, Components::Constant<1, Scalar> >;
};

// Set the spatial parameters
SET_TYPE_PROP(TissueTypeTag, SpatialParams,
              TissueSpatialParams<typename GET_PROP_TYPE(TypeTag, FVGridGeometry),
                                  typename GET_PROP_TYPE(TypeTag, Scalar)>);
} // end namespace Properties


/*!
 * \ingroup OnePTests
 */
template <class TypeTag>
class TissueProblem : public PorousMediumFlowProblem<TypeTag>
{
    using ParentType = PorousMediumFlowProblem<TypeTag>;
    using Scalar = typename GET_PROP_TYPE(TypeTag, Scalar);
    using FVGridGeometry = typename GET_PROP_TYPE(TypeTag, FVGridGeometry);
    using FVElementGeometry = typename FVGridGeometry::LocalView;
    using SubControlVolume = typename FVGridGeometry::SubControlVolume;
    using GridView = typename FVGridGeometry::GridView;
    using PrimaryVariables = typename GET_PROP_TYPE(TypeTag, PrimaryVariables);
    using SolutionVector = typename GET_PROP_TYPE(TypeTag, SolutionVector);
    using GridVariables = typename GET_PROP_TYPE(TypeTag, GridVariables);
    using BoundaryTypes = typename GET_PROP_TYPE(TypeTag, BoundaryTypes);
    using PointSource = typename GET_PROP_TYPE(TypeTag, PointSource);
    using Indices = typename GET_PROP_TYPE(TypeTag, ModelTraits)::Indices;
    using Element = typename GridView::template Codim<0>::Entity;
    using GlobalPosition = typename FVGridGeometry::GlobalCoordinate;

    using CouplingManager = typename GET_PROP_TYPE(TypeTag, CouplingManager);

public:
    TissueProblem(std::shared_ptr<const FVGridGeometry> fvGridGeometry,
                  std::shared_ptr<CouplingManager> couplingManager)
    : ParentType(fvGridGeometry)
    , couplingManager_(couplingManager)
    {
        //read parameters from input file
        name_ = getParam<std::string>("Problem.Name") + "_3d";

        exactPressure_.resize(this->fvGridGeometry().numDofs());
        for (const auto& element : elements(this->fvGridGeometry().gridView()))
        {
            auto fvGeometry = localView(this->fvGridGeometry());
            fvGeometry.bindElement(element);

            for (auto&& scv : scvs(fvGeometry))
                exactPressure_[scv.dofIndex()] = exactSolution(scv.dofPosition());
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

        if (globalPos[2] > this->fvGridGeometry().bBoxMax()[2] - eps_
            || globalPos[2] < this->fvGridGeometry().bBoxMin()[2] + eps_)
            values.setAllNeumann();
        return values;
    }

    /*!
     * \brief Evaluate the boundary conditions for a dirichlet
     *        control volume face.
     *
     * \param element The finite element
     * \param scvf the sub control volume face
     * \note used for cell-centered discretization schemes
     *
     * The method returns the boundary types information.
     */
    template<class SubControlVolumeFace>
    PrimaryVariables dirichlet(const Element &element, const SubControlVolumeFace &scvf) const
    {
        const auto& globalPos = scvf.ipGlobal();
        PrimaryVariables values(0.0);
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
    template<class ElementVolumeVariables>
    void pointSource(PointSource& source,
                     const Element &element,
                     const FVElementGeometry& fvGeometry,
                     const ElementVolumeVariables& elemVolVars,
                     const SubControlVolume &scv) const
    {
        // compute source at every integration point
        static const Scalar pressure1D = getParam<Scalar>("BoundaryConditions1D.PressureInput");
        static const Scalar lp = getParam<Scalar>("Problem.Lp");
        static const Scalar rho = getParam<Scalar>("Component.LiquidDensity");

        // calculate the source
        const Scalar radius = this->couplingManager().radius(source.id());
        const Scalar pressure3D = this->couplingManager().bulkPriVars(source.id())[Indices::pressureIdx];// + exactSolution(GlobalPosition({radius, 0, 0}));
        const Scalar beta = 2*M_PI*radius*lp*rho;
        const Scalar sourceValue = beta*(pressure1D - pressure3D);
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
        static const Scalar pressure1D = getParam<Scalar>("BoundaryConditions1D.PressureInput");
        static const Scalar lp = getParam<Scalar>("Problem.Lp");
        static const Scalar k = getParam<Scalar>("SpatialParams.PermeabilityTissue");
        static const Scalar mu = getParam<Scalar>("Component.LiquidKinematicViscosity")*getParam<Scalar>("Component.LiquidDensity");
        static const Scalar r = getParam<Scalar>("SpatialParams.Radius");
        static const Scalar d = lp*mu/k;

        Dune::FieldVector<double, 2> xy({globalPos[0], globalPos[1]});
        const auto p = -d*r*pressure1D/(1.0 - d*r*std::log(r))*std::log(xy.two_norm());
        return std::isinf(p) ? 1e20 : p;

    }

    Scalar exactGradient(const GlobalPosition& globalPos) const
    {
        static const Scalar pressure1D = getParam<Scalar>("BoundaryConditions1D.PressureInput");
        static const Scalar lp = getParam<Scalar>("Problem.Lp");
        static const Scalar k = getParam<Scalar>("SpatialParams.PermeabilityTissue");
        static const Scalar mu = getParam<Scalar>("Component.LiquidKinematicViscosity")*getParam<Scalar>("Component.LiquidDensity");
        static const Scalar r = getParam<Scalar>("SpatialParams.Radius");
        static const Scalar d = lp*mu/k;

        Dune::FieldVector<double, 2> xy({globalPos[0], globalPos[1]});
        const auto gradp =  -d*r*pressure1D/(1.0 - d*r*std::log(r))/xy.two_norm();
        return std::isinf(gradp) ? 1e20 : gradp;
    }

    //! Called after every time step
    //! Output the total global exchange term
    void computeSourceIntegral(const SolutionVector& sol, const GridVariables& gridVars)
    {
        PrimaryVariables source(0.0);
        for (const auto& element : elements(this->fvGridGeometry().gridView()))
        {
            auto fvGeometry = localView(this->fvGridGeometry());
            fvGeometry.bindElement(element);

            auto elemVolVars = localView(gridVars.curGridVolVars());
            elemVolVars.bindElement(element, fvGeometry, sol);

            for (auto&& scv : scvs(fvGeometry))
            {
                auto pointSources = this->scvPointSources(element, fvGeometry, elemVolVars, scv);
                pointSources *= scv.volume()*elemVolVars[scv].extrusionFactor();
                source += pointSources;
            }
        }

        std::cout << "Global integrated source (3D): " << source << '\n';
    }

    /*!
     * \brief Adds additional VTK output data to the VTKWriter. Function is called by the output module on every write.
     */
    template<class VtkOutputModule>
    void addVtkOutputFields(VtkOutputModule& vtk) const
    {
        vtk.addField(exactPressure_, "exact pressure");
        vtk.addField(corPressure_, "corrected pressure");
    }

    template<class SolutionVector, class GridVariables>
    void updatePressure(const GridVariables& gridVariables, const SolutionVector& curSol)
    {
        corPressure_.resize(this->fvGridGeometry().numDofs());
        for (const auto& element : elements(this->fvGridGeometry().gridView()))
        {
            auto fvGeometry = localView(this->fvGridGeometry());
            fvGeometry.bindElement(element);

            auto elemVolVars = localView(gridVariables.curGridVolVars());
            elemVolVars.bindElement(element, fvGeometry, curSol);

            for (auto&& scv : scvs(fvGeometry))
            {
                auto key = std::make_pair(this->fvGridGeometry().elementMapper().index(element), 0);
                if (this->pointSourceMap().count(key))
                    corPressure_[scv.dofIndex()] = elemVolVars[scv].pressure() + exactSolution(scv.dofPosition());
                else
                    corPressure_[scv.dofIndex()] = elemVolVars[scv].pressure();
            }
        }
    }

    //! Set the coupling manager
    void setCouplingManager(std::shared_ptr<CouplingManager> cm)
    { couplingManager_ = cm; }

    //! Get the coupling manager
    const CouplingManager& couplingManager() const
    { return *couplingManager_; }
private:
    std::vector<Scalar> exactPressure_;
    std::vector<Scalar> corPressure_;

    static constexpr Scalar eps_ = 1.5e-7;
    std::string name_;

    std::shared_ptr<CouplingManager> couplingManager_;
};

} //end namespace Dumux

#endif
