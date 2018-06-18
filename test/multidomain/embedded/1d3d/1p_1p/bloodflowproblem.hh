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
 * \brief A test problem for the one-phase blood flow model:
 * Blood is flowing through a 1d network grid.
 */
#ifndef DUMUX_BLOOD_FLOW_PROBLEM_HH
#define DUMUX_BLOOD_FLOW_PROBLEM_HH

#include <dumux/common/parameters.hh>
#include <dumux/common/properties.hh>
#include <dumux/discretization/cellcentered/tpfa/properties.hh>

#include <dumux/porousmediumflow/1p/model.hh>
#include <dumux/porousmediumflow/problem.hh>
#include <dumux/porousmediumflow/1p/incompressiblelocalresidual.hh>

#include <dumux/material/components/constant.hh>
#include <dumux/material/fluidsystems/1pliquid.hh>

#include "bloodflowspatialparams.hh"

namespace Dumux {
// forward declaration
template <class TypeTag> class BloodFlowProblem;

namespace Properties {

NEW_TYPE_TAG(BloodFlowTypeTag, INHERITS_FROM(OneP));
NEW_TYPE_TAG(BloodFlowCCTypeTag, INHERITS_FROM(CCTpfaModel, BloodFlowTypeTag));

// Set the grid type
SET_TYPE_PROP(BloodFlowTypeTag, Grid, Dune::FoamGrid<1, 3>);

SET_BOOL_PROP(BloodFlowTypeTag, EnableFVGridGeometryCache, true);
SET_BOOL_PROP(BloodFlowTypeTag, EnableGridVolumeVariablesCache, true);
SET_BOOL_PROP(BloodFlowTypeTag, EnableGridFluxVariablesCache, true);
SET_BOOL_PROP(BloodFlowTypeTag, SolutionDependentAdvection, false);
SET_BOOL_PROP(BloodFlowTypeTag, SolutionDependentMolecularDiffusion, false);
SET_BOOL_PROP(BloodFlowTypeTag, SolutionDependentHeatConduction, false);

// Set the problem property
SET_TYPE_PROP(BloodFlowTypeTag, Problem, BloodFlowProblem<TypeTag>);

// the fluid system
SET_PROP(BloodFlowTypeTag, FluidSystem)
{
    using Scalar = typename GET_PROP_TYPE(TypeTag, Scalar);
    using type = FluidSystems::OnePLiquid<Scalar, Components::Constant<1, Scalar> >;
};

// Set the problem property
SET_TYPE_PROP(BloodFlowTypeTag, LocalResidual, OnePIncompressibleLocalResidual<TypeTag>);

// Set the spatial parameters
SET_TYPE_PROP(BloodFlowTypeTag, SpatialParams,
              BloodFlowSpatialParams<typename GET_PROP_TYPE(TypeTag, FVGridGeometry),
                                     typename GET_PROP_TYPE(TypeTag, Scalar)>);
} // end namespace Properties

/*!
 * \ingroup OnePTests
 * \brief Exact solution 1D-3D
 */
template <class TypeTag>
class BloodFlowProblem : public PorousMediumFlowProblem<TypeTag>
{
    using ParentType = PorousMediumFlowProblem<TypeTag>;
    using Scalar = typename GET_PROP_TYPE(TypeTag, Scalar);
    using PointSource = typename GET_PROP_TYPE(TypeTag, PointSource);
    using Indices = typename GET_PROP_TYPE(TypeTag, ModelTraits)::Indices;
    using PrimaryVariables = typename GET_PROP_TYPE(TypeTag, PrimaryVariables);
    using BoundaryTypes = typename GET_PROP_TYPE(TypeTag, BoundaryTypes);
    using FVGridGeometry = typename GET_PROP_TYPE(TypeTag, FVGridGeometry);
    using GridView = typename FVGridGeometry::GridView;
    using FVElementGeometry = typename FVGridGeometry::LocalView;
    using SubControlVolume = typename FVGridGeometry::SubControlVolume;
    using SolutionVector = typename GET_PROP_TYPE(TypeTag, SolutionVector);
    using GridVariables = typename GET_PROP_TYPE(TypeTag, GridVariables);
    using Element = typename GridView::template Codim<0>::Entity;
    using GlobalPosition = typename FVGridGeometry::GlobalCoordinate;
    using CouplingManager = typename GET_PROP_TYPE(TypeTag, CouplingManager);

public:
    BloodFlowProblem(std::shared_ptr<const FVGridGeometry> fvGridGeometry,
                     std::shared_ptr<CouplingManager> couplingManager)
    : ParentType(fvGridGeometry)
    , couplingManager_(couplingManager)
    {
        //read parameters from input file
        name_ = getParam<std::string>("Problem.Name") + "_1d";
        p_in_ = getParam<Scalar>("BoundaryConditions1D.PressureInput");
        delta_p_ = getParam<Scalar>("BoundaryConditions1D.DeltaPressure");
        sourceFactor_ = getParam<Scalar>("MixedDimension.SourceFactor", 1.0);
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
     * \brief Return how much the domain is extruded at a given sub-control volume.
     *
     * The extrusion factor here makes extrudes the 1d line to a circular tube with
     * cross-section area pi*r^2.
     */
    template<class ElementSolution>
    Scalar extrusionFactor(const Element &element,
                           const SubControlVolume &scv,
                           const ElementSolution& elemSol) const
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
     *
     */
    Scalar temperature() const
    { return 273.15 + 37.0; } // Body temperature

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
        values.setAllDirichlet();
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
    {
        PrimaryVariables values;
        if(globalPos[2] > 0.5)
            values[Indices::pressureIdx] = p_in_;
        else
            values[Indices::pressureIdx] = p_in_-delta_p_;
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
    template<class ElementVolumeVariables>
    void pointSource(PointSource& source,
                     const Element &element,
                     const FVElementGeometry& fvGeometry,
                     const ElementVolumeVariables& elemVolVars,
                     const SubControlVolume &scv) const
    {
        // compute source at every integration point
        // const auto& bulkVolVars = this->couplingManager().bulkVolVars(source.id());
        const Scalar pressure1D = this->couplingManager().lowDimPriVars(source.id())[Indices::pressureIdx];
        const Scalar pressure3D = this->couplingManager().bulkPriVars(source.id())[Indices::pressureIdx];

        // calculate the source
        const Scalar radius = this->couplingManager().radius(source.id());
        const Scalar beta = 2*M_PI/(2*M_PI + std::log(radius));
        const Scalar sourceValue = beta*(pressure3D - pressure1D)*sourceFactor_;//*bulkVolVars.density();

        source = sourceValue*source.quadratureWeight()*source.integrationElement();
    }

    //! evaluate coupling residual for the derivative bulk DOF with respect to low dim DOF
    //! we only need to evaluate the part of the residual that will be influence by the low dim DOF
    template<class MatrixBlock, class VolumeVariables>
    void addSourceDerivatives(MatrixBlock& block,
                              const Element& element,
                              const FVElementGeometry& fvGeometry,
                              const VolumeVariables& curElemVolVars,
                              const SubControlVolume& scv) const
    {
        const auto eIdx = this->fvGridGeometry().elementMapper().index(element);

        auto key = std::make_pair(eIdx, 0);
        if (this->pointSourceMap().count(key))
        {
            // call the solDependent function. Herein the user might fill/add values to the point sources
            // we make a copy of the local point sources here
            auto pointSources = this->pointSourceMap().at(key);

            // add the point source values to the local residual (negative sign is convention for source term)
            for (const auto& source : pointSources)
                block[0][0] -= this->couplingManager().pointSourceDerivative(source, Dune::index_constant<1>{}, Dune::index_constant<1>{});
        }
    }

    /*!
     * \brief Evaluate the initial value for a control volume.
     *
     * For this method, the \a priVars parameter stores primary
     * variables.
     */
    PrimaryVariables initialAtPos(const GlobalPosition &globalPos) const
    { return PrimaryVariables(0.0); }

    // \}

    //! The exact pressure solution
    Scalar exactSolution(const GlobalPosition &globalPos) const
    { return 1 + globalPos[2]; }

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

        std::cout << "Global integrated source (1D): " << source << '\n';
    }

    /*!
     * \brief Adds additional VTK output data to the VTKWriter. Function is called by the output module on every write.
     */
    template<class VtkOutputModule>
    void addVtkOutputFields(VtkOutputModule& vtk) const
    {
        vtk.addField(exactPressure_, "exact pressure");
    }

    //! Set the coupling manager
    void setCouplingManager(std::shared_ptr<CouplingManager> cm)
    { couplingManager_ = cm; }

    //! Get the coupling manager
    const CouplingManager& couplingManager() const
    { return *couplingManager_; }

private:
    std::vector<Scalar> exactPressure_;

    static constexpr Scalar eps_ = 1.5e-7;
    std::string name_;

    Scalar p_in_;
    Scalar delta_p_;
    Scalar sourceFactor_;

    std::shared_ptr<CouplingManager> couplingManager_;
};

} //end namespace Dumux

#endif
