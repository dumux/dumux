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
 * \brief Channel flow test for the staggered grid (Navier-)Stokes model
 */
#ifndef DUMUX_DENSITY_FLOW_NC_TEST_PROBLEM_HH
#define DUMUX_DENSITY_FLOW_NC_TEST_PROBLEM_HH

#include <dumux/implicit/staggered/properties.hh>
#include <dumux/freeflow/staggerednc/model.hh>
#include <dumux/implicit/problem.hh>
#include <dumux/material/components/simpleh2o.hh>

#include <dumux/material/fluidsystems/h2oair.hh>

namespace Dumux
{
template <class TypeTag>
class DensityDrivenFlowProblem;

namespace Capabilities
{
    template<class TypeTag>
    struct isStationary<DensityDrivenFlowProblem<TypeTag>>
    { static const bool value = false; };
}

namespace Properties
{
NEW_TYPE_TAG(DensityDrivenFlowProblem, INHERITS_FROM(StaggeredModel, NavierStokesNC));

NEW_PROP_TAG(FluidSystem);

// Select the fluid system
SET_TYPE_PROP(DensityDrivenFlowProblem, FluidSystem,
              FluidSystems::H2OAir<typename GET_PROP_TYPE(TypeTag, Scalar)/*, SimpleH2O<typename GET_PROP_TYPE(TypeTag, Scalar)>, false*/>);

SET_PROP(DensityDrivenFlowProblem, PhaseIdx)
{
private:
    using FluidSystem = typename GET_PROP_TYPE(TypeTag, FluidSystem);
public:
    static constexpr int value = FluidSystem::wPhaseIdx;
};

SET_INT_PROP(DensityDrivenFlowProblem, ReplaceCompEqIdx, 0);

// Set the grid type
SET_TYPE_PROP(DensityDrivenFlowProblem, Grid, Dune::YaspGrid<2>);

// Set the problem property
SET_TYPE_PROP(DensityDrivenFlowProblem, Problem, Dumux::DensityDrivenFlowProblem<TypeTag> );

SET_BOOL_PROP(DensityDrivenFlowProblem, EnableGlobalFVGeometryCache, true);

SET_BOOL_PROP(DensityDrivenFlowProblem, EnableGlobalFluxVariablesCache, true);
SET_BOOL_PROP(DensityDrivenFlowProblem, EnableGlobalVolumeVariablesCache, true);

// Enable gravity
SET_BOOL_PROP(DensityDrivenFlowProblem, ProblemEnableGravity, true);
SET_BOOL_PROP(DensityDrivenFlowProblem, UseMoles, true);

#if ENABLE_NAVIERSTOKES
SET_BOOL_PROP(DensityDrivenFlowProblem, EnableInertiaTerms, true);
#else
SET_BOOL_PROP(DensityDrivenFlowProblem, EnableInertiaTerms, false);
#endif
}

/*!
 * \brief  Test problem for the one-phase model:
   \todo doc me!
 */
template <class TypeTag>
class DensityDrivenFlowProblem : public NavierStokesProblem<TypeTag>
{
    typedef NavierStokesProblem<TypeTag> ParentType;

    typedef typename GET_PROP_TYPE(TypeTag, GridView) GridView;
    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
    using FluidSystem = typename GET_PROP_TYPE(TypeTag, FluidSystem);

    // copy some indices for convenience
    typedef typename GET_PROP_TYPE(TypeTag, Indices) Indices;
    enum {
        // Grid and world dimension
        dim = GridView::dimension,
        dimWorld = GridView::dimensionworld
    };
    enum {
        massBalanceIdx = Indices::massBalanceIdx,
        transportEqIdx = 1,
        momentumBalanceIdx = Indices::momentumBalanceIdx,
        momentumXBalanceIdx = Indices::momentumXBalanceIdx,
        momentumYBalanceIdx = Indices::momentumYBalanceIdx,
        pressureIdx = Indices::pressureIdx,
        velocityXIdx = Indices::velocityXIdx,
        velocityYIdx = Indices::velocityYIdx,
        transportCompIdx = 1/*FluidSystem::wCompIdx*/
    };

    typedef typename GET_PROP_TYPE(TypeTag, BoundaryTypes) BoundaryTypes;
    typedef typename GET_PROP_TYPE(TypeTag, TimeManager) TimeManager;

    typedef typename GridView::template Codim<0>::Entity Element;
    typedef typename GridView::Intersection Intersection;

    typedef typename GET_PROP_TYPE(TypeTag, FVElementGeometry) FVElementGeometry;
    typedef typename GET_PROP_TYPE(TypeTag, SubControlVolume) SubControlVolume;

    typedef Dune::FieldVector<Scalar, dimWorld> GlobalPosition;


    using CellCenterPrimaryVariables = typename GET_PROP_TYPE(TypeTag, CellCenterPrimaryVariables);
    using FacePrimaryVariables = typename GET_PROP_TYPE(TypeTag, FacePrimaryVariables);

    using BoundaryValues = typename GET_PROP_TYPE(TypeTag, BoundaryValues);
    using InitialValues = typename GET_PROP_TYPE(TypeTag, BoundaryValues);
    using SourceValues = typename GET_PROP_TYPE(TypeTag, BoundaryValues);

public:
    DensityDrivenFlowProblem(TimeManager &timeManager, const GridView &gridView)
    : ParentType(timeManager, gridView)
    {
        name_ = GET_RUNTIME_PARAM_FROM_GROUP(TypeTag,
                                             std::string,
                                             Problem,
                                             Name);

        useWholeLength_ = GET_RUNTIME_PARAM_FROM_GROUP(TypeTag,
                                             bool,
                                             Problem,
                                             UseWholeLength);

        FluidSystem::init();
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
    std::string name() const
    {
        return name_;
    }

    bool shouldWriteRestartFile() const
    {
        return false;
    }

    /*!
     * \brief Return the temperature within the domain in [K].
     *
     * This problem assumes a temperature of 10 degrees Celsius.
     */
    Scalar temperature() const
    { return 273.15 + 10; } // 10C

    /*!
     * \brief Return the sources within the domain.
     *
     * \param globalPos The global position
     */
    SourceValues sourceAtPos(const GlobalPosition &globalPos) const
    {
        return SourceValues(0.0);
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
        BoundaryTypes values;

        // set Dirichlet values for the velocity everywhere
        values.setDirichlet(momentumBalanceIdx);
        values.setOutflow(transportEqIdx);
        values.setOutflow(massBalanceIdx);

        if(globalPos[1] <  eps_)
            values.setDirichlet(massBalanceIdx);

        if(globalPos[1] > this->bBoxMax()[1] - eps_)
        {
            if(useWholeLength_)
                values.setDirichlet(transportEqIdx);
            else
                if(globalPos[0] > 0.4 && globalPos[0] < 0.6)
                    values.setDirichlet(transportEqIdx);
        }

        return values;
    }

    /*!
     * \brief Evaluate the boundary conditions for a dirichlet
     *        control volume.
     *
     * \param globalPos The center of the finite volume which ought to be set.
     */
    BoundaryValues dirichletAtPos(const GlobalPosition &globalPos) const
    {
        BoundaryValues values;

        values[pressureIdx] = 1.1e+5;
        values[transportCompIdx] = 1e-3;
        values[velocityXIdx] = 0.0;
        values[velocityYIdx] = 0.0;

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
     * \param globalPos The global position
     */
    InitialValues initialAtPos(const GlobalPosition &globalPos) const
    {
        InitialValues values;
        values[pressureIdx] = 1.1e+5;
        values[transportCompIdx] = 0.0;
        values[velocityXIdx] = 0.0;
        values[velocityYIdx] = 0.0;

        return values;
    }

    /*!
     * \brief Adds additional VTK output data to the VTKWriter. Function is called by the output module on every write.
     */
    template<class VtkOutputModule>
    void addVtkOutputFields(VtkOutputModule& outputModule) const
    {
        auto& delRho = outputModule.createScalarField("delRho", 0);

        for (const auto& element : elements(this->gridView()))
        {
            auto fvGeometry = localView(this->model().globalFvGeometry());
            fvGeometry.bindElement(element);
            for (auto&& scv : scvs(fvGeometry))
            {
                auto ccDofIdx = scv.dofIndex();

                auto elemVolVars = localView(this->model().curGlobalVolVars());
                elemVolVars.bind(element, fvGeometry, this->model().curSol());

                delRho[ccDofIdx] = elemVolVars[scv].density() - 999.694;
            }
        }
    }

    // \}

private:
    const Scalar eps_{1e-6};
    std::string name_;
    bool useWholeLength_;
};
} //end namespace

#endif
