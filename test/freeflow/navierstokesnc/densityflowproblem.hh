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
 * \ingroup NavierStokesNCTests
 * \brief Channel flow test for the staggered grid (Navier-)Stokes model
 */
#ifndef DUMUX_DENSITY_FLOW_NC_TEST_PROBLEM_HH
#define DUMUX_DENSITY_FLOW_NC_TEST_PROBLEM_HH

#include <dune/grid/yaspgrid.hh>

#include <dumux/material/components/simpleh2o.hh>
#include <dumux/material/fluidsystems/h2oair.hh>
#include <dumux/material/fluidsystems/1padapter.hh>

#include <dumux/discretization/staggered/freeflow/properties.hh>
#include <dumux/freeflow/compositional/navierstokesncmodel.hh>
#include <dumux/freeflow/navierstokes/problem.hh>


namespace Dumux {

template <class TypeTag>
class DensityDrivenFlowProblem;

namespace Properties {

NEW_TYPE_TAG(DensityDrivenFlowTypeTag, INHERITS_FROM(StaggeredFreeFlowModel, NavierStokesNC));

// Select the fluid system
SET_PROP(DensityDrivenFlowTypeTag, FluidSystem)
{
    using H2OAir = FluidSystems::H2OAir<typename GET_PROP_TYPE(TypeTag, Scalar)>;
    static constexpr int phaseIdx = H2OAir::liquidPhaseIdx;
    using type = FluidSystems::OnePAdapter<H2OAir, phaseIdx>;
};

SET_INT_PROP(DensityDrivenFlowTypeTag, ReplaceCompEqIdx, 0);

// Set the grid type
SET_TYPE_PROP(DensityDrivenFlowTypeTag, Grid, Dune::YaspGrid<2>);

// Set the problem property
SET_TYPE_PROP(DensityDrivenFlowTypeTag, Problem, Dumux::DensityDrivenFlowProblem<TypeTag> );

SET_BOOL_PROP(DensityDrivenFlowTypeTag, EnableFVGridGeometryCache, true);

SET_BOOL_PROP(DensityDrivenFlowTypeTag, EnableGridFluxVariablesCache, true);
SET_BOOL_PROP(DensityDrivenFlowTypeTag, EnableGridVolumeVariablesCache, true);

SET_BOOL_PROP(DensityDrivenFlowTypeTag, UseMoles, true);

#if ENABLE_NAVIERSTOKES
SET_BOOL_PROP(DensityDrivenFlowTypeTag, EnableInertiaTerms, true);
#else
SET_BOOL_PROP(DensityDrivenFlowTypeTag, EnableInertiaTerms, false);
#endif
}

/*!
 * \ingroup NavierStokesNCTests
 * \brief  Test problem for the one-phase model.
 * \todo doc me!
 */
template <class TypeTag>
class DensityDrivenFlowProblem : public NavierStokesProblem<TypeTag>
{
    using ParentType = NavierStokesProblem<TypeTag>;

    using BoundaryTypes = typename GET_PROP_TYPE(TypeTag, BoundaryTypes);
    using FluidSystem = typename GET_PROP_TYPE(TypeTag, FluidSystem);
    using FVGridGeometry = typename GET_PROP_TYPE(TypeTag, FVGridGeometry);
    using Indices = typename GET_PROP_TYPE(TypeTag, ModelTraits)::Indices;
    using NumEqVector = typename GET_PROP_TYPE(TypeTag, NumEqVector);
    using PrimaryVariables = typename GET_PROP_TYPE(TypeTag, PrimaryVariables);
    using Scalar = typename GET_PROP_TYPE(TypeTag, Scalar);

    using Element = typename FVGridGeometry::GridView::template Codim<0>::Entity;
    using GlobalPosition = typename Element::Geometry::GlobalCoordinate;

    using TimeLoopPtr = std::shared_ptr<CheckPointTimeLoop<Scalar>>;

    static constexpr auto dimWorld = GET_PROP_TYPE(TypeTag, GridView)::dimensionworld;

    static constexpr auto transportCompIdx = Indices::conti0EqIdx + 1;
    static constexpr auto transportEqIdx = Indices::conti0EqIdx + 1;

public:
    DensityDrivenFlowProblem(std::shared_ptr<const FVGridGeometry> fvGridGeometry)
    : ParentType(fvGridGeometry), eps_(1e-6)
    {
        useWholeLength_ = getParam<bool>("Problem.UseWholeLength");
        FluidSystem::init();
        deltaRho_.resize(this->fvGridGeometry().numCellCenterDofs());

        using CellArray = std::array<unsigned int, dimWorld>;
        const CellArray numCells = getParam<CellArray>("Grid.Cells");
        cellSizeX_ = this->fvGridGeometry().bBoxMax()[0] / numCells[0];
    }

   /*!
     * \name Problem parameters
     */
    // \{

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
    NumEqVector sourceAtPos(const GlobalPosition &globalPos) const
    {
        return NumEqVector(0.0);
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
        values.setDirichlet(Indices::velocityXIdx);
        values.setDirichlet(Indices::velocityYIdx);
        values.setNeumann(Indices::conti0EqIdx);
        values.setNeumann(transportEqIdx);

        if(isLowerLeftCell_(globalPos))
            values.setDirichletCell(Indices::pressureIdx);

        if(globalPos[1] > this->fvGridGeometry().bBoxMax()[1] - eps_)
        {
            if(useWholeLength_)
                values.setDirichlet(transportCompIdx);
            else
                if(globalPos[0] > 0.4 && globalPos[0] < 0.6)
                    values.setDirichlet(transportCompIdx);
        }

        return values;
    }

   /*!
     * \brief Evaluate the boundary conditions for a dirichlet
     *        control volume.
     *
     * \param globalPos The center of the finite volume which ought to be set.
     */
    PrimaryVariables dirichletAtPos(const GlobalPosition &globalPos) const
    {
        PrimaryVariables values;

        values[Indices::pressureIdx] = 1.1e+5;
        values[transportCompIdx] = 1e-3;
        values[Indices::velocityXIdx] = 0.0;
        values[Indices::velocityYIdx] = 0.0;

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
    PrimaryVariables initialAtPos(const GlobalPosition &globalPos) const
    {
        PrimaryVariables values;
        values[Indices::pressureIdx] = 1.1e+5;
        values[transportCompIdx] = 0.0;
        values[Indices::velocityXIdx] = 0.0;
        values[Indices::velocityYIdx] = 0.0;

        return values;
    }

   /*!
     * \brief Adds additional VTK output data to the VTKWriter. Function is called by the output module on every write.
     *
     * \param gridVariables The grid variables
     * \param sol The solution vector
     */
    template<class GridVariables, class SolutionVector>
    void calculateDeltaRho(const GridVariables& gridVariables, const SolutionVector& sol)
    {
        for (const auto& element : elements(this->fvGridGeometry().gridView()))
        {
            auto fvGeometry = localView(this->fvGridGeometry());
            fvGeometry.bindElement(element);
            for (auto&& scv : scvs(fvGeometry))
            {
                auto ccDofIdx = scv.dofIndex();

                auto elemVolVars = localView(gridVariables.curGridVolVars());
                elemVolVars.bind(element, fvGeometry, sol);

                deltaRho_[ccDofIdx] = elemVolVars[scv].density() - 999.694;
            }
        }
    }

    auto& getDeltaRho() const
    { return deltaRho_; }


    bool isLowerLeftCell_(const GlobalPosition& globalPos) const
    {
        return globalPos[0] < (0.5*cellSizeX_ + eps_) && globalPos[1] < eps_;
    }

    // \}

private:
    const Scalar eps_;
    bool useWholeLength_;
    std::vector<Scalar> deltaRho_;
    Scalar cellSizeX_;
};
} //end namespace

#endif
