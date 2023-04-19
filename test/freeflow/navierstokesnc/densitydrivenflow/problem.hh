// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup NavierStokesNCTests
 * \brief Test for the compositional staggered grid (Navier-)Stokes model.
 */

#ifndef DUMUX_DENSITY_FLOW_NC_TEST_PROBLEM_HH
#define DUMUX_DENSITY_FLOW_NC_TEST_PROBLEM_HH

#include <dumux/common/parameters.hh>
#include <dumux/common/properties.hh>
#include <dumux/common/numeqvector.hh>

#include <dumux/freeflow/navierstokes/boundarytypes.hh>
#include <dumux/freeflow/navierstokes/staggered/problem.hh>

namespace Dumux {

/*!
 * \ingroup NavierStokesNCTests
 * \brief  Test problem for the one-phase model.
 *
 * Here, a quadratic two-dimensional domain with closed and non-moving walls at
 * all sides is considered. Initially, the domain is filled with pure water.
 * At the top, a fixed concentration of the air component is set.
 * The air slowly dissolves in the water which leads to a local increase of density.
 * Due to the influence of gravity and
 * small numerical instabilities, fingers of denser water will form and sink downwards.
 */
template <class TypeTag>
class DensityDrivenFlowProblem : public NavierStokesStaggeredProblem<TypeTag>
{
    using ParentType = NavierStokesStaggeredProblem<TypeTag>;

    using BoundaryTypes = Dumux::NavierStokesBoundaryTypes<GetPropType<TypeTag, Properties::ModelTraits>::numEq()>;
    using FluidSystem = GetPropType<TypeTag, Properties::FluidSystem>;
    using GridGeometry = GetPropType<TypeTag, Properties::GridGeometry>;
    using Indices = typename GetPropType<TypeTag, Properties::ModelTraits>::Indices;
    using PrimaryVariables = GetPropType<TypeTag, Properties::PrimaryVariables>;
    using NumEqVector = Dumux::NumEqVector<PrimaryVariables>;
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;

    using Element = typename GridGeometry::GridView::template Codim<0>::Entity;
    using GlobalPosition = typename Element::Geometry::GlobalCoordinate;

    static constexpr auto dimWorld = GridGeometry::GridView::dimensionworld;

    static constexpr auto transportCompIdx = Indices::conti0EqIdx + 1;
    static constexpr auto transportEqIdx = Indices::conti0EqIdx + 1;

public:
    DensityDrivenFlowProblem(std::shared_ptr<const GridGeometry> gridGeometry)
    : ParentType(gridGeometry), eps_(1e-6)
    {
        useWholeLength_ = getParam<bool>("Problem.UseWholeLength");
        FluidSystem::init();
        deltaRho_.resize(this->gridGeometry().numCellCenterDofs());
    }

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

        if(globalPos[1] > this->gridGeometry().bBoxMax()[1] - eps_)
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
     * \brief Returns whether a fixed Dirichlet value shall be used at a given cell.
     *
     * \param element The finite element
     * \param fvGeometry The finite-volume geometry
     * \param scv The sub control volume
     * \param pvIdx The primary variable index in the solution vector
     */
    template<class FVElementGeometry, class SubControlVolume>
    bool isDirichletCell(const Element& element,
                         const FVElementGeometry& fvGeometry,
                         const SubControlVolume& scv,
                         int pvIdx) const
    {
        // set a fixed pressure in one cell
        return (isLowerLeftCell_(scv) && pvIdx == Indices::pressureIdx);
    }

   /*!
     * \brief Evaluates the boundary conditions for a Dirichlet control volume.
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
     * \brief Evaluates the initial value for a control volume.
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
     * \brief Adds additional VTK output data to the VTKWriter.
     *
     * Function is called by the output module on every write.
     *
     * \param gridVariables The grid variables
     * \param sol The solution vector
     */
    template<class GridVariables, class SolutionVector>
    void calculateDeltaRho(const GridVariables& gridVariables, const SolutionVector& sol)
    {
        auto fvGeometry = localView(this->gridGeometry());
        auto elemVolVars = localView(gridVariables.curGridVolVars());
        for (const auto& element : elements(this->gridGeometry().gridView()))
        {
            fvGeometry.bindElement(element);
            elemVolVars.bind(element, fvGeometry, sol);
            for (auto&& scv : scvs(fvGeometry))
            {
                const auto ccDofIdx = scv.dofIndex();
                deltaRho_[ccDofIdx] = elemVolVars[scv].density() - 999.694;
            }
        }
    }

    const std::vector<Scalar>& getDeltaRho() const
    { return deltaRho_; }



    // \}

private:

    template<class SubControlVolume>
    bool isLowerLeftCell_(const SubControlVolume& scv) const
    { return scv.dofIndex() == 0; }

    const Scalar eps_;
    bool useWholeLength_;
    std::vector<Scalar> deltaRho_;
};
} // end namespace Dumux

#endif
