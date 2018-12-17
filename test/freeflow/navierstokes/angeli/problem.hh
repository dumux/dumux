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
 * \ingroup NavierStokesTests
 * \brief Test for the instationary staggered grid Navier-Stokes model with analytical solution (Angeli et al., 2017)
 */
#ifndef DUMUX_ANGELI_TEST_PROBLEM_HH
#define DUMUX_ANGELI_TEST_PROBLEM_HH

#include <dune/grid/yaspgrid.hh>

#include <dumux/material/fluidsystems/1pliquid.hh>
#include <dumux/material/components/constant.hh>

#include <dumux/freeflow/navierstokes/problem.hh>
#include <dumux/discretization/staggered/freeflow/properties.hh>
#include <dumux/freeflow/navierstokes/model.hh>
#include "../l2error.hh"

namespace Dumux
{
template <class TypeTag>
class AngeliTestProblem;

namespace Properties
{
// Create new type tags
namespace TTag {
struct AngeliTest { using InheritsFrom = std::tuple<NavierStokes, StaggeredFreeFlowModel>; };
} // end namespace TTag

// the fluid system
template<class TypeTag>
struct FluidSystem<TypeTag, TTag::AngeliTest>
{
private:
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
public:
    using type = FluidSystems::OnePLiquid<Scalar, Components::Constant<1, Scalar> >;
};

// Set the grid type
template<class TypeTag>
struct Grid<TypeTag, TTag::AngeliTest> { using type = Dune::YaspGrid<2, Dune::EquidistantOffsetCoordinates<GetPropType<TypeTag, Properties::Scalar>, 2> >; };

// Set the problem property
template<class TypeTag>
struct Problem<TypeTag, TTag::AngeliTest> { using type = Dumux::AngeliTestProblem<TypeTag> ; };

template<class TypeTag>
struct EnableFVGridGeometryCache<TypeTag, TTag::AngeliTest> { static constexpr bool value = true; };

template<class TypeTag>
struct EnableGridFluxVariablesCache<TypeTag, TTag::AngeliTest> { static constexpr bool value = true; };
template<class TypeTag>
struct EnableGridVolumeVariablesCache<TypeTag, TTag::AngeliTest> { static constexpr bool value = true; };

}

/*!
 * \ingroup NavierStokesTests
 * \brief  Test problem for the staggered grid (Angeli 1947)
 * \todo doc me!
 */
template <class TypeTag>
class AngeliTestProblem : public NavierStokesProblem<TypeTag>
{
    using ParentType = NavierStokesProblem<TypeTag>;

    using BoundaryTypes = GetPropType<TypeTag, Properties::BoundaryTypes>;
    using FVGridGeometry = GetPropType<TypeTag, Properties::FVGridGeometry>;
    using Indices = typename GetPropType<TypeTag, Properties::ModelTraits>::Indices;
    using ModelTraits = GetPropType<TypeTag, Properties::ModelTraits>;
    using NumEqVector = GetPropType<TypeTag, Properties::NumEqVector>;
    using PrimaryVariables = GetPropType<TypeTag, Properties::PrimaryVariables>;
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using SolutionVector = GetPropType<TypeTag, Properties::SolutionVector>;
    using TimeLoopPtr = std::shared_ptr<TimeLoop<Scalar>>;

    static constexpr auto dimWorld = GetPropType<TypeTag, Properties::GridView>::dimensionworld;
    using Element = typename FVGridGeometry::GridView::template Codim<0>::Entity;
    using GlobalPosition = typename Element::Geometry::GlobalCoordinate;
    using VelocityVector = Dune::FieldVector<Scalar, dimWorld>;

public:
    AngeliTestProblem(std::shared_ptr<const FVGridGeometry> fvGridGeometry)
    : ParentType(fvGridGeometry), eps_(1e-6)
    {
        printL2Error_ = getParam<bool>("Problem.PrintL2Error");
        kinematicViscosity_ = getParam<Scalar>("Component.LiquidKinematicViscosity", 1.0);
    }

   /*!
     * \name Problem parameters
     */
    // \{

    bool shouldWriteRestartFile() const
    {
        return false;
    }

    void postTimeStep(const SolutionVector& curSol) const
    {
        if(printL2Error_)
        {
            using L2Error = NavierStokesTestL2Error<Scalar, ModelTraits, PrimaryVariables>;
            const auto l2error = L2Error::calculateL2Error(*this, curSol);
            const int numCellCenterDofs = this->fvGridGeometry().numCellCenterDofs();
            const int numFaceDofs = this->fvGridGeometry().numFaceDofs();
            std::cout << std::setprecision(8) << "** L2 error (abs/rel) for "
                    << std::setw(6) << numCellCenterDofs << " cc dofs and " << numFaceDofs << " face dofs (total: " << numCellCenterDofs + numFaceDofs << "): "
                    << std::scientific
                    << "L2(p) = " << l2error.first[Indices::pressureIdx] << " / " << l2error.second[Indices::pressureIdx]
                    << ", L2(vx) = " << l2error.first[Indices::velocityXIdx] << " / " << l2error.second[Indices::velocityXIdx]
                    << ", L2(vy) = " << l2error.first[Indices::velocityYIdx] << " / " << l2error.second[Indices::velocityYIdx]
                    << std::endl;
        }
    }

   /*!
     * \brief Return the temperature within the domain in [K].
     *
     * This problem assumes a temperature of 10 degrees Celsius.
     */
    Scalar temperature() const
    { return 298.0; }


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
    bool isDirichletCell(const Element& element,
                         const typename FVGridGeometry::LocalView& fvGeometry,
                         const typename FVGridGeometry::SubControlVolume& scv,
                         int pvIdx) const
    {
        // set a fixed pressure in all cells at the boundary
        for (const auto& scvf : scvfs(fvGeometry))
            if (scvf.boundary())
                return true;

        return false;
    }

   /*!
     * \brief Return dirichlet boundary values at a given position
     *
     * \param globalPos The global position
     */
    PrimaryVariables dirichletAtPos(const GlobalPosition & globalPos) const
    {
        // use the values of the analytical solution
        return analyticalSolution(globalPos, time());
    }

    /*!
     * \brief Return the analytical solution of the problem at a given position
     *
     * \param globalPos The global position
     * \param time The current simulation time
     */
    PrimaryVariables analyticalSolution(const GlobalPosition& globalPos, const Scalar time) const
    {
        const Scalar x = globalPos[0];
        const Scalar y = globalPos[1];

        const Scalar t = time + timeLoop_->timeStepSize();

        PrimaryVariables values;

        values[Indices::pressureIdx] = - 0.25 * std::exp(-10.0 * kinematicViscosity_ * M_PI * M_PI * t) * M_PI * M_PI * (4.0 * std::cos(2.0 * M_PI * x) + std::cos(4.0 * M_PI * y));
        values[Indices::velocityXIdx] = - 2.0 * M_PI * std::exp(- 5.0 * kinematicViscosity_ * M_PI * M_PI * t) * std::cos(M_PI * x) * std::sin(2.0 * M_PI * y);
        values[Indices::velocityYIdx] = M_PI * std::exp(- 5.0 * kinematicViscosity_ * M_PI * M_PI * t) * std::sin(M_PI * x) * std::cos(2.0 * M_PI * y);

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
        return analyticalSolution(globalPos, -timeLoop_->timeStepSize());
    }

   /*!
     * \brief Returns the analytical solution for the pressure
     */
    auto& getAnalyticalPressureSolution() const
    {
        return analyticalPressure_;
    }

   /*!
     * \brief Returns the analytical solution for the velocity
     */
    auto& getAnalyticalVelocitySolution() const
    {
        return analyticalVelocity_;
    }

   /*!
     * \brief Returns the analytical solution for the velocity at the faces
     */
    auto& getAnalyticalVelocitySolutionOnFace() const
    {
        return analyticalVelocityOnFace_;
    }

    void setTimeLoop(TimeLoopPtr timeLoop)
    {
        timeLoop_ = timeLoop;
        createAnalyticalSolution();
    }

    Scalar time() const
    {
        return timeLoop_->time();
    }

   /*!
     * \brief Adds additional VTK output data to the VTKWriter. Function is called by the output module on every write.
     */
    void createAnalyticalSolution()
    {
        analyticalPressure_.resize(this->fvGridGeometry().numCellCenterDofs());
        analyticalVelocity_.resize(this->fvGridGeometry().numCellCenterDofs());
        analyticalVelocityOnFace_.resize(this->fvGridGeometry().numFaceDofs());

        for (const auto& element : elements(this->fvGridGeometry().gridView()))
        {
            auto fvGeometry = localView(this->fvGridGeometry());
            fvGeometry.bindElement(element);
            for (auto&& scv : scvs(fvGeometry))
            {
                auto ccDofIdx = scv.dofIndex();
                auto ccDofPosition = scv.dofPosition();
                auto analyticalSolutionAtCc = analyticalSolution(ccDofPosition, time());

                // velocities on faces
                for (auto&& scvf : scvfs(fvGeometry))
                {
                    const auto faceDofIdx = scvf.dofIndex();
                    const auto faceDofPosition = scvf.center();
                    const auto dirIdx = scvf.directionIndex();
                    const auto analyticalSolutionAtFace = analyticalSolution(faceDofPosition, time());
                    analyticalVelocityOnFace_[faceDofIdx][dirIdx] = analyticalSolutionAtFace[Indices::velocity(dirIdx)];
                }

                analyticalPressure_[ccDofIdx] = analyticalSolutionAtCc[Indices::pressureIdx];

                for(int dirIdx = 0; dirIdx < ModelTraits::dim(); ++dirIdx)
                    analyticalVelocity_[ccDofIdx][dirIdx] = analyticalSolutionAtCc[Indices::velocity(dirIdx)];
            }
        }
    }

private:

    Scalar eps_;

    Scalar kinematicViscosity_;
    TimeLoopPtr timeLoop_;
    bool printL2Error_;
    std::vector<Scalar> analyticalPressure_;
    std::vector<VelocityVector> analyticalVelocity_;
    std::vector<VelocityVector> analyticalVelocityOnFace_;
};
} //end namespace

#endif
