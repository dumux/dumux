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
 * \brief Test for the 1-D Navier-Stokes model with an analytical solution.
 *
 * \copydoc Dumux::NavierStokesAnalyticProblem
 */
#ifndef DUMUX_DONEA_TEST_PROBLEM_HH
#define DUMUX_DONEA_TEST_PROBLEM_HH

#include <dune/common/fmatrix.hh>
#include <dune/grid/yaspgrid.hh>

#include <dumux/discretization/staggered/freeflow/properties.hh>

#include <dumux/freeflow/navierstokes/boundarytypes.hh>
#include <dumux/freeflow/navierstokes/model.hh>
#include <dumux/freeflow/navierstokes/problem.hh>

#include <dumux/material/components/constant.hh>
#include <dumux/material/fluidsystems/1pliquid.hh>

#include "../../l2error.hh"

namespace Dumux {
template <class TypeTag>
class NavierStokesAnalyticProblem;

namespace Properties {
// Create new type tags
namespace TTag {
struct NavierStokesAnalytic { using InheritsFrom = std::tuple<NavierStokes, StaggeredFreeFlowModel>; };
} // end namespace TTag

// the fluid system
template<class TypeTag>
struct FluidSystem<TypeTag, TTag::NavierStokesAnalytic>
{
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using type = FluidSystems::OnePLiquid<Scalar, Components::Constant<1, Scalar> >;
};

// Set the grid type
template<class TypeTag>
struct Grid<TypeTag, TTag::NavierStokesAnalytic> { using type = Dune::YaspGrid<1>; };

// Set the problem property
template<class TypeTag>
struct Problem<TypeTag, TTag::NavierStokesAnalytic> { using type = Dumux::NavierStokesAnalyticProblem<TypeTag> ; };

template<class TypeTag>
struct EnableGridGeometryCache<TypeTag, TTag::NavierStokesAnalytic> { static constexpr bool value = true; };

template<class TypeTag>
struct EnableGridFluxVariablesCache<TypeTag, TTag::NavierStokesAnalytic> { static constexpr bool value = true; };
template<class TypeTag>
struct EnableGridVolumeVariablesCache<TypeTag, TTag::NavierStokesAnalytic> { static constexpr bool value = true; };

template<class TypeTag>
struct NormalizePressure<TypeTag, TTag::NavierStokesAnalytic> { static constexpr bool value = false; };
} // end namespace Properties

/*!
 * \ingroup NavierStokesTests
 * \brief Test for the 1-D Navier-Stokes model with an analytical solution.
 *
 * The 1-D analytic solution is given by
 * \f[ p = 2 - 2 \cdot x \f]
 * \f[ v_\text{x} = 2 \cdot x^3 \f].
 */
template <class TypeTag>
class NavierStokesAnalyticProblem : public NavierStokesProblem<TypeTag>
{
    using ParentType = NavierStokesProblem<TypeTag>;

    using BoundaryTypes = Dumux::NavierStokesBoundaryTypes<GetPropType<TypeTag, Properties::ModelTraits>::numEq()>;
    using GridGeometry = GetPropType<TypeTag, Properties::GridGeometry>;
    using Indices = typename GetPropType<TypeTag, Properties::ModelTraits>::Indices;
    using NumEqVector = GetPropType<TypeTag, Properties::NumEqVector>;
    using ModelTraits = GetPropType<TypeTag, Properties::ModelTraits>;
    using PrimaryVariables = GetPropType<TypeTag, Properties::PrimaryVariables>;
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using SolutionVector = GetPropType<TypeTag, Properties::SolutionVector>;

    static constexpr auto dimWorld = GridGeometry::GridView::dimensionworld;
    using Element = typename GridGeometry::GridView::template Codim<0>::Entity;
    using GlobalPosition = typename Element::Geometry::GlobalCoordinate;
    using DimVector = GlobalPosition;
    using DimMatrix = Dune::FieldMatrix<Scalar, dimWorld, dimWorld>;

public:
    NavierStokesAnalyticProblem(std::shared_ptr<const GridGeometry> gridGeometry)
    : ParentType(gridGeometry)
    {
        printL2Error_ = getParam<bool>("Problem.PrintL2Error");
        density_ = getParam<Scalar>("Component.LiquidDensity");
        kinematicViscosity_ = getParam<Scalar>("Component.LiquidKinematicViscosity");
        createAnalyticalSolution_();
    }

   /*!
     * \name Problem parameters
     */
    // \{

    void printL2Error(const SolutionVector& curSol) const
    {
        if(printL2Error_)
        {
            using L2Error = NavierStokesTestL2Error<Scalar, ModelTraits, PrimaryVariables>;
            const auto l2error = L2Error::calculateL2Error(*this, curSol);
            const int numCellCenterDofs = this->gridGeometry().numCellCenterDofs();
            const int numFaceDofs = this->gridGeometry().numFaceDofs();
            std::cout << std::setprecision(8) << "** L2 error (abs/rel) for "
                    << std::setw(6) << numCellCenterDofs << " cc dofs and " << numFaceDofs << " face dofs (total: " << numCellCenterDofs + numFaceDofs << "): "
                    << std::scientific
                    << "L2(p) = " << l2error.first[Indices::pressureIdx] << " / " << l2error.second[Indices::pressureIdx]
                    << " , L2(vx) = " << l2error.first[Indices::velocityXIdx] << " / " << l2error.second[Indices::velocityXIdx]
                    << std::endl;
        }
    }

   /*!
     * \brief Returns the temperature within the domain in [K].
     *
     * This problem assumes a temperature of 10 degrees Celsius.
     */
    Scalar temperature() const
    { return 298.0; }

   /*!
     * \brief Returns the sources within the domain.
     *
     * \param globalPos The global position
     */
    NumEqVector sourceAtPos(const GlobalPosition &globalPos) const
    {
        NumEqVector source(0.0);

        // mass balance - term div(rho*v)
        for (unsigned int dimIdx = 0; dimIdx < dimWorld; ++dimIdx)
        {
            source[Indices::conti0EqIdx] += dvdx(globalPos)[dimIdx][dimIdx];
        }
        source[Indices::conti0EqIdx] *= density_;

        // momentum balance
        for (unsigned int velIdx = 0; velIdx < dimWorld; ++velIdx)
        {
            for (unsigned int dimIdx = 0; dimIdx < dimWorld; ++dimIdx)
            {
                // inertia term
                if (this->enableInertiaTerms())
                  source[Indices::velocity(velIdx)] += density_ * dv2dx(globalPos)[velIdx][dimIdx];

                // viscous term (molecular)
                source[Indices::velocity(velIdx)] -= density_ * kinematicViscosity_* dvdx2(globalPos)[velIdx][dimIdx];
                static const bool enableUnsymmetrizedVelocityGradient = getParam<bool>("FreeFlow.EnableUnsymmetrizedVelocityGradient", false);
                if (!enableUnsymmetrizedVelocityGradient)
                    source[Indices::velocity(velIdx)] -= density_ * kinematicViscosity_* dvdx2(globalPos)[dimIdx][velIdx];
            }
            // pressure term
            source[Indices::velocity(velIdx)] += dpdx(globalPos)[velIdx];

            // gravity term
            static const bool enableGravity = getParam<bool>("Problem.EnableGravity");
            if (enableGravity)
            {
                source[Indices::velocity(velIdx)] -= density_ * this->gravity()[velIdx];
            }
        }

        return source;
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
    BoundaryTypes boundaryTypesAtPos(const GlobalPosition& globalPos) const
    {
        BoundaryTypes values;

        // set Dirichlet values for the velocity everywhere
        values.setDirichlet(Indices::momentumXBalanceIdx);

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
                         const typename GridGeometry::LocalView& fvGeometry,
                         const typename GridGeometry::SubControlVolume& scv,
                         int pvIdx) const
    {
        // set a fixed pressure in all cells at the boundary
        for (const auto& scvf : scvfs(fvGeometry))
            if (scvf.boundary())
                return true;

        return false;
    }

   /*!
     * \brief Returns Dirichlet boundary values at a given position
     *
     * \param globalPos The global position
     */
    PrimaryVariables dirichletAtPos(const GlobalPosition& globalPos) const
    {
        // use the values of the analytical solution
        return analyticalSolution(globalPos);
    }

    /*!
     * \brief Returns the analytical solution of the problem at a given position
     *
     * \param globalPos The global position
     */
    PrimaryVariables analyticalSolution(const GlobalPosition& globalPos) const
    {
        PrimaryVariables values;
        values[Indices::pressureIdx] = p(globalPos);
        values[Indices::velocityXIdx] = v(globalPos);
        return values;
    }

    //! \brief The velocity
    const DimVector v(const DimVector& globalPos) const
    {
        DimVector v(0.0);
        v[0] = 2.0 * globalPos[0] * globalPos[0] * globalPos[0];
        return v;
    }

    //! \brief The velocity gradient
    const DimMatrix dvdx(const DimVector& globalPos) const
    {
        DimMatrix dvdx(0.0);
        dvdx[0][0] = 6.0 * globalPos[0] * globalPos[0];
        return dvdx;
    }

    //! \brief The gradient of the velocity squared (using product rule -> nothing to do here)
    const DimMatrix dv2dx(const DimVector& globalPos) const
    {
        DimMatrix dv2dx;
        for (unsigned int velIdx = 0; velIdx < dimWorld; ++velIdx)
        {
            for (unsigned int dimIdx = 0; dimIdx < dimWorld; ++dimIdx)
            {
                dv2dx[velIdx][dimIdx] = dvdx(globalPos)[velIdx][dimIdx] * v(globalPos)[dimIdx]
                                        + dvdx(globalPos)[dimIdx][dimIdx] * v(globalPos)[velIdx];
            }
        }
        return dv2dx;
    }

    //! \brief The gradient of the velocity gradient
    const DimMatrix dvdx2(const DimVector& globalPos) const
    {
        DimMatrix dvdx2(0.0);
        dvdx2[0][0] = 12.0 * globalPos[0];
        return dvdx2;
    }

    //! \brief The pressure
    const Scalar p(const DimVector& globalPos) const
    { return 2.0 - 2.0 * globalPos[0]; }

    //! \brief The pressure gradient
    const DimVector dpdx(const DimVector& globalPos) const
    {
        DimVector dpdx(0.0);
        dpdx[0] = -2.0;
        return dpdx;
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
    PrimaryVariables initialAtPos(const GlobalPosition& globalPos) const
    {
        return analyticalSolution(globalPos);
    }

   /*!
     * \brief Returns the analytical solution for the pressure
     */
    auto& getAnalyticalPressureSolution() const
    {
        return analyticalPressure_;
    }

   /*!
     * \brief Returns the analytical solution for the velocity.
     */
    auto& getAnalyticalVelocitySolution() const
    {
        return analyticalVelocity_;
    }

   /*!
     * \brief Returns the analytical solution for the velocity at the faces.
     */
    auto& getAnalyticalVelocitySolutionOnFace() const
    {
        return analyticalVelocityOnFace_;
    }

private:

   /*!
     * \brief Adds additional VTK output data to the VTKWriter.
     *
     * Function is called by the output module on every write.
     */
    void createAnalyticalSolution_()
    {
        analyticalPressure_.resize(this->gridGeometry().numCellCenterDofs());
        analyticalVelocity_.resize(this->gridGeometry().numCellCenterDofs());
        analyticalVelocityOnFace_.resize(this->gridGeometry().numFaceDofs());

        for (const auto& element : elements(this->gridGeometry().gridView()))
        {
            auto fvGeometry = localView(this->gridGeometry());
            fvGeometry.bindElement(element);
            for (auto&& scv : scvs(fvGeometry))
            {
                auto ccDofIdx = scv.dofIndex();
                auto ccDofPosition = scv.dofPosition();
                auto analyticalSolutionAtCc = analyticalSolution(ccDofPosition);

                // velocities on faces
                for (auto&& scvf : scvfs(fvGeometry))
                {
                    const auto faceDofIdx = scvf.dofIndex();
                    const auto faceDofPosition = scvf.center();
                    const auto dirIdx = scvf.directionIndex();
                    const auto analyticalSolutionAtFace = analyticalSolution(faceDofPosition);
                    analyticalVelocityOnFace_[faceDofIdx][dirIdx] = analyticalSolutionAtFace[Indices::velocity(dirIdx)];
                }

                analyticalPressure_[ccDofIdx] = analyticalSolutionAtCc[Indices::pressureIdx];

                for(int dirIdx = 0; dirIdx < ModelTraits::dim(); ++dirIdx)
                    analyticalVelocity_[ccDofIdx][dirIdx] = analyticalSolutionAtCc[Indices::velocity(dirIdx)];
            }
        }
     }

    bool printL2Error_;
    Scalar density_;
    Scalar kinematicViscosity_;
    std::vector<Scalar> analyticalPressure_;
    std::vector<GlobalPosition> analyticalVelocity_;
    std::vector<GlobalPosition> analyticalVelocityOnFace_;
};
} // end namespace Dumux

#endif
