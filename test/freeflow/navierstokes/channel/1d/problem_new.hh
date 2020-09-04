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

#include <dumux/freeflow/navierstokes/momentum/model.hh>
#include <dumux/freeflow/navierstokes/problem.hh>
#include <dumux/freeflow/navierstokes/mass/1p/model.hh>
#include <dumux/discretization/fcstaggered.hh>
#include <dumux/discretization/cctpfa.hh>

#include <dumux/material/components/constant.hh>
#include <dumux/material/fluidsystems/1pliquid.hh>

#include "../../l2error.hh"

namespace Dumux {
template <class TypeTag>
class NavierStokesAnalyticProblem;

namespace Properties {
// Create new type tags
namespace TTag {
struct NavierStokesAnalytic {};
struct NavierStokesAnalyticMomentum { using InheritsFrom = std::tuple<NavierStokesAnalytic, NavierStokesMomentum, FaceCenteredStaggeredModel>; };
struct NavierStokesAnalyticMass { using InheritsFrom = std::tuple<NavierStokesAnalytic, NavierStokesMassOneP, CCTpfaModel>; };
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

    using BoundaryTypes = typename ParentType::BoundaryTypes;
    using GridGeometry = GetPropType<TypeTag, Properties::GridGeometry>;
    using FVElementGeometry = typename GridGeometry::LocalView;
    using SubControlVolume = typename FVElementGeometry::SubControlVolume;
    using SubControlVolumeFace = typename FVElementGeometry::SubControlVolumeFace;
    using Indices = typename GetPropType<TypeTag, Properties::ModelTraits>::Indices;
    using NumEqVector = typename ParentType::NumEqVector;
    using ModelTraits = GetPropType<TypeTag, Properties::ModelTraits>;
    using PrimaryVariables = typename ParentType::PrimaryVariables;
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using SolutionVector = GetPropType<TypeTag, Properties::SolutionVector>;

    static constexpr auto dimWorld = GridGeometry::GridView::dimensionworld;
    using Element = typename GridGeometry::GridView::template Codim<0>::Entity;
    using GlobalPosition = typename Element::Geometry::GlobalCoordinate;
    using VelocityVector = Dune::FieldVector<Scalar, dimWorld>;

    using CouplingManager = GetPropType<TypeTag, Properties::CouplingManager>;
    using DimVector = GlobalPosition;
    using DimMatrix = Dune::FieldMatrix<Scalar, dimWorld, dimWorld>;

public:
    NavierStokesAnalyticProblem(std::shared_ptr<const GridGeometry> gridGeometry, std::shared_ptr<CouplingManager> couplingManager)
    : ParentType(gridGeometry, couplingManager)
    {
        printL2Error_ = getParam<bool>("Problem.PrintL2Error");
        density_ = getParam<Scalar>("Component.LiquidDensity");
        kinematicViscosity_ = getParam<Scalar>("Component.LiquidKinematicViscosity");
    }

   /*!
     * \name Problem parameters
     */
    // \{


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

        if constexpr (!ParentType::isMomentumProblem())
        {
            // mass balance - term div(rho*v)
            for (unsigned int dimIdx = 0; dimIdx < dimWorld; ++dimIdx)
            {
                source[Indices::conti0EqIdx] += dvdx(globalPos)[dimIdx][dimIdx];
            }
            source[Indices::conti0EqIdx] *= density_;
        }
        else
        {
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
        values.setAllDirichlet();
        return values;
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
     * \brief Return the analytical solution of the problem at a given position
     *
     * \param globalPos The global position
     */
    PrimaryVariables analyticalSolution(const GlobalPosition& globalPos) const
    {
        PrimaryVariables values;

        if constexpr (ParentType::isMomentumProblem())
            values[Indices::velocityXIdx] = v(globalPos);
        else
            values[Indices::pressureIdx] = p(globalPos);

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

    //! Enable internal Dirichlet constraints
    static constexpr bool enableInternalDirichletConstraints()
    { return !ParentType::isMomentumProblem(); }

    /*!
     * \brief Tag a degree of freedom to carry internal Dirichlet constraints.
     *        If true is returned for a dof, the equation for this dof is replaced
     *        by the constraint that its primary variable values must match the
     *        user-defined values obtained from the function internalDirichlet(),
     *        which must be defined in the problem.
     *
     * \param element The finite element
     * \param scv The sub-control volume
     */
    std::bitset<PrimaryVariables::dimension> hasInternalDirichletConstraint(const Element& element, const SubControlVolume& scv) const
    {
        std::bitset<PrimaryVariables::dimension> values;

        auto fvGeometry = localView(this->gridGeometry());
        fvGeometry.bindElement(element);

        bool onBoundary = false;
        for (const auto& scvf : scvfs(fvGeometry))
            onBoundary = std::max(onBoundary, scvf.boundary());

        if (onBoundary)
            values.set(0);

        // TODO: only use one cell or pass fvGeometry to hasInternalDirichletConstraint

        // if (scv.dofIndex() == 0)
        //     values.set(0);
        // the pure Neumann problem is only defined up to a constant
        // we create a well-posed problem by fixing the pressure at one dof
        return values;
    }

    /*!
     * \brief Define the values of internal Dirichlet constraints for a degree of freedom.
     * \param element The finite element
     * \param scv The sub-control volume
     */
    PrimaryVariables internalDirichlet(const Element& element, const SubControlVolume& scv) const
    { return PrimaryVariables(analyticalSolution(scv.center())[Indices::pressureIdx]); }

    // \}

   /*!
     * \name Volume terms
     */
    // \{


private:

    bool printL2Error_;
    Scalar density_;
    Scalar kinematicViscosity_;

};
} // end namespace Dumux

#endif
