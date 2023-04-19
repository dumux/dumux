// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
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

#include <dumux/common/properties.hh>
#include <dumux/common/parameters.hh>
#include <dumux/freeflow/navierstokes/boundarytypes.hh>

namespace Dumux {

/*!
 * \ingroup NavierStokesTests
 * \brief Test for the 1-D Navier-Stokes model with an analytical solution.
 */
template <class TypeTag, class BaseProblem>
class NavierStokesAnalyticProblem : public BaseProblem
{
    using ParentType = BaseProblem;

    using BoundaryTypes = Dumux::NavierStokesBoundaryTypes<GetPropType<TypeTag, Properties::ModelTraits>::numEq()>;
    using GridGeometry = GetPropType<TypeTag, Properties::GridGeometry>;
    using InitialValues = typename ParentType::InitialValues;
    using Sources = typename ParentType::Sources;
    using DirichletValues = typename ParentType::DirichletValues;
    using BoundaryFluxes = typename ParentType::BoundaryFluxes;
    using ModelTraits = GetPropType<TypeTag, Properties::ModelTraits>;
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using SolutionVector = GetPropType<TypeTag, Properties::SolutionVector>;
    using SubControlVolume = typename GridGeometry::SubControlVolume;
    using SubControlVolumeFace = typename GridGeometry::SubControlVolumeFace;
    using FVElementGeometry = typename GridGeometry::LocalView;

    static constexpr auto dimWorld = GridGeometry::GridView::dimensionworld;
    using Element = typename GridGeometry::GridView::template Codim<0>::Entity;
    using GlobalPosition = typename Element::Geometry::GlobalCoordinate;
    using CouplingManager = GetPropType<TypeTag, Properties::CouplingManager>;
    using DimVector = GlobalPosition;
    using DimMatrix = Dune::FieldMatrix<Scalar, dimWorld, dimWorld>;

public:
    using Indices = typename GetPropType<TypeTag, Properties::ModelTraits>::Indices;

    NavierStokesAnalyticProblem(std::shared_ptr<const GridGeometry> gridGeometry, std::shared_ptr<CouplingManager> couplingManager)
    : ParentType(gridGeometry, couplingManager)
    {
        density_ = getParam<Scalar>("Component.LiquidDensity");
        kinematicViscosity_ = getParam<Scalar>("Component.LiquidKinematicViscosity");
    }

    /*!
     * \brief Returns the sources within the domain.
     *
     * \param globalPos The global position
     */
    Sources sourceAtPos(const GlobalPosition &globalPos) const
    {
        Sources source(0.0);

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

        if constexpr (ParentType::isMomentumProblem())
        {
            // set Dirichlet values for the velocity everywhere
            values.setDirichlet(Indices::momentumXBalanceIdx);
        }
        else
            values.setNeumann(Indices::pressureIdx);

        return values;
    }

    /*!
     * \brief Returns Dirichlet boundary values at a given position
     *
     * \param globalPos The global position
     */
    DirichletValues dirichletAtPos(const GlobalPosition& globalPos) const
    {
        // use the values of the analytical solution
        return analyticalSolution(globalPos);
    }

    /*!
     * \brief Returns the analytical solution of the problem at a given position
     *
     * \param globalPos The global position
     * \param time A parameter for consistent signatures. It is ignored here as this is a stationary test.
     */
    DirichletValues analyticalSolution(const GlobalPosition& globalPos, Scalar time = 0.0) const
    {
        DirichletValues values;

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
    std::bitset<DirichletValues::dimension> hasInternalDirichletConstraint(const Element& element, const SubControlVolume& scv) const
    {
        std::bitset<DirichletValues::dimension> values;

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
    DirichletValues internalDirichlet(const Element& element, const SubControlVolume& scv) const
    { return DirichletValues(analyticalSolution(scv.center())[Indices::pressureIdx]); }

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
    InitialValues initialAtPos(const GlobalPosition& globalPos) const
    {
        return analyticalSolution(globalPos);
    }

private:
    Scalar density_;
    Scalar kinematicViscosity_;
};
} // end namespace Dumux

#endif
