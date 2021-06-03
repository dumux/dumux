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
 * \brief Test for the diamond (Navier-)Stokes model with analytical solution (Donea 2003, \cite Donea2003).
 */
#ifndef DUMUX_DONEA_TEST_PROBLEM_DIAMOND_HH_NEW
#define DUMUX_DONEA_TEST_PROBLEM_DIAMOND_HH_NEW

#ifndef ENABLECACHING
#define ENABLECACHING 0
#endif

#ifndef GRIDTYPE
#define GRIDTYPE Dune::YaspGrid<2>
#endif

#include <dune/grid/yaspgrid.hh>
#include <dune/alugrid/grid.hh>

#include <dumux/material/fluidsystems/1pliquid.hh>
#include <dumux/material/components/constant.hh>

#include <dumux/freeflow/navierstokes/momentum/diamond/model.hh>
#include <dumux/freeflow/navierstokes/problem.hh>
#include <dumux/freeflow/navierstokes/mass/1p/model.hh>
#include <dumux/discretization/fcdiamond.hh>
#include <dumux/discretization/cctpfa.hh>
#include <dumux/discretization/box.hh>

namespace Dumux
{
template <class TypeTag>
class DoneaTestProblemNew;


namespace Properties
{
// Create new type tags
namespace TTag {
struct DoneaTestNew {};
struct DoneaTestNewMomentum { using InheritsFrom = std::tuple<DoneaTestNew, NavierStokesMomentumDiamond, FaceCenteredDiamondModel>; };
struct DoneaTestNewMass { using InheritsFrom = std::tuple<DoneaTestNew, NavierStokesMassOneP, BoxModel>; };
} // end namespace TTag

// Set the problem property
template<class TypeTag>
struct Problem<TypeTag, TTag::DoneaTestNew>
{
    using type = Dumux::DoneaTestProblemNew<TypeTag> ;
};

// the fluid system
template<class TypeTag>
struct FluidSystem<TypeTag, TTag::DoneaTestNew>
{
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using type = FluidSystems::OnePLiquid<Scalar, Components::Constant<1, Scalar> >;
};

// Set the grid type
template<class TypeTag>
struct Grid<TypeTag, TTag::DoneaTestNew> { using type = GRIDTYPE; };

template<class TypeTag>
struct EnableGridGeometryCache<TypeTag, TTag::DoneaTestNew> { static constexpr bool value = true; };
template<class TypeTag>
struct EnableGridFluxVariablesCache<TypeTag, TTag::DoneaTestNew> { static constexpr bool value = true; };
template<class TypeTag>
struct EnableGridVolumeVariablesCache<TypeTag, TTag::DoneaTestNew> { static constexpr bool value = true; };

}

/*!
 * \ingroup NavierStokesTests
 * \brief  Test problem for the staggered grid (Donea 2003, \cite Donea2003).
 *
 * A two-dimensional Stokes flow in a square domain is considered.
 * With the source terms as given in Donea 2003 \cite Donea2003, an analytical solution
 * is available and can be compared to the numerical solution.
 */
template <class TypeTag>
class DoneaTestProblemNew : public NavierStokesProblem<TypeTag>
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

public:
    DoneaTestProblemNew(std::shared_ptr<const GridGeometry> gridGeometry, std::shared_ptr<CouplingManager> couplingManager)
    : ParentType(gridGeometry, couplingManager)
    {
        useNeumann_ = getParam<bool>("Problem.UseNeumann", false);
    }

    DoneaTestProblemNew(std::shared_ptr<const GridGeometry> gridGeometry)
    : ParentType(gridGeometry)
    {
        useNeumann_ = getParam<bool>("Problem.UseNeumann", false);
    }

   /*!
     * \name Problem parameters
     */
    // \{
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
        if constexpr (ParentType::isMomentumProblem())
        {
            NumEqVector source;
            const Scalar x = globalPos[0];
            const Scalar y = globalPos[1];

            source[Indices::momentumXBalanceIdx] = (12.0-24.0*y) * x*x*x*x + (-24.0 + 48.0*y)* x*x*x
                                                 + (-48.0*y + 72.0*y*y - 48.0*y*y*y + 12.0)* x*x
                                                 + (-2.0 + 24.0*y - 72.0*y*y + 48.0*y*y*y)*x
                                                 + 1.0 - 4.0*y + 12.0*y*y - 8.0*y*y*y;
            source[Indices::momentumYBalanceIdx] = (8.0 - 48.0*y + 48.0*y*y)*x*x*x + (-12.0 + 72.0*y - 72.0*y*y)*x*x
                                                 + (4.0 - 24.0*y + 48.0*y*y - 48.0*y*y*y + 24.0*y*y*y*y)*x - 12.0*y*y
                                                 + 24.0*y*y*y - 12.0*y*y*y*y;
            return source;
        }
        else
        {
            return NumEqVector(0.0);
        }
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

        // set Dirichlet values for the velocity and pressure everywhere
        if constexpr (ParentType::isMomentumProblem())
        {
            if (useNeumann_)
            {
                static constexpr Scalar eps = 1e-8;
                if ((globalPos[0] > this->gridGeometry().bBoxMax()[0] - eps) || (globalPos[1] > this->gridGeometry().bBoxMax()[1] - eps))
                    values.setAllNeumann();
                else
                    values.setAllDirichlet();
            }
            else
            {
                values.setDirichlet(Indices::velocityXIdx);
                values.setDirichlet(Indices::velocityYIdx);
            }
        }
        else
            values.setNeumann(Indices::conti0EqIdx);

        return values;
    }

   /*!
     * \brief Return dirichlet boundary values at a given position
     *
     * \param globalPos The global position
     */
    PrimaryVariables dirichletAtPos(const GlobalPosition& globalPos) const
    {
        return analyticalSolution(globalPos);
    }

    /*!
     * \brief Evaluates the boundary conditions for a Neumann control volume.
     *
     * \param element The element for which the Neumann boundary condition is set
     * \param fvGeometry The fvGeometry
     * \param elemVolVars The element volume variables
     * \param elemFaceVars The element face variables
     * \param scvf The boundary sub control volume face
     */
    template<class ElementVolumeVariables, class ElementFluxVariablesCache>
    NumEqVector neumann(const Element& element,
                        const FVElementGeometry& fvGeometry,
                        const ElementVolumeVariables& elemVolVars,
                        const ElementFluxVariablesCache& elemFluxVarsCache,
                        const SubControlVolumeFace& scvf) const
    {
        NumEqVector values(0.0);

        if constexpr (ParentType::isMomentumProblem())
        {
            const auto x = scvf.ipGlobal()[0];
            const auto y = scvf.ipGlobal()[1];

            Dune::FieldMatrix<Scalar, dimWorld, dimWorld> momentumFlux(0.0);
            momentumFlux[0][0] = x*(-2*x*y*(2*x - 2.0)*(4.0*y*y - 6.0*y + 2.0) - x - 4*y*(x - 1.0)*(x - 1.0)*(4.0*y*y - 6.0*y + 2.0) + 1.0);
            momentumFlux[0][1] = x*x*(x - 1.0)*(x - 1.0)*(-12.0*y*y + 12.0*y - 2.0) + y*y*(y - 1.0)*(y - 1.0)*(12.0*x*x - 12.0*x + 2.0);
            momentumFlux[1][0] = momentumFlux[0][1];
            momentumFlux[1][1] = x*(-x + 2.0*y*y*(2*y - 2.0)*(4.0*x*x - 6.0*x + 2.0) + 4.0*y*(y - 1.0)*(y - 1.0)*(4.0*x*x - 6.0*x + 2.0) + 1.0);

            const auto normal = scvf.unitOuterNormal();
            momentumFlux.mv(normal, values);
        }
        else
        {
            const auto insideDensity = elemVolVars[scvf.insideScvIdx()].density();
            values[Indices::conti0EqIdx] = this->faceVelocity(element, fvGeometry, scvf) * insideDensity * scvf.unitOuterNormal();
        }

        return values;
    }

    /*!
     * \brief Return the analytical solution of the problem at a given position
     *
     * \param globalPos The global position
     */
    PrimaryVariables analyticalSolution(const GlobalPosition& globalPos) const
    {
        PrimaryVariables values;
        const Scalar x = globalPos[0];
        [[maybe_unused]] const Scalar y = globalPos[1];

        if constexpr (ParentType::isMomentumProblem())
        {
            values[Indices::velocityXIdx] = x*x * (1.0 - x)*(1.0 - x) * (2.0*y - 6.0*y*y + 4.0*y*y*y);
            values[Indices::velocityYIdx] = -1.0*y*y * (1.0 - y)*(1.0 - y) * (2.0*x - 6.0*x*x + 4.0*x*x*x);
        }
        else
            values[Indices::pressureIdx] = x * (1.0-x);

        return values;
    }

    // \}

   /*!
     * \name Volume terms
     */
    // \{

    Scalar pressureAtPos(const GlobalPosition& globalPos) const
    {
        const Scalar x = globalPos[0];
        return x * (1.0-x);
    }

    Scalar densityAtPos(const GlobalPosition& globalPos) const
    { return 1; }

    Scalar effectiveViscosityAtPos(const GlobalPosition& globalPos) const
    { return 1; }

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

        if (!useNeumann_)
        {
            if (onBoundary_(scv.dofPosition()))
                values.set(0);
        }

        return values;
    }

    /*!
     * \brief Define the values of internal Dirichlet constraints for a degree of freedom.
     * \param element The finite element
     * \param scv The sub-control volume
     */
    // PrimaryVariables internalDirichlet(const Element& element, const SubControlVolume& scv) const
    // { return PrimaryVariables(analyticalSolution(scv.dofPosition())[Indices::pressureIdx]); }

    template<class SolutionVector>
    auto calculateL2Error(const SolutionVector& curSol)
    {
        Scalar l2Error = 0.0;

        for (const auto& element : elements(this->gridGeometry().gridView()))
        {
            auto fvGeometry = localView(this->gridGeometry());
            fvGeometry.bindElement(element);

            for (auto&& scv : scvs(fvGeometry))
            {
                PrimaryVariables v(analyticalSolution(scv.dofPosition()));
                v -= curSol[scv.dofIndex()];
                l2Error += scv.volume()*v.two_norm2();
            }
        }

        return std::sqrt(l2Error);
    }

    /*!
     * \brief Define the values of internal Dirichlet constraints for a degree of freedom.
     * \param element The finite element
     * \param scv The sub-control volume
     */
    PrimaryVariables internalDirichlet(const Element& element, const SubControlVolume& scv) const
    { return PrimaryVariables(analyticalSolution(scv.dofPosition())[Indices::pressureIdx]); }

private:
    bool onLeftBoundary_(const GlobalPosition &globalPos) const
    {
        return globalPos[0] < this->gridGeometry().bBoxMin()[0] + eps_;
    }

    bool onRightBoundary_(const GlobalPosition &globalPos) const
    {
        return globalPos[0] > this->gridGeometry().bBoxMax()[0] - eps_;
    }

    bool onLowerBoundary_(const GlobalPosition &globalPos) const
    {
        return globalPos[1] < this->gridGeometry().bBoxMin()[1] + eps_;
    }

    bool onUpperBoundary_(const GlobalPosition &globalPos) const
    {
        return globalPos[1] > this->gridGeometry().bBoxMax()[1] - eps_;
    }


    bool onBoundary_(const GlobalPosition& globalPos) const
    {
        return onLeftBoundary_(globalPos) || onRightBoundary_(globalPos)
               || onLowerBoundary_(globalPos) || onUpperBoundary_(globalPos);
    }

    bool useNeumann_;
    static constexpr Scalar eps_ = 1.5e-7;
};
} // end namespace Dumux

#endif