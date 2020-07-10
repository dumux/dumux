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
 * \brief Test for the staggered grid Navier-Stokes model with analytical solution (Kovasznay 1948, \cite Kovasznay1948)
 */

#ifndef DUMUX_KOVASZNAY_TEST_PROBLEM_NEW_HH
#define DUMUX_KOVASZNAY_TEST_PROBLEM_NEW_HH

#ifndef UPWINDSCHEMEORDER
#define UPWINDSCHEMEORDER 0
#endif

#include <dune/grid/yaspgrid.hh>

#include <dumux/material/fluidsystems/1pliquid.hh>
#include <dumux/material/components/constant.hh>

#include <dumux/common/boundarytypes.hh>
#include <dumux/freeflow/navierstokes/momentum/model.hh>
#include <dumux/freeflow/navierstokes/momentum/problem.hh>
#include <dumux/freeflow/navierstokes/massandenergy/model.hh>
#include <dumux/freeflow/navierstokes/massandenergy/problem.hh>
#include <dumux/discretization/fcstaggered.hh>
#include <dumux/discretization/cctpfa.hh>
#include "../l2error.hh"

#include <dumux/multidomain/traits.hh>
#include <dumux/multidomain/staggeredfreeflow/couplingmanager.hh>

// #include "../l2error.hh"

namespace Dumux {

template <class TypeTag>
class KovasznayTestProblem;

namespace Properties {
// Create new type tags
namespace TTag {
struct KovasznayTest {};
struct KovasznayTestMomentum { using InheritsFrom = std::tuple<KovasznayTest, NavierStokesMomentum, FaceCenteredStaggeredModel>; };
struct KovasznayTestMass { using InheritsFrom = std::tuple<KovasznayTest, NavierStokesMassAndEnergy, CCTpfaModel>; };
} // end namespace TTag


// the fluid system
template<class TypeTag>
struct FluidSystem<TypeTag, TTag::KovasznayTest>
{
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using type = FluidSystems::OnePLiquid<Scalar, Components::Constant<1, Scalar> >;
};

// Set the grid type
template<class TypeTag>
struct Grid<TypeTag, TTag::KovasznayTest> { using type = Dune::YaspGrid<2, Dune::EquidistantOffsetCoordinates<GetPropType<TypeTag, Properties::Scalar>, 2> >; };

// Set the problem property
template<class TypeTag>
struct Problem<TypeTag, TTag::KovasznayTest> { using type = Dumux::KovasznayTestProblem<TypeTag> ; };

template<class TypeTag>
struct EnableGridGeometryCache<TypeTag, TTag::KovasznayTest> { static constexpr bool value = true; };

template<class TypeTag>
struct EnableGridFluxVariablesCache<TypeTag, TTag::KovasznayTest> { static constexpr bool value = true; };
template<class TypeTag>
struct EnableGridVolumeVariablesCache<TypeTag, TTag::KovasznayTest> { static constexpr bool value = true; };

template<class TypeTag>
struct UpwindSchemeOrder<TypeTag, TTag::KovasznayTest> { static constexpr int value = UPWINDSCHEMEORDER; };
} // end namespace Properties

namespace Impl {

template<class TypeTag>
constexpr bool isMomentumProblem()
{
    return GetPropType<TypeTag, Properties::GridGeometry>::discMethod == DiscretizationMethod::fcstaggered;
};

template<class TypeTag>
using BaseProblem = std::conditional_t<isMomentumProblem<TypeTag>(),
                                       NavierStokesMomentumProblem<TypeTag>,
                                       NavierStokesMassAndEnergyProblem<TypeTag>>;

template<class TypeTag>
using BoundaryTypes = std::conditional_t<isMomentumProblem<TypeTag>(),
                                         Dumux::BoundaryTypes<GetPropType<TypeTag, Properties::ModelTraits>::dim()>,
                                         Dumux::BoundaryTypes<GetPropType<TypeTag, Properties::ModelTraits>::numEq()>>;
template<class TypeTag>
using NumEqVector = std::conditional_t<isMomentumProblem<TypeTag>(),
                                       Dune::FieldVector<GetPropType<TypeTag, Properties::Scalar>, GetPropType<TypeTag, Properties::ModelTraits>::dim()>,
                                       GetPropType<TypeTag, Properties::NumEqVector>>;

template<class TypeTag>
using PrimaryVariables = NumEqVector<TypeTag>;

}

/*!
 * \ingroup NavierStokesTests
 * \brief  Test problem for the staggered grid (Kovasznay 1948, \cite Kovasznay1948)
 *
 * A two-dimensional Navier-Stokes flow with a periodicity in one direction
 * is considered. The set-up represents a wake behind a two-dimensional grid
 * and is chosen in a way such that an exact solution is available.
 */
template <class TypeTag>
class KovasznayTestProblem : public Impl::BaseProblem<TypeTag>
{
    using ParentType = Impl::BaseProblem<TypeTag>;

    using BoundaryTypes = Impl::BoundaryTypes<TypeTag>;
    using GridGeometry = GetPropType<TypeTag, Properties::GridGeometry>;
    using FVElementGeometry = typename GridGeometry::LocalView;
    using SubControlVolume = typename FVElementGeometry::SubControlVolume;
    using SubControlVolumeFace = typename FVElementGeometry::SubControlVolumeFace;
    using Indices = typename GetPropType<TypeTag, Properties::ModelTraits>::Indices;
    using NumEqVector = Impl::NumEqVector<TypeTag>;
    using ModelTraits = GetPropType<TypeTag, Properties::ModelTraits>;
    using PrimaryVariables = Impl::PrimaryVariables<TypeTag>;
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using SolutionVector = GetPropType<TypeTag, Properties::SolutionVector>;

    static constexpr auto dimWorld = GridGeometry::GridView::dimensionworld;
    using Element = typename GridGeometry::GridView::template Codim<0>::Entity;
    using GlobalPosition = typename Element::Geometry::GlobalCoordinate;
    using VelocityVector = Dune::FieldVector<Scalar, dimWorld>;

    using CouplingManager = GetPropType<TypeTag, Properties::CouplingManager>;

    static constexpr auto upwindSchemeOrder = getPropValue<TypeTag, Properties::UpwindSchemeOrder>();

public:
    KovasznayTestProblem(std::shared_ptr<const GridGeometry> gridGeometry, std::shared_ptr<CouplingManager> couplingManager)
    : ParentType(gridGeometry, couplingManager)
    , couplingManager_(couplingManager)
    {
        // printL2Error_ = getParam<bool>("Problem.PrintL2Error");
        // std::cout<< "upwindSchemeOrder is: " << GridGeometry::upwindStencilOrder() << "\n";
        kinematicViscosity_ = getParam<Scalar>("Component.LiquidKinematicViscosity", 1.0);
        Scalar reynoldsNumber = 1.0 / kinematicViscosity_;
        lambda_ = 0.5 * reynoldsNumber
                        - std::sqrt(reynoldsNumber * reynoldsNumber * 0.25 + 4.0 * M_PI * M_PI);

        // createAnalyticalSolution_();
    }

   /*!
     * \name Problem parameters
     */
    // \{

    // void printL2Error(const SolutionVector& curSol) const
    // {
    //     if(printL2Error_)
    //     {
    //         using L2Error = NavierStokesTestL2Error<Scalar, ModelTraits, PrimaryVariables>;
    //         const auto l2error = L2Error::calculateL2Error(*this, curSol);
    //         const int numCellCenterDofs = this->gridGeometry().numCellCenterDofs();
    //         const int numFaceDofs = this->gridGeometry().numFaceDofs();
    //         std::cout << std::setprecision(8) << "** L2 error (abs/rel) for "
    //                 << std::setw(6) << numCellCenterDofs << " cc dofs and " << numFaceDofs << " face dofs (total: " << numCellCenterDofs + numFaceDofs << "): "
    //                 << std::scientific
    //                 << "L2(p) = " << l2error.first[Indices::pressureIdx] << " / " << l2error.second[Indices::pressureIdx]
    //                 << " , L2(vx) = " << l2error.first[Indices::velocityXIdx] << " / " << l2error.second[Indices::velocityXIdx]
    //                 << " , L2(vy) = " << l2error.first[Indices::velocityYIdx] << " / " << l2error.second[Indices::velocityYIdx]
    //                 << std::endl;
    //     }
    // }

   /*!
     * \brief Returns the temperature within the domain in [K].
     *
     * This problem assumes a temperature of 10 degrees Celsius.
     */
    Scalar temperature() const
    { return 298.0; }

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

        values.setAllDirichlet(); // does not really make sense for the mass balance

        // set Dirichlet values for the velocity everywhere
        // values.setDirichlet(Indices::velocityXIdx);
        // values.setDirichlet(Indices::velocityYIdx);

        return values;
    }

   /*!
     * \brief Returns Dirichlet boundary values at a given position.
     *
     * \param globalPos The global position
     */
    PrimaryVariables dirichletAtPos(const GlobalPosition & globalPos) const
    {
        PrimaryVariables values(0.0);
        const auto sol = analyticalSolution(globalPos);
        // use the values of the analytical solution
        if constexpr (Impl::isMomentumProblem<TypeTag>())
        {
            values[Indices::velocityXIdx] = sol[0];
            values[Indices::velocityYIdx] = sol[1];
        }
        else
            values[0] = sol[2];

        return values;
    }

   /*!
     * \brief Returns the analytical solution of the problem at a given position.
     *
     * \param globalPos The global position
     */
    std::array<Scalar, 3> analyticalSolution(const GlobalPosition& globalPos) const
    {
        Scalar x = globalPos[0];
        Scalar y = globalPos[1];

        std::array<Scalar, 3> values;
        values[2] = 0.5 * (1.0 - std::exp(2.0 * lambda_ * x));
        values[0] = 1.0 - std::exp(lambda_ * x) * std::cos(2.0 * M_PI * y);
        values[1] = 0.5 * lambda_ / M_PI * std::exp(lambda_ * x) * std::sin(2.0 * M_PI * y);

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
        return PrimaryVariables(0.0);
    }

//    /*!
//      * \brief Returns the analytical solution for the pressure.
//      */
//     auto& getAnalyticalPressureSolution() const
//     {
//         return analyticalPressure_;
//     }

//    /*!
//      * \brief Returns the analytical solution for the velocity
//      */
//     auto& getAnalyticalVelocitySolution() const
//     {
//         return analyticalVelocity_;
//     }

//    /*!
//      * \brief Returns the analytical solution for the velocity at the faces.
//      */
//     auto& getAnalyticalVelocitySolutionOnFace() const
//     {
//         return analyticalVelocityOnFace_;
//     }

        //! Enable internal Dirichlet constraints
    static constexpr bool enableInternalDirichletConstraints()
    { return !Impl::isMomentumProblem<TypeTag>(); }

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

        auto isAtLeftBoundary = [&](const FVElementGeometry& fvGeometry)
        {
            if (fvGeometry.hasBoundaryScvf())
            {
                for (const auto& scvf : scvfs(fvGeometry))
                    if (scvf.boundary() && scvf.center()[0] < this->gridGeometry().bBoxMin()[0] + eps_)
                        return true;
            }
            return false;
        };

        if (isAtLeftBoundary(fvGeometry))
            values.set(0);

        // TODO: only use one cell or pass fvGeometry to hasInternalDirichletConstraint
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
    { return PrimaryVariables(analyticalSolution(scv.center())[2]); }

private:


    static constexpr Scalar eps_=1e-6;

    Scalar kinematicViscosity_;
    Scalar lambda_;
    std::shared_ptr<CouplingManager> couplingManager_;
};
} // end namespace Dumux

#endif
