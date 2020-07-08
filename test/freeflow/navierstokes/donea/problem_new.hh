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
 * \brief Test for the staggered grid (Navier-)Stokes model with analytical solution (Donea 2003, \cite Donea2003).
 */
#ifndef DUMUX_DONEA_TEST_PROBLEM_HH_NEW
#define DUMUX_DONEA_TEST_PROBLEM_HH_NEW

#ifndef ENABLECACHING
#define ENABLECACHING 0
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


namespace Dumux
{
template <class TypeTag>
class DoneaTestProblemNew;


namespace Properties
{
// Create new type tags
namespace TTag {
struct DoneaTestNew {};
struct DoneaTestNewMomentum { using InheritsFrom = std::tuple<DoneaTestNew, NavierStokesMomentum, FaceCenteredStaggeredModel>; };
struct DoneaTestNewMass { using InheritsFrom = std::tuple<DoneaTestNew, NavierStokesMassAndEnergy, CCTpfaModel>; };
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
struct Grid<TypeTag, TTag::DoneaTestNew> { using type = Dune::YaspGrid<2>; };

template<class TypeTag>
struct EnableGridGeometryCache<TypeTag, TTag::DoneaTestNew> { static constexpr bool value = true; };
template<class TypeTag>
struct EnableGridFluxVariablesCache<TypeTag, TTag::DoneaTestNew> { static constexpr bool value = true; };
template<class TypeTag>
struct EnableGridVolumeVariablesCache<TypeTag, TTag::DoneaTestNew> { static constexpr bool value = true; };

}

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
 * \brief  Test problem for the staggered grid (Donea 2003, \cite Donea2003).
 *
 * A two-dimensional Stokes flow in a square domain is considered.
 * With the source terms as given in Donea 2003 \cite Donea2003, an analytical solution
 * is available and can be compared to the numerical solution.
 */
template <class TypeTag>
class DoneaTestProblemNew : public Impl::BaseProblem<TypeTag>
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

public:
    DoneaTestProblemNew(std::shared_ptr<const GridGeometry> gridGeometry, std::shared_ptr<CouplingManager> couplingManager)
    : ParentType(gridGeometry, couplingManager)
    , couplingManager_(couplingManager)
    {
        printL2Error_ = getParam<bool>("Problem.PrintL2Error");
    }

    DoneaTestProblemNew(std::shared_ptr<const GridGeometry> gridGeometry)
    : ParentType(gridGeometry)
    {
        printL2Error_ = getParam<bool>("Problem.PrintL2Error");
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
    //         const int numCellCenterDofs = this->gridGeometry().gridView().size(0);
    //         const int numFaceDofs = this->gridGeometry().numDofs();
    //         std::cout << std::setprecision(8) << "** L2 error (abs/rel) for "
    //                 << std::setw(6) << numCellCenterDofs << " cc dofs and " << numFaceDofs << " face dofs (total: " << numCellCenterDofs + numFaceDofs << "): "
    //                 << std::scientific
    //                 << "L2(p) = " << l2error.first[/*Indices::pressureIdx*/0] << " / " << l2error.second[/*Indices::pressureIdx*/0]
    //                 << " , L2(vx) = " << l2error.first[Indices::velocityXIdx] << " / " << l2error.second[Indices::velocityXIdx]
    //                 << " , L2(vy) = " << l2error.first[Indices::velocityYIdx] << " / " << l2error.second[Indices::velocityYIdx]
    //                 << std::endl;
    //     }
    // }

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
        if constexpr (Impl::isMomentumProblem<TypeTag>())
        {
            NumEqVector source;
            Scalar x = globalPos[0];
            Scalar y = globalPos[1];

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
        if constexpr (Impl::isMomentumProblem<TypeTag>())
        {
            values.setDirichlet(Indices::velocityXIdx);
            values.setDirichlet(Indices::velocityYIdx);
        }
        else
        {
            values.setDirichlet(0/* Indices::pressureIdx */);
        }

        return values;
    }

   /*!
     * \brief Return dirichlet boundary values at a given position
     *
     * \param globalPos The global position
     */
    PrimaryVariables dirichletAtPos(const GlobalPosition& globalPos) const
    {
        const auto sol = analyticalSolution(globalPos);
        PrimaryVariables values;

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
     * \brief Return the analytical solution of the problem at a given position
     *
     * \param globalPos The global position
     */
    auto analyticalSolution(const GlobalPosition& globalPos) const
    {
        Scalar x = globalPos[0];
        Scalar y = globalPos[1];
        std::array<Scalar, 3> values;

        values[0] = x*x * (1.0 - x)*(1.0 - x) * (2.0*y - 6.0*y*y + 4.0*y*y*y);        // vx
        values[1] = -1.0*y*y * (1.0 - y)*(1.0 - y) * (2.0*x - 6.0*x*x + 4.0*x*x*x);   // vy
        values[2] = x * (1.0-x);                                                      //p
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
    PrimaryVariables initialAtPos(const GlobalPosition& globalPos) const
    {
        return PrimaryVariables(0.0);
    }

    Scalar pressureAtPos(const GlobalPosition& globalPos) const
    {
        return analyticalSolution(globalPos)[2];
    }

    Scalar densityAtPos(const GlobalPosition& globalPos) const
    {
        return 1;
    }

    Scalar effectiveViscosityAtPos(const GlobalPosition& globalPos) const
    { return 1; }


   /*!
     * \brief Returns the analytical solution for the pressure
     */
    template<bool enable = !Impl::isMomentumProblem<TypeTag>(), std::enable_if_t<enable, int> = 0>
    const std::vector<Scalar> getAnalyticalPressureSolution() const
    {
        std::vector<Scalar> analyticalPressure(this->gridGeometry().gridView().size(0));
        for (const auto& element : elements(this->gridGeometry().gridView()))
        {
            auto fvGeometry = localView(this->gridGeometry());
            fvGeometry.bindElement(element);
            for (auto&& scv : scvs(fvGeometry))
            {
                const auto ccDofIdx = scv.dofIndex();
                const auto ccDofPosition = scv.dofPosition();
                const auto analyticalSolutionAtCc = analyticalSolution(ccDofPosition);
                analyticalPressure[ccDofIdx] = analyticalSolutionAtCc[/*Indices::pressureIdx*/2];
            }
        }

        return analyticalPressure;
    }

   /*!
     * \brief Returns the analytical solution for the velocity
     */
    template<bool enable = Impl::isMomentumProblem<TypeTag>(), std::enable_if_t<enable, int> = 0>
    const std::vector<VelocityVector> getAnalyticalVelocitySolution() const
    {
        std::vector<VelocityVector> analyticalVelocity(this->gridGeometry().gridView().size(0));
        for (const auto& element : elements(this->gridGeometry().gridView()))
        {
            const auto& pos = element.geometry().center();
            const auto eIdx = this->gridGeometry().elementMapper().index(element);
            for(int dirIdx = 0; dirIdx < ModelTraits::dim(); ++dirIdx)
                analyticalVelocity[eIdx][dirIdx] = analyticalSolution(pos)[Indices::velocity(dirIdx)];
        }

        return analyticalVelocity;
    }

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
    { return PrimaryVariables(analyticalSolution(scv.center())[2]); }


private:

    bool printL2Error_;
    std::shared_ptr<CouplingManager> couplingManager_;
};
} // end namespace Dumux

#endif
