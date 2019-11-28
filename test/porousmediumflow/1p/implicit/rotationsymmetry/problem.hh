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
 * \ingroup OnePTests
 * \brief The properties and problem setup for rotation symmetry test
 */

#ifndef DUMUX_ONEP_ROTATIONSYM_TEST_PROBLEM_HH
#define DUMUX_ONEP_ROTATIONSYM_TEST_PROBLEM_HH

#include <dune/grid/yaspgrid.hh>

#include <dumux/discretization/cctpfa.hh>
#include <dumux/discretization/box.hh>
#include <dumux/discretization/rotationsymmetricgridgeometrytraits.hh>

#include <dumux/porousmediumflow/problem.hh>
#include <dumux/porousmediumflow/1p/model.hh>
#include <dumux/porousmediumflow/1p/incompressiblelocalresidual.hh>

#include <dumux/material/components/constant.hh>
#include <dumux/material/fluidsystems/1pliquid.hh>

#include "spatialparams.hh"

namespace Dumux {
// forward declarations
template<class TypeTag> class OnePTestProblem;

namespace Properties {
// Create new type tags
namespace TTag {
struct OnePRotSym { using InheritsFrom = std::tuple<OneP>; };
struct OnePRotSymTpfa { using InheritsFrom = std::tuple<OnePRotSym, CCTpfaModel>; };
struct OnePRotSymBox { using InheritsFrom = std::tuple<OnePRotSym, BoxModel>; };
} // end namespace TTag

// Set the grid type
template<class TypeTag>
struct Grid<TypeTag, TTag::OnePRotSym> { using type = Dune::YaspGrid<1, Dune::EquidistantOffsetCoordinates<GetPropType<TypeTag, Properties::Scalar>, 1> >; };

template<class TypeTag>
struct GridGeometry<TypeTag, TTag::OnePRotSymTpfa>
{
private:
    static constexpr bool enableCache = getPropValue<TypeTag, Properties::EnableGridGeometryCache>();
    using GridView = GetPropType<TypeTag, Properties::GridView>;
    using GGTraits = RotationSymmetricGridGeometryTraits<CCTpfaDefaultGridGeometryTraits<GridView>, RotationPolicy::disc>;
public:
    using type = CCTpfaFVGridGeometry<GridView, enableCache, GGTraits>;
};

template<class TypeTag>
struct GridGeometry<TypeTag, TTag::OnePRotSymBox>
{
private:
    static constexpr bool enableCache = getPropValue<TypeTag, Properties::EnableGridGeometryCache>();
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using GridView = GetPropType<TypeTag, Properties::GridView>;
    using GGTraits = RotationSymmetricGridGeometryTraits<BoxDefaultGridGeometryTraits<GridView>, RotationPolicy::disc>;
public:
    using type = BoxFVGridGeometry<Scalar, GridView, enableCache, GGTraits>;
};

// Set the problem type
template<class TypeTag>
struct Problem<TypeTag, TTag::OnePRotSym> { using type = OnePTestProblem<TypeTag>; };

// set the spatial params
template<class TypeTag>
struct SpatialParams<TypeTag, TTag::OnePRotSym>
{
    using GridGeometry = GetPropType<TypeTag, Properties::GridGeometry>;
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using type = OnePTestSpatialParams<GridGeometry, Scalar>;
};

// use the incompressible local residual (provides analytic jacobian)
template<class TypeTag>
struct LocalResidual<TypeTag, TTag::OnePRotSym> { using type = OnePIncompressibleLocalResidual<TypeTag>; };

// the fluid system
template<class TypeTag>
struct FluidSystem<TypeTag, TTag::OnePRotSym>
{
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using type = FluidSystems::OnePLiquid<Scalar, Components::Constant<1, Scalar> >;
};

// Enable caching
template<class TypeTag>
struct EnableGridVolumeVariablesCache<TypeTag, TTag::OnePRotSym> { static constexpr bool value = false; };
template<class TypeTag>
struct EnableGridFluxVariablesCache<TypeTag, TTag::OnePRotSym> { static constexpr bool value = false; };
template<class TypeTag>
struct EnableGridGeometryCache<TypeTag, TTag::OnePRotSym> { static constexpr bool value = false; };
} // end namespace Properties

/*!
 * \ingroup OnePTests
 * \brief  Test problem for the incompressible one-phase model
 */
template<class TypeTag>
class OnePTestProblem : public PorousMediumFlowProblem<TypeTag>
{
    using ParentType = PorousMediumFlowProblem<TypeTag>;
    using GridView = GetPropType<TypeTag, Properties::GridView>;
    using Element = typename GridView::template Codim<0>::Entity;
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using PrimaryVariables = GetPropType<TypeTag, Properties::PrimaryVariables>;
    using GridGeometry = GetPropType<TypeTag, Properties::GridGeometry>;
    using BoundaryTypes = GetPropType<TypeTag, Properties::BoundaryTypes>;
    static constexpr int dimWorld = GridView::dimensionworld;
    using GlobalPosition = Dune::FieldVector<Scalar, dimWorld>;
    using NumEqVector = GetPropType<TypeTag, Properties::NumEqVector>;
public:
    OnePTestProblem(std::shared_ptr<const GridGeometry> gridGeometry)
    : ParentType(gridGeometry), q_(0.0)
    {
        k_ = getParam<Scalar>("SpatialParams.Permeability");
        nu_ = getParam<Scalar>("Component.LiquidKinematicViscosity");
        q_ = getParam<Scalar>("Problem.Source");
        pW_ = getParam<Scalar>("Problem.WellPressure");
        rW_ = gridGeometry->bBoxMin()[0];

        pExact_.resize(gridGeometry->numDofs());

        for (const auto& element : elements(gridGeometry->gridView()))
        {
            auto fvGeometry = localView(*gridGeometry);
            fvGeometry.bindElement(element);
            for (auto&& scv : scvs(fvGeometry))
            {
                pExact_[scv.dofIndex()] = exactSolution(scv.dofPosition());
            }
        }
    }

    /*!
     * \brief Specifies which kind of boundary condition should be
     *        used for which equation on a given boundary control volume.
     *
     * \param globalPos The position of the center of the finite volume
     */
    BoundaryTypes boundaryTypesAtPos(const GlobalPosition &globalPos) const
    {
        BoundaryTypes values;
        values.setAllDirichlet();

        return values;
    }

    /*!
     * \brief Evaluates the boundary conditions for a Dirichlet control volume.
     *
     * \param globalPos The center of the finite volume which ought to be set.
     *
     * For this method, the \a values parameter stores primary variables.
     */
    PrimaryVariables dirichletAtPos(const GlobalPosition &globalPos) const
    {
        return exactSolution(globalPos);
    }

    /*!
     * \brief Returns the temperature \f$\mathrm{[K]}\f$ for an isothermal problem.
     *
     * This is not specific to the discretization. By default it just
     * throws an exception so it must be overloaded by the problem if
     * no energy equation is used.
     */
    Scalar temperature() const
    {
        return 283.15; // 10°C
    }

    /*!
     * \brief Adds additional VTK output data to the VTKWriter. Function is called by the output module on every write.
     */
    template<class VTKWriter>
    void addVtkFields(VTKWriter& vtk) const
    {
        vtk.addField(pExact_, "pExact");
    }

    /*!
     * \brief Calculate the L2 error between the exact and numerical solution
     * For the tpfa scheme the L2 error is calculated for the primary cells.
     * For the box scheme it is calculated on the dual mesh using a discrete L2-norm.
     * A second-order convergence rate is expected with respect to this norm
     *
     * \param curSol The current solution vector
     *
     */
    template<class SolutionVector>
    Scalar calculateL2Error(const SolutionVector& curSol)
    {
        Scalar l2error = 0.0;

        for (const auto& element :elements(this->gridGeometry().gridView()))
        {
            auto fvGeometry = localView(this->gridGeometry());
            fvGeometry.bindElement(element);
            for (auto&& scv : scvs(fvGeometry))
            {
                Scalar diff = (curSol[scv.dofIndex()] - pExact_[scv.dofIndex()]);
                l2error += scv.volume()*diff*diff;
            }
        }
        using std::sqrt;
        return sqrt(l2error);
    }

    /*!
     * \brief Writes to l2 error
     */
    template<class SolutionVector>
    void writeOutput(const SolutionVector& curSol)
    {

        Scalar l2error = calculateL2Error(curSol);

        // compute L2 error
        std::cout.precision(8);
        std::cout << "L2 error for "
                  << std::setw(6) << this->gridGeometry().numDofs()
                  << " dofs: "
                  << std::scientific
                  << l2error
                  << std::endl;
    }

private:
    /*!
     * \brief The exact solution
     * The exact solution is calculated such that the mass flux over the surface of circular disc with radius rW is q
     * and the pressure at this surface is given by pW.
     */
    PrimaryVariables exactSolution(const GlobalPosition &globalPos) const
    {
        PrimaryVariables priVars(0.0);
        const auto r = globalPos[0];
        priVars[0] = pW_ - 1.0/(2*M_PI)*nu_/k_*q_*std::log(r/rW_);
        return priVars;
    }

    Scalar q_, k_, nu_, rW_;
    GlobalPosition pW_;
    static constexpr Scalar eps_ = 1.5e-7;
    std::vector<Scalar> pExact_;
};

} // end namespace Dumux

#endif
