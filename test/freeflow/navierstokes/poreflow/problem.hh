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
 * \brief Pore flow test for the staggered grid (Navier-)Stokes model.
 */

#ifndef DUMUX_POREFLOW_TEST_PROBLEM_HH
#define DUMUX_POREFLOW_TEST_PROBLEM_HH

#include <dumux/common/parameters.hh>
#include <dumux/common/properties.hh>

#include <dumux/freeflow/navierstokes/boundarytypes.hh>
#include <dumux/freeflow/navierstokes/problem.hh>

#include <dumux/freeflow/navierstokes/momentum/fluxhelper.hh>
#include <dumux/freeflow/navierstokes/scalarfluxhelper.hh>
#include <dumux/freeflow/navierstokes/mass/1p/advectiveflux.hh>

namespace Dumux {

/*!
 * \ingroup NavierStokesTests
 * \brief  Test problem for the one-phase (Navier-) Stokes problem in a channel.
 */
template <class TypeTag>
class PoreFlowTestProblem : public NavierStokesProblem<TypeTag>
{
    using ParentType = NavierStokesProblem<TypeTag>;

    using BoundaryTypes = typename ParentType::BoundaryTypes;
    using GridGeometry = GetPropType<TypeTag, Properties::GridGeometry>;
    using FVElementGeometry = typename GridGeometry::LocalView;
    using SubControlVolumeFace = typename GridGeometry::SubControlVolumeFace;
    using InitialValues = typename ParentType::InitialValues;
    using Sources = typename ParentType::Sources;
    using DirichletValues = typename ParentType::DirichletValues;
    using BoundaryFluxes = typename ParentType::BoundaryFluxes;
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using ModelTraits = GetPropType<TypeTag, Properties::ModelTraits>;
    using SolutionVector = GetPropType<TypeTag, Properties::SolutionVector>;

    static constexpr auto dimWorld = GridGeometry::GridView::dimensionworld;
    using Element = typename FVElementGeometry::Element;
    using GlobalPosition = typename Element::Geometry::GlobalCoordinate;

    using CouplingManager = GetPropType<TypeTag, Properties::CouplingManager>;

    // the types of outlet boundary conditions
    enum class OutletCondition
    {
        outflow, unconstrainedOutflow, neumannXdirichletY, neumannXneumannY
    };

public:
    using Indices = typename GetPropType<TypeTag, Properties::ModelTraits>::Indices;

    PoreFlowTestProblem(std::shared_ptr<const GridGeometry> gridGeometry,
                        std::shared_ptr<CouplingManager> couplingManager)
    : ParentType(gridGeometry, couplingManager)
    {
        dynamicViscosity_ = getParam<Scalar>("Component.LiquidKinematicViscosity", 1.0)
                            * getParam<Scalar>("Component.LiquidDensity", 1.0);
        outletPressure_ = getParam<Scalar>("Problem.OutletPressure", 1.1e5);
        pressureDifference_ = getParam<Scalar>("Problem.PressureDifference", 50.0);
    }

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
            if (isInlet_(globalPos) || isOutlet_(globalPos))
                values.setAllNeumann();
            else
                values.setAllDirichlet();
        }
        else
            values.setAllNeumann();

        return values;
    }


    /*!
     * \brief Evaluates the boundary conditions for a Dirichlet control volume.
     *
     * \param globalPos The center of the finite volume which ought to be set.
     */
    DirichletValues dirichlet(const Element& element, const SubControlVolumeFace& scvf) const
    { return DirichletValues(0.0); }

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
    BoundaryFluxes neumann(const Element& element,
                           const FVElementGeometry& fvGeometry,
                           const ElementVolumeVariables& elemVolVars,
                           const ElementFluxVariablesCache& elemFluxVarsCache,
                           const SubControlVolumeFace& scvf) const
    {
        BoundaryFluxes values(0.0);
        const auto& globalPos = scvf.ipGlobal();

        if constexpr (ParentType::isMomentumProblem())
        {
            using MomentumFluxHelper = NavierStokesMomentumBoundaryFluxHelper;
            if (isOutlet_(globalPos))
            {
                values = MomentumFluxHelper::fixedPressureMomentumFlux(*this, fvGeometry, scvf,
                                                                       elemVolVars, elemFluxVarsCache,
                                                                       outletPressure_, true);
            }
            else if (isInlet_(globalPos))
            {
                values = MomentumFluxHelper::fixedPressureMomentumFlux(*this, fvGeometry, scvf,
                                                                       elemVolVars, elemFluxVarsCache,
                                                                       outletPressure_ + pressureDifference_, true);
            }
        }
        else
        {
            values = NavierStokesScalarBoundaryFluxHelper<AdvectiveFlux<ModelTraits>>::scalarOutflowFlux(
                *this, element, fvGeometry, scvf, elemVolVars);
        }

        return values;
    }

    /*!
     * \brief Returns a reference pressure at a given sub control volume face.
     *        This pressure is substracted from the actual pressure for the momentum balance
     *        which potentially helps to improve numerical accuracy by avoiding issues related do floating point arithmetic.
     */
    Scalar referencePressure(const Element& element,
                             const FVElementGeometry& fvGeometry,
                             const SubControlVolumeFace& scvf) const
    { return outletPressure_; }

private:

    bool isInlet_(const GlobalPosition& globalPos) const
    { return isLeftThroat_(globalPos) || isTopThroat_(globalPos);}

    bool isOutlet_(const GlobalPosition& globalPos) const
    { return isRightThroat_(globalPos) || isBottomThroat_(globalPos);}

    bool isLeftThroat_(const GlobalPosition& globalPos) const
    { return globalPos[0] < this->gridGeometry().bBoxMin()[0] + eps_;}

    bool isRightThroat_(const GlobalPosition& globalPos) const
    { return globalPos[0] > this->gridGeometry().bBoxMax()[0] - eps_; }

    bool isBottomThroat_(const GlobalPosition& globalPos) const
    { return globalPos[1] < this->gridGeometry().bBoxMin()[1] + eps_; }

    bool isTopThroat_(const GlobalPosition& globalPos) const
    { return globalPos[1] > this->gridGeometry().bBoxMax()[1] - eps_; }

    static constexpr Scalar eps_=1e-6;
    Scalar dynamicViscosity_;
    Scalar outletPressure_;
    Scalar pressureDifference_;
};
} // end namespace Dumux

#endif
