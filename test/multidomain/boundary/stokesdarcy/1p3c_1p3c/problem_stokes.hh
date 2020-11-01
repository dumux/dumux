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
 * \ingroup BoundaryTests
 * \brief A Navier-Stokes test problem using Maxwell-Stefan diffusion.
 */

#ifndef DUMUX_STOKES_SUBPROBLEM_ONEPTHREEC_HH
#define DUMUX_STOKES_SUBPROBLEM_ONEPTHREEC_HH

#include <dune/grid/yaspgrid.hh>

#include <dumux/discretization/staggered/freeflow/properties.hh>

#include <dumux/freeflow/compositional/navierstokesncmodel.hh>
#include <dumux/freeflow/navierstokes/boundarytypes.hh>
#include <dumux/freeflow/navierstokes/problem.hh>

#include <dumux/multidomain/boundary/stokesdarcy/couplingdata.hh>

#include <dumux/flux/maxwellstefanslaw.hh>

#include "h2n2co2fluidsystem.hh"

namespace Dumux {
template <class TypeTag>
class StokesSubProblem;

namespace Properties {
// Create new type tags
namespace TTag {
struct StokesOnePThreeC { using InheritsFrom = std::tuple<NavierStokesNC, StaggeredFreeFlowModel>; };
} // end namespace TTag

// Set the fluid system
template<class TypeTag>
struct FluidSystem<TypeTag, TTag::StokesOnePThreeC> { using type = FluidSystems::H2N2CO2FluidSystem<GetPropType<TypeTag, Properties::Scalar>>; };

// Set the grid type
template<class TypeTag>
struct Grid<TypeTag, TTag::StokesOnePThreeC> { using type = Dune::YaspGrid<2, Dune::EquidistantOffsetCoordinates<GetPropType<TypeTag, Properties::Scalar>, 2> >; };

// Set the problem property
template<class TypeTag>
struct Problem<TypeTag, TTag::StokesOnePThreeC> { using type = Dumux::StokesSubProblem<TypeTag> ; };

template<class TypeTag>
struct EnableGridGeometryCache<TypeTag, TTag::StokesOnePThreeC> { static constexpr bool value = true; };
template<class TypeTag>
struct EnableGridFluxVariablesCache<TypeTag, TTag::StokesOnePThreeC> { static constexpr bool value = true; };
template<class TypeTag>
struct EnableGridVolumeVariablesCache<TypeTag, TTag::StokesOnePThreeC> { static constexpr bool value = true; };

// Use moles
template<class TypeTag>
struct UseMoles<TypeTag, TTag::StokesOnePThreeC> { static constexpr bool value = true; };

// Set the grid type
template<class TypeTag>
struct MolecularDiffusionType<TypeTag, TTag::StokesOnePThreeC> { using type = MaxwellStefansLaw<TypeTag>; };

// Do not replace one equation with a total mass balance
template<class TypeTag>
struct ReplaceCompEqIdx<TypeTag, TTag::StokesOnePThreeC> { static constexpr int value = 3; };
} // end namespace Properties

/*!
 * \ingroup BoundaryTests
 * \brief Test problem for the 1pnc (Navier-) Stokes problem.
 *
 * Horizontal flow from left to right with a parabolic velocity profile.
 */
template <class TypeTag>
class StokesSubProblem : public NavierStokesProblem<TypeTag>
{
    using ParentType = NavierStokesProblem<TypeTag>;
    using GridView = typename GetPropType<TypeTag, Properties::GridGeometry>::GridView;
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using Indices = typename GetPropType<TypeTag, Properties::ModelTraits>::Indices;
    using FluidSystem = GetPropType<TypeTag, Properties::FluidSystem>;
    using BoundaryTypes = Dumux::NavierStokesBoundaryTypes<GetPropType<TypeTag, Properties::ModelTraits>::numEq()>;
    using GridGeometry = GetPropType<TypeTag, Properties::GridGeometry>;
    using FVElementGeometry = typename GridGeometry::LocalView;
    using SubControlVolumeFace = typename FVElementGeometry::SubControlVolumeFace;
    using PrimaryVariables = GetPropType<TypeTag, Properties::PrimaryVariables>;
    using NumEqVector = GetPropType<TypeTag, Properties::NumEqVector>;
    using DiffusionCoefficientAveragingType = typename StokesDarcyCouplingOptions::DiffusionCoefficientAveragingType;

    using Element = typename GridView::template Codim<0>::Entity;
    using GlobalPosition = typename Element::Geometry::GlobalCoordinate;

    using CouplingManager = GetPropType<TypeTag, Properties::CouplingManager>;

public:
    StokesSubProblem(std::shared_ptr<const GridGeometry> gridGeometry, std::shared_ptr<CouplingManager> couplingManager)
    : ParentType(gridGeometry, "Stokes"), eps_(1e-6),
    couplingManager_(couplingManager)
    {
        inletVelocity_ = getParamFromGroup<Scalar>(this->paramGroup(), "Problem.Velocity");
        pressure_ = getParamFromGroup<Scalar>(this->paramGroup(), "Problem.Pressure");
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
    { return 273.15 + 10; } // 10°C

    // \}

   /*!
     * \name Boundary conditions
     */
    // \{

    /*!
     * \brief Specifies which kind of boundary condition should be
     *        used for which equation on a given boundary segment.
     *
     * \param element The finite element
     * \param scvf The sub control volume face
     */
    BoundaryTypes boundaryTypes(const Element& element,
                                const SubControlVolumeFace& scvf) const
    {
        BoundaryTypes values;

        const auto& globalPos = scvf.dofPosition();

        if(onLeftBoundary_(globalPos))
        {
            values.setDirichlet(Indices::conti0EqIdx + 2);
            values.setDirichlet(Indices::conti0EqIdx + 1);
            values.setDirichlet(Indices::velocityXIdx);
            values.setDirichlet(Indices::velocityYIdx);
        }

        else if(onRightBoundary_(globalPos))
        {
            values.setDirichlet(Indices::pressureIdx);
            values.setOutflow(Indices::conti0EqIdx + 1);
            values.setOutflow(Indices::conti0EqIdx + 2);

        }
        else
        {
            values.setDirichlet(Indices::velocityXIdx);
            values.setDirichlet(Indices::velocityYIdx);
            values.setNeumann(Indices::conti0EqIdx);
            values.setNeumann(Indices::conti0EqIdx + 1);
            values.setNeumann(Indices::conti0EqIdx + 2);
        }

        if(couplingManager().isCoupledEntity(couplingManager().stokesIdx, scvf))
        {
            values.setNeumann(Indices::conti0EqIdx);
            values.setNeumann(Indices::conti0EqIdx+1);
            values.setNeumann(Indices::conti0EqIdx+2);
            values.setNeumann(Indices::momentumYBalanceIdx);
            values.setBeaversJoseph(Indices::velocityXIdx);
        }

        return values;
    }

    /*!
     * \brief Evaluates the boundary conditions for a Dirichlet control volume.
     */
    PrimaryVariables dirichletAtPos(const GlobalPosition& globalPos) const
    {
        PrimaryVariables values(0.0);
        values = initialAtPos(globalPos);
        return values;
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
    template<class ElementVolumeVariables, class ElementFaceVariables>
    NumEqVector neumann(const Element& element,
                        const FVElementGeometry& fvGeometry,
                        const ElementVolumeVariables& elemVolVars,
                        const ElementFaceVariables& elemFaceVars,
                        const SubControlVolumeFace& scvf) const
    {
        NumEqVector values(0.0);

        if(couplingManager().isCoupledEntity(couplingManager().stokesIdx, scvf))
        {
            values[Indices::momentumYBalanceIdx] = couplingManager().couplingData().momentumCouplingCondition(element, fvGeometry, elemVolVars, elemFaceVars, scvf);

            const auto tmp = couplingManager().couplingData().massCouplingCondition(element, fvGeometry, elemVolVars, elemFaceVars, scvf, DiffusionCoefficientAveragingType::harmonic);
            values[Indices::conti0EqIdx] = tmp[0];
            values[Indices::conti0EqIdx + 1] = tmp[1];
            values[Indices::conti0EqIdx + 2] = tmp[2];
        }
        return values;
    }

    // \}

    //! Get the coupling manager
    const CouplingManager& couplingManager() const
    { return *couplingManager_; }

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
        PrimaryVariables values(0.0);
        values[Indices::pressureIdx] = pressure_;
        values[Indices::velocityXIdx] = inletVelocity_ * (globalPos[1] - this->gridGeometry().bBoxMin()[1])
                                              * (this->gridGeometry().bBoxMax()[1] - globalPos[1])
                                              / (0.25 * (this->gridGeometry().bBoxMax()[1] - this->gridGeometry().bBoxMin()[1])
                                              * (this->gridGeometry().bBoxMax()[1] - this->gridGeometry().bBoxMin()[1]));
        return values;
    }

    /*!
     * \brief Returns the intrinsic permeability of required as input parameter
     *        for the Beavers-Joseph-Saffman boundary condition.
     */
    Scalar permeability(const Element& element, const SubControlVolumeFace& scvf) const
    {
        return couplingManager().couplingData().darcyPermeability(element, scvf);
    }

    /*!
     * \brief Returns the alpha value required as input parameter for the
     *        Beavers-Joseph-Saffman boundary condition.
     */
    Scalar alphaBJ(const SubControlVolumeFace& scvf) const
    {
        return 1.0;
    }

   // \}

private:
    bool onLeftBoundary_(const GlobalPosition &globalPos) const
    { return globalPos[0] < this->gridGeometry().bBoxMin()[0] + eps_; }

    bool onRightBoundary_(const GlobalPosition &globalPos) const
    { return globalPos[0] > this->gridGeometry().bBoxMax()[0] - eps_; }

    bool onLowerBoundary_(const GlobalPosition &globalPos) const
    { return globalPos[1] < this->gridGeometry().bBoxMin()[1] + eps_; }

    bool onUpperBoundary_(const GlobalPosition &globalPos) const
    { return globalPos[1] > this->gridGeometry().bBoxMax()[1] - eps_; }

    Scalar eps_;
    Scalar inletVelocity_;
    Scalar pressure_;
    std::shared_ptr<CouplingManager> couplingManager_;
};
} // end namespace Dumux

#endif
