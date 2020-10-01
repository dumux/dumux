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
 * \brief A simple Darcy test problem (cell-centered finite volume method).
 */

#ifndef DUMUX_DARCY_SUBPROBLEM_HH
#define DUMUX_DARCY_SUBPROBLEM_HH

#include <dune/grid/yaspgrid.hh>

#include <dumux/discretization/cctpfa.hh>

#include <dumux/common/boundarytypes.hh>

#include <dumux/material/fluidmatrixinteractions/diffusivityconstanttortuosity.hh>
#include <dumux/material/fluidsystems/1padapter.hh>
#include <dumux/material/fluidsystems/h2oair.hh>

#include <dumux/porousmediumflow/1pnc/model.hh>
#include <dumux/porousmediumflow/problem.hh>

#include "spatialparams.hh"

namespace Dumux {
template <class TypeTag>
class DarcySubProblem;

namespace Properties {
// Create new type tags
namespace TTag {
struct DarcyOnePTwoC { using InheritsFrom = std::tuple<OnePNC, CCTpfaModel>; };
} // end namespace TTag

// Set the problem property
template<class TypeTag>
struct Problem<TypeTag, TTag::DarcyOnePTwoC> { using type = Dumux::DarcySubProblem<TypeTag>; };

// The fluid system
template<class TypeTag>
struct FluidSystem<TypeTag, TTag::DarcyOnePTwoC>
{
  using H2OAir = FluidSystems::H2OAir<GetPropType<TypeTag, Properties::Scalar>>;
  static constexpr auto phaseIdx = H2OAir::liquidPhaseIdx; // simulate the water phase
  using type = FluidSystems::OnePAdapter<H2OAir, phaseIdx>;
};

// Use moles
template<class TypeTag>
struct UseMoles<TypeTag, TTag::DarcyOnePTwoC> { static constexpr bool value = true; };

// Do not replace one equation with a total mass balance
template<class TypeTag>
struct ReplaceCompEqIdx<TypeTag, TTag::DarcyOnePTwoC> { static constexpr int value = 3; };

//! Use a model with constant tortuosity for the effective diffusivity
template<class TypeTag>
struct EffectiveDiffusivityModel<TypeTag, TTag::DarcyOnePTwoC>
{ using type = DiffusivityConstantTortuosity<GetPropType<TypeTag, Properties::Scalar>>; };

// Set the grid type
template<class TypeTag>
struct Grid<TypeTag, TTag::DarcyOnePTwoC> { using type = Dune::YaspGrid<2>; };

// Set the spatial paramaters type
template<class TypeTag>
struct SpatialParams<TypeTag, TTag::DarcyOnePTwoC>
{
    using GridGeometry = GetPropType<TypeTag, Properties::GridGeometry>;
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using type = OnePSpatialParams<GridGeometry, Scalar>;
};
} // end namespace Dumux

template <class TypeTag>
class DarcySubProblem : public PorousMediumFlowProblem<TypeTag>
{
    using ParentType = PorousMediumFlowProblem<TypeTag>;
    using GridGeometry = GetPropType<TypeTag, Properties::GridGeometry>;
    using GridView = typename GridGeometry::GridView;
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using PrimaryVariables = GetPropType<TypeTag, Properties::PrimaryVariables>;
    using NumEqVector = GetPropType<TypeTag, Properties::NumEqVector>;
    using BoundaryTypes = Dumux::BoundaryTypes<GetPropType<TypeTag, Properties::ModelTraits>::numEq()>;
    using FVElementGeometry = typename GetPropType<TypeTag, Properties::GridGeometry>::LocalView;
    using SubControlVolume = typename FVElementGeometry::SubControlVolume;
    using SubControlVolumeFace = typename FVElementGeometry::SubControlVolumeFace;
    using Element = typename GridView::template Codim<0>::Entity;
    using Indices = typename GetPropType<TypeTag, Properties::ModelTraits>::Indices;

    using GlobalPosition = typename Element::Geometry::GlobalCoordinate;

    using CouplingManager = GetPropType<TypeTag, Properties::CouplingManager>;
    using TimeLoopPtr = std::shared_ptr<TimeLoop<Scalar>>;

public:
    DarcySubProblem(std::shared_ptr<const GridGeometry> gridGeometry,
                   std::shared_ptr<CouplingManager> couplingManager)
    : ParentType(gridGeometry, "Darcy"), eps_(1e-7), couplingManager_(couplingManager)
    {
        problemName_  =  getParam<std::string>("Vtk.OutputName") + "_" + getParamFromGroup<std::string>(this->paramGroup(), "Problem.Name");

        // determine whether to simulate a vertical or horizontal flow configuration
        verticalFlow_ = problemName_.find("vertical") != std::string::npos;
    }

    /*!
     * \brief The problem name.
     */
    const std::string& name() const
    {
        return problemName_;
    }

    /*!
     * \brief Returns the temperature within the domain in [K].
     *
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
      *        used for which equation on a given boundary control volume.
      *
      * \param element The element
      * \param scvf The boundary sub control volume face
      */
    BoundaryTypes boundaryTypes(const Element& element, const SubControlVolumeFace& scvf) const
    {
        BoundaryTypes values;
        values.setAllNeumann();

        if (couplingManager().isCoupledEntity(CouplingManager::darcyIdx, scvf))
            values.setAllCouplingNeumann();

        if (verticalFlow_)
        {
            if (onLowerBoundary_(scvf.center()))
                values.setAllDirichlet();
        }

        return values;
    }

    /*!
     * \brief Evaluates the boundary conditions for a Dirichlet control volume.
     *
     * \param element The element for which the Dirichlet boundary condition is set
     * \param scvf The boundary subcontrolvolumeface
     *
     */
    PrimaryVariables dirichlet(const Element &element, const SubControlVolumeFace &scvf) const
    {
        PrimaryVariables values(0.0);
        values = initial(element);

        if (verticalFlow_)
        {
            // Check if this a pure diffusion problem.
            static const bool isDiffusionProblem = problemName_.find("diffusion") != std::string::npos;

            Scalar bottomMoleFraction = 0.0;

            if (isDiffusionProblem)
            {
                // For the diffusion problem, change the top mole fraction after some time
                // in order to revert the concentration gradient.
                if (time() >= 1e10)
                    bottomMoleFraction = 1e-3;
            }

            if(onLowerBoundary_(scvf.center()))
                values[Indices::conti0EqIdx + 1] = bottomMoleFraction;
        }

        return values;
    }

    /*!
     * \brief Evaluates the boundary conditions for a Neumann control volume.
     *
     * \param element The element for which the Neumann boundary condition is set
     * \param fvGeometry The fvGeometry
     * \param elemVolVars The element volume variables
     * \param elemFluxVarsCache Flux variables caches for all faces in stencil
     * \param scvf The boundary sub control volume face
     *
     * For this method, the \a values variable stores primary variables.
     */
    template<class ElementVolumeVariables, class ElementFluxVarsCache>
    NumEqVector neumann(const Element& element,
                        const FVElementGeometry& fvGeometry,
                        const ElementVolumeVariables& elemVolVars,
                        const ElementFluxVarsCache& elemFluxVarsCache,
                        const SubControlVolumeFace& scvf) const
    {
        NumEqVector values(0.0);

        if (couplingManager().isCoupledEntity(CouplingManager::darcyIdx, scvf))
            values = couplingManager().couplingData().massCouplingCondition(element, fvGeometry, elemVolVars, scvf);

        return values;
    }

    // \}

    /*!
     * \name Volume terms
     */
    // \{
    /*!
     * \brief Evaluates the source term for all phases within a given
     *        sub control volume.
     *
     * \param element The element for which the source term is set
     * \param fvGeometry The fvGeometry
     * \param elemVolVars The element volume variables
     * \param scv The sub control volume
     */
    template<class ElementVolumeVariables>
    NumEqVector source(const Element &element,
                       const FVElementGeometry& fvGeometry,
                       const ElementVolumeVariables& elemVolVars,
                       const SubControlVolume &scv) const
    { return NumEqVector(0.0); }

    // \}

    /*!
     * \brief Evaluates the initial value for a control volume.
     *
     * \param element The element
     *
     * For this method, the \a priVars parameter stores primary
     * variables.
     */
    PrimaryVariables initial(const Element &element) const
    {
        PrimaryVariables values(0.0);
        values[Indices::pressureIdx] = 1e5;

        return values;
    }

    // \}

    //! Get the coupling manager
    const CouplingManager& couplingManager() const
    { return *couplingManager_; }

    /*!
     * \brief Sets the time loop pointer.
     */
    void setTimeLoop(TimeLoopPtr timeLoop)
    { timeLoop_ = timeLoop; }

    /*!
     * \brief Returns the time.
     */
    Scalar time() const
    { return timeLoop_->time(); }

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
    std::string problemName_;
    bool verticalFlow_;
    std::shared_ptr<CouplingManager> couplingManager_;
    TimeLoopPtr timeLoop_;
};
} // end namespace Dumux

#endif
