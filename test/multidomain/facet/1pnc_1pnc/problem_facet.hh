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
 * \ingroup FacetTests
 * \brief The problem for the lower-dimensional domain in the 1pnc facet coupling test.
 */
#ifndef DUMUX_TEST_TPFAFACETCOUPLING_ONEPNC_FACETPROBLEM_HH
#define DUMUX_TEST_TPFAFACETCOUPLING_ONEPNC_FACETPROBLEM_HH

#include <dumux/common/boundarytypes.hh>
#include <dumux/porousmediumflow/problem.hh>

namespace Dumux {

/*!
 * \ingroup FacetTests
 * \brief The lower-dimensional test problem for the
 *        1pnc model with coupling across the bulk grid facets.
 */
template<class TypeTag>
class OnePNCLowDimProblem : public PorousMediumFlowProblem<TypeTag>
{
    using ParentType = PorousMediumFlowProblem<TypeTag>;

    using GridVariables = GetPropType<TypeTag, Properties::GridVariables>;
    using ElementVolumeVariables = typename GridVariables::GridVolumeVariables::LocalView;
    using PrimaryVariables = typename GridVariables::PrimaryVariables;
    using Scalar = typename GridVariables::Scalar;

    using GridGeometry = typename GridVariables::GridGeometry;
    using FVElementGeometry = typename GridGeometry::LocalView;
    using SubControlVolume = typename GridGeometry::SubControlVolume;
    using GridView = typename GridGeometry::GridView;
    using Element = typename GridView::template Codim<0>::Entity;
    using GlobalPosition = typename Element::Geometry::GlobalCoordinate;

    using BoundaryTypes = Dumux::BoundaryTypes<GetPropType<TypeTag, Properties::ModelTraits>::numEq()>;
    using CouplingManager = GetPropType<TypeTag, Properties::CouplingManager>;
    using NumEqVector = GetPropType<TypeTag, Properties::NumEqVector>;
    using Indices = typename GetPropType<TypeTag, Properties::ModelTraits>::Indices;
    using FluidSystem = GetPropType<TypeTag, Properties::FluidSystem>;

    using ModelTraits = GetPropType<TypeTag, Properties::ModelTraits>;
    static constexpr bool enableHeatConduction = ModelTraits::enableEnergyBalance();

    enum
    {
        // indices of the primary variables
        pressureIdx = Indices::pressureIdx,
        N2Idx = FluidSystem::compIdx(FluidSystem::MultiPhaseFluidSystem::N2Idx)
    };

public:
    OnePNCLowDimProblem(std::shared_ptr<const GridGeometry> gridGeometry,
                        std::shared_ptr<typename ParentType::SpatialParams> spatialParams,
                        std::shared_ptr<CouplingManager> couplingManager,
                        const std::string& paramGroup = "")
    : ParentType(gridGeometry, spatialParams, paramGroup)
    , couplingManagerPtr_(couplingManager)
    , extrusion_(getParam<Scalar>("Problem.FacetExtrusion"))
    {
        problemName_  =  getParam<std::string>("Vtk.OutputName") + "_" +
                         getParamFromGroup<std::string>(this->paramGroup(), "Problem.Name");
    }

    //! The problem name.
    const std::string& name() const
    { return problemName_; }

    //! Specifies the type of boundary condition at a given position.
    BoundaryTypes boundaryTypesAtPos(const GlobalPosition& globalPos) const
    {
        BoundaryTypes values;
        values.setAllNeumann();
        return values;
    }

    //! Evaluates the source term at a given position.
    NumEqVector source(const Element& element,
                       const FVElementGeometry& fvGeometry,
                       const ElementVolumeVariables& elemVolVars,
                       const SubControlVolume& scv) const
    {
        // evaluate sources from bulk domain
        auto source = couplingManagerPtr_->evalSourcesFromBulk(element, fvGeometry, elemVolVars, scv);
        source /= scv.volume()*elemVolVars[scv].extrusionFactor();
        return source;
    }

    //! Evaluates the Dirichlet boundary condition for a given position.
    PrimaryVariables dirichletAtPos(const GlobalPosition& globalPos) const
    { return initialAtPos(globalPos); }

    //! Returns the extrusion factor.
    Scalar extrusionFactorAtPos(const GlobalPosition& globalPos) const
    { return extrusion_; }

    //! Evaluates the initial conditions.
    PrimaryVariables initialAtPos(const GlobalPosition& globalPos) const
    {
        PrimaryVariables values;
        values[pressureIdx] = 1.0e5;
        values[N2Idx] = 0.0;

        if constexpr (enableHeatConduction)
        {
            values[Indices::temperatureIdx] = temperature();
        }

        return values;
    }

    //! Returns the temperature in \f$\mathrm{[K]}\f$ in the domain.
    Scalar temperature() const
    { return 283.15; /*10°*/ }

    //! Returns reference to the coupling manager.
    const CouplingManager& couplingManager() const
    { return *couplingManagerPtr_; }

private:
    std::shared_ptr<CouplingManager> couplingManagerPtr_;
    Scalar extrusion_;
    std::string problemName_;
};

} // end namespace Dumux

#endif
