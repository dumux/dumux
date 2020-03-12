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
 * \brief The problem for the bulk domain in the elastic
 *        single-phase facet coupling test.
 */
#ifndef DUMUX_NETWORK_2D_BULK_FLOW_PROBLEM_HH
#define DUMUX_NETWORK_2D_BULK_FLOW_PROBLEM_HH

#include <dune/alugrid/grid.hh>

#include <dumux/material/components/constant.hh>
#include <dumux/material/fluidsystems/1pliquid.hh>

#include <dumux/multidomain/facet/cellcentered/tpfa/properties.hh>
#include <dumux/multidomain/facet/cellcentered/mpfa/properties.hh>
#include <dumux/multidomain/facet/box/properties.hh>

#include <dumux/porousmediumflow/problem.hh>
#include <dumux/porousmediumflow/1p/model.hh>

#include "spatialparams_bulk_onep.hh"

namespace Dumux {
// forward declarations
template<class TypeTag> class OnePBulkProblem;

namespace Properties {
// create the type tag nodes
namespace TTag {
struct OnePBulk { using InheritsFrom = std::tuple<OneP>; };
struct OnePBulkTpfa { using InheritsFrom = std::tuple<CCTpfaFacetCouplingModel, OnePBulk>; };
struct OnePBulkMpfa { using InheritsFrom = std::tuple<CCMpfaFacetCouplingModel, OnePBulk>; };
struct OnePBulkBox { using InheritsFrom = std::tuple<BoxFacetCouplingModel, OnePBulk>; };
} // end namespace TTag

// Set the grid type
template<class TypeTag>
struct Grid<TypeTag, TTag::OnePBulk> { using type = Dune::ALUGrid<2, 2, Dune::simplex, Dune::nonconforming>; };

// Set the problem type
template<class TypeTag>
struct Problem<TypeTag, TTag::OnePBulk> { using type = OnePBulkProblem<TypeTag>; };

// set the spatial params
template<class TypeTag>
struct SpatialParams<TypeTag, TTag::OnePBulk>
{
private:
    using FVGridGeometry = GetPropType<TypeTag, Properties::FVGridGeometry>;
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using CouplingManager = GetPropType<TypeTag, Properties::CouplingManager>;

public:
    using type = OnePBulkSpatialParams<FVGridGeometry, Scalar, CouplingManager>;
};

// the fluid system
template<class TypeTag>
struct FluidSystem<TypeTag, TTag::OnePBulk>
{
private:
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
public:
    using type = FluidSystems::OnePLiquid< Scalar, Components::Constant<0, Scalar> >;
};

// solution-independent permeabilities
template<class TypeTag>
struct SolutionDependentAdvection<TypeTag, TTag::OnePBulk>
{ static constexpr bool value = false; };

} // end namespace Properties

/*!
 * \brief The problem for the bulk domain in the elastic
 *        single-phase facet coupling test.
 */
template<class TypeTag>
class OnePBulkProblem : public PorousMediumFlowProblem<TypeTag>
{
    using ParentType = PorousMediumFlowProblem<TypeTag>;

    using NumEqVector = GetPropType<TypeTag, Properties::NumEqVector>;
    using GridVariables = GetPropType<TypeTag, Properties::GridVariables>;
    using PrimaryVariables = typename GridVariables::PrimaryVariables;
    using Scalar = typename GridVariables::Scalar;

    using FVGridGeometry = typename GridVariables::GridGeometry;
    using FVElementGeometry = typename FVGridGeometry::LocalView;
    using SubControlVolumeFace = typename FVGridGeometry::SubControlVolumeFace;
    using GridView = typename FVGridGeometry::GridView;
    using Element = typename GridView::template Codim<0>::Entity;
    using GlobalPosition = typename Element::Geometry::GlobalCoordinate;

    using BoundaryTypes = GetPropType<TypeTag, Properties::BoundaryTypes>;
    using CouplingManager = GetPropType<TypeTag, Properties::CouplingManager>;

    static constexpr bool isBox = FVGridGeometry::discMethod == DiscretizationMethod::box;

public:
    //! The constructor
    OnePBulkProblem(std::shared_ptr<const FVGridGeometry> fvGridGeometry,
                    std::shared_ptr<typename ParentType::SpatialParams> spatialParams,
                    std::shared_ptr<CouplingManager> couplingManagerPtr,
                    const std::string& paramGroup = "")
    : ParentType(fvGridGeometry, spatialParams, paramGroup)
    , couplingManagerPtr_(couplingManagerPtr)
    {
        problemName_  =  getParam<std::string>("Vtk.OutputName") + "_" + getParamFromGroup<std::string>(this->paramGroup(), "Problem.Name");
        extractionPressure_ = getParam<Scalar>("Problem.ExtractionPressure");
        injectionOverPressure_ = getParam<Scalar>("Problem.InjectionOverPressure");
        porosity_.resize(fvGridGeometry->gridView().size(0), getParam<Scalar>("SpatialParams.InitialPorosity"));
    }

    //! The problem name.
    const std::string& name() const
    { return problemName_; }

    //! Specifies the kind of boundary condition on a given boundary position.
    BoundaryTypes boundaryTypesAtPos(const GlobalPosition& globalPos) const
    {
        BoundaryTypes values;
        values.setAllNeumann();

        if ( isOnInlet(globalPos) || isOnOutlet(globalPos) )
            values.setAllDirichlet();

        return values;
    }

    //! Specifies which kind of interior boundary condition should be
    //! used for which equation on a given sub-control volume face
    //! that couples to a facet element.
    BoundaryTypes interiorBoundaryTypes(const Element& element, const SubControlVolumeFace& scvf) const
    {
        BoundaryTypes values;
        static const bool useInteriorDirichlet = getParam<bool>("Problem.UseInteriorDirichlet");
        if (useInteriorDirichlet)
            values.setAllDirichlet();
        else
            values.setAllNeumann();
        return values;
    }

    //! Evaluates the Dirichlet boundary conditions at a given position.
    PrimaryVariables dirichletAtPos(const GlobalPosition& globalPos) const
    {
        auto values = initialAtPos(globalPos);
        if (isOnInlet(globalPos))
            values[0] += injectionOverPressure_;
        return values;
    }

    //! evaluate the Neumann boundary conditions
    template<class ElementVolumeVariables>
    NumEqVector neumann(const Element& element,
                        const FVElementGeometry& fvGeometry,
                        const ElementVolumeVariables& elemVolVars,
                        const SubControlVolumeFace& scvf) const
    { return NumEqVector(0.0); }

    //! Evaluates the initial conditions.
    PrimaryVariables initialAtPos(const GlobalPosition& globalPos) const
    { return PrimaryVariables(extractionPressure_); }

    //! Returns the temperature in \f$\mathrm{[K]}\f$ in the domain.
    Scalar temperature() const
    { return 283.15; /*10°C*/ }

    //! Returns const reference to the coupling manager.
    const CouplingManager& couplingManager() const
    { return *couplingManagerPtr_; }

    //! The inlet is on the left side of the domain
    bool isOnInlet(const GlobalPosition& globalPos) const
    { return globalPos[0] < this->fvGridGeometry().bBoxMin()[0] + 1e-6; }

    //! The inlet is on the left side of the domain
    bool isOnOutlet(const GlobalPosition& globalPos) const
    { return globalPos[0] > this->fvGridGeometry().bBoxMax()[0] - 1e-6; }

    //! adds additional output fields to a vtk writer
    template<class OutputModule>
    void addOutputFields(OutputModule& outputModule)
    {
        outputModule.addField(porosity_, "porosity", OutputModule::FieldType::element);
    }

    //! updates the output fields for a given solution
    template<class SolutionVector>
    void updateOutputFields(const SolutionVector& x)
    {
        for (const auto& element : elements(this->fvGridGeometry().gridView()))
        {
            const auto elemSol = elementSolution(element, x, this->fvGridGeometry());
            auto fvGeometry = localView(this->fvGridGeometry());
            fvGeometry.bindElement(element);

            Scalar porosity = 0.0;
            for (const auto& scv : scvs(fvGeometry))
                porosity += this->spatialParams().porosity(element, scv, elemSol);
            porosity /= fvGeometry.numScv();

            porosity_[this->fvGridGeometry().elementMapper().index(element)] = porosity;
        }
    }

private:
    std::string problemName_;
    std::shared_ptr<CouplingManager> couplingManagerPtr_;
    Scalar extractionPressure_;
    Scalar injectionOverPressure_;
    std::vector<Scalar> porosity_;
};

} // end namespace Dumux

#endif
