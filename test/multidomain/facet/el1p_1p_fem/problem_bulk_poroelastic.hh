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
 * \brief The problem for the poromechanical domain in the elastic
 *        single-phase facet coupling test.
 */
#ifndef DUMUX_ANALYTIC_CRACK_BULK_POROELASTIC_PROBLEM_HH
#define DUMUX_ANALYTIC_CRACK_BULK_POROELASTIC_PROBLEM_HH

#include <dune/alugrid/grid.hh>
#include <dune/functions/functionspacebases/lagrangebasis.hh>

#include <dumux/geomechanics/poroelastic/model.hh>
#include <dumux/geomechanics/fvproblem.hh>

#include <dumux/discretization/fem.hh>
#include <dumux/discretization/fem/ipdata.hh>
#include <dumux/discretization/fem/elementsolution.hh>
#include <dumux/discretization/fem/fegridgeometry.hh>
#include <dumux/multidomain/facet/box/properties.hh>

#include <dumux/material/fluidsystems/1pliquid.hh>
#include <dumux/material/components/constant.hh>

#include "spatialparams_bulk_poroelastic.hh"

// order of the ansatz space in FEM
#ifndef FEMORDER
#define FEMORDER 1
#endif

namespace Dumux {

// forward declaration of the problem class
template <class TypeTag>
class PoroElasticSubProblem;

namespace Properties {

// Create new type tags
namespace TTag {
struct PoroElasticBulk { using InheritsFrom = std::tuple<PoroElastic>; };
struct PoroElasticBulkFem { using InheritsFrom = std::tuple<PoroElasticBulk, FiniteElementModel>; };
struct PoroElasticBulkBox { using InheritsFrom = std::tuple<BoxFacetCouplingModel, PoroElasticBulk>; };

} // end namespace TTag
// Set the grid type
template<class TypeTag>
struct Grid<TypeTag, TTag::PoroElasticBulk> { using type = Dune::ALUGrid<2, 2, Dune::simplex, Dune::nonconforming>; };

// Set the problem property
template<class TypeTag>
struct Problem<TypeTag, TTag::PoroElasticBulk> { using type = Dumux::PoroElasticSubProblem<TypeTag>; };

// The fluid phase consists of one constant component
template<class TypeTag>
struct FluidSystem<TypeTag, TTag::PoroElasticBulk>
{
    using type = Dumux::FluidSystems::OnePLiquid< GetPropType<TypeTag, Properties::Scalar>,
                                                  Dumux::Components::Constant<0, GetPropType<TypeTag, Properties::Scalar>> >;
};
// The spatial parameters property
template<class TypeTag>
struct SpatialParams<TypeTag, TTag::PoroElasticBulk>
{
    using type = PoroElasticSpatialParams< GetPropType<TypeTag, Properties::Scalar>,
                                           GetPropType<TypeTag, Properties::GridGeometry> >;
};

// We use a lagrange basis of first order here
template<class TypeTag>
struct GridGeometry<TypeTag, TTag::PoroElasticBulkFem>
{
private:
    using GridView = GetPropType<TypeTag, Properties::GridView>;
    using FEBasis = Dune::Functions::LagrangeBasis<GridView, FEMORDER>;
public:
    using type = FEGridGeometry<FEBasis>;
};

// TODO: Remove this once fix is finished in dumux
template<class TypeTag>
struct FVGridGeometry<TypeTag, TTag::PoroElasticBulkFem>
{ using type = GetPropType<TypeTag, Properties::GridGeometry>; };

} // end namespace Properties

/*!
 * \brief The problem for the poromechanical domain in the elastic
 *        single-phase facet coupling test.
 */
template<class TypeTag>
class PoroElasticSubProblem : public GeomechanicsFVProblem<TypeTag>
{
    using ParentType = GeomechanicsFVProblem<TypeTag>;

    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using Indices = typename GetPropType<TypeTag, Properties::ModelTraits>::Indices;
    using BoundaryTypes = GetPropType<TypeTag, Properties::BoundaryTypes>;
    using CouplingManager = GetPropType<TypeTag, Properties::CouplingManager>;

    using GridVariables = GetPropType<TypeTag, Properties::GridVariables>;
    using PrimaryVariables = GetPropType<TypeTag, Properties::PrimaryVariables>;
    using NumEqVector = GetPropType<TypeTag, Properties::NumEqVector>;

    using GridGeometry = GetPropType<TypeTag, Properties::GridGeometry>;
    using GGLocalView = typename GridGeometry::LocalView;

    using GridView = GetPropType<TypeTag, Properties::GridView>;
    using Element = typename GridView::template Codim<0>::Entity;
    using Intersection = typename GridView::Intersection;
    using GlobalPosition = typename Element::Geometry::GlobalCoordinate;

    using StressType = GetPropType<TypeTag, Properties::StressType>;
    using Force = typename StressType::ForceVector;

    static constexpr bool isFem = GridGeometry::discMethod == DiscretizationMethod::fem;

public:
    PoroElasticSubProblem(std::shared_ptr<const GridGeometry> gridGeometry,
                          std::shared_ptr<typename ParentType::SpatialParams> spatialParams,
                          std::shared_ptr<CouplingManager> couplingManagerPtr,
                          const std::string& paramGroup = "")
    : ParentType(gridGeometry, spatialParams, paramGroup)
    , couplingManagerPtr_(couplingManagerPtr)
    {
        problemName_  =  getParam<std::string>("Vtk.OutputName") + "_" + getParamFromGroup<std::string>(this->paramGroup(), "Problem.Name");
        boundaryStress_ = getParamFromGroup<Scalar>(paramGroup, "Problem.InjectionOverPressure");

        sigma_x_.resize(gridGeometry->gridView().size(0), Force(0.0));
        sigma_y_.resize(gridGeometry->gridView().size(0), Force(0.0));
    }

    //! The problem name.
    const std::string& name() const
    { return problemName_; }

    //! Returns the temperature in the domain.
    static constexpr Scalar temperature()
    { return 273.15; }

    //! Evaluates the initial conditions for a given position.
    PrimaryVariables initialAtPos(const GlobalPosition& globalPos) const
    { return PrimaryVariables(0.0); }

    //! Evaluates the boundary conditions for a given position.
    PrimaryVariables dirichletAtPos(const GlobalPosition& globalPos) const
    { return PrimaryVariables(0.0); }

    //! Evaluate the boundary conditions for a neumann boundary segment (fem implementation).
    template<class ElementSolution, class IpData, class SecondaryVariables>
    NumEqVector neumann(const Element& element,
                        const Intersection& is,
                        const ElementSolution& elemSol,
                        const IpData& ipData,
                        const SecondaryVariables& secVars) const
    {
        PrimaryVariables values(0.0);

        // apply contact traction on interior boundaries
        // neglect facet pressure here!!
        if (couplingManager().isOnInteriorBoundary(element, is))
        {
            // Contribution of the fluid pressure
            const auto& facetVolVars = couplingManager().getLowDimVolVars(element, is, ipData.ipGlobal());
            Force force = is.centerUnitOuterNormal();
            force *= facetVolVars.pressure();
            force *= -1.0;

            // Contribution of the contact mechanics
            force += couplingManager().getContactTraction(element, is, ipData.ipGlobal());

            values = force;
        }

        // apply boundary pressures on inlet/outlet
        else if (onStressBoundary_(ipData.ipGlobal()))
        {
            values = is.centerUnitOuterNormal();
            values *= boundaryStress_;
            values *= -1.0;
        }

        return values;
    }

    //! Evaluates the boundary conditions for a Neumann boundary segment (box implementation).
    template<class ElementVolumeVariables, class SubControlVolumeFace>
    NumEqVector neumann(const Element& element,
                        const GGLocalView& fvGeometry,
                        const ElementVolumeVariables& elemVolVars,
                        const SubControlVolumeFace& scvf) const
    {
        PrimaryVariables values(0.0);

        // apply boundary pressures on inlet/outlet
        if (onStressBoundary_(scvf.ipGlobal()))
        {
            values = scvf.unitOuterNormal();
            values *= boundaryStress_;
            values *= -1.0;
        }

        return values;
    }

    //! Specifies which kind of interior boundary condition should be
    //! used for which equation on a given sub-control volume face
    //! that couples to a facet element. (Only relevant for box)
    template<class SubControlVolumeFace>
    BoundaryTypes interiorBoundaryTypes(const Element& element,
                                        const SubControlVolumeFace& scvf) const
    {
        BoundaryTypes values;
        values.setAllNeumann();
        return values;
    }

    //! Applies the forces stemming from pressure and contact forces (box implementation).
    template<class ElementVolumeVariables, class SubControlVolumeFace>
    Force interiorBoundaryForce(const Element& element,
                                const GGLocalView& fvGeometry,
                                const ElementVolumeVariables& elemVolVars,
                                const SubControlVolumeFace& scvf) const
    {
        const auto& facetVolVars = couplingManager().getLowDimVolVars(element, scvf);

        // Contribution of the fluid pressure
        Force force = scvf.unitOuterNormal();
        force *= facetVolVars.pressure();
        force *= scvf.area();
        force *= -1.0;

        // Contribution of the contact mechanics
        force += couplingManager().getContactForce(element, scvf);

        return force;
    }

    //! Returns the effective fluid density (fem implementation).
    template<class IpData, class ElemSol>
    Scalar effectiveFluidDensity(const Element& element,
                                 const IpData& ipData,
                                 const ElemSol& elemSol) const
    { return couplingManager().getPMFlowVolVars(element, ipData.ipGlobal()).density(); }

    //! Returns the effective pore pressure (fem implementation).
    template< class IpData, class ElemSol, class SecondaryVariables >
    Scalar effectivePorePressure(const Element& element,
                                 const GGLocalView& feGeometry,
                                 const ElemSol& elemSol,
                                 const IpData& ipData,
                                 const SecondaryVariables& secVars) const
    { return couplingManager().getPMFlowVolVars(element, ipData.ipGlobal()).pressure(); }

    //! Returns the effective fluid density in an scv (box implementation).
    template<class SubControlVolume>
    Scalar effectiveFluidDensity(const Element& element,
                                 const SubControlVolume& scv) const
    { return couplingManager().getPMFlowVolVars(element, scv.center()).density(); }

    //! Returns the effective pore pressure at a flux integration point (box implementation).
    template< class ElementVolumeVariables, class FluxVarsCache >
    Scalar effectivePorePressure(const Element& element,
                                 const GGLocalView& fvGeometry,
                                 const ElementVolumeVariables& elemVolVars,
                                 const FluxVarsCache& fluxVarsCache) const
    { return couplingManager().getPMFlowVolVars(element, fluxVarsCache.ipGlobal()).pressure(); }

    //! Specifies which kind of boundary condition should be
    //! used for which equation on a given boundary segment.
    BoundaryTypes boundaryTypesAtPos(const GlobalPosition& globalPos) const
    {
        BoundaryTypes values;
        values.setAllNeumann();

        static const auto xMax = this->gridGeometry().bBoxMax()[0];
        static const auto yMin = this->gridGeometry().bBoxMin()[1];
        static const auto yMax = this->gridGeometry().bBoxMax()[1];

        if (globalPos[0] > xMax - 1e-6)
            values.setDirichlet(Indices::uxIdx);
        if (globalPos[1] < yMin + 1e-6 || globalPos[1] > yMax - 1e-6)
            values.setDirichlet(Indices::uyIdx);

        return values;
    }

    //! Returns reference to the coupling manager.
    const CouplingManager& couplingManager() const
    { return *couplingManagerPtr_; }

    //! update the output fields (fem implementation)
    template<class GridVariables, class Assembler, class SolutionVector, class Id,
             bool isFE = isFem, std::enable_if_t<isFE, int> = 0>
    void updateOutputFields(const GridVariables& gridVariables,
                            const Assembler& assembler,
                            const SolutionVector& x,
                            Id id)
    {
        using AnsatzBasis = typename GridGeometry::AnsatzSpaceBasis;
        using FiniteElement = typename AnsatzBasis::LocalView::Tree::FiniteElement;
        using LocalBasis = typename FiniteElement::Traits::LocalBasisType;

        for (const auto& element : elements(this->gridGeometry().gridView()))
        {
            const auto eIdx = this->gridGeometry().elementMapper().index(element);

            auto feGeometry = localView(this->gridGeometry());
            feGeometry.bind(element);

            const auto& eg = element.geometry();
            const auto elemSol = elementSolution(element, x, this->gridGeometry());
            const auto& localBasis = feGeometry.feBasisLocalView().tree().finiteElement().localBasis();

            FEIntegrationPointData<GlobalPosition, LocalBasis> ipData(eg, eg.local(eg.center()), localBasis);
            typename GridVariables::SecondaryVariables secVars;
            secVars.update(elemSol, *this, element, ipData);

            const auto sigma = StressType::effectiveStressTensor(*this, element, feGeometry, elemSol, ipData, secVars);
            sigma_x_[eIdx] = sigma[Indices::uxIdx];
            sigma_y_[eIdx] = sigma[Indices::uyIdx];
        }
    }

    //! update the output fields (box implementation)
    template<class GridVariables, class Assembler, class SolutionVector, class Id,
             bool isFE = isFem, std::enable_if_t<!isFE, int> = 0>
    void updateOutputFields(const GridVariables& gridVariables,
                            const Assembler& assembler,
                            const SolutionVector& x,
                            Id id)
    {
        for (const auto& element : elements(this->fvGridGeometry().gridView()))
        {
            const auto eIdx = this->fvGridGeometry().elementMapper().index(element);

            auto fvGeometry = localView(this->fvGridGeometry());
            auto elemVolVars = localView(gridVariables.curGridVolVars());

            couplingManagerPtr_->bindCouplingContext(id, element, assembler);
            fvGeometry.bindElement(element);
            elemVolVars.bindElement(element, fvGeometry, x);

            using StressVariablesCache = GetPropType<TypeTag, Properties::FluxVariablesCache>;
            StressVariablesCache cache;
            cache.update(*this, element, fvGeometry, elemVolVars, element.geometry().center());

            const auto sigma = StressType::effectiveStressTensor(*this, element, fvGeometry, elemVolVars, cache);
            sigma_x_[eIdx] = sigma[Indices::uxIdx];
            sigma_y_[eIdx] = sigma[Indices::uyIdx];
        }
    }

    //! return references to the stress tensors
    const std::vector<Force>& sigma_x() const { return sigma_x_; }
    const std::vector<Force>& sigma_y() const { return sigma_y_; }

private:
    bool onStressBoundary_(const GlobalPosition& globalPos) const
    { return globalPos[0] < this->gridGeometry().bBoxMin()[0] + 1e-6; }

    std::shared_ptr<CouplingManager> couplingManagerPtr_;
    static constexpr Scalar eps_ = 3e-6;
    std::string problemName_;

    Scalar boundaryStress_;
    std::vector<Force> sigma_x_;
    std::vector<Force> sigma_y_;
};

} // end namespace Dumux

#endif
