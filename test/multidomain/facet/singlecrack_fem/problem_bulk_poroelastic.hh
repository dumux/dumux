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
 * \brief The problem for the bulk domain in the elastic
 *        single-phase facet coupling test.
 */
#ifndef DUMUX_TEST_FACETCOUPLING_ELONEP_ONEP_BULK_POROELASTIC_PROBLEM_HH
#define DUMUX_TEST_FACETCOUPLING_ELONEP_ONEP_BULK_POROELASTIC_PROBLEM_HH

#include <dune/alugrid/grid.hh>
#include <dune/functions/functionspacebases/lagrangebasis.hh>

#include <dumux/geomechanics/poroelastic/model.hh>
#include <dumux/geomechanics/fvproblem.hh>

#include <dumux/discretization/fem.hh>
#include <dumux/discretization/fem/ipdata.hh>
#include <dumux/discretization/fem/elementsolution.hh>
#include <dumux/discretization/fem/fegridgeometry.hh>
#include <dumux/material/fluidsystems/1pliquid.hh>
#include <dumux/material/components/constant.hh>

#include "spatialparams_bulk_poroelastic.hh"

namespace Dumux {

// forward declaration of the problem class
template <class TypeTag>
class PoroElasticSubProblem;

namespace Properties {

// Create new type tags
namespace TTag {
struct PoroElasticBulk { using InheritsFrom = std::tuple<PoroElastic, FiniteElementModel>; };
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
struct GridGeometry<TypeTag, TTag::PoroElasticBulk>
{
private:
    using GridView = GetPropType<TypeTag, Properties::GridView>;
    using FEBasis = Dune::Functions::LagrangeBasis<GridView, 2>;
public:
    using type = FEGridGeometry<FEBasis>;
};

// We use a lagrange basis of zero-th order here
template<class TypeTag>
struct FVGridGeometry<TypeTag, TTag::PoroElasticBulk>
{ using type = GetPropType<TypeTag, Properties::GridGeometry>; };

} // end namespace Properties

/*!
 * \ingroup FacetTests
 * \brief The problem for the bulk domain in the elastic
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
    using SecondaryVariables = typename GridVariables::SecondaryVariables;

    using PrimaryVariables = GetPropType<TypeTag, Properties::PrimaryVariables>;
    using NumEqVector = GetPropType<TypeTag, Properties::NumEqVector>;

    using GridGeometry = GetPropType<TypeTag, Properties::GridGeometry>;
    using FEElementGeometry = typename GridGeometry::LocalView;
    using AnsatzBasis = typename GridGeometry::AnsatzSpaceBasis;
    using FiniteElement = typename AnsatzBasis::LocalView::Tree::FiniteElement;
    using LocalBasis = typename FiniteElement::Traits::LocalBasisType;

    using GridView = GetPropType<TypeTag, Properties::GridView>;
    using Element = typename GridView::template Codim<0>::Entity;
    using Intersection = typename GridView::Intersection;
    using GlobalPosition = typename Element::Geometry::GlobalCoordinate;

    using StressType = GetPropType<TypeTag, Properties::StressType>;
    using Force = typename StressType::ForceVector;

public:
    PoroElasticSubProblem(std::shared_ptr<const GridGeometry> gridGeometry,
                          std::shared_ptr<typename ParentType::SpatialParams> spatialParams,
                          std::shared_ptr<CouplingManager> couplingManagerPtr,
                          const std::string& paramGroup = "")
    : ParentType(gridGeometry, spatialParams, paramGroup)
    , couplingManagerPtr_(couplingManagerPtr)
    {
        problemName_  =  getParam<std::string>("Vtk.OutputName") + "_" + getParamFromGroup<std::string>(this->paramGroup(), "Problem.Name");
        boundaryStress_ = getParamFromGroup<Scalar>(paramGroup, "Problem.BoundaryStress");

        sigma_x_.resize(gridGeometry->gridView().size(0), Force(0.0));
        sigma_y_.resize(gridGeometry->gridView().size(0), Force(0.0));
    }

    /*!
     * \brief The problem name.
     */
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

    /*!
     * \brief Evaluate the boundary conditions for a neumann
     *        boundary segment.
     *
     * \param element The grid element
     * \param is The intersection
     * \param elemSol The element solution vector
     * \param ipData Shape function values/gradients evaluated at an integration point
     * \param secVars The primary/secondary variables evaluated at an integration point
     *
     * Negative values mean influx.
     * E.g. for a mass balance that would the mass flux in \f$ [ kg / (m^2 \cdot s)] \f$.
     */
    template<class ElementSolution, class IpData>
    NumEqVector neumann(const Element& element,
                        const Intersection& is,
                        const ElementSolution& elemSol,
                        const IpData& ipData,
                        const SecondaryVariables& secVars) const
    {
        PrimaryVariables values(0.0);

        // apply facet domain pressure and/or contact traction on interior boundaries
        if (couplingManager().isOnInteriorBoundary(element, is))
        {
            // Contribution of the fluid pressure
            Force force = is.centerUnitOuterNormal();
            force *= couplingManager().getLowDimVolVars(element, is, ipData.ipGlobal()).pressure();

            // Contribution of the contact mechanics
            force += couplingManager().getContactTraction(element, is, ipData.ipGlobal());

            values = force;
        }

        // apply boundary pressures on inlet/outlet
        const auto& ip = ipData.ipGlobal();
        if (ip[0] < this->gridGeometry().bBoxMin()[0] + 1e-6
            || ip[0] > this->gridGeometry().bBoxMax()[0] - 1e-6)
        {
            values = is.centerUnitOuterNormal();
            values *= boundaryStress_;
        }

        return values;
    }

    /*!
     * \brief Returns the effective fluid density (box implementation).
     */
    template<class IpData, class ElemSol>
    Scalar effectiveFluidDensity(const Element& element,
                                 const IpData& ipData,
                                 const ElemSol& elemSol) const
    { return couplingManager().getPMFlowVolVars(element).density(); }

    /*!
     * \brief Returns the effective pore pressure (fem implementation).
     */
    template< class IpData, class ElemSol, class SecondaryVariables >
    Scalar effectivePorePressure(const Element& element,
                                 const FEElementGeometry& feGeometry,
                                 const ElemSol& elemSol,
                                 const IpData& ipData,
                                 const SecondaryVariables& secVars) const
    { return couplingManager().getPMFlowVolVars(element).pressure(); }

    /*!
     * \brief Specifies which kind of boundary condition should be
     *        used for which equation on a given boundary segment.
     *
     * \param globalPos The global position
     */
    BoundaryTypes boundaryTypesAtPos(const GlobalPosition& globalPos) const
    {
        BoundaryTypes values;
        values.setAllNeumann();

        static const auto xMax = this->gridGeometry().bBoxMax()[0];
        static const auto xMin = this->gridGeometry().bBoxMax()[0];
        static const auto dx = xMax - xMin;
        static const auto c = xMin + 0.5*dx;
        static const auto minThres = c - 0.1*dx;
        static const auto maxThres = c + 0.1*dx;

        if (globalPos[1] < this->gridGeometry().bBoxMin()[1] + 1e-6)
        {
            values.setDirichlet(Indices::uyIdx);
            const auto x = globalPos[0];
            if (x > minThres && x < maxThres)
                values.setDirichlet(Indices::uxIdx);
        }

        return values;
    }

    //! Returns reference to the coupling manager.
    const CouplingManager& couplingManager() const
    { return *couplingManagerPtr_; }

    //! update the output fields
    template<class GridVariables, class Assembler, class SolutionVector, class Id>
    void updateOutputFields(const GridVariables& gridVariables,
                            const Assembler& assembler,
                            const SolutionVector& x,
                            Id id)
    {
        for (const auto& element : elements(this->gridGeometry().gridView()))
        {
            const auto eIdx = this->gridGeometry().elementMapper().index(element);

            auto feGeometry = localView(this->gridGeometry());
            feGeometry.bind(element);

            const auto& eg = element.geometry();
            const auto elemSol = elementSolution(element, x, this->gridGeometry());
            const auto& localBasis = feGeometry.feBasisLocalView().tree().finiteElement().localBasis();

            FEIntegrationPointData<GlobalPosition, LocalBasis> ipData(eg, eg.local(eg.center()), localBasis);
            SecondaryVariables secVars;
            secVars.update(elemSol, *this, element, ipData);

            const auto sigma = StressType::effectiveStressTensor(*this, element, feGeometry, elemSol, ipData, secVars);
            sigma_x_[eIdx] = sigma[Indices::uxIdx];
            sigma_y_[eIdx] = sigma[Indices::uyIdx];
        }
    }

    //! return references to the stress tensors
    const std::vector<Force>& sigma_x() const { return sigma_x_; }
    const std::vector<Force>& sigma_y() const { return sigma_y_; }

private:
    std::shared_ptr<CouplingManager> couplingManagerPtr_;
    static constexpr Scalar eps_ = 3e-6;
    std::string problemName_;

    Scalar boundaryStress_;
    std::vector<Force> sigma_x_;
    std::vector<Force> sigma_y_;
};

} // end namespace Dumux

#endif // DUMUX_TEST_FACETCOUPLING_ELONEP_ONEP_BULK_POROELASTIC_PROBLEM_HH
