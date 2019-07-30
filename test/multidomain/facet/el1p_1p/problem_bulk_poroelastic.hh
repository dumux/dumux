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

#include <dumux/geomechanics/poroelastic/model.hh>
#include <dumux/geomechanics/fvproblem.hh>

#include <dumux/multidomain/facet/box/properties.hh>
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
struct PoroElasticBulk { using InheritsFrom = std::tuple<BoxFacetCouplingModel, PoroElastic>; };
} // end namespace TTag
// Set the grid type
template<class TypeTag>
struct Grid<TypeTag, TTag::PoroElasticBulk> { using type = Dune::ALUGrid<2, 2, Dune::cube, Dune::nonconforming>; };
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
                                           GetPropType<TypeTag, Properties::FVGridGeometry> >;
};
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
    using PrimaryVariables = GetPropType<TypeTag, Properties::PrimaryVariables>;
    using ElementVolumeVariables = typename GetPropType<TypeTag, Properties::GridVolumeVariables>::LocalView;

    using FVGridGeometry = GetPropType<TypeTag, Properties::FVGridGeometry>;
    using FVElementGeometry = typename FVGridGeometry::LocalView;
    using SubControlVolume = typename FVGridGeometry::SubControlVolume;
    using SubControlVolumeFace = typename FVGridGeometry::SubControlVolumeFace;

    using GridView = GetPropType<TypeTag, Properties::GridView>;
    using Element = typename GridView::template Codim<0>::Entity;
    using GlobalPosition = typename Element::Geometry::GlobalCoordinate;

    using StressType = GetPropType<TypeTag, Properties::StressType>;
    using Force = typename StressType::ForceVector;

public:
    PoroElasticSubProblem(std::shared_ptr<const FVGridGeometry> fvGridGeometry,
                          std::shared_ptr<typename ParentType::SpatialParams> spatialParams,
                          std::shared_ptr<CouplingManager> couplingManagerPtr,
                          const std::string& paramGroup = "")
    : ParentType(fvGridGeometry, spatialParams, paramGroup)
    , couplingManagerPtr_(couplingManagerPtr)
    {
        problemName_  =  getParam<std::string>("Vtk.OutputName") + "_" + getParamFromGroup<std::string>(this->paramGroup(), "Problem.Name");
        injectionPressure_ = getParam<Scalar>("Problem.InjectionPressure");

        sigma_x_.resize(fvGridGeometry->gridView().size(0), Force(0.0));
        sigma_y_.resize(fvGridGeometry->gridView().size(0), Force(0.0));
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

    //! Evaluates the boundary conditions for a Neumann boundary segment.
    PrimaryVariables neumann(const Element& element,
                             const FVElementGeometry& fvGeometry,
                             const ElementVolumeVariables& elemVolVars,
                             const SubControlVolumeFace& scvf) const
    {
        PrimaryVariables values(0.0);

        // apply boundary pressures on inlet/outlet
        if (isOnInlet_(scvf.ipGlobal()))
        {
            values = scvf.unitOuterNormal();
            values *= injectionPressure_;
            values *= -1.0;
        }
        // else if (isOnOutlet_(scvf.ipGlobal()))
        // {
        //     values = scvf.unitOuterNormal();
        //     values *= extractionPressure_;
        //     values *= -1.0;
        // }
        // // there might be neumann faces on the upper or lower boundary
        // // if there are fractures extending to them. In this case, compute
        // // the boundary force from the urrent stress such that displacement
        // // in y results in 0
        else if (isOnUpperBoundary_(scvf.ipGlobal()) || isOnLowerBoundary_(scvf.ipGlobal()))
        {
            // make element solution and set displacements on boundary to zero
            auto elemSol = elementSolution(element, elemVolVars, fvGeometry);

            for (const auto& curScvf : scvfs(fvGeometry))
                if ( curScvf.boundary() &&
                     (isOnUpperBoundary_(curScvf.ipGlobal()) || isOnLowerBoundary_(curScvf.ipGlobal())) )
                { elemSol[fvGeometry.scv(curScvf.insideScvIdx()).localDofIndex()][Indices::uyIdx] = 0.0; }

            // with this element sol, make new elemVolVars
            auto modifiedElemVolVars = elemVolVars;
            for (const auto& scv : scvs(fvGeometry))
                modifiedElemVolVars[scv].update(elemSol, *this, element, scv);

            using StressVariablesCache = GetPropType<TypeTag, Properties::FluxVariablesCache>;
            StressVariablesCache cache;
            cache.update(*this, element, fvGeometry, modifiedElemVolVars, scvf);

            // compute the force f = -sigma*n
            const auto sigma = StressType::stressTensor(*this, element, fvGeometry, modifiedElemVolVars, cache);
            sigma.mv(scvf.unitOuterNormal(), values);
        }

        return values;
    }

    /*!
     * \brief Returns the effective fluid density in an scv.
     */
    Scalar effectiveFluidDensity(const Element& element,
                                 const SubControlVolume& scv) const
    { return couplingManager().getPMFlowVolVars(element).density(); }

    /*!
     * \brief Returns the effective pore pressure at a flux integration point.
     */
    template< class FluxVarsCache >
    Scalar effectivePorePressure(const Element& element,
                                 const FVElementGeometry& fvGeometry,
                                 const ElementVolumeVariables& elemVolVars,
                                 const FluxVarsCache& fluxVarsCache) const
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

        const auto xMax = this->fvGridGeometry().bBoxMax()[0];
        const auto yMin = this->fvGridGeometry().bBoxMin()[1];
        const auto yMax = this->fvGridGeometry().bBoxMax()[1];

        if (globalPos[0] > xMax - 1e-6)
            values.setDirichlet(Indices::uxIdx);
        else if (globalPos[1] < yMin + 1e-6 || globalPos[1] > yMax - 1e-6)
            values.setDirichlet(Indices::uyIdx);

        return values;
    }

    /*!
     * \brief Specifies which kind of interior boundary condition should be
     *        used for which equation on a given sub-control volume face
     *        that couples to a facet element.
     *
     * \param element The finite element the scvf is embedded in
     * \param scvf The sub-control volume face
     */
    BoundaryTypes interiorBoundaryTypes(const Element& element, const SubControlVolumeFace& scvf) const
    {
        BoundaryTypes values;
        values.setAllNeumann();
        return values;
    }

    /*!
     * \brief Specifies which kind of interior boundary condition should be
     *        used for which equation on a given sub-control volume face
     *        that couples to a facet element.
     *
     * \param element The finite element the scvf is embedded in
     * \param scvf The sub-control volume face
     */
    template<class ElementVolumeVariables>
    Force interiorBoundaryForce(const Element& element,
                                const FVElementGeometry& fvGeometry,
                                const ElementVolumeVariables& elemVolVars,
                                const SubControlVolumeFace& scvf) const
    {
        const auto& facetVolVars = couplingManager().getLowDimVolVars(element, scvf);

        // Contribution of the fluid pressure
        Force force = scvf.unitOuterNormal();
        force *= -1.0;
        force *= facetVolVars.pressure();
        force *= scvf.area();

        // Contribution of the contact mechanics
        force += couplingManager().getContactForce(element, scvf);
        std::cout << "ContactForce = " << couplingManager().getContactForce(element, scvf) << std::endl;

        return force;
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
    //! The inlet is on the left side of the domain
    bool isOnInlet_(const GlobalPosition& globalPos) const
    { return globalPos[0] < this->fvGridGeometry().bBoxMin()[0] + 1e-6; }

    //! The outlet is on the right side of the domain
    bool isOnOutlet_(const GlobalPosition& globalPos) const
    { return globalPos[0] > this->fvGridGeometry().bBoxMax()[0] - 1e-6; }

    //! Returns true if position is on upper boundary
    bool isOnUpperBoundary_(const GlobalPosition& globalPos) const
    { return globalPos[1] > this->fvGridGeometry().bBoxMax()[1] - 1e-6; }

    //! Returns true if position is on lower boundary
    bool isOnLowerBoundary_(const GlobalPosition& globalPos) const
    { return globalPos[1] < this->fvGridGeometry().bBoxMin()[1] + 1e-6; }

    std::shared_ptr<CouplingManager> couplingManagerPtr_;
    static constexpr Scalar eps_ = 3e-6;
    std::string problemName_;

    Scalar injectionPressure_;
    std::vector<Force> sigma_x_;
    std::vector<Force> sigma_y_;
};

} // end namespace Dumux

#endif // DUMUX_TEST_FACETCOUPLING_ELONEP_ONEP_BULK_POROELASTIC_PROBLEM_HH
