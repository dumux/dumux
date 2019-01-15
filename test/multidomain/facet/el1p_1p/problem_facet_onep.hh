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
 * \brief The problem for the (d-1)-dimensional facet domain in the elastic
 *        single-phase facet coupling test.
 */
#ifndef DUMUX_TEST_FACETCOUPLING_ELONEP_ONEP_FACET_FLOW_PROBLEM_HH
#define DUMUX_TEST_FACETCOUPLING_ELONEP_ONEP_FACET_FLOW_PROBLEM_HH

#include <dune/foamgrid/foamgrid.hh>

#include <dumux/material/components/simpleh2o.hh>
#include <dumux/material/fluidsystems/1pliquid.hh>

#include <dumux/discretization/cctpfa.hh>
#include <dumux/porousmediumflow/problem.hh>
#include <dumux/porousmediumflow/1p/model.hh>

#include "spatialparams_facet_onep.hh"

namespace Dumux {
// forward declarations
template<class TypeTag> class OnePFacetProblem;

namespace Properties {
// Create new type tags
namespace TTag {
struct OnePFacet { using InheritsFrom = std::tuple<OneP>; };
struct OnePFacetTpfa { using InheritsFrom = std::tuple<OnePFacet, CCTpfaModel>; };
} // end namespace TTag

// Set the grid type
template<class TypeTag>
struct Grid<TypeTag, TTag::OnePFacet> { using type = Dune::FoamGrid<1, 2>; };
// Set the problem type
template<class TypeTag>
struct Problem<TypeTag, TTag::OnePFacet> { using type = OnePFacetProblem<TypeTag>; };
// set the spatial params
template<class TypeTag>
struct SpatialParams<TypeTag, TTag::OnePFacet>
{
private:
    using FVGridGeometry = GetPropType<TypeTag, Properties::FVGridGeometry>;
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using CouplingManager = GetPropType<TypeTag, Properties::CouplingManager>;

public:
    using type = OnePFacetSpatialParams<FVGridGeometry, Scalar, CouplingManager>;
};

// the fluid system
template<class TypeTag>
struct FluidSystem<TypeTag, TTag::OnePFacet>
{
private:
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
public:
    using type = FluidSystems::OnePLiquid< Scalar, Components::SimpleH2O<Scalar> >;
};

} // end namespace Properties

/*!
 * \ingroup FacetTests
 * \brief The problem for the (d-1)-dimensional facet domain in the elastic
 *        single-phase facet coupling test.
 */
template<class TypeTag>
class OnePFacetProblem : public PorousMediumFlowProblem<TypeTag>
{
    using ParentType = PorousMediumFlowProblem<TypeTag>;

    using NumEqVector = GetPropType<TypeTag, Properties::NumEqVector>;
    using GridVariables = GetPropType<TypeTag, Properties::GridVariables>;
    using ElementVolumeVariables = typename GridVariables::GridVolumeVariables::LocalView;
    using PrimaryVariables = typename GridVariables::PrimaryVariables;
    using Scalar = typename GridVariables::Scalar;

    using FVGridGeometry = typename GridVariables::GridGeometry;
    using FVElementGeometry = typename FVGridGeometry::LocalView;
    using SubControlVolume = typename FVGridGeometry::SubControlVolume;
    using SubControlVolumeFace = typename FVGridGeometry::SubControlVolumeFace;
    using GridView = typename FVGridGeometry::GridView;
    using Element = typename GridView::template Codim<0>::Entity;
    using GlobalPosition = typename Element::Geometry::GlobalCoordinate;

    using BoundaryTypes = GetPropType<TypeTag, Properties::BoundaryTypes>;
    using CouplingManager = GetPropType<TypeTag, Properties::CouplingManager>;

public:
    OnePFacetProblem(std::shared_ptr<const FVGridGeometry> fvGridGeometry,
                     std::shared_ptr<typename ParentType::SpatialParams> spatialParams,
                     std::shared_ptr<CouplingManager> couplingManagerPtr,
                     const std::string& paramGroup = "")
    : ParentType(fvGridGeometry, spatialParams, paramGroup)
    , couplingManagerPtr_(couplingManagerPtr)
    , initialAperture_(getParamFromGroup<Scalar>(paramGroup, "SpatialParams.InitialAperture"))
    , injectionPressure_(getParamFromGroup<Scalar>(paramGroup, "Problem.InjectionPressure"))
    {
        problemName_  =  getParam<std::string>("Vtk.OutputName") + "_" + getParamFromGroup<std::string>(this->paramGroup(), "Problem.Name");

        apertures_.resize(fvGridGeometry->gridView().size(0), initialAperture_);
        permeabilities_.resize(fvGridGeometry->gridView().size(0), initialAperture_*initialAperture_/12.0);
    }

    /*!
     * \brief The problem name.
     */
    const std::string& name() const
    { return problemName_; }

    //! Specifies the kind of boundary condition at a boundary position.
    BoundaryTypes boundaryTypesAtPos(const GlobalPosition& globalPos) const
    {
        BoundaryTypes values;
        values.setAllNeumann();
        if (globalPos[0] < this->fvGridGeometry().bBoxMin()[0] + 1e-6)
            values.setAllDirichlet();
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
     * \brief Evaluates the source term for all phases within a given
     *        sub-control volume.
     */
    NumEqVector source(const Element& element,
                       const FVElementGeometry& fvGeometry,
                       const ElementVolumeVariables& elemVolVars,
                       const SubControlVolume& scv) const
    {
        // obtain the sources stemming from the bulk flow domain
        auto source = couplingManagerPtr_->evalSourcesFromBulk(element, fvGeometry, elemVolVars, scv);
        source /= scv.volume()*elemVolVars[scv].extrusionFactor();
        return source;
    }

    //! Sets the aperture as extrusion factor.
    template<class ElementSolution>
    Scalar extrusionFactor(const Element& element,
                           const SubControlVolume& scv,
                           const ElementSolution& elemSol) const
    { return couplingManager().computeAperture(element, scv, initialAperture_); }

    //! Evaluates the initial conditions.
    PrimaryVariables initialAtPos(const GlobalPosition& globalPos) const
    { return PrimaryVariables(1.0e5); }

    //! Evaluates the Dirichlet boundary conditions at a given position.
    PrimaryVariables dirichletAtPos(const GlobalPosition& globalPos) const
    {
        if (globalPos[0] < this->fvGridGeometry().bBoxMin()[0] + 1e-6)
            return PrimaryVariables(injectionPressure_);
        else
            return initialAtPos(globalPos);
    }

    //! Returns the temperature \f$\mathrm{[K]}\f$ for an isothermal problem.
    Scalar temperature() const
    { return 283.15; /*10°*/ }

    //! Returns const reference to the coupling manager.
    const CouplingManager& couplingManager() const
    { return *couplingManagerPtr_; }

    //! returns the vector of apertures
    const std::vector<Scalar>& apertures() const
    { return apertures_; }

    //! returns the vector of permeabilities
    const std::vector<Scalar>& permeabilities() const
    { return permeabilities_; }

    //! update the output fields
    template<class SolutionVector>
    void updateOutputFields(SolutionVector x)
    {
        for (const auto& element : elements(this->fvGridGeometry().gridView()))
        {
            const auto elemSol = elementSolution(element, x, this->fvGridGeometry());
            auto fvGeometry = localView(this->fvGridGeometry());
            fvGeometry.bindElement(element);

            for (const auto& scv : scvs(fvGeometry))
            {
                apertures_[scv.elementIndex()] = extrusionFactor(element, scv, elemSol);
                permeabilities_[scv.elementIndex()] = this->spatialParams().permeability(element, scv, elemSol);
            }
        }
    }

private:
    std::shared_ptr<CouplingManager> couplingManagerPtr_;
    std::string problemName_;
    Scalar initialAperture_;
    Scalar injectionPressure_;

    // fields to be added to output
    std::vector<Scalar> apertures_;
    std::vector<Scalar> permeabilities_;
};

} // end namespace Dumux

#endif
