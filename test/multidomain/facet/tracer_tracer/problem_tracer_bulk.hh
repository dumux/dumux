// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*****************************************************************************
 *   See the file COPYING for full copying permissions.                      *
 *                                                                           *
 *   This program is free software: you can redistribute it and/or modify    *
 *   it under the terms of the GNU General Public License as published by    *
 *   the Free Software Foundation, either version 2 of the License, or       *
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
/**
 * \file
 * \ingroup TracerTests
 * \brief The problem for the bulk domain of the tracer facet coupling test
 */
#ifndef DUMUX_TEST_TPFAFACETCOUPLING_TRACER_BULK_PROBLEM_HH
#define DUMUX_TEST_TPFAFACETCOUPLING_TRACER_BULK_PROBLEM_HH

#include <dune/alugrid/grid.hh>

#include <dumux/multidomain/facet/box/properties.hh>
#include <dumux/multidomain/facet/cellcentered/tpfa/properties.hh>

#include <dumux/porousmediumflow/tracer/model.hh>
#include <dumux/porousmediumflow/problem.hh>

#include <dumux/material/fluidsystems/base.hh>

#include "spatialparams_tracer.hh"
#include "tracerfluidsystem.hh"
#include "tracermodeltraits.hh"

namespace Dumux {

template <class TypeTag>
class TracerBulkProblem;

namespace Properties {
NEW_TYPE_TAG(TracerTestBulk, INHERITS_FROM(Tracer));

// define the type tags
NEW_TYPE_TAG(TracerBulkTpfa, INHERITS_FROM(TracerTestBulk, CCTpfaFacetCouplingModel));
NEW_TYPE_TAG(TracerBulkBox, INHERITS_FROM(TracerTestBulk, BoxFacetCouplingModel));

// Set the grid type
SET_TYPE_PROP(TracerTestBulk, Grid, Dune::ALUGrid<2, 2, Dune::simplex, Dune::conforming>);

//! Overwrite the advection type property
SET_TYPE_PROP(TracerBulkTpfa, AdvectionType, StationaryVelocityField<typename GET_PROP_TYPE(TypeTag, Scalar)>);
SET_TYPE_PROP(TracerBulkBox, AdvectionType, StationaryVelocityField<typename GET_PROP_TYPE(TypeTag, Scalar)>);

// Set the problem property
SET_TYPE_PROP(TracerTestBulk, Problem, TracerBulkProblem<TypeTag>);

// Set the spatial parameters
SET_PROP(TracerTestBulk, SpatialParams)
{
    using FVGridGeometry = typename GET_PROP_TYPE(TypeTag, FVGridGeometry);
    using Scalar = typename GET_PROP_TYPE(TypeTag, Scalar);
    using type = TracerSpatialParams<FVGridGeometry, Scalar>;
};

// Define whether mole(true) or mass (false) fractions are used
SET_BOOL_PROP(TracerTestBulk, UseMoles, false);

//! set the model traits (with disabled diffusion)
SET_PROP(TracerTestBulk, ModelTraits)
{
private:
    using FluidSystem = typename GET_PROP_TYPE(TypeTag, FluidSystem);
public:
    using type = TracerTestModelTraits<FluidSystem::numComponents, GET_PROP_VALUE(TypeTag, UseMoles)>;
};

// use the test-specific fluid system
SET_TYPE_PROP(TracerTestBulk, FluidSystem, TracerFluidSystem<TypeTag>);
} // end namespace Properties


/*!
 * \ingroup MultiDomain
 * \ingroup TracerModel
 * \brief The problem for the bulk domain of the tracer facet coupling test
 */
template <class TypeTag>
class TracerBulkProblem : public PorousMediumFlowProblem<TypeTag>
{
    using ParentType = PorousMediumFlowProblem<TypeTag>;

    using Scalar = typename GET_PROP_TYPE(TypeTag, Scalar);
    using Indices = typename GET_PROP_TYPE(TypeTag, ModelTraits)::Indices;
    using GridView = typename GET_PROP_TYPE(TypeTag, GridView);
    using FVGridGeometry = typename GET_PROP_TYPE(TypeTag, FVGridGeometry);
    using FVElementGeometry = typename FVGridGeometry::LocalView;
    using BoundaryTypes = typename GET_PROP_TYPE(TypeTag, BoundaryTypes);
    using PrimaryVariables = typename GET_PROP_TYPE(TypeTag, PrimaryVariables);
    using NumEqVector = typename GET_PROP_TYPE(TypeTag, NumEqVector);
    using FluidSystem = typename GET_PROP_TYPE(TypeTag, FluidSystem);
    using CouplingManager = typename GET_PROP_TYPE(TypeTag, CouplingManager);

    //! property that defines whether mole or mass fractions are used
    static constexpr bool useMoles = GET_PROP_VALUE(TypeTag, UseMoles);

    using Element = typename FVGridGeometry::GridView::template Codim<0>::Entity;
    using GlobalPosition = typename Element::Geometry::GlobalCoordinate;

public:
    using typename ParentType::SpatialParams;

    TracerBulkProblem(std::shared_ptr<const FVGridGeometry> fvGridGeom,
                      std::shared_ptr<SpatialParams> spatialParams,
                      std::shared_ptr<CouplingManager> couplingManager,
                      const std::string& paramGroup = "")
    : ParentType(fvGridGeom, spatialParams, paramGroup)
    , couplingManagerPtr_(couplingManager)
    , initialMassFraction_(getParamFromGroup<Scalar>(paramGroup, "Problem.ContaminationMassFraction"))
    {
        // stating in the console whether mole or mass fractions are used
        const auto problemName = getParamFromGroup<std::string>(this->paramGroup(), "Problem.Name");
        std::cout<< "problem " << problemName << " uses " << (useMoles ? "mole" : "mass") << " fractions" << '\n';
        problemName_  =  getParamFromGroup<std::string>(this->paramGroup(), "Vtk.OutputName") + "_" + problemName;
    }

    //! The problem name.
    const std::string& name() const
    { return problemName_; }

    /*!
     * \brief Specifies which kind of boundary condition should be
     *        used for which equation on a given boundary segment.
     * \param globalPos The position for which the bc type should be evaluated
     */
    BoundaryTypes boundaryTypesAtPos(const GlobalPosition &globalPos) const
    {
        BoundaryTypes values;
        values.setAllNeumann();
        return values;
    }

    /*!
     * \brief Evaluate the initial value for a control volume.
     * \param globalPos The position for which the initial condition should be evaluated
     */
    PrimaryVariables initialAtPos(const GlobalPosition &globalPos) const
    {
        if (isInContaminatedRegion_(globalPos))
            return PrimaryVariables(initialMassFraction_);
        return PrimaryVariables(0.0);
    }

    /*!
     * \brief Evaluate the boundary conditions for a neumann
     *        boundary segment.
     *
     * \param element The finite element
     * \param fvGeometry The finite-volume geometry
     * \param elemVolVars All volume variables for the element
     * \param scvf The sub control volume face
     *
     * Negative values mean influx.
     * E.g. for the mass balance that would the mass flux in \f$ [ kg / (m^2 \cdot s)] \f$.
     */
    template<class ElementVolumeVariables, class SubControlVolumeFace>
    NumEqVector neumann(const Element& element,
                        const FVElementGeometry& fvGeometry,
                        const ElementVolumeVariables& elemVolVars,
                        const SubControlVolumeFace& scvf) const
    {
        // get the volume flux on this segment
        const auto flux = this->spatialParams().volumeFlux(element, fvGeometry, elemVolVars, scvf);
        if (flux > 0.0)
        {
            const auto& insideVolVars = elemVolVars[fvGeometry.scv(scvf.insideScvIdx())];
            const auto tracerFlux = insideVolVars.massFraction(/*phaseIdx*/0, /*compIdx*/0)*flux;
            return NumEqVector(tracerFlux);
        }

        return NumEqVector(0.0);
    }

    //! returns reference to the coupling manager.
    const CouplingManager& couplingManager() const
    { return *couplingManagerPtr_; }

private:
    bool isInContaminatedRegion_(const GlobalPosition& globalPos) const
    {
        return globalPos[0] > 0.75 && globalPos[0] < 1.3
               && globalPos[1] > 5.3 && globalPos[1] < 5.8;
    }

    std::shared_ptr<CouplingManager> couplingManagerPtr_;
    Scalar initialMassFraction_;
    std::string problemName_;
};

} //end namespace Dumux

#endif
