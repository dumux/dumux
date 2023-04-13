// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup FacetTests
 * \brief The problem for the bulk domain of the single-phase problem
 *        within the tracer facet coupling test.
 */

#ifndef DUMUX_TEST_TPFAFACETCOUPLING_TRACER_ONEP_BULKPROBLEM_HH
#define DUMUX_TEST_TPFAFACETCOUPLING_TRACER_ONEP_BULKPROBLEM_HH

#include <dumux/common/boundarytypes.hh>
#include <dumux/common/properties.hh>
#include <dumux/common/parameters.hh>
#include <dumux/common/numeqvector.hh>

#include <dumux/porousmediumflow/problem.hh>

namespace Dumux {

/*!
 * \ingroup FacetTests
 * \brief The problem for the bulk domain of the single-phase problem
 *        within the tracer facet coupling test.
 */
template<class TypeTag>
class OnePBulkProblem : public PorousMediumFlowProblem<TypeTag>
{
    using ParentType = PorousMediumFlowProblem<TypeTag>;

    using GridVariables = GetPropType<TypeTag, Properties::GridVariables>;
    using PrimaryVariables = typename GridVariables::PrimaryVariables;
    using Scalar = typename GridVariables::Scalar;

    using GridGeometry = typename GridVariables::GridGeometry;
    using FVElementGeometry = typename GridGeometry::LocalView;
    using SubControlVolume = typename GridGeometry::SubControlVolume;
    using SubControlVolumeFace = typename GridGeometry::SubControlVolumeFace;
    using GridView = typename GridGeometry::GridView;
    using Element = typename GridView::template Codim<0>::Entity;
    using GlobalPosition = typename Element::Geometry::GlobalCoordinate;

    using BoundaryTypes = Dumux::BoundaryTypes<GetPropType<TypeTag, Properties::ModelTraits>::numEq()>;
    using CouplingManager = GetPropType<TypeTag, Properties::CouplingManager>;
    using NumEqVector = Dumux::NumEqVector<PrimaryVariables>;
    using Indices = typename GetPropType<TypeTag, Properties::ModelTraits>::Indices;

public:
    OnePBulkProblem(std::shared_ptr<const GridGeometry> gridGeometry,
                    std::shared_ptr<typename ParentType::SpatialParams> spatialParams,
                    std::shared_ptr<CouplingManager> couplingManager,
                    const std::string& paramGroup = "")
    : ParentType(gridGeometry, spatialParams, paramGroup)
    , couplingManagerPtr_(couplingManager)
    , overPressure_(getParamFromGroup<Scalar>(paramGroup, "Problem.BoundaryOverpressure"))
    {
        problemName_  =  getParamFromGroup<std::string>(this->paramGroup(), "Vtk.OutputName") + "_" +
                         getParamFromGroup<std::string>(this->paramGroup(), "Problem.Name");
    }

    //! The problem name.
    const std::string& name() const
    { return problemName_; }

    //! Specifies the type of boundary condition at a given position.
    BoundaryTypes boundaryTypesAtPos(const GlobalPosition &globalPos) const
    {
        BoundaryTypes values;
        values.setAllNeumann();
        if ( (globalPos[0] < 1e-6 && globalPos[1] < 4.0) || globalPos[0] > this->gridGeometry().bBoxMax()[0] - 1e-6)
            values.setAllDirichlet();
        return values;
    }

    //! Specifies the type of interior boundary condition at a given position.
    BoundaryTypes interiorBoundaryTypes(const Element& element, const SubControlVolumeFace& scvf) const
    {
        BoundaryTypes values;
        values.setAllNeumann();
        return values;
    }

    //! Evaluates the source term at a given position.
    NumEqVector sourceAtPos(const GlobalPosition& globalPos) const
    { return NumEqVector(0.0); }

    //! Evaluates the Dirichlet boundary condition for a given position.
    PrimaryVariables dirichletAtPos(const GlobalPosition& globalPos) const
    {
        auto values = initialAtPos(globalPos);
        if (globalPos[0] < 1e-6)
            values[Indices::pressureIdx] += overPressure_;
        return values;
    }

    //! Evaluates the Neumann boundary condition for a boundary segment.
    NumEqVector neumannAtPos(const GlobalPosition& globalPos) const
    { return NumEqVector(0.0); }

    //! Evaluates the initial conditions.
    PrimaryVariables initialAtPos(const GlobalPosition& globalPos) const
    { return PrimaryVariables(1.0e5); }

    //! Returns reference to the coupling manager.
    const CouplingManager& couplingManager() const
    { return *couplingManagerPtr_; }

private:
    std::shared_ptr<CouplingManager> couplingManagerPtr_;
    Scalar overPressure_;
    std::string problemName_;
};

} // end namespace Dumux

#endif
