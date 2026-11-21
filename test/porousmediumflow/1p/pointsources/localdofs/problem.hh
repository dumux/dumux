// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup OnePTests
 * \brief Single-phase flow with point sources and zero pressure on the boundary
 *
 * Darcy flow with unit permeability and unit viscosity on the unit domain, \f$ -\Delta p = \sum_k q_k \delta(\mathbf{x} - \mathbf{x}_k) \f$.
 * The point sources are set after construction.
 */
#ifndef DUMUX_ONEP_POINTSOURCE_LOCALDOFS_PROBLEM_HH
#define DUMUX_ONEP_POINTSOURCE_LOCALDOFS_PROBLEM_HH

#include <cstddef>
#include <memory>
#include <utility>
#include <vector>

#include <dune/common/fvector.hh>

#include <dumux/common/boundarytypes_.hh>
#include <dumux/common/constraintinfo.hh>
#include <dumux/common/indextraits.hh>
#include <dumux/common/numeqvector.hh>
#include <dumux/common/pointsources.hh>
#include <dumux/common/problemwithspatialparams.hh>
#include <dumux/common/properties.hh>
#include <dumux/discretization/dirichletconstraints.hh>

namespace Dumux {

template<class TypeTag>
class OnePPointSourceProblem : public Experimental::ProblemWithSpatialParams<TypeTag>
{
    using ParentType = Experimental::ProblemWithSpatialParams<TypeTag>;
    using GridDiscretization = GetPropType<TypeTag, Properties::GridGeometry>;
    using GridView = typename GridDiscretization::GridView;
    using Element = typename GridView::template Codim<0>::Entity;
    using GlobalPosition = typename Element::Geometry::GlobalCoordinate;
    using LocalPosition = typename Element::Geometry::LocalCoordinate;
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using PrimaryVariables = GetPropType<TypeTag, Properties::PrimaryVariables>;
    using NumEqVector = Dumux::NumEqVector<PrimaryVariables>;
    using ModelTraits = GetPropType<TypeTag, Properties::ModelTraits>;
    using BoundaryTypes = Dumux::Experimental::BoundaryTypes<ModelTraits::numEq()>;
    using ConstraintInfo = Dumux::DirichletConstraintInfo<ModelTraits::numEq()>;
    using ConstraintValues = Dune::FieldVector<Scalar, ModelTraits::numEq()>;
    using GridIndexType = typename IndexTraits<GridView>::GridIndex;
    using DirichletConstraintData = Dumux::DirichletConstraintData<ConstraintInfo, ConstraintValues, GridIndexType>;

    //! the strength of a point source is given by the index of its context
    struct Strengths
    {
        std::vector<Scalar> values;

        template<class ElementDiscretization, class ElementVariables, class Context>
        NumEqVector operator()(const ElementDiscretization&, const ElementVariables&, const Context& context) const
        { return NumEqVector(values[context.id()]); }
    };

    using PointSourceContext = Experimental::IdPointSourceContext<Experimental::PointSourceContext<GlobalPosition>, std::size_t>;
    using PointSources = Experimental::PointSources<PointSourceContext, Strengths, LocalPosition>;

public:
    OnePPointSourceProblem(std::shared_ptr<const GridDiscretization> gridDiscretization)
    : ParentType(gridDiscretization)
    {
        appendBoundaryDirichletConstraints(
            constraints_, this->gridDiscretization(),
            [&] (const auto& elemDisc, const auto& boundaryFace, const auto& localDof)
            {
                ConstraintInfo info;
                info.setAll();
                return info;
            },
            [&] (const auto& elemDisc, const auto& localDof)
            { return ConstraintValues(0.0); }
        );
    }

    //! Replace the point sources by sources of the given strengths at the given positions
    void setPointSources(const std::vector<GlobalPosition>& positions, std::vector<Scalar> strengths)
    {
        std::vector<PointSourceContext> contexts;
        for (std::size_t i = 0; i < positions.size(); ++i)
            contexts.emplace_back(Experimental::PointSourceContext<GlobalPosition>(positions[i]), i);
        pointSources_ = Experimental::makePointSources(this->gridDiscretization(), contexts, Strengths{std::move(strengths)});
    }

    const PointSources& pointSources() const
    { return pointSources_; }

    BoundaryTypes boundaryTypesAtPos(const GlobalPosition& globalPos) const
    { return BoundaryTypes(); }

    NumEqVector sourceAtPos(const GlobalPosition& globalPos) const
    { return NumEqVector(0.0); }

    PrimaryVariables initialAtPos(const GlobalPosition& globalPos) const
    { return PrimaryVariables(0.0); }

    const std::vector<DirichletConstraintData>& constraints() const
    { return constraints_; }

private:
    std::vector<DirichletConstraintData> constraints_;
    PointSources pointSources_;
};

} // end namespace Dumux

#endif
