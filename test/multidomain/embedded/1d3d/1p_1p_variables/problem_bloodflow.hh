// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup EmbeddedTests
 * \brief The vessel problem of the embedded 1d-3d test stated with the variables of the local dofs
 */
#ifndef DUMUX_EMBEDDED_1P1P_VARIABLES_PROBLEM_BLOODFLOW_HH
#define DUMUX_EMBEDDED_1P1P_VARIABLES_PROBLEM_BLOODFLOW_HH

#include <cmath>
#include <memory>
#include <utility>
#include <vector>

#include <dune/common/fvector.hh>

#include <dumux/common/boundarytypes_.hh>
#include <dumux/common/constraintinfo.hh>
#include <dumux/common/indextraits.hh>
#include <dumux/common/numeqvector.hh>
#include <dumux/common/parameters.hh>
#include <dumux/common/problemwithspatialparams.hh>
#include <dumux/common/properties.hh>
#include <dumux/discretization/dirichletconstraints.hh>

#include <dumux/common/pointsources.hh>
#include <dumux/multidomain/embedded/couplingmanager1d3d.hh> // for coupling mode
#include <dumux/multidomain/embedded/pointsourcecontext.hh>

namespace Dumux {

template<class TypeTag>
class BloodFlowProblem : public Experimental::ProblemWithSpatialParams<TypeTag>
{
    using ParentType = Experimental::ProblemWithSpatialParams<TypeTag>;
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using GridDiscretization = GetPropType<TypeTag, Properties::GridGeometry>;
    using GridView = typename GridDiscretization::GridView;
    using GlobalPosition = typename GridView::template Codim<0>::Entity::Geometry::GlobalCoordinate;
    using PrimaryVariables = GetPropType<TypeTag, Properties::PrimaryVariables>;
    using ModelTraits = GetPropType<TypeTag, Properties::ModelTraits>;
    using NumEqVector = Dumux::NumEqVector<PrimaryVariables>;
    using BoundaryTypes = Dumux::Experimental::BoundaryTypes<ModelTraits::numEq()>;
    using ConstraintInfo = Dumux::DirichletConstraintInfo<ModelTraits::numEq()>;
    using ConstraintValues = Dune::FieldVector<Scalar, ModelTraits::numEq()>;
    using GridIndexType = typename IndexTraits<GridView>::GridIndex;
    using DirichletConstraintData = Dumux::DirichletConstraintData<ConstraintInfo, ConstraintValues, GridIndexType>;
    using CouplingManager = GetPropType<TypeTag, Properties::CouplingManager>;

    struct Exchange
    {
        const CouplingManager* couplingManager;
        Scalar beta;

        template<class ElementDiscretization, class ElementVariables, class Context>
        NumEqVector operator()(const ElementDiscretization&, const ElementVariables&, const Context& context) const
        {
            const auto id = context.id();
            const Scalar pressure3D = couplingManager->bulkPriVars(id)[0];
            const Scalar pressure1D = couplingManager->lowDimPriVars(id)[0];
            Scalar value = beta*(pressure3D - pressure1D);
            if constexpr (CouplingManager::couplingMode == Embedded1d3dCouplingMode::kernel)
                value *= couplingManager->fluxScalingFactor(id);
            return NumEqVector(value*context.weight());
        }
    };

    using LocalPosition = typename GridView::template Codim<0>::Entity::Geometry::LocalCoordinate;
    using PointSources = Experimental::PointSources<EmbeddedCoupling::PointSourceContext<GlobalPosition>, Exchange, LocalPosition>;

public:
    BloodFlowProblem(std::shared_ptr<const GridDiscretization> gridDiscretization,
                     std::shared_ptr<CouplingManager> couplingManager)
    : ParentType(gridDiscretization, "Vessel")
    , couplingManager_(couplingManager)
    {
        const auto pIn = getParam<Scalar>("BoundaryConditions1D.PressureInput");
        const auto deltaP = getParam<Scalar>("BoundaryConditions1D.DeltaPressure");
        appendBoundaryDirichletConstraints(
            constraints_, this->gridDiscretization(),
            [&] (const auto&, const auto&, const auto&)
            {
                ConstraintInfo info;
                info.setAll();
                return info;
            },
            [&] (const auto& elemDisc, const auto& localDof)
            { return ConstraintValues(ipData(elemDisc, localDof).global()[2] > 0.5 ? pIn : pIn - deltaP); }
        );
    }

    //! Set the point sources after the coupling manager has computed them
    void setPointSources()
    {
        static const auto radius = getParam<Scalar>("SpatialParams.Radius");
        pointSources_ = couplingManager_->makePointSources(
            Dune::index_constant<1>{}, Exchange{couplingManager_.get(), 2*M_PI/(2*M_PI + std::log(radius))}
        );
    }

    const PointSources& pointSources() const
    { return pointSources_; }

    BoundaryTypes boundaryTypesAtPos(const GlobalPosition&) const
    { return BoundaryTypes(); }

    NumEqVector sourceAtPos(const GlobalPosition&) const
    { return NumEqVector(0.0); }

    PrimaryVariables initialAtPos(const GlobalPosition&) const
    { return PrimaryVariables(0.0); }

    const std::vector<DirichletConstraintData>& constraints() const
    { return constraints_; }

    Scalar exactSolution(const GlobalPosition& globalPos) const
    { return 1 + globalPos[2]; }

    const CouplingManager& couplingManager() const
    { return *couplingManager_; }

private:
    std::shared_ptr<CouplingManager> couplingManager_;
    std::vector<DirichletConstraintData> constraints_;
    PointSources pointSources_;
};

} // end namespace Dumux

#endif
