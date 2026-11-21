// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup EmbeddedTests
 * \brief The tissue problem of the embedded 1d-3d test stated with the variables of the local dofs
 *
 * The Dirichlet values are imposed as constraints.
 */
#ifndef DUMUX_EMBEDDED_1P1P_VARIABLES_PROBLEM_TISSUE_HH
#define DUMUX_EMBEDDED_1P1P_VARIABLES_PROBLEM_TISSUE_HH

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
class TissueProblem : public Experimental::ProblemWithSpatialParams<TypeTag>
{
    using ParentType = Experimental::ProblemWithSpatialParams<TypeTag>;
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using GridDiscretization = GetPropType<TypeTag, Properties::GridGeometry>;
    using ElementDiscretization = typename GridDiscretization::LocalView;
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

    //! The exchange term of a point source, the kernel method scales it with the flux scaling factor
    struct Exchange
    {
        const CouplingManager* couplingManager;
        Scalar beta;

        template<class ElementDiscretizationType, class ElementVariables, class Context>
        NumEqVector operator()(const ElementDiscretizationType&, const ElementVariables&, const Context& context) const
        {
            const auto id = context.id();
            const Scalar pressure3D = couplingManager->bulkPriVars(id)[0];
            const Scalar pressure1D = couplingManager->lowDimPriVars(id)[0];
            Scalar value = beta*(pressure1D - pressure3D)*context.weight();
            if constexpr (CouplingManager::couplingMode == Embedded1d3dCouplingMode::kernel)
                value *= couplingManager->fluxScalingFactor(id);
            return NumEqVector(value);
        }
    };

    using LocalPosition = typename GridView::template Codim<0>::Entity::Geometry::LocalCoordinate;
    using PointSources = Experimental::PointSources<EmbeddedCoupling::PointSourceContext<GlobalPosition>, Exchange, LocalPosition>;

public:
    TissueProblem(std::shared_ptr<const GridDiscretization> gridDiscretization,
                  std::shared_ptr<CouplingManager> couplingManager)
    : ParentType(gridDiscretization, "Tissue")
    , couplingManager_(couplingManager)
    , dirichletOnly_(getParam<bool>("Problem.DirichletOnly", false))
    {
        appendBoundaryDirichletConstraints(
            constraints_, this->gridDiscretization(),
            [&] (const auto& elemDisc, const auto& boundaryFace, const auto& localDof)
            {
                ConstraintInfo info;
                if (dirichletOnly_ || isLateralPosition_(ipData(elemDisc, localDof).global()))
                    info.setAll();
                return info;
            },
            [&] (const auto& elemDisc, const auto& localDof)
            { return ConstraintValues(exactSolution(ipData(elemDisc, localDof).global())); }
        );
    }

    //! Set the point sources after the coupling manager has computed them
    void setPointSources()
    {
        static const auto radius = getParam<Scalar>("SpatialParams.Radius");
        pointSources_ = couplingManager_->makePointSources(
            Dune::index_constant<0>{}, Exchange{couplingManager_.get(), 2*M_PI/(2*M_PI + std::log(radius))}
        );
    }

    const PointSources& pointSources() const
    { return pointSources_; }

    BoundaryTypes boundaryTypes(const ElementDiscretization&, const typename ElementDiscretization::BoundaryFace& boundaryFace) const
    {
        BoundaryTypes values;
        if (!dirichletOnly_ && isFluxPosition_(boundaryFace.center()))
            values.setAllFluxBoundary();
        return values;
    }

    //! The outflow through the top and bottom face per area, the derivative in z-direction is the exact solution evaluated with z=0
    template<class ElementVariables, class FaceIpData>
    NumEqVector boundaryFlux(const ElementDiscretization&, const ElementVariables&, const FaceIpData& faceIpData) const
    {
        auto pos = faceIpData.global();
        const bool top = pos[2] > this->gridDiscretization().bBoxMax()[2] - eps_;
        pos[2] = 0.0;
        const auto derivative = exactSolution(pos);
        return NumEqVector(top ? -derivative : derivative);
    }

    NumEqVector sourceAtPos(const GlobalPosition&) const
    { return NumEqVector(0.0); }

    PrimaryVariables initialAtPos(const GlobalPosition&) const
    { return PrimaryVariables(0.0); }

    const std::vector<DirichletConstraintData>& constraints() const
    { return constraints_; }

    Scalar p1DExact(const GlobalPosition& globalPos) const
    { return 1.0 + globalPos[2]; }

    Scalar exactSolution(const GlobalPosition& globalPos) const
    {
        const auto r = std::hypot(globalPos[0], globalPos[1]);
        static const auto R = getParam<Scalar>("SpatialParams.Radius");

        if constexpr (CouplingManager::couplingMode == Embedded1d3dCouplingMode::kernel)
        {
            static const auto rho = getParam<Scalar>("MixedDimension.KernelWidthFactor")*R;
            return -1.0*p1DExact(globalPos)/(2*M_PI)*radialProfile_(r, rho);
        }
        else if constexpr (CouplingManager::couplingMode == Embedded1d3dCouplingMode::surface)
        {
            if (r > R)
                return -1.0*p1DExact(globalPos)/(2*M_PI)*std::log(r);
            else
                return -1.0*p1DExact(globalPos)/(2*M_PI)*std::log(R);
        }
        else
            return -1.0*p1DExact(globalPos)/(2*M_PI)*std::log(r);
    }

    //! The flux scaling factor (xi) for a distance r for a vessel with radius R and kernel width rho
    Scalar fluxScalingFactor(const Scalar r, const Scalar R, const Scalar rho) const
    {
        using std::log;
        static const Scalar beta = 2*M_PI/(2*M_PI + log(R));
        return 1.0/(1.0 + beta/(2*M_PI)*(radialProfile_(r, rho) - log(R)));
    }

    const CouplingManager& couplingManager() const
    { return *couplingManager_; }

private:
    /*!
     * \brief The radial profile g of the pressure -q/(2 pi) g(r) around a line source of strength q
     *        distributed with a kernel of width rho
     *
     * Outside of the kernel, g = ln(r). Inside, with u = r^2/rho^2, the uniform kernel gives
     * g = ln(rho) - (1 - u)/2 and the kernel 3/(pi rho^2) (1 - u)^2 (the second-moment kernel
     * of the coupling manager with second moment rho^2/4) gives
     * g = ln(rho) - (3(1 - u) - 3(1 - u^2)/2 + (1 - u^3)/3)/2.
     */
    Scalar radialProfile_(const Scalar r, const Scalar rho) const
    {
        using std::log;
        if (r > rho)
            return log(r);

        const auto u = r*r/(rho*rho);
        static const bool smoothKernel = hasParam("MixedDimension.KernelSecondMoment");
        if (!smoothKernel)
            return log(rho) - 0.5*(1.0 - u);

        static const auto secondMoment = getParam<Scalar>("MixedDimension.KernelSecondMoment");
        if (std::abs(secondMoment - 0.25*rho*rho) > 1e-12*rho*rho)
            DUNE_THROW(Dune::NotImplemented, "Exact solution for a kernel second moment other than rho^2/4");

        return log(rho) - 0.5*(3.0*(1.0 - u) - 1.5*(1.0 - u*u) + (1.0 - u*u*u)/3.0);
    }

    bool isFluxPosition_(const GlobalPosition& globalPos) const
    {
        return globalPos[2] > this->gridDiscretization().bBoxMax()[2] - eps_
               || globalPos[2] < this->gridDiscretization().bBoxMin()[2] + eps_;
    }

    bool isLateralPosition_(const GlobalPosition& globalPos) const
    {
        const auto& bBoxMin = this->gridDiscretization().bBoxMin();
        const auto& bBoxMax = this->gridDiscretization().bBoxMax();
        return globalPos[0] < bBoxMin[0] + eps_ || globalPos[0] > bBoxMax[0] - eps_
               || globalPos[1] < bBoxMin[1] + eps_ || globalPos[1] > bBoxMax[1] - eps_;
    }

    static constexpr Scalar eps_ = 1.5e-7;
    std::shared_ptr<CouplingManager> couplingManager_;
    bool dirichletOnly_;
    std::vector<DirichletConstraintData> constraints_;
    PointSources pointSources_;
};

} // end namespace Dumux

#endif
