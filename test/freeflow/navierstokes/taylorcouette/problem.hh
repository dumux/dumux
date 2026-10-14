// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup NavierStokesTests
 * \brief Taylor-Couette flow between two concentric rotating cylinders
 */
#ifndef DUMUX_TEST_FREEFLOW_NAVIERSTOKES_TAYLORCOUETTE_PROBLEM_HH
#define DUMUX_TEST_FREEFLOW_NAVIERSTOKES_TAYLORCOUETTE_PROBLEM_HH

#include <cmath>
#include <vector>

#include <dumux/common/constraintinfo.hh>
#include <dumux/common/indextraits.hh>
#include <dumux/common/parameters.hh>
#include <dumux/common/properties.hh>
#include <dumux/discretization/cvfe/localdof.hh>
#include <dumux/discretization/dirichletconstraints.hh>

namespace Dumux {

/*!
 * \ingroup NavierStokesTests
 * \brief Taylor-Couette flow between two concentric rotating cylinders
 *
 * The inner cylinder of radius \f$ R_1 \f$ rotates with the angular velocity \f$ \Omega_1 \f$,
 * the outer cylinder of radius \f$ R_2 \f$ with \f$ \Omega_2 \f$. The circular Couette flow
 * \f$ u_\theta = A r + B/r \f$ with
 * \f$ A = (\Omega_2 R_2^2 - \Omega_1 R_1^2)/(R_2^2 - R_1^2) \f$ and
 * \f$ B = (\Omega_1 - \Omega_2) R_1^2 R_2^2/(R_2^2 - R_1^2) \f$
 * is a steady solution of the Navier-Stokes equations, and the pressure follows from
 * \f$ \partial p/\partial r = \rho u_\theta^2/r \f$ \cite taylor1923.
 * The velocity is prescribed on both cylinders.
 */
template <class TypeTag, class BaseProblem>
class TaylorCouetteProblem : public BaseProblem
{
    using ParentType = BaseProblem;
    using GridGeometry = GetPropType<TypeTag, Properties::GridGeometry>;
    using FVElementGeometry = typename GridGeometry::LocalView;
    using SubControlVolumeFace = typename FVElementGeometry::SubControlVolumeFace;
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using ModelTraits = GetPropType<TypeTag, Properties::ModelTraits>;
    using Element = typename GridGeometry::GridView::template Codim<0>::Entity;
    using GlobalPosition = typename Element::Geometry::GlobalCoordinate;
    using DirichletValues = typename ParentType::DirichletValues;
    using BoundaryFluxes = typename ParentType::BoundaryFluxes;
    using BoundaryTypes = typename ParentType::BoundaryTypes;
    using CouplingManager = GetPropType<TypeTag, Properties::CouplingManager>;
    using ConstraintInfo = DirichletConstraintInfo<ModelTraits::numEq()>;
    using GridIndexType = typename IndexTraits<typename GridGeometry::GridView>::GridIndex;
    using Constraint = DirichletConstraintData<ConstraintInfo, DirichletValues, GridIndexType>;

public:
    using Indices = typename ModelTraits::Indices;

    TaylorCouetteProblem(std::shared_ptr<const GridGeometry> gridGeometry, std::shared_ptr<CouplingManager> couplingManager)
    : ParentType(gridGeometry, couplingManager)
    {
        const auto radii = getParam<std::vector<Scalar>>("Grid.Radial0");
        const auto r1 = radii.front();
        const auto r2 = radii.back();
        const auto omega1 = getParam<Scalar>("Problem.Omega1");
        const auto omega2 = getParam<Scalar>("Problem.Omega2");
        a_ = (omega2*r2*r2 - omega1*r1*r1)/(r2*r2 - r1*r1);
        b_ = (omega1 - omega2)*r1*r1*r2*r2/(r2*r2 - r1*r1);
        density_ = getParam<Scalar>("Component.LiquidDensity");

        if constexpr (ParentType::isMomentumProblem())
            setVelocityOnCylinders_(*gridGeometry);
        else
            fixPressureInOneCell_(*gridGeometry);
    }

    BoundaryTypes boundaryTypesAtPos(const GlobalPosition& globalPos) const
    {
        BoundaryTypes values;
        if constexpr (!ParentType::isMomentumProblem())
            values.setAllNeumann();
        return values;
    }

    template<class ElementVolumeVariables, class ElementFluxVariablesCache>
    BoundaryFluxes neumann(const Element& element,
                           const FVElementGeometry& fvGeometry,
                           const ElementVolumeVariables& elemVolVars,
                           const ElementFluxVariablesCache& elemFluxVarsCache,
                           const SubControlVolumeFace& scvf) const
    {
        BoundaryFluxes values(0.0);

        if constexpr (!ParentType::isMomentumProblem())
            values[Indices::conti0EqIdx] = this->faceVelocity(element, fvGeometry, scvf)
                                           *density_*scvf.unitOuterNormal();

        return values;
    }

    const std::vector<Constraint>& constraints() const
    { return constraints_; }

    DirichletValues analyticalSolution(const GlobalPosition& globalPos, Scalar time = 0.0) const
    {
        DirichletValues values;
        const auto r2 = globalPos.two_norm2();

        if constexpr (ParentType::isMomentumProblem())
        {
            const auto angularVelocity = a_ + b_/r2;
            values[Indices::velocityXIdx] = -globalPos[1]*angularVelocity;
            values[Indices::velocityYIdx] = globalPos[0]*angularVelocity;
        }
        else
        {
            using std::log;
            values[Indices::pressureIdx] = density_*(0.5*a_*a_*r2 + a_*b_*log(r2) - 0.5*b_*b_/r2);
        }

        return values;
    }

private:
    void setVelocityOnCylinders_(const GridGeometry& gridGeometry)
    {
        auto fvGeometry = localView(gridGeometry);
        for (const auto& element : elements(gridGeometry.gridView()))
        {
            fvGeometry.bind(element);
            for (const auto& boundaryFace : boundaryFaces(fvGeometry))
            {
                for (const auto& localDof : localDofs(fvGeometry, boundaryFace))
                {
                    ConstraintInfo info;
                    info.setAll();
                    constraints_.push_back(Constraint{std::move(info),
                                                      analyticalSolution(ipData(fvGeometry, localDof).global()),
                                                      localDof.dofIndex()});
                }
            }
        }
    }

    // the pressure is only defined up to a constant
    void fixPressureInOneCell_(const GridGeometry& gridGeometry)
    {
        auto fvGeometry = localView(gridGeometry);
        fvGeometry.bindElement(*elements(gridGeometry.gridView()).begin());
        for (const auto& scv : scvs(fvGeometry))
        {
            ConstraintInfo info;
            info.setAll();
            constraints_.push_back(Constraint{std::move(info), analyticalSolution(scv.dofPosition()), scv.dofIndex()});
        }
    }

    Scalar a_, b_;
    Scalar density_;
    std::vector<Constraint> constraints_;
};

} // end namespace Dumux

#endif
