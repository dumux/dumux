// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \brief Momentum problem for the Cryer sphere consolidation benchmark.
 *
 * Geometry: 1/8 sphere (R = 1 m), symmetry planes at x=0, y=0, z=0.
 *
 * Boundary conditions:
 * - Symmetry planes (x=0, y=0, z=0): normal displacement = 0 (Dirichlet per component)
 * - Sphere surface: Neumann pressure load P = p_0 · n (inward)
 *
 * Uses the new experimental boundary condition interface:
 *   `boundaryTypes(FVElementGeometry, BoundaryFace)` returning `Experimental::BoundaryTypes<dim>`.
 * Dirichlet DOFs are built globally via `CVFE::appendDirichletConstraints`.
 *
 * The 1st Piola-Kirchhoff effective stress:
 *   P_eff = P_iso(F_bar) + (p_s - alpha_B * p_f) * J * F^{-T}
 * p_s and p_f are retrieved from the coupling manager.
 */
#ifndef DUMUX_CRYER_PROBLEM_MOMENTUM_HH
#define DUMUX_CRYER_PROBLEM_MOMENTUM_HH

#include <dumux/common/boundarytypes_.hh>
#include <dumux/common/fvproblemwithspatialparams.hh>
#include <dumux/common/numeqvector.hh>
#include <dumux/common/parameters.hh>
#include <dumux/common/properties.hh>
#include <dumux/common/constraintinfo.hh>
#include <dumux/common/indextraits.hh>
#include <dumux/discretization/dirichletconstraints.hh>
#include <dumux/discretization/cvfe/appenddirichletconstraints_.hh>
#include <dumux/discretization/cvfe/quadraturerules.hh>

namespace Dumux {

template<class TypeTag>
class CryerMomentumProblem : public FVProblemWithSpatialParams<TypeTag>
{
    using ParentType = FVProblemWithSpatialParams<TypeTag>;
    using GridGeometry = GetPropType<TypeTag, Properties::GridGeometry>;
    using FVElementGeometry = typename GridGeometry::LocalView;
    using SubControlVolumeFace = typename GridGeometry::SubControlVolumeFace;
    using BoundaryFace = typename GridGeometry::BoundaryFace;
    using GridView = typename GridGeometry::GridView;
    using Element = typename GridView::template Codim<0>::Entity;
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using PrimaryVariables = GetPropType<TypeTag, Properties::PrimaryVariables>;
    using NumEqVector = Dumux::NumEqVector<PrimaryVariables>;
    using BoundaryTypes = Dumux::Experimental::BoundaryTypes<GetPropType<TypeTag, Properties::ModelTraits>::numEq()>;
    using GlobalPosition = typename Element::Geometry::GlobalCoordinate;
    using GridIndexType = typename IndexTraits<GridView>::GridIndex;
    using CouplingManager = GetPropType<TypeTag, Properties::CouplingManager>;
    using Tensor = Dune::FieldMatrix<Scalar, GridView::dimension, GridView::dimension>;

    static constexpr int dim = GridView::dimension;
    static constexpr int numEq = GetPropType<TypeTag, Properties::ModelTraits>::numEq();
    using ConstraintInfo = Dumux::DirichletConstraintInfo<numEq>;
    using ConstraintValues = Dune::FieldVector<Scalar, numEq>;
    using DirichletConstraintData = Dumux::DirichletConstraintData<ConstraintInfo, ConstraintValues, GridIndexType>;

public:
    CryerMomentumProblem(std::shared_ptr<const GridGeometry> gridGeometry,
                          std::shared_ptr<CouplingManager> couplingManager)
    : ParentType(gridGeometry)
    , couplingManager_(couplingManager)
    , p0_(getParam<Scalar>("Problem.PressureLoad"))
    {
        CVFE::appendDirichletConstraints(*this,
            [&](const auto& fvGeometry, const auto& /*bf*/, const auto& localDof) {
                return PrimaryVariables(0.0);
            },
            constraints_);
    }

    const CouplingManager& couplingManager() const { return *couplingManager_; }

    //! New boundary types interface: FVElementGeometry + BoundaryFace
    BoundaryTypes boundaryTypes(const FVElementGeometry&, const BoundaryFace& bf) const
    { return boundaryTypesAtPos(bf.center()); }

    BoundaryTypes boundaryTypesAtPos(const GlobalPosition& pos) const
    {
        BoundaryTypes bct;
        // Default: all Dirichlet (isFluxBoundary = false for all)
        if (!isOnSymmetryPlane_(pos))
        {
            // Sphere surface: all Neumann (traction load)
            bct.setAllFluxBoundary();
        }
        else
        {
            // Symmetry plane: all Neumann initially, then constrain the normal direction
            bct.setAllFluxBoundary();
            for (int d = 0; d < dim; ++d)
                if (std::abs(pos[d]) < eps_)
                    bct.resetEq(d);   // normal component: NOT flux → Dirichlet (u_n = 0)
        }
        return bct;
    }

    PrimaryVariables dirichletAtPos(const GlobalPosition&) const
    { return PrimaryVariables(0.0); }

    //! Neumann traction per quadrature point (sphere surface: inward pressure)
    template<class ElementVariables, class IpData>
    NumEqVector boundaryFlux(const FVElementGeometry&, const ElementVariables&,
                             const IpData& ipd) const
    {
        NumEqVector flux(0.0);
        const auto& pos = ipd.global();
        const Scalar r = pos.two_norm();
        if (!isOnSymmetryPlane_(pos))
        {
            // sphere surface: apply pressure load along outward radial direction
            const Scalar rSafe = std::max(r, eps_);
            for (int d = 0; d < dim; ++d)
                flux[d] = p0_ * pos[d] / rSafe;
        }
        return flux;
    }

    template<class ElementVariables>
    NumEqVector boundaryFluxIntegral(const FVElementGeometry& fvGeometry,
                                      const ElementVariables& elemVars,
                                      const SubControlVolumeFace& scvf) const
    {
        NumEqVector flux(0.0);
        for (const auto& qpData : CVFE::quadratureRule(fvGeometry, scvf))
            flux += qpData.weight() * this->boundaryFlux(fvGeometry, elemVars, qpData.ipData());
        return flux;
    }

    const auto& constraints() const { return constraints_; }

    /*!
     * \brief Neo-Hookean effective 1st PK stress with Biot coupling:
     *   P_eff = P_iso(F_bar) + (p_s - alpha_B*p_f)*J*F^{-T}
     */
    Tensor firstPiolaKirchhoffStressTensor(const Tensor& F,
                                            const FVElementGeometry& fvGeometry,
                                            const GlobalPosition& ip) const
    {
        const Scalar G = this->spatialParams().shearModulus();
        const Scalar alphaB = this->spatialParams().biotCoefficient(
            fvGeometry.element(), fvGeometry.scv(0));

        const Scalar J = F.determinant();
        const Scalar Jm23 = std::pow(J, -2.0/3.0);

        // Right Cauchy-Green C = F^T F and its first invariant
        Tensor C(0.0);
        for (int i = 0; i < dim; ++i)
            for (int j = 0; j < dim; ++j)
                for (int k = 0; k < dim; ++k)
                    C[i][j] += F[k][i] * F[k][j];
        Scalar I1 = 0.0;
        for (int i = 0; i < dim; ++i) I1 += C[i][i];

        Tensor Finv(F); Finv.invert();

        // P_iso = G * J^{-2/3} * (F - I1/3 * F^{-T})
        Tensor P(0.0);
        for (int i = 0; i < dim; ++i)
            for (int j = 0; j < dim; ++j)
                P[i][j] = G * Jm23 * (F[i][j] - I1/3.0 * Finv[j][i]);

        // Effective pressure contribution
        const Scalar ps = couplingManager_->solidPressureAtPoint(fvGeometry, ip);
        const Scalar pf = couplingManager_->fluidPressureAtPoint(fvGeometry, ip);
        const Scalar pEff = ps - alphaB * pf;

        for (int i = 0; i < dim; ++i)
            for (int j = 0; j < dim; ++j)
                P[i][j] += pEff * J * Finv[j][i];

        return P;
    }

private:
    bool isOnSymmetryPlane_(const GlobalPosition& pos) const
    {
        for (int d = 0; d < dim; ++d)
            if (std::abs(pos[d]) < eps_) return true;
        return false;
    }

    static constexpr Scalar eps_ = 1e-8;
    std::shared_ptr<CouplingManager> couplingManager_;
    Scalar p0_;
    std::vector<DirichletConstraintData> constraints_;
};

} // end namespace Dumux
#endif
