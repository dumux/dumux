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
 * - Symmetry planes (x=0, y=0, z=0): zero normal displacement (constraint per component)
 * - Sphere surface: pressure load p_0 on the outer surface
 *
 * The problem supplies the isochoric Neo-Hookean part of the skeleton stress;
 * the model adds the pressure contribution (p_s - alpha_B p_f) J F^{-T}.
 */
#ifndef DUMUX_CRYER_PROBLEM_MOMENTUM_HH
#define DUMUX_CRYER_PROBLEM_MOMENTUM_HH

#include <dumux/common/boundarytypes_.hh>
#include <dumux/common/problemwithspatialparams.hh>
#include <dumux/common/numeqvector.hh>
#include <dumux/common/parameters.hh>
#include <dumux/common/properties.hh>
#include <dumux/common/constraintinfo.hh>
#include <dumux/common/indextraits.hh>
#include <dumux/discretization/dirichletconstraints.hh>
#include <dumux/discretization/cvfe/quadraturerules.hh>

namespace Dumux {

template<class TypeTag>
class CryerMomentumProblem : public Experimental::ProblemWithSpatialParams<TypeTag>
{
    using ParentType = Experimental::ProblemWithSpatialParams<TypeTag>;
    using GridGeometry = GetPropType<TypeTag, Properties::GridGeometry>;
    using ElementDiscretization = typename GridGeometry::LocalView;
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
    { appendDirichletConstraints_(); }

    const CouplingManager& couplingManager() const { return *couplingManager_; }

    BoundaryTypes boundaryTypesAtPos(const GlobalPosition& pos) const
    {
        BoundaryTypes bct;
        bct.setAllFluxBoundary();
        for (int d = 0; d < dim; ++d)
            if (std::abs(pos[d]) < eps_)
                bct.resetEq(d);
        return bct;
    }

    //! The pressure load on the sphere surface, \f$ \mathbf{g} = -\mathbf{P}\mathbf{N} = p_0 \mathbf{n} \f$
    template<class ElementVariables, class IpData>
    NumEqVector boundaryFlux(const ElementDiscretization&, const ElementVariables&,
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

    const auto& constraints() const { return constraints_; }

    /*!
     * \brief Isochoric part of the Neo-Hookean first Piola-Kirchhoff stress
     *   \f$ \mathbf{P}_\text{iso} = G J^{-2/3} (\mathbf{F} - \tfrac{1}{3} I_1 \mathbf{F}^{-\top}) \f$
     */
    Tensor isochoricFirstPiolaKirchhoffStressTensor(const Tensor& F,
                                                    const ElementDiscretization&,
                                                    const GlobalPosition&) const
    {
        const Scalar G = this->spatialParams().shearModulus();
        const Scalar Jm23 = std::pow(F.determinant(), -2.0/3.0);

        Scalar I1 = 0.0;
        for (int i = 0; i < dim; ++i)
            for (int k = 0; k < dim; ++k)
                I1 += F[k][i]*F[k][i];

        Tensor Finv(F);
        Finv.invert();

        Tensor P(0.0);
        for (int i = 0; i < dim; ++i)
            for (int j = 0; j < dim; ++j)
                P[i][j] = G*Jm23*(F[i][j] - I1/3.0*Finv[j][i]);

        return P;
    }

private:
    //! The displacement components without traction condition are zero
    void appendDirichletConstraints_()
    {
        auto elemDisc = localView(this->gridDiscretization());
        for (const auto& element : elements(this->gridDiscretization().gridView()))
        {
            elemDisc.bind(element);
            for (const auto& boundaryFace : boundaryFaces(elemDisc))
            {
                const auto bcTypes = this->boundaryTypes(elemDisc, boundaryFace);
                if (bcTypes.hasOnlyFluxBoundary())
                    continue;

                ConstraintInfo info;
                for (int eqIdx = 0; eqIdx < numEq; ++eqIdx)
                    if (!bcTypes.isFluxBoundary(eqIdx))
                        info.set(eqIdx);

                for (const auto& localDof : localDofs(elemDisc, boundaryFace))
                    constraints_.push_back(DirichletConstraintData{info, ConstraintValues(0.0), localDof.dofIndex()});
            }
        }
    }

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
