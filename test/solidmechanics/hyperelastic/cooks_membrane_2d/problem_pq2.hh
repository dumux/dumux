// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
#ifndef DUMUX_HYPERELASTIC_COOKS_MEMBRANE_PQ2_PROBLEM_HH
#define DUMUX_HYPERELASTIC_COOKS_MEMBRANE_PQ2_PROBLEM_HH

#include <vector>

#include <dune/common/fmatrix.hh>

#include <dumux/common/boundarytypes.hh>
#include <dumux/common/constraintinfo.hh>
#include <dumux/common/indextraits.hh>
#include <dumux/common/numeqvector.hh>
#include <dumux/common/parameters.hh>
#include <dumux/common/properties.hh>
#include <dumux/common/math.hh>
#include <dumux/common/problemwithspatialparams.hh>
#include <dumux/discretization/cvfe/localdof.hh>
#include <dumux/discretization/cvfe/quadraturerules.hh>
#include <dumux/discretization/dirichletconstraints.hh>

namespace Dumux {

/*!
 * \brief 2D Cook's membrane — SPP-1748 benchmark with hyperelastic material ψ₁.
 *        Uses the current Experimental::Assembler CVFE interface.
 *
 * Geometry: tapered panel (0,0)–(48,44)–(48,60)–(0,44) mm.
 * BCs:
 *   x=0:  clamped (u=0), applied via per-dof Dirichlet constraints
 *   x=48: shear traction T=(0, shearTraction)
 *   rest: traction-free
 *
 * Material ψ₁ (SPP-1748 Eq. 2.1):
 *   P = µ(F − F⁻ᵀ) + λ/2·(J²−1)·F⁻ᵀ
 */
template<class TypeTag>
class HyperelasticCooksMembraneP2Problem : public Dumux::Experimental::ProblemWithSpatialParams<TypeTag>
{
    using ParentType = Dumux::Experimental::ProblemWithSpatialParams<TypeTag>;
    using GridGeometry = GetPropType<TypeTag, Properties::GridGeometry>;
    using FVElementGeometry = typename GridGeometry::LocalView;
    using GridView = typename GridGeometry::GridView;
    using Element = typename GridView::template Codim<0>::Entity;
    using GridIndexType = typename IndexTraits<GridView>::GridIndex;
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using PrimaryVariables = GetPropType<TypeTag, Properties::PrimaryVariables>;
    using NumEqVector = Dumux::NumEqVector<PrimaryVariables>;
    using GlobalPosition = typename Element::Geometry::GlobalCoordinate;
    static constexpr int dim = GridView::dimension;
    static constexpr int numEq = GetPropType<TypeTag, Properties::ModelTraits>::numEq();
    using Tensor = Dune::FieldMatrix<Scalar, dim, dim>;
    using BoundaryTypes = Dumux::Experimental::BoundaryTypes<numEq>;
    using ConstraintInfo = Dumux::DirichletConstraintInfo<numEq>;
    using ConstraintValues = Dune::FieldVector<Scalar, numEq>;
    using DirichletConstraintData = Dumux::DirichletConstraintData<ConstraintInfo, ConstraintValues, GridIndexType>;

public:
    HyperelasticCooksMembraneP2Problem(std::shared_ptr<const GridGeometry> gridGeometry)
    : ParentType(gridGeometry)
    , shearTraction_(getParam<Scalar>("Problem.ShearTraction"))
    {
        appendDirichletConstraints_();
    }

    //! Flux (traction) boundaries: all faces except the clamped x=0 face.
    BoundaryTypes boundaryTypesAtPos(const GlobalPosition& globalPos) const
    {
        BoundaryTypes values{};
        if (globalPos[0] >= eps_)
            values.setAllFluxBoundary();
        return values;
    }

    //! Per-QP traction: shear T=(0,-shearTraction) on the x=48 face (box flux sign).
    template<class EV, class IpData>
    NumEqVector boundaryFlux(const FVElementGeometry&, const EV&, const IpData& ipData) const
    {
        NumEqVector traction(0.0);
        if (ipData.global()[0] > 48.0 - eps_)
            traction[1] = -shearTraction_;
        return traction;
    }

    const auto& constraints() const { return constraints_; }

    // SPP-1748 ψ₁: P = µ(F − F⁻ᵀ) + λ/2·(J²−1)·F⁻ᵀ
    Tensor firstPiolaKirchhoffStressTensor(Tensor F) const
    {
        const auto J = F.determinant();
        auto invFT = F; invFT.invert(); invFT = transpose(invFT);

        const auto mu     = this->spatialParams().shearModulus();
        const auto lambda = this->spatialParams().firstLameParameter();

        Tensor P = F;
        P.axpy(-1.0, invFT);
        P *= mu;
        P.axpy(0.5*lambda*(J*J - 1.0), invFT);
        return P;
    }

private:
    void appendDirichletConstraints_()
    {
        auto fvGeometry = localView(this->gridGeometry());
        for (const auto& element : elements(this->gridGeometry().gridView()))
        {
            fvGeometry.bind(element);
            for (const auto& boundaryFace : boundaryFaces(fvGeometry))
            {
                if (boundaryFace.center()[0] >= eps_) // only the clamped x=0 face
                    continue;
                for (const auto& localDof : localDofs(fvGeometry, boundaryFace))
                {
                    ConstraintInfo info;
                    info.setAll(); // clamp all components
                    ConstraintValues values(0.0);
                    constraints_.push_back(DirichletConstraintData{std::move(info), std::move(values), localDof.dofIndex()});
                }
            }
        }
    }

    static constexpr Scalar eps_ = 1e-6;
    Scalar shearTraction_;
    std::vector<DirichletConstraintData> constraints_;
};

} // end namespace Dumux
#endif
