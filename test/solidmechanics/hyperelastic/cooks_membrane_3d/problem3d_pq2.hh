// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
#ifndef DUMUX_HYPERELASTIC_COOKS_MEMBRANE_3D_PQ2_PROBLEM_HH
#define DUMUX_HYPERELASTIC_COOKS_MEMBRANE_3D_PQ2_PROBLEM_HH

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
 * \brief 3D Cook's membrane SPP-1748 benchmark (PQ2 single-field, current
 *        Experimental::Assembler interface).
 *
 * Cross-section (0,0)-(48,44)-(48,60)-(0,44) mm extruded by depth=1 mm in z.
 *
 * BCs:
 *   x=0 face:  clamped (u=0 in all directions)
 *   x=48 face: shear traction T=(0, P₀, 0)  (P₀ = 20 MPa)
 *   z=0 face:  u_z=0  (plane-strain symmetry → results match 2D and SPP-1748 3D p-FEM)
 *   all other: traction-free
 *
 * Measurement: tip displacement u_y at (48, 60, 0).
 * Reference:   SPP-1748 Table 2.9 (3D p-FEM values).
 *
 * Boundary conditions use the current CVFE interface: flux (traction) boundaries
 * via boundaryTypesAtPos()/boundaryFlux(), Dirichlet boundaries via per-dof
 * constraints assembled once in the constructor (constraints()).
 */
template<class TypeTag>
class CooksMembrane3DProblemPQ2 : public Dumux::Experimental::ProblemWithSpatialParams<TypeTag>
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
    CooksMembrane3DProblemPQ2(std::shared_ptr<const GridGeometry> gridGeometry)
    : ParentType(gridGeometry)
    , shearTraction_(getParam<Scalar>("Problem.ShearTraction"))
    {
        appendDirichletConstraints_();
    }

    //! Flux (traction) boundaries: the loaded/free faces. x=0 is fully Dirichlet
    //! (clamped) and the z=0 face is Dirichlet in u_z only (plane strain), both via
    //! constraints; so x=0 carries no flux and z=0 carries flux on u_x,u_y only.
    BoundaryTypes boundaryTypesAtPos(const GlobalPosition& globalPos) const
    {
        BoundaryTypes values{};
        if (globalPos[0] < eps_)
            return values; // clamped: all Dirichlet, no flux
        if (globalPos[2] < eps_)
        {
            values.setFluxBoundary(0); // plane strain: u_x,u_y traction free
            values.setFluxBoundary(1);
            return values;
        }
        values.setAllFluxBoundary();
        return values;
    }

    //! Per-QP boundary traction: shear T=(0,-P₀,0) on the x=48 face, traction free
    //! elsewhere.
    template<class ElementVariables, class IpData>
    NumEqVector boundaryFlux(const FVElementGeometry&, const ElementVariables&, const IpData& ipData) const
    {
        NumEqVector traction(0.0);
        if (ipData.global()[0] > 48.0 - eps_)
            traction[1] = -shearTraction_;
        return traction;
    }

    const auto& constraints() const { return constraints_; }

    // SPP-1748 benchmark material ψ₁:
    //   P = µ(F - F⁻ᵀ) + λ/2·(J²-1)·F⁻ᵀ
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
    //! per-component Dirichlet map (old BoundaryTypes) used to build constraints:
    //! x=0 fully clamped; z=0 fixes only u_z (plane strain).
    Dumux::BoundaryTypes<numEq> dirichletMapAtPos_(const GlobalPosition& globalPos) const
    {
        Dumux::BoundaryTypes<numEq> values;
        values.setAllNeumann();
        if (globalPos[0] < eps_)
        { values.setDirichlet(0); values.setDirichlet(1); values.setDirichlet(2); }
        else if (globalPos[2] < eps_)
            values.setDirichlet(2);
        return values;
    }

    void appendDirichletConstraints_()
    {
        auto fvGeometry = localView(this->gridGeometry());
        for (const auto& element : elements(this->gridGeometry().gridView()))
        {
            fvGeometry.bind(element);
            for (const auto& boundaryFace : boundaryFaces(fvGeometry))
            {
                const auto bcTypes = dirichletMapAtPos_(boundaryFace.center());
                if (!bcTypes.hasDirichlet())
                    continue;
                for (const auto& localDof : localDofs(fvGeometry, boundaryFace))
                {
                    ConstraintInfo info;
                    for (int eqIdx = 0; eqIdx < numEq; ++eqIdx)
                        if (bcTypes.isDirichlet(eqIdx))
                            info.set(bcTypes.eqToDirichletIndex(eqIdx), eqIdx);
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
