// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
#ifndef DUMUX_HYPERELASTIC_COOKS_MEMBRANE_PROBLEM_HH
#define DUMUX_HYPERELASTIC_COOKS_MEMBRANE_PROBLEM_HH

#include <dumux/common/boundarytypes.hh>
#include <dumux/common/fvproblemwithspatialparams.hh>
#include <dumux/common/numeqvector.hh>
#include <dumux/common/parameters.hh>
#include <dumux/common/properties.hh>
#include <dumux/common/math.hh>

namespace Dumux {

// Cook's membrane benchmark with large deformations using a hyperelastic material.
// Reference: R.D. Cook, "Improved two-dimensional finite element",
// ASCE J. Struct. Div. 100(9):1851-1863, 1974
//
// Geometry: tapered panel with vertices (0,0), (48,44), (48,60), (0,44)
// Boundary conditions:
//   - Left face (x=0): clamped (zero displacement)
//   - Right face (x=48): uniform shear traction in y-direction
//   - Top/bottom edges: traction-free
template<class TypeTag>
class HyperelasticCooksMembraneProblem : public FVProblemWithSpatialParams<TypeTag>
{
    using ParentType = FVProblemWithSpatialParams<TypeTag>;
    using GridGeometry = GetPropType<TypeTag, Properties::GridGeometry>;
    using FVElementGeometry = typename GridGeometry::LocalView;
    using SubControlVolumeFace = typename GridGeometry::SubControlVolumeFace;
    using GridView = typename GridGeometry::GridView;
    using Element = typename GridView::template Codim<0>::Entity;
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using PrimaryVariables = GetPropType<TypeTag, Properties::PrimaryVariables>;
    using NumEqVector = Dumux::NumEqVector<PrimaryVariables>;
    using BoundaryTypes = Dumux::BoundaryTypes<GetPropType<TypeTag, Properties::ModelTraits>::numEq()>;
    using GlobalPosition = typename Element::Geometry::GlobalCoordinate;
    static constexpr int dim = GridView::dimension;
    using Tensor = Dune::FieldMatrix<Scalar, dim, dim>;

public:
    HyperelasticCooksMembraneProblem(std::shared_ptr<const GridGeometry> gridGeometry)
    : ParentType(gridGeometry)
    , shearTraction_(getParam<Scalar>("Problem.ShearTraction"))
    {}

    BoundaryTypes boundaryTypesAtPos(const GlobalPosition& globalPos) const
    {
        BoundaryTypes values;
        if (globalPos[0] < eps_)
            values.setAllDirichlet();
        else
            values.setAllNeumann();
        return values;
    }

    PrimaryVariables dirichletAtPos(const GlobalPosition&) const
    { return PrimaryVariables(0.0); }

    template<class ElementVolumeVariables, class ElementFluxVariablesCache>
    NumEqVector neumann(const Element& element,
                        const FVElementGeometry& fvGeometry,
                        const ElementVolumeVariables& elemVolVars,
                        const ElementFluxVariablesCache& elemFluxVarsCache,
                        const SubControlVolumeFace& scvf) const
    {
        // In DuMux's hyperelastic model, flux = -P·N·area, so the Neumann return
        // value equals -T_applied (negative of the physical applied traction).
        // Upward shear traction T_y = +shearTraction_ → neumann_y = -shearTraction_.
        NumEqVector traction(0.0);
        if (scvf.ipGlobal()[0] > 48.0 - eps_)
            traction[1] = -shearTraction_;
        return traction;
    }

    // SPP-1748 benchmark material ψ₁ (Schröder et al., Eq. 2.1):
    // ψ₁ = µ/2*(I_C - 3) + λ/4*(J² - 1) - (λ/2 + µ)*ln J
    // P  = µ*(F - F^{-T}) + λ/2*(J² - 1)*F^{-T}          (Eq. 2.2 in Piola form)
    // Stress-free at F=I; linearises to σ = λ tr(ε)I + 2µε.
    Tensor firstPiolaKirchhoffStressTensor(Tensor F) const
    {
        const auto J = F.determinant();

        auto invFT = F;
        invFT.invert();
        invFT = transpose(invFT);

        const auto mu     = this->spatialParams().shearModulus();
        const auto lambda = this->spatialParams().firstLameParameter();

        // P = µ*(F - F^-T) + λ/2*(J²-1)*F^-T
        Tensor P = F;
        P.axpy(-1.0, invFT);
        P *= mu;
        P.axpy(0.5*lambda*(J*J - 1.0), invFT);
        return P;
    }

private:
    static constexpr Scalar eps_ = 1e-6;
    Scalar shearTraction_;
};

} // end namespace Dumux

#endif
