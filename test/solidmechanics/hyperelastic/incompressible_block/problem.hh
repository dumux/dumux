// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
#ifndef DUMUX_HYPERELASTIC_INCOMPRESSIBLE_BLOCK_PROBLEM_HH
#define DUMUX_HYPERELASTIC_INCOMPRESSIBLE_BLOCK_PROBLEM_HH

#include <dumux/common/boundarytypes.hh>
#include <dumux/common/fvproblemwithspatialparams.hh>
#include <dumux/common/numeqvector.hh>
#include <dumux/common/parameters.hh>
#include <dumux/common/properties.hh>
#include <dumux/common/math.hh>

namespace Dumux {

// SPP-1748 benchmark: quasi-incompressible cube under constant partial load
// (Kollmannsberger, Bog, D'Angella, Rank, Wriggers — Section 3 of the PDF)
//
// Geometry: h=50, w=100, l=100 mm; quarter model exploiting symmetry
//   Domain (quarter): [0, w/2] x [0, l/2] x [0, h] = [0,50]x[0,50]x[0,50]
//
// BCs (quarter model):
//   z=0 : fixed in z (bottom)
//   z=h : fixed in x and y (top surface sliding constraint per PDF)
//   x=w/2: symmetry → fixed in x
//   y=l/2: symmetry → fixed in y
//   Load q on area a×b=[25×25] at z=h (x∈[0,a], y∈[0,b])
//
// Reference: |u_z| at point P (center of load, x=y=0, z=h) ≈ 20 mm (Table 3.3)
template<class TypeTag>
class HyperelasticIncompressibleBlockProblem : public FVProblemWithSpatialParams<TypeTag>
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

    // quarter model geometry
    static constexpr Scalar h_ = 50.0;   // height (z)
    static constexpr Scalar w2_ = 50.0;  // half-width (x) = w/2
    static constexpr Scalar l2_ = 50.0;  // half-length (y) = l/2
    static constexpr Scalar a_ = 25.0;   // load patch x-extent
    static constexpr Scalar b_ = 25.0;   // load patch y-extent

public:
    HyperelasticIncompressibleBlockProblem(std::shared_ptr<const GridGeometry> gridGeometry)
    : ParentType(gridGeometry)
    , load_(getParam<Scalar>("Problem.Load"))
    , loadFactor_(1.0)
    {}

    // Scale the applied load — used for pseudo-static load stepping in main.cc
    void setLoadFactor(Scalar f) { loadFactor_ = f; }

    BoundaryTypes boundaryTypesAtPos(const GlobalPosition& globalPos) const
    {
        BoundaryTypes values;
        values.setAllNeumann(); // default: traction-free

        // z=0: fixed in z-direction only (bottom)
        if (globalPos[2] < eps_)
        {
            values.setDirichlet(2); // u_z = 0
        }
        // z=h: fixed in x and y (top surface — sliding per PDF Section 3.1)
        else if (globalPos[2] > h_ - eps_)
        {
            values.setDirichlet(0); // u_x = 0
            values.setDirichlet(1); // u_y = 0
        }
        // x=w/2: symmetry plane, fixed in x
        else if (globalPos[0] > w2_ - eps_)
        {
            values.setDirichlet(0); // u_x = 0
        }
        // y=l/2: symmetry plane, fixed in y
        else if (globalPos[1] > l2_ - eps_)
        {
            values.setDirichlet(1); // u_y = 0
        }

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
        // Downward compressive load q on the loaded patch at z=h.
        // neumann = -T_applied; T_applied_z = -q (down) → neumann_z = +q
        const auto& pos = scvf.ipGlobal();
        NumEqVector traction(0.0);
        // a and b are the half-widths of the full (symmetric) load patch centred at (w/2, l/2).
        // Quarter patch: x∈[w/2-a, w/2] × y∈[l/2-b, l/2] = [25,50]×[25,50]
        if (pos[2] > h_ - eps_
            && pos[0] > w2_ - a_ - eps_
            && pos[1] > l2_ - b_ - eps_)
            traction[2] = loadFactor_ * load_;
        return traction;
    }

    // SPP-1748 material (Kollmannsberger et al., Eq. 3.2):
    // W = µ/2*(tr C - 3 - 2 ln J) + λ/2*(J-1)²
    // P = µ*(F - F^{-T}) + λ*(J-1)*J*F^{-T}
    Tensor firstPiolaKirchhoffStressTensor(Tensor F) const
    {
        const auto J = F.determinant();

        auto invFT = F;
        invFT.invert();
        invFT = transpose(invFT);

        const auto mu     = this->spatialParams().shearModulus();
        const auto lambda = this->spatialParams().firstLameParameter();

        // P = µ*(F - F^-T) + λ*(J-1)*J*F^-T
        Tensor P = F;
        P.axpy(-1.0, invFT);
        P *= mu;
        P.axpy(lambda*(J - 1.0)*J, invFT);
        return P;
    }

private:
    static constexpr Scalar eps_ = 1e-6;
    Scalar load_, loadFactor_;
};

} // end namespace Dumux

#endif
