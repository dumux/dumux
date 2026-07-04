// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
#ifndef DUMUX_HYPERELASTIC_COOKS_MEMBRANE_3D_PROBLEM_HH
#define DUMUX_HYPERELASTIC_COOKS_MEMBRANE_3D_PROBLEM_HH

#include <dune/common/fmatrix.hh>

#include <dumux/common/boundarytypes.hh>
#include <dumux/common/fvproblemwithspatialparams.hh>
#include <dumux/common/numeqvector.hh>
#include <dumux/common/parameters.hh>
#include <dumux/common/properties.hh>
#include <dumux/common/math.hh>

namespace Dumux {

/*!
 * \brief 3D Cook's membrane SPP-1748 benchmark (box / P1 displacement, classic
 *        FVAssembler interface).
 *
 * Cross-section (0,0)-(48,44)-(48,60)-(0,44) mm extruded by depth=1 mm in z.
 *
 * BCs:
 *   x=0 face:  clamped (u=0 in all directions)
 *   x=48 face: shear traction T=(0, P₀, 0)  (P₀ = 20 MPa)
 *   z=0 face:  u_z=0  (plane-strain symmetry)
 *   all other: traction-free
 */
template<class TypeTag>
class CooksMembrane3DProblem : public FVProblemWithSpatialParams<TypeTag>
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
    CooksMembrane3DProblem(std::shared_ptr<const GridGeometry> gridGeometry)
    : ParentType(gridGeometry)
    , shearTraction_(getParam<Scalar>("Problem.ShearTraction"))
    {}

    BoundaryTypes boundaryTypesAtPos(const GlobalPosition& globalPos) const
    {
        BoundaryTypes values;
        if (globalPos[0] < eps_)
            values.setAllDirichlet(); // x=0 clamped
        else
        {
            values.setAllNeumann();
            if (globalPos[2] < eps_)
                values.setDirichlet(2); // z=0 plane strain: u_z = 0
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
        NumEqVector traction(0.0);
        if (scvf.ipGlobal()[0] > 48.0 - eps_)
            traction[1] = -shearTraction_;
        return traction;
    }

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
    static constexpr Scalar eps_ = 1e-6;
    Scalar shearTraction_;
};

} // end namespace Dumux
#endif
