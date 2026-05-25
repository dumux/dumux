// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
#ifndef DUMUX_HYPERELASTIC_COOKS_MEMBRANE_3D_PROBLEM_HH
#define DUMUX_HYPERELASTIC_COOKS_MEMBRANE_3D_PROBLEM_HH

#include <dumux/common/boundarytypes.hh>
#include <dumux/common/constraintinfo.hh>
#include <dumux/common/indextraits.hh>
#include <dumux/common/fvproblemwithspatialparams.hh>
#include <dumux/common/numeqvector.hh>
#include <dumux/common/parameters.hh>
#include <dumux/common/properties.hh>
#include <dumux/common/math.hh>
#include <dumux/discretization/cvfe/quadraturerules.hh>
#include <dumux/discretization/dirichletconstraints.hh>

namespace Dumux {

/*!
 * \brief 3D Cook's membrane SPP-1748 benchmark.
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
    using Intersection = typename GridView::Intersection;
    using GridIndexType = typename IndexTraits<GridView>::GridIndex;
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using PrimaryVariables = GetPropType<TypeTag, Properties::PrimaryVariables>;
    using NumEqVector = Dumux::NumEqVector<PrimaryVariables>;
    using BoundaryTypes = Dumux::BoundaryTypes<GetPropType<TypeTag, Properties::ModelTraits>::numEq()>;
    using GlobalPosition = typename Element::Geometry::GlobalCoordinate;
    static constexpr int dim = GridView::dimension;
    static constexpr int numEq = GetPropType<TypeTag, Properties::ModelTraits>::numEq();
    using Tensor = Dune::FieldMatrix<Scalar, dim, dim>;
    using ConstraintInfo = Dumux::DirichletConstraintInfo<numEq>;
    using ConstraintValues = Dune::FieldVector<Scalar, numEq>;
    using DirichletConstraintData = Dumux::DirichletConstraintData<ConstraintInfo, ConstraintValues, GridIndexType>;

public:
    CooksMembrane3DProblem(std::shared_ptr<const GridGeometry> gridGeometry)
    : ParentType(gridGeometry)
    , shearTraction_(getParam<Scalar>("Problem.ShearTraction"))
    {
        CVFE::appendDirichletConstraints(*this,
            [&](const auto& fvGeometry, const auto&, const auto& localDof) {
                return this->dirichletAtPos(ipData(fvGeometry, localDof).global());
            },
            constraints_);
    }

    // New intersection-based BC type (for appendDirichletConstraints + Experimental::Assembler)
    BoundaryTypes boundaryTypes(const FVElementGeometry&, const Intersection& is) const
    { return boundaryTypesAtPos(is.geometry().center()); }

    //! Dirichlet constraints applied globally by Experimental::Assembler.
    const auto& constraints() const { return constraints_; }

    BoundaryTypes boundaryTypesAtPos(const GlobalPosition& globalPos) const
    {
        BoundaryTypes values;
        if (globalPos[0] < eps_ || globalPos[2] < eps_)
            values.setAllDirichlet();
        else
            values.setAllNeumann();

        // Only z-component is Dirichlet on z=0 face (plane-strain symmetry)
        if (globalPos[2] < eps_ && globalPos[0] > eps_)
        {
            values.setAllNeumann();
            values.setDirichlet(2); // u_z = 0
        }
        return values;
    }

    PrimaryVariables dirichletAtPos(const GlobalPosition& globalPos) const
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

    //! Per-QP traction for the new Experimental interface.
    template<class EV, class IpData>
    NumEqVector boundaryFlux(const FVElementGeometry&, const EV&, const IpData& ipData) const
    {
        NumEqVector traction(0.0);
        if (ipData.global()[0] > 48.0 - eps_)
            traction[1] = -shearTraction_;
        return traction;
    }

    //! Integrated Neumann traction over boundary scvf (new Experimental interface).
    template<class EV>
    NumEqVector boundaryFluxIntegral(const FVElementGeometry& fvGeometry,
                                     const EV& elemVars,
                                     const SubControlVolumeFace& scvf) const
    {
        NumEqVector flux(0.0);
        for (const auto& qpData : CVFE::quadratureRule(fvGeometry, scvf))
            flux += qpData.weight() * this->boundaryFlux(fvGeometry, elemVars, qpData.ipData());
        return flux;
    }

    //! Boundary flux integral over intersection (for FE edge-midpoint contributions).
    template<class EV, class Intersection, class BoundaryTypes>
    void addBoundaryFluxIntegrals(auto& residual,
                                  const FVElementGeometry& fvGeometry,
                                  const EV& elemVars,
                                  const Intersection& intersection,
                                  const BoundaryTypes& bcTypes) const
    {
        for (const auto& qpData : CVFE::quadratureRule(fvGeometry, intersection))
        {
            const auto flux = qpData.weight() * this->boundaryFlux(fvGeometry, elemVars, qpData.ipData());
            const auto& ipCache = cache(elemVars, qpData.ipData());
            const auto& shapeValues = ipCache.shapeValues();
            for (const auto& nonCVdof : nonCVLocalDofs(fvGeometry))
            {
                const auto idx = nonCVdof.index();
                for (int eqIdx = 0; eqIdx < NumEqVector::dimension; ++eqIdx)
                    if (bcTypes.isNeumann(eqIdx))
                        residual[idx][eqIdx] += double(shapeValues[idx]) * flux[eqIdx];
            }
        }
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
    std::vector<DirichletConstraintData> constraints_;
};

} // end namespace Dumux
#endif
