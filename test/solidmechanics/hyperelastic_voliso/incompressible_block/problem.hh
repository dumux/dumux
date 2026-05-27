// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
#ifndef DUMUX_HYPERELASTIC_INCOMPRESSIBLE_BLOCK_PROBLEM_HH
#define DUMUX_HYPERELASTIC_INCOMPRESSIBLE_BLOCK_PROBLEM_HH

#include <dumux/common/boundarytypes.hh>
#include <dumux/common/problemwithspatialparams.hh>
#include <dumux/common/numeqvector.hh>
#include <dumux/common/parameters.hh>
#include <dumux/common/properties.hh>
#include <dumux/common/constraintinfo.hh>
#include <dumux/common/indextraits.hh>
#include <dumux/discretization/dirichletconstraints.hh>
#include <dumux/discretization/cvfe/quadraturerules.hh>

namespace Dumux {

// SPP-1748 Benchmark 3: quasi-incompressible cube under constant partial load
// (Kollmannsberger et al., Section 3, SPP-1748 Benchmark Collection, 2020).
//
// Full model: w=100mm, l=100mm, h=50mm (Table 3.1). Axial symmetry → quarter model
// [0, w/2]×[0, l/2]×[0, h] = [0,50]³ mm is simulated.
// Strain energy (Eq. 3.2): W = µ/2·(trC − 3 − 2·ln J) + λ/2·(J−1)²
// Material (Table 3.2): λ=499.92568 MPa, µ=1.61148 MPa (ν≈0.4983)
//
// BCs (Section 3.1):
//   z=0   : u_z=0 (bottom fixed in z)
//   z=h   : u_x=0, u_y=0 (top surface sliding constraint)
//   x=w/2 : u_x=0 (symmetry plane)
//   y=l/2 : u_y=0 (symmetry plane)
//   Neumann load q=3 MPa at z=h on quarter patch x∈[w/2−a, w/2]×y∈[l/2−b, l/2] = [25,50]²
//   with a=25mm, b=25mm (Table 3.1) the load half-widths in x and y
//
// Probe: |u_z(P)| at P=(w/2, l/2, h) = (50,50,50) mm (Tables 3.3 & 3.4, ref ≈ 20 mm)
template<class TypeTag>
class IncompressibleBlockMomentumProblem : public Dumux::Experimental::ProblemWithSpatialParams<TypeTag>
{
    using ParentType = Dumux::Experimental::ProblemWithSpatialParams<TypeTag>;
    using GridGeometry = GetPropType<TypeTag, Properties::GridGeometry>;
    using FVElementGeometry = typename GridGeometry::LocalView;
    using GridView = typename GridGeometry::GridView;
    using Intersection = typename GridView::Intersection;
    using Element = typename GridView::template Codim<0>::Entity;
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using PrimaryVariables = GetPropType<TypeTag, Properties::PrimaryVariables>;
    using NumEqVector = Dumux::NumEqVector<PrimaryVariables>;
    using BoundaryTypes = Dumux::Experimental::BoundaryTypes<PrimaryVariables::size()>;
    using GlobalPosition = typename Element::Geometry::GlobalCoordinate;
    using GridIndexType = typename IndexTraits<GridView>::GridIndex;
    using CouplingManager = GetPropType<TypeTag, Properties::CouplingManager>;

    static constexpr int numEq = GetPropType<TypeTag, Properties::ModelTraits>::numEq();
    using ConstraintInfo = Dumux::DirichletConstraintInfo<numEq>;
    using ConstraintValues = Dune::FieldVector<Scalar, numEq>;
    using DirichletConstraintData = Dumux::DirichletConstraintData<
        ConstraintInfo, ConstraintValues, GridIndexType>;

    static constexpr Scalar h_  = 50.0;   // height z = h (Table 3.1)
    static constexpr Scalar w2_ = 50.0;   // half-width x = w/2 (Table 3.1: w=100mm)
    static constexpr Scalar l2_ = 50.0;   // half-length y = l/2 (Table 3.1: l=100mm)
    static constexpr Scalar a_  = 25.0;   // load half-width in x (Table 3.1: a=25mm)
    static constexpr Scalar b_  = 25.0;   // load half-width in y (Table 3.1: b=25mm)

public:
    IncompressibleBlockMomentumProblem(std::shared_ptr<const GridGeometry> gridGeometry,
                                       std::shared_ptr<CouplingManager> couplingManager)
    : ParentType(gridGeometry)
    , couplingManager_(couplingManager)
    , load_(getParam<Scalar>("Problem.Load"))
    , loadFactor_(1.0)
    {
        appendDirichletConstraints_();
    }

    void setLoadFactor(Scalar f) { loadFactor_ = f; }

    const CouplingManager& couplingManager() const { return *couplingManager_; }

    // First Piola-Kirchhoff stress P = ∂W/∂F from the vol-iso split of Eq. 3.2:
    //   W = µ/2·(trC − 3 − 2·ln J) + λ/2·(J−1)²
    // Isochoric part gives µ·(F − F^{−T}), volumetric part gives J·p_s·F^{−T}
    // where p_s is the mixed pressure variable (= λ·(J−1) at equilibrium).
    template<class Tensor, class FVElementGeometry_, class GlobalPosition_>
    Tensor firstPiolaKirchhoffStressTensor(const Tensor& F,
                                           const FVElementGeometry_& fvGeometry,
                                           const GlobalPosition_& globalPos) const
    {
        const Scalar mu = this->spatialParams().shearModulus();
        const Scalar ps = couplingManager().pressureAtPoint(fvGeometry, globalPos);

        Tensor Finv = F;
        Finv.invert();
        const Tensor FinvT = transpose(Finv);

        Tensor P(0.0);
        const Scalar J = F.determinant();
        for (int i = 0; i < P.N(); ++i)
            for (int j = 0; j < P.M(); ++j)
                P[i][j] = mu*(F[i][j] - FinvT[i][j]) + J*ps*FinvT[i][j];

        return P;
    }

    BoundaryTypes boundaryTypesAtPos(const GlobalPosition& pos) const
    {
        BoundaryTypes values{};

        if (pos[2] > h_ - eps_ && pos[0] > w2_ - a_ - eps_ && pos[1] > l2_ - b_ - eps_)
            values.setAllFluxBoundary();

        return values;
    }

    template<class ElementVariables, class IpData>
    NumEqVector boundaryFlux(const FVElementGeometry&, const ElementVariables&,
                             const IpData& ipData) const
    {
        NumEqVector flux(0.0);
        const auto& pos = ipData.global();
        // Compressive load q (Table 3.1) at z=h on patch x∈[w/2−a, w/2]×y∈[l/2−b, l/2]
        // (Section 3.1). Applied traction T = −q·e_z → boundaryFlux[z] = +q.
        if (pos[2] > h_ - eps_ && pos[0] > w2_ - a_ - eps_ && pos[1] > l2_ - b_ - eps_)
            flux[2] = loadFactor_ * load_;

        return flux;
    }

    const auto& constraints() const { return constraints_; }

private:
    auto boundaryTypesAtPos_(const GlobalPosition& globalPos) const
    {
        Dumux::BoundaryTypes<numEq> values;
        values.setAllNeumann();

        if (globalPos[2] < eps_)
            values.setDirichlet(2);                     // bottom: u_z = 0
        else if (globalPos[2] > h_ - eps_)
        {
            values.setDirichlet(0);                     // top: u_x = 0
            values.setDirichlet(1);                     // top: u_y = 0
        }
        else if (globalPos[0] > w2_ - eps_)
            values.setDirichlet(0);                     // symmetry x: u_x = 0
        else if (globalPos[1] > l2_ - eps_)
            values.setDirichlet(1);                     // symmetry y: u_y = 0

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
                const auto& bcTypes = this->boundaryTypesAtPos_(boundaryFace.center());
                if (bcTypes.hasDirichlet())
                {
                    for (const auto& localDof : localDofs(fvGeometry, boundaryFace))
                    {
                        ConstraintInfo info;
                        // set the Dirichlet constraints
                        for (int eqIdx = 0; eqIdx < numEq; ++eqIdx)
                            if (bcTypes.isDirichlet(eqIdx))
                                info.set(bcTypes.eqToDirichletIndex(eqIdx), eqIdx);

                        auto dirichletValues = PrimaryVariables(0.0);
                        constraints_.push_back(DirichletConstraintData{std::move(info), std::move(dirichletValues), localDof.dofIndex()});
                    }
                }
            }
        }

    }

    static constexpr Scalar eps_ = 1e-6;
    std::shared_ptr<CouplingManager> couplingManager_;
    Scalar load_, loadFactor_;
    std::vector<DirichletConstraintData> constraints_;
};


// Pressure subdomain problem for the incompressible block P1BP1 (mini/mixed) formulation.
// The equation is the local bulk pressure constraint — no spatial BCs needed.
template<class TypeTag>
class IncompressibleBlockPressureProblem : public Dumux::Experimental::ProblemWithSpatialParams<TypeTag>
{
    using ParentType = Dumux::Experimental::ProblemWithSpatialParams<TypeTag>;
    using GridGeometry = GetPropType<TypeTag, Properties::GridGeometry>;
    using GridView = typename GridGeometry::GridView;
    using Element = typename GridView::template Codim<0>::Entity;
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using PrimaryVariables = GetPropType<TypeTag, Properties::PrimaryVariables>;
    using NumEqVector = Dumux::NumEqVector<PrimaryVariables>;
    using BoundaryTypes = Dumux::Experimental::BoundaryTypes<PrimaryVariables::size()>;
    using GlobalPosition = typename Element::Geometry::GlobalCoordinate;
    using CouplingManager = GetPropType<TypeTag, Properties::CouplingManager>;

public:
    IncompressibleBlockPressureProblem(std::shared_ptr<const GridGeometry> gridGeometry,
                                       std::shared_ptr<CouplingManager> couplingManager)
    : ParentType(gridGeometry)
    , couplingManager_(couplingManager)
    {}

    const CouplingManager& couplingManager() const { return *couplingManager_; }

    // Equilibrium pressure from W_vol = λ/2·(J−1)² (Eq. 3.2): p_s = dW_vol/dJ = λ·(J−1)
    Scalar volumetricPressure(Scalar J) const
    { return this->spatialParams().firstLameParameter() * (J - 1.0); }

    BoundaryTypes boundaryTypesAtPos(const GlobalPosition&) const
    {
        BoundaryTypes values;
        values.setAllFluxBoundary();
        return values;
    }

    PrimaryVariables initialAtPos(const GlobalPosition&) const
    { return PrimaryVariables(0.0); }

private:
    std::shared_ptr<CouplingManager> couplingManager_;
};

} // end namespace Dumux
#endif
