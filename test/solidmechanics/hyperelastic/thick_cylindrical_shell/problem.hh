// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
#ifndef DUMUX_HYPERELASTIC_THICK_CYLINDRICAL_SHELL_PROBLEM_HH
#define DUMUX_HYPERELASTIC_THICK_CYLINDRICAL_SHELL_PROBLEM_HH

#include <vector>

#include <dumux/common/boundarytypes_.hh>
#include <dumux/common/constraintinfo.hh>
#include <dumux/common/indextraits.hh>
#include <dumux/common/problemwithspatialparams.hh>
#include <dumux/common/numeqvector.hh>
#include <dumux/common/parameters.hh>
#include <dumux/common/properties.hh>
#include <dumux/common/math.hh>
#include <dumux/discretization/cvfe/appenddirichletconstraints_.hh>
#include <dumux/discretization/dirichletconstraints.hh>

namespace Dumux {

/*!
 * \brief Thick cylindrical shell under pressure — Elguedj et al. (2008) F-bar locking benchmark.
 *   DOI: 10.1016/j.cma.2008.01.012, Section 5.4.
 *
 * Geometry: half of a thick cylindrical shell (the x <= 0 half of an annulus),
 *   inner radius 0.008 m, outer radius 0.010 m (thickness 0.002 m), height 0.015 m.
 *   The cylinder axis is z. The crown of the arch is at (-0.01, 0, z).
 *
 * Boundary conditions (cf. Reese 2000, Fig. 1; Elguedj 2008, Fig. 22):
 *   x = 0 plane (the two radial cut sections):  u_x = 0   (symmetry, "X1 = 0")
 *   z = 0 plane:                                 u_z = 0   (symmetry, "X3 = 0")
 *   bottom section (x = 0, y < 0):               u_y = 0   (support)
 *   top section   (x = 0, y > 0):                downward line load (load)
 *
 * Load: Reese 2000 thick shell (t = 2 mm) uses load factor ν = p/p0 = 7500 with
 *   p0 = 1/15 N/mm, i.e. line load p = 500 N/mm => 7500 N total on the quarter,
 *   applied here as a traction over the top cut face (see params.input).
 *
 * Material: compressible neo-Hookean (Reese 2000 §4.1 / Elguedj 2008, Eq. 109),
 *   Ψ = µ/2(tr[C] − 3) − µ ln J + Λ/2(ln J)²,
 *   P = µ F + (Λ ln J − µ) F⁻ᵀ,
 * where Λ = λ (Lamé first parameter) in the small-strain limit.
 * Parameters: µ = 6000 N/mm², Λ = 240000 N/mm²  (=> ν ≈ 0.488, E ≈ 17.85 GPa).
 *
 * Uses the new experimental boundary-condition interface:
 *   boundaryTypes(FVElementGeometry, BoundaryFace) -> Experimental::BoundaryTypes<dim>
 * and builds Dirichlet DOFs globally via CVFE::appendDirichletConstraints.
 */
template<class TypeTag>
class ThickCylindricalShellProblem : public Experimental::ProblemWithSpatialParams<TypeTag>
{
    using ParentType = Experimental::ProblemWithSpatialParams<TypeTag>;
    using GridGeometry = GetPropType<TypeTag, Properties::GridGeometry>;
    using FVElementGeometry = typename GridGeometry::LocalView;
    using BoundaryFace = typename GridGeometry::BoundaryFace;
    using GridView = typename GridGeometry::GridView;
    using Element = typename GridView::template Codim<0>::Entity;
    using GridIndexType = typename IndexTraits<GridView>::GridIndex;
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using PrimaryVariables = GetPropType<TypeTag, Properties::PrimaryVariables>;
    using NumEqVector = Dumux::NumEqVector<PrimaryVariables>;
    using BoundaryTypes = Dumux::Experimental::BoundaryTypes<GetPropType<TypeTag, Properties::ModelTraits>::numEq()>;
    using GlobalPosition = typename Element::Geometry::GlobalCoordinate;
    static constexpr int dim = GridView::dimension;
    static constexpr int numEq = GetPropType<TypeTag, Properties::ModelTraits>::numEq();
    using Tensor = Dune::FieldMatrix<Scalar, dim, dim>;
    using ConstraintInfo = Dumux::DirichletConstraintInfo<numEq>;
    using ConstraintValues = Dune::FieldVector<Scalar, numEq>;
    using DirichletConstraintData = Dumux::DirichletConstraintData<ConstraintInfo, ConstraintValues, GridIndexType>;

public:
    ThickCylindricalShellProblem(std::shared_ptr<const GridGeometry> gridGeometry)
    : ParentType(gridGeometry)
    , load_(getParam<Scalar>("Problem.Load"))
    {
        CVFE::appendDirichletConstraints(*this,
            [&](const auto& fvGeometry, const auto& /*bf*/, const auto& localDof) {
                return PrimaryVariables(0.0);
            },
            constraints_);
    }

    //! New boundary-types interface: FVElementGeometry + BoundaryFace.
    BoundaryTypes boundaryTypes(const FVElementGeometry&, const BoundaryFace& bf) const
    { return boundaryTypesAtPos(bf.center()); }

    const auto& constraints() const { return constraints_; }

    BoundaryTypes boundaryTypesAtPos(const GlobalPosition& globalPos) const
    {
        BoundaryTypes values;
        values.setAllFluxBoundary();                 // traction-free (Neumann) by default

        if (globalPos[0] > -eps_)                    // x = 0 plane (both cut sections)
            values.resetEq(0);                       // u_x = 0 (symmetry)
        if (globalPos[2] < eps_)                     // z = 0 plane
            values.resetEq(2);                       // u_z = 0 (symmetry)
        if (globalPos[0] > -eps_ && globalPos[1] < -eps_)
            values.resetEq(1);                       // bottom section: u_y = 0 (support)

        return values;
    }

    PrimaryVariables dirichletAtPos(const GlobalPosition&) const
    { return PrimaryVariables(0.0); }

    //! Set the load ramp factor (load stepping in main).
    void setLoadFactor(Scalar f) { loadFactor_ = f; }

    //! Downward traction on the full top cut section (x = 0, y > 0).
    //! Sign convention: flux = -T_applied, so traction[1] = +load_ encodes a
    //! downward (-y) applied force.
    template<class EV, class IpData>
    NumEqVector boundaryFlux(const FVElementGeometry&, const EV&, const IpData& ipData) const
    {
        NumEqVector traction(0.0);
        const auto& gp = ipData.global();
        if (gp[0] > -eps_ && gp[1] > eps_)
            traction[1] = loadFactor_*load_;
        return traction;
    }

    // Compressible neo-Hookean (Reese 2000 §4.1 / Elguedj 2008, Eq. 109):
    //   Ψ = µ/2 (tr[C] − 3) − µ ln J + Λ/2 (ln J)²
    //   P = µ F + (Λ ln J − µ) F⁻ᵀ
    //
    // Λ = λ (Lamé first parameter) in the small-strain limit.
    // Reduces to σ = 2µ ε + λ tr(ε) I at small strain.
    Tensor firstPiolaKirchhoffStressTensor(Tensor F) const
    {
        using std::log;
        const auto J = F.determinant();
        auto invFT = F; invFT.invert(); invFT = transpose(invFT);

        const auto mu     = this->spatialParams().shearModulus();
        const auto lambda = this->spatialParams().firstLameParameter();

        Tensor P = F;
        P *= mu;                                            // µ F
        P.axpy(lambda * log(J) - mu, invFT);               // + (Λ ln J − µ) F⁻ᵀ
        return P;
    }

private:
    static constexpr Scalar eps_ = 1e-7;
    Scalar load_;
    Scalar loadFactor_ = 1.0;
    std::vector<DirichletConstraintData> constraints_;
};

} // end namespace Dumux
#endif
