// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
#ifndef DUMUX_COOKS_MEMBRANE_MOMENTUM_MIXED_PROBLEM_HH
#define DUMUX_COOKS_MEMBRANE_MOMENTUM_MIXED_PROBLEM_HH

#include <dumux/common/boundarytypes.hh>
#include <dumux/common/fvproblemwithspatialparams.hh>
#include <dumux/common/numeqvector.hh>
#include <dumux/common/parameters.hh>
#include <dumux/common/properties.hh>
#include <dumux/common/constraintinfo.hh>
#include <dumux/common/indextraits.hh>
#include <dumux/discretization/dirichletconstraints.hh>
#include <dumux/discretization/cvfe/quadraturerules.hh>

namespace Dumux {

/*!
 * \brief Momentum problem for the Cook's membrane MINI mixed formulation.
 *
 * Uses the new Donea-style problem interface:
 *  - `constraints()`: global Dirichlet constraints built via
 *    CVFE::appendDirichletConstraints — handles ALL DOFs incl. edge midpoints.
 *  - `boundaryFlux(fvGeometry, elemVars, ipCache)`: called per boundary QP
 *    from the local residual's evalFluxAndSource.
 *  - `boundaryTypes(fvGeometry, intersection)`: new intersection-based interface
 *    used by appendDirichletConstraints.
 */
template<class TypeTag>
class CooksMembraneMomentumMixedProblem : public FVProblemWithSpatialParams<TypeTag>
{
    using ParentType = FVProblemWithSpatialParams<TypeTag>;
    using GridGeometry = GetPropType<TypeTag, Properties::GridGeometry>;
    using FVElementGeometry = typename GridGeometry::LocalView;
    using SubControlVolumeFace = typename GridGeometry::SubControlVolumeFace;
    using GridView = typename GridGeometry::GridView;
    using Intersection = typename GridView::Intersection;
    using Element = typename GridView::template Codim<0>::Entity;
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using PrimaryVariables = GetPropType<TypeTag, Properties::PrimaryVariables>;
    using NumEqVector = Dumux::NumEqVector<PrimaryVariables>;
    using BoundaryTypes = Dumux::BoundaryTypes<GetPropType<TypeTag, Properties::ModelTraits>::numEq()>;
    using GlobalPosition = typename Element::Geometry::GlobalCoordinate;
    using GridIndexType = typename IndexTraits<GridView>::GridIndex;
    using CouplingManager = GetPropType<TypeTag, Properties::CouplingManager>;

    static constexpr int numEq = GetPropType<TypeTag, Properties::ModelTraits>::numEq();
    using ConstraintInfo = Dumux::DirichletConstraintInfo<numEq>;
    using ConstraintValues = Dune::FieldVector<Scalar, numEq>;
    using DirichletConstraintData = Dumux::DirichletConstraintData<
        ConstraintInfo, ConstraintValues, GridIndexType>;

public:
    CooksMembraneMomentumMixedProblem(std::shared_ptr<const GridGeometry> gridGeometry,
                                      std::shared_ptr<CouplingManager> couplingManager)
    : ParentType(gridGeometry)
    , couplingManager_(couplingManager)
    , shearTraction_(getParam<Scalar>("Problem.ShearTraction"))
    {
        // Build Dirichlet constraints for ALL DOFs on the left face (x=0).
        // Works for PQ1Bubble (vertex + bubble) and PQ2 (vertex + edge midpoints).
        CVFE::appendDirichletConstraints(*this,
            [&](const auto& fvGeometry, const auto&, const auto& localDof) {
                return this->dirichletAtPos(ipData(fvGeometry, localDof).global());
            },
            constraints_);
    }

    const CouplingManager& couplingManager() const { return *couplingManager_; }

    // New intersection-based BC type (used by appendDirichletConstraints)
    BoundaryTypes boundaryTypes(const FVElementGeometry&, const Intersection& is) const
    { return boundaryTypesAtPos(is.geometry().center()); }

    // Position-based BC type.
    // 2D: left face clamped, others free.
    // 3D: additionally z=0 face has u_z=0 (plane-strain symmetry).
    BoundaryTypes boundaryTypesAtPos(const GlobalPosition& globalPos) const
    {
        BoundaryTypes values;
        if (globalPos[0] < eps_)
        {
            values.setAllDirichlet(); // left face: fully clamped
        }
        else if constexpr (GlobalPosition::dimension == 3)
        {
            if (globalPos[2] < eps_) // z=0 plane-strain symmetry face
            {
                values.setAllNeumann();
                values.setDirichlet(2); // u_z = 0 only
            }
            else
                values.setAllNeumann();
        }
        else
            values.setAllNeumann();
        return values;
    }

    PrimaryVariables dirichletAtPos(const GlobalPosition&) const
    { return PrimaryVariables(0.0); }

    //! Neumann traction per quadrature point. Used by both CV (boundaryFluxIntegral)
    //! and FE (addToElementFluxAndSourceResidual) Neumann integration paths.
    template<class ElementVariables, class IpData>
    NumEqVector boundaryFlux(const FVElementGeometry&, const ElementVariables&,
                             const IpData& ipData) const
    {
        NumEqVector flux(0.0);
        if (ipData.global()[0] > 48.0 - eps_)
            flux[1] = -shearTraction_;
        return flux;
    }

    //! Integrated Neumann traction over a boundary scvf (CV path).
    //! Called by Experimental::CVFELocalResidual::evalFlux for Neumann scvfs.
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

    //! Dirichlet constraints — applied globally by the assembler after assembly.
    const auto& constraints() const { return constraints_; }

private:
    static constexpr Scalar eps_ = 1e-6;
    std::shared_ptr<CouplingManager> couplingManager_;
    Scalar shearTraction_;
    std::vector<DirichletConstraintData> constraints_;
};

} // end namespace Dumux
#endif
