// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
#ifndef DUMUX_COOKS_MEMBRANE_PROBLEM_HH
#define DUMUX_COOKS_MEMBRANE_PROBLEM_HH

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

// Cook's membrane benchmark for volumetric locking
// Reference: R.D. Cook, "Improved two-dimensional finite element",
// ASCE J. Struct. Div. 100(9):1851-1863, 1974
//
// Geometry: tapered panel with vertices (0,0), (48,44), (48,60), (0,44)
// Boundary conditions:
//   - Left face (x=0): clamped (zero displacement)
//   - Right face (x=48): uniform shear traction in y-direction
//   - Top/bottom edges: traction-free
//
// To study volumetric locking, use nearly incompressible material (nu -> 0.5).
template<class TypeTag>
class CooksMembraneProblem : public FVProblemWithSpatialParams<TypeTag>
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

    static constexpr int numEq = GetPropType<TypeTag, Properties::ModelTraits>::numEq();
    using ConstraintInfo = Dumux::DirichletConstraintInfo<numEq>;
    using ConstraintValues = Dune::FieldVector<Scalar, numEq>;
    using DirichletConstraintData = Dumux::DirichletConstraintData<ConstraintInfo, ConstraintValues, GridIndexType>;

public:
    CooksMembraneProblem(std::shared_ptr<const GridGeometry> gridGeometry)
    : ParentType(gridGeometry)
    , shearTraction_(getParam<Scalar>("Problem.ShearTraction"))
    {
        // build Dirichlet constraints for ALL DOFs on the left face (x=0)
        // this works for both CV DOFs (vertices) and non-CV DOFs (edge midpoints)
        CVFE::appendDirichletConstraints(*this,
            [&](const auto& fvGeometry, const auto&, const auto& localDof) {
                return this->dirichletAtPos(ipData(fvGeometry, localDof).global());
            },
            constraints_
        );
    }

    // boundary types at a given position (old-style, used for vertex SCVs)
    BoundaryTypes boundaryTypesAtPos(const GlobalPosition& globalPos) const
    {
        BoundaryTypes values;
        if (globalPos[0] < eps_)
            values.setAllDirichlet();
        else
            values.setAllNeumann();
        return values;
    }

    // boundary types for an intersection (new-style, used by appendDirichletConstraints)
    BoundaryTypes boundaryTypes(const FVElementGeometry&, const Intersection& intersection) const
    { return boundaryTypesAtPos(intersection.geometry().center()); }

    PrimaryVariables dirichletAtPos(const GlobalPosition& globalPos) const
    {
        // left face (x=0): clamped — zero displacement in all directions
        assert(globalPos[0] < eps_);
        return PrimaryVariables(0.0);
    }

    //! Old Neumann interface (FVAssembler path): traction vector with positive sign convention.
    template<class ElementVolumeVariables, class ElementFluxVariablesCache>
    NumEqVector neumann(const Element& element,
                        const FVElementGeometry& fvGeometry,
                        const ElementVolumeVariables& elemVolVars,
                        const ElementFluxVariablesCache& elemFluxVarsCache,
                        const SubControlVolumeFace& scvf) const
    {
        NumEqVector traction(0.0);
        if (scvf.ipGlobal()[0] > 48.0 - eps_)
            traction[1] = shearTraction_;
        return traction;
    }

    //! Per-QP traction for the new Experimental interface (sign: -T_applied).
    template<class EV, class IpData>
    NumEqVector boundaryFlux(const FVElementGeometry&, const EV&,
                             const IpData& ipData) const
    {
        NumEqVector traction(0.0);
        if (ipData.global()[0] > 48.0 - eps_)
            traction[1] = -shearTraction_;
        return traction;
    }

    //! Integrated Neumann traction over boundary scvf (Experimental assembler path).
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

    //! Dirichlet constraints applied globally by FVAssembler / Experimental::Assembler.
    const auto& constraints() const
    { return constraints_; }

private:
    static constexpr Scalar eps_ = 1e-6;
    Scalar shearTraction_;
    std::vector<DirichletConstraintData> constraints_;
};

} // end namespace Dumux

#endif
