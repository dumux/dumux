// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
#ifndef DUMUX_KIRCHHOFF_LOVE_PLATE_FREE_CORNER_TEST_PROBLEM_HH
#define DUMUX_KIRCHHOFF_LOVE_PLATE_FREE_CORNER_TEST_PROBLEM_HH

#include <cmath>

#include <dumux/common/boundarytypes.hh>
#include <dumux/common/fvproblem.hh>
#include <dumux/common/numeqvector.hh>
#include <dumux/common/parameters.hh>
#include <dumux/common/properties.hh>
#include <dumux/common/math.hh>

namespace Dumux {

/*!
 * \brief Rotation sub-problem of a plate with clamped and free edges
 *
 * The free-edge traction imposes M_nn = 0 and psi = M_ns. Continuity of psi
 * enforces matching twisting moments at corners where both incident edges are free.
 */
template<class TypeTag>
class FreeCornerProblemRotation : public FVProblem<TypeTag>
{
    using ParentType = FVProblem<TypeTag>;
    using GridGeometry = GetPropType<TypeTag, Properties::GridGeometry>;
    using FVElementGeometry = typename GridGeometry::LocalView;
    using SubControlVolume = typename GridGeometry::SubControlVolume;
    using SubControlVolumeFace = typename GridGeometry::SubControlVolumeFace;
    using GridView = typename GridGeometry::GridView;
    using Element = typename GridView::template Codim<0>::Entity;
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using PrimaryVariables = GetPropType<TypeTag, Properties::PrimaryVariables>;
    using NumEqVector = Dumux::NumEqVector<PrimaryVariables>;
    using BoundaryTypes = Dumux::BoundaryTypes<GetPropType<TypeTag, Properties::ModelTraits>::numEq()>;
    using GlobalPosition = typename Element::Geometry::GlobalCoordinate;

    using CouplingManager = GetPropType<TypeTag, Properties::CouplingManager>;

public:
    FreeCornerProblemRotation(std::shared_ptr<const GridGeometry> gridGeometry,
                              std::shared_ptr<CouplingManager> couplingManager)
    : ParentType(gridGeometry)
    , couplingManager_(couplingManager)
    , poissonRatio_(getParam<Scalar>("Problem.PoissonRatio"))
    , length_(getParam<Scalar>("Problem.Length", 1.0))
    , eps_(1e-7*getParam<Scalar>("Problem.Length", 1.0))
    {
        const auto E = getParam<Scalar>("Problem.E");
        const auto t = getParam<Scalar>("Problem.Thickness");
        stiffness_ = E*t*t*t/(12.0*(1.0 - poissonRatio_*poissonRatio_));
    }

    Scalar D(const GlobalPosition& globalPos) const { return stiffness_; }
    Scalar poissonRatio(const GlobalPosition& globalPos) const { return poissonRatio_; }

    BoundaryTypes boundaryTypes(const Element& element, const SubControlVolume& scv) const
    {
        BoundaryTypes values;
        if (onClampedEdge(scv.dofPosition()))
            values.setAllDirichlet();
        else
            values.setAllNeumann();
        return values;
    }

    template<class ElementVolumeVariables, class ElementFluxVariablesCache>
    NumEqVector neumann(const Element& element,
                        const FVElementGeometry& fvGeometry,
                        const ElementVolumeVariables& elemVolVars,
                        const ElementFluxVariablesCache& elemFluxVarsCache,
                        const SubControlVolumeFace& scvf) const
    {
        NumEqVector values(0.0);
        const auto vars = this->couplingManager().deformationAndPotentials(fvGeometry, scvf);
        const auto phi = vars[this->couplingManager().shearGradPotentialIdx()];
        values.axpy(phi, scvf.unitOuterNormal());
        return values;
    }

    PrimaryVariables dirichletAtPos(const GlobalPosition& globalPos) const
    { return PrimaryVariables(0.0); }

    PrimaryVariables initialAtPos(const GlobalPosition& globalPos) const
    { return PrimaryVariables(0.0); }

    bool onClampedEdge(const GlobalPosition& globalPos) const
    {
        return (clampX0_ && globalPos[0] < eps_) || (clampY0_ && globalPos[1] < eps_);
    }

    const CouplingManager& couplingManager() const { return *couplingManager_; }

private:
    std::shared_ptr<CouplingManager> couplingManager_;
    Scalar poissonRatio_, length_, eps_, stiffness_;
    bool clampX0_ = getParam<bool>("Problem.ClampX0", true);
    bool clampY0_ = getParam<bool>("Problem.ClampY0", true);
};

/*!
 * \brief Deformation sub-problem of a plate with clamped and free edges
 */
template<class TypeTag>
class FreeCornerProblemDeformation : public FVProblem<TypeTag>
{
    using ParentType = FVProblem<TypeTag>;
    using GridGeometry = GetPropType<TypeTag, Properties::GridGeometry>;
    using FVElementGeometry = typename GridGeometry::LocalView;
    using SubControlVolume = typename GridGeometry::SubControlVolume;
    using SubControlVolumeFace = typename GridGeometry::SubControlVolumeFace;
    using GridView = typename GridGeometry::GridView;
    using Element = typename GridView::template Codim<0>::Entity;
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using PrimaryVariables = GetPropType<TypeTag, Properties::PrimaryVariables>;
    using NumEqVector = Dumux::NumEqVector<PrimaryVariables>;
    using BoundaryTypes = Dumux::BoundaryTypes<GetPropType<TypeTag, Properties::ModelTraits>::numEq()>;
    using GlobalPosition = typename Element::Geometry::GlobalCoordinate;
    using Indices = typename GetPropType<TypeTag, Properties::ModelTraits>::Indices;

    using CouplingManager = GetPropType<TypeTag, Properties::CouplingManager>;

public:
    FreeCornerProblemDeformation(std::shared_ptr<const GridGeometry> gridGeometry,
                                 std::shared_ptr<CouplingManager> couplingManager)
    : ParentType(gridGeometry)
    , couplingManager_(couplingManager)
    , force_(getParam<Scalar>("Problem.Force"))
    , eps_(1e-7*getParam<Scalar>("Problem.Length", 1.0))
    {}

    BoundaryTypes boundaryTypes(const Element& element, const SubControlVolume& scv) const
    {
        BoundaryTypes values;
        if (onClampedEdge(scv.dofPosition()))
        {
            values.setDirichlet(Indices::shearGradPotentialIdx);
            values.setDirichlet(Indices::verticalDeformationIdx);
            values.setNeumann(Indices::shearCurlPotentialEqIdx);
        }
        else
            values.setAllNeumann();
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
        NumEqVector values(0.0);
        const auto rotation = this->couplingManager().rotation(fvGeometry, scvf);
        const auto tangent = [&](){
            auto tangent = scvf.unitOuterNormal();
            std::swap(tangent[0], tangent[1]);
            tangent[1] = -tangent[1];
            return tangent;
        }();
        values[Indices::shearCurlPotentialEqIdx] = -vtmv(tangent, 1.0, rotation);
        return values;
    }

    NumEqVector sourceAtPos(const GlobalPosition& globalPos) const
    {
        NumEqVector source(0.0);
        source[Indices::shearGradPotentialEqIdx] = force_;
        return source;
    }

    PrimaryVariables initialAtPos(const GlobalPosition& globalPos) const
    { return PrimaryVariables(0.0); }

    bool onClampedEdge(const GlobalPosition& globalPos) const
    {
        return (clampX0_ && globalPos[0] < eps_) || (clampY0_ && globalPos[1] < eps_);
    }

    const CouplingManager& couplingManager() const { return *couplingManager_; }

private:
    std::shared_ptr<CouplingManager> couplingManager_;
    Scalar force_, eps_;
    bool clampX0_ = getParam<bool>("Problem.ClampX0", true);
    bool clampY0_ = getParam<bool>("Problem.ClampY0", true);
};

} // end namespace Dumux

#endif
