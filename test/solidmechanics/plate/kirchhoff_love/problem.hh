// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
#ifndef DUMUX_KIRCHHOFF_LOVE_PLATE_TEST_PROBLEM_HH
#define DUMUX_KIRCHHOFF_LOVE_PLATE_TEST_PROBLEM_HH

#include <cmath>

#include <dumux/common/boundarytypes.hh>
#include <dumux/common/fvproblem.hh>
#include <dumux/common/numeqvector.hh>
#include <dumux/common/parameters.hh>
#include <dumux/common/properties.hh>
#include <dumux/common/math.hh>

namespace Dumux {

template<class TypeTag>
class KirchhoffLovePlateTestProblemRotation : public FVProblem<TypeTag>
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
    static constexpr int dimWorld = GridView::dimensionworld;
    using GlobalPosition = typename Element::Geometry::GlobalCoordinate;
    using Indices = typename GetPropType<TypeTag, Properties::ModelTraits>::Indices;

    using CouplingManager = GetPropType<TypeTag, Properties::CouplingManager>;

public:
    KirchhoffLovePlateTestProblemRotation(std::shared_ptr<const GridGeometry> gridGeometry,
                                         std::shared_ptr<CouplingManager> couplingManager)
    : ParentType(gridGeometry)
    , couplingManager_(couplingManager)
    , youngsModulus_(getParam<Scalar>("Problem.E"))
    , poissonRatio_(getParam<Scalar>("Problem.PoissonRatio"))
    , thickness_(getParam<Scalar>("Problem.Thickness"))
    {
        const auto E = youngsModulus_;
        const auto t = thickness_;
        const auto nu = poissonRatio_;
        stiffness_ = E*t*t*t/(12*(1-nu*nu));
    }

    Scalar D(const GlobalPosition& globalPos) const { return stiffness_; }
    Scalar poissonRatio(const GlobalPosition& globalPos) const { return poissonRatio_; }

    BoundaryTypes boundaryTypes(const Element& element, const SubControlVolume& scv) const
    {
        BoundaryTypes values;
        values.setAllDirichlet();
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

    const CouplingManager& couplingManager() const
    { return *couplingManager_; }

private:
    std::shared_ptr<CouplingManager> couplingManager_;
    Scalar youngsModulus_, stiffness_, poissonRatio_, thickness_;
};

template<class TypeTag>
class KirchhoffLovePlateTestProblemDeformation : public FVProblem<TypeTag>
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
    static constexpr int dimWorld = GridView::dimensionworld;
    using GlobalPosition = typename Element::Geometry::GlobalCoordinate;
    using Indices = typename GetPropType<TypeTag, Properties::ModelTraits>::Indices;

    using CouplingManager = GetPropType<TypeTag, Properties::CouplingManager>;

public:
    KirchhoffLovePlateTestProblemDeformation(std::shared_ptr<const GridGeometry> gridGeometry,
                                            std::shared_ptr<CouplingManager> couplingManager)
    : ParentType(gridGeometry)
    , couplingManager_(couplingManager)
    , force_(getParam<Scalar>("Problem.Force", 0.0))
    , youngsModulus_(getParam<Scalar>("Problem.E"))
    , radius_(getParam<Scalar>("Problem.Radius", 1.0))
    , thickness_(getParam<Scalar>("Problem.Thickness"))
    , poissonRatio_(getParam<Scalar>("Problem.PoissonRatio"))
    {
        const auto E = youngsModulus_;
        const auto t = thickness_;
        const auto nu = poissonRatio_;
        stiffness_ = E*t*t*t/(12*(1-nu*nu));
    }

    BoundaryTypes boundaryTypes(const Element& element, const SubControlVolume& scv) const
    {
        BoundaryTypes values;
        values.setDirichlet(Indices::shearGradPotentialIdx);
        values.setDirichlet(Indices::verticalDeformationIdx);
        values.setNeumann(Indices::shearCurlPotentialEqIdx);
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

    static constexpr bool enableInternalDirichletConstraints()
    { return true; }

    std::bitset<NumEqVector::dimension>
    hasInternalDirichletConstraint(const Element& element, const SubControlVolume& scv) const
    {
        std::bitset<NumEqVector::dimension> values;

        // Pin the shear curl potential at one interior DOF to fix its gauge freedom
        // (it is only defined up to a constant, so the system is otherwise singular).
        // Here we arbitrarily choose a DOF near the center of the plate.
        if (scv.dofIndex() == static_cast<std::size_t>(this->gridGeometry().numDofs()/2))
            values.set(Indices::shearCurlPotentialIdx);

        return values;
    }

    PrimaryVariables internalDirichlet(const Element& element, const SubControlVolume& scv) const
    { return PrimaryVariables(0.0); }

    PrimaryVariables initialAtPos(const GlobalPosition& globalPos) const
    { return PrimaryVariables(0.0); }

    const CouplingManager& couplingManager() const
    { return *couplingManager_; }

    // Analytic solution: clamped circular plate under uniform load (Kirchhoff-Love)
    Scalar analyticDeformation(const GlobalPosition& globalPos) const
    {
        const auto R = radius_;
        const auto r = std::hypot(globalPos[0], globalPos[1]);
        return -force_/(64*stiffness_)*(R*R - r*r)*(R*R - r*r);
    }

private:
    std::shared_ptr<CouplingManager> couplingManager_;
    Scalar force_, youngsModulus_, radius_, thickness_, poissonRatio_, stiffness_;
};

} // end namespace Dumux

#endif
