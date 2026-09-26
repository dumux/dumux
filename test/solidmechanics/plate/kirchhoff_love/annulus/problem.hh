// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
#ifndef DUMUX_KIRCHHOFF_LOVE_PLATE_ANNULUS_TEST_PROBLEM_HH
#define DUMUX_KIRCHHOFF_LOVE_PLATE_ANNULUS_TEST_PROBLEM_HH

#include <cmath>

#include <dune/common/fmatrix.hh>
#include <dune/common/fvector.hh>

#include <dumux/common/boundarytypes.hh>
#include <dumux/common/fvproblem.hh>
#include <dumux/common/numeqvector.hh>
#include <dumux/common/parameters.hh>
#include <dumux/common/properties.hh>
#include <dumux/common/math.hh>

#include "../../annulusplate.hh"

namespace Dumux {

/*!
 * \brief Rotation sub-problem of the annular plate clamped at the inner and free at the outer edge
 *
 * On the free edge the rotations are not prescribed. Instead the boundary traction of
 * the divergence-free tensor \f$ \mathbf{T} = \mathbf{M} - \varphi\mathbf{I} - \psi\mathbf{J} \f$
 * is set to \f$ \mathbf{T}\mathbf{n} = -\varphi\mathbf{n} \f$. Component-wise this reads
 * \f$ M_{nn} = 0 \f$, the moment condition of a free edge, and \f$ \psi = M_{ns} \f$, which
 * ties the curl potential to the twisting moment and thereby turns the natural flux
 * \f$ \partial_n\varphi \f$ of the deformation sub-problem into the Kirchhoff effective
 * shear \f$ V_n = q_n + \partial_s M_{ns} \f$.
 */
template<class TypeTag>
class KirchhoffLovePlateAnnulusProblemRotation : public FVProblem<TypeTag>
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
    KirchhoffLovePlateAnnulusProblemRotation(std::shared_ptr<const GridGeometry> gridGeometry,
                                             std::shared_ptr<CouplingManager> couplingManager)
    : ParentType(gridGeometry)
    , couplingManager_(couplingManager)
    , poissonRatio_(getParam<Scalar>("Problem.PoissonRatio"))
    , innerRadius_(getParam<Scalar>("Problem.InnerRadius"))
    , outerRadius_(getParam<Scalar>("Problem.OuterRadius"))
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
    { return globalPos.two_norm() < 0.5*(innerRadius_ + outerRadius_); }

    const CouplingManager& couplingManager() const
    { return *couplingManager_; }

private:
    std::shared_ptr<CouplingManager> couplingManager_;
    Scalar poissonRatio_, innerRadius_, outerRadius_, stiffness_;
};

/*!
 * \brief Deformation sub-problem of the annular plate clamped at the inner and free at the outer edge
 *
 * The clamped edge prescribes \f$ w = 0 \f$ and fixes the gauge of the Helmholtz
 * potentials with \f$ \varphi = 0 \f$. On the free edge all three equations are natural:
 * the flux of the \f$ \varphi \f$ equation is the effective shear and vanishes, the flux
 * of the \f$ w \f$ equation vanishes because \f$ \boldsymbol{\theta} = \nabla w \f$, and
 * the flux of the \f$ \psi \f$ equation is prescribed with its own value, which leaves the
 * curl constraint unconstrained on the boundary.
 */
template<class TypeTag>
class KirchhoffLovePlateAnnulusProblemDeformation : public FVProblem<TypeTag>
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
    KirchhoffLovePlateAnnulusProblemDeformation(std::shared_ptr<const GridGeometry> gridGeometry,
                                                std::shared_ptr<CouplingManager> couplingManager)
    : ParentType(gridGeometry)
    , couplingManager_(couplingManager)
    , force_(getParam<Scalar>("Problem.Force"))
    , poissonRatio_(getParam<Scalar>("Problem.PoissonRatio"))
    , innerRadius_(getParam<Scalar>("Problem.InnerRadius"))
    , outerRadius_(getParam<Scalar>("Problem.OuterRadius"))
    {
        const auto E = getParam<Scalar>("Problem.E");
        const auto t = getParam<Scalar>("Problem.Thickness");
        stiffness_ = E*t*t*t/(12.0*(1.0 - poissonRatio_*poissonRatio_));
        coefficients_ = annulusPlateCoefficients(innerRadius_, outerRadius_, force_, stiffness_, poissonRatio_);
    }

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
    { return globalPos.two_norm() < 0.5*(innerRadius_ + outerRadius_); }

    const CouplingManager& couplingManager() const
    { return *couplingManager_; }

    Scalar analyticDeformation(const GlobalPosition& globalPos) const
    { return annulusPlateDeflection(globalPos.two_norm(), force_, stiffness_, coefficients_); }

    /*!
     * \brief Analytic shear gradient potential
     *
     * The free edge makes \f$ \varphi \f$ the potential of the effective shear: it solves
     * \f$ \Delta\varphi = F \f$ with \f$ \partial_n\varphi = 0 \f$ at the free edge and the
     * gauge \f$ \varphi = 0 \f$ at the clamped one, so axisymmetry gives
     * \f$ \varphi' = Fr/2 - Fb^2/(2r) \f$. That derivative is identically the shear
     * resultant \f$ q_r = -D\partial_r\Delta w \f$ of the closed form above, which is the
     * statement \f$ \partial_n\varphi = V_n \f$ that the free edge rests on.
     */
    Scalar analyticShearGradPotential(const GlobalPosition& globalPos) const
    {
        using std::log;
        const auto r = globalPos.two_norm();
        const auto a = innerRadius_, b = outerRadius_;
        return 0.25*force_*(r*r - a*a) - 0.5*force_*b*b*log(r/a);
    }

private:
    std::shared_ptr<CouplingManager> couplingManager_;
    Scalar force_, poissonRatio_, innerRadius_, outerRadius_, stiffness_;
    Dune::FieldVector<Scalar, 4> coefficients_;
};

} // end namespace Dumux

#endif
