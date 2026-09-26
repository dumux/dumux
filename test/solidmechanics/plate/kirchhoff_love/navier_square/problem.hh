// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
#ifndef DUMUX_KIRCHHOFF_LOVE_PLATE_NAVIER_SQUARE_TEST_PROBLEM_HH
#define DUMUX_KIRCHHOFF_LOVE_PLATE_NAVIER_SQUARE_TEST_PROBLEM_HH

#include <bitset>
#include <cmath>
#include <limits>
#include <string>

#include <dune/common/exceptions.hh>

#include <dumux/common/boundarytypes.hh>
#include <dumux/common/fvproblem.hh>
#include <dumux/common/numeqvector.hh>
#include <dumux/common/parameters.hh>
#include <dumux/common/properties.hh>
#include <dumux/common/math.hh>

namespace Dumux {

enum class SimplySupportedMapping { traction, tangential };

inline SimplySupportedMapping simplySupportedMapping()
{
    const auto name = getParam<std::string>("Problem.Mapping");
    if (name == "traction")
        return SimplySupportedMapping::traction;
    if (name == "tangential")
        return SimplySupportedMapping::tangential;
    DUNE_THROW(Dune::InvalidStateException, "Unknown mapping " << name);
}

/*!
 * \brief Rotation sub-problem of the simply supported square plate
 */
template<class TypeTag>
class NavierSquareProblemRotation : public FVProblem<TypeTag>
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
    NavierSquareProblemRotation(std::shared_ptr<const GridGeometry> gridGeometry,
                                std::shared_ptr<CouplingManager> couplingManager)
    : ParentType(gridGeometry)
    , couplingManager_(couplingManager)
    , mapping_(simplySupportedMapping())
    , poissonRatio_(getParam<Scalar>("Problem.PoissonRatio"))
    , length_(getParam<Scalar>("Problem.Length", 1.0))
    , eps_(1e-7*length_)
    {
        const auto E = getParam<Scalar>("Problem.E");
        const auto t = getParam<Scalar>("Problem.Thickness");
        stiffness_ = E*t*t*t/(12.0*(1.0 - poissonRatio_*poissonRatio_));
    }

    Scalar D(const GlobalPosition& globalPos) const { return stiffness_; }
    Scalar poissonRatio(const GlobalPosition& globalPos) const { return poissonRatio_; }

    //! On the edges of the square the tangential rotation is one Cartesian component
    BoundaryTypes boundaryTypes(const Element& element, const SubControlVolume& scv) const
    {
        BoundaryTypes values;
        values.setAllNeumann();
        if (mapping_ == SimplySupportedMapping::tangential)
        {
            const auto& pos = scv.dofPosition();
            if (pos[1] < eps_ || pos[1] > length_ - eps_)
                values.setDirichlet(Indices::rotation0Idx);
            if (pos[0] < eps_ || pos[0] > length_ - eps_)
                values.setDirichlet(Indices::rotation1Idx);
        }
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

    const CouplingManager& couplingManager() const { return *couplingManager_; }

private:
    std::shared_ptr<CouplingManager> couplingManager_;
    SimplySupportedMapping mapping_;
    Scalar poissonRatio_, length_, eps_, stiffness_;
};

/*!
 * \brief Deformation sub-problem of the simply supported square plate
 */
template<class TypeTag>
class NavierSquareProblemDeformation : public FVProblem<TypeTag>
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
    NavierSquareProblemDeformation(std::shared_ptr<const GridGeometry> gridGeometry,
                                   std::shared_ptr<CouplingManager> couplingManager)
    : ParentType(gridGeometry)
    , couplingManager_(couplingManager)
    , mapping_(simplySupportedMapping())
    , force_(getParam<Scalar>("Problem.Force"))
    , length_(getParam<Scalar>("Problem.Length", 1.0))
    , eps_(1e-7*length_)
    {
        GlobalPosition edgeMidpoint(0.0), centre(0.5*length_);
        edgeMidpoint[0] = 0.5*length_;
        Scalar edgeDistance = std::numeric_limits<Scalar>::max();
        Scalar centreDistance = std::numeric_limits<Scalar>::max();
        for (const auto& vertex : vertices(this->gridGeometry().gridView()))
        {
            const auto pos = vertex.geometry().center();
            const auto idx = this->gridGeometry().vertexMapper().index(vertex);
            if ((pos - edgeMidpoint).two_norm() < edgeDistance)
            { edgeDistance = (pos - edgeMidpoint).two_norm(); edgePinDof_ = idx; }
            if ((pos - centre).two_norm() < centreDistance)
            { centreDistance = (pos - centre).two_norm(); centrePinDof_ = idx; }
        }
    }

    /*!
     * \brief The support prescribes \f$ w = 0 \f$ in place of the transverse balance
     *
     * The deflection takes the row of the transverse balance, whose boundary flux becomes
     * the reaction, and the potential \f$ \varphi \f$ the row of the compatibility equation.
     */
    BoundaryTypes boundaryTypes(const Element& element, const SubControlVolume& scv) const
    {
        BoundaryTypes values;
        values.setAllNeumann();
        values.setDirichlet(Indices::verticalDeformationIdx, Indices::shearGradPotentialEqIdx);
        if (mapping_ == SimplySupportedMapping::tangential || scv.dofIndex() == edgePinDof_)
            values.setDirichlet(Indices::shearGradPotentialIdx, Indices::deformationEqIdx);
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

    //! With the tangential rotation prescribed on every edge nothing fixes the constant of \f$ \psi \f$
    std::bitset<NumEqVector::dimension>
    hasInternalDirichletConstraint(const Element& element, const SubControlVolume& scv) const
    {
        std::bitset<NumEqVector::dimension> values;
        if (mapping_ == SimplySupportedMapping::tangential && scv.dofIndex() == centrePinDof_)
            values.set(Indices::shearCurlPotentialIdx);
        return values;
    }

    PrimaryVariables internalDirichlet(const Element& element, const SubControlVolume& scv) const
    { return PrimaryVariables(0.0); }

    PrimaryVariables initialAtPos(const GlobalPosition& globalPos) const
    { return PrimaryVariables(0.0); }

    const CouplingManager& couplingManager() const { return *couplingManager_; }

private:
    std::shared_ptr<CouplingManager> couplingManager_;
    SimplySupportedMapping mapping_;
    Scalar force_, length_, eps_;
    std::size_t edgePinDof_ = 0, centrePinDof_ = 0;
};

} // end namespace Dumux

#endif
