// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \brief The Laplace equation \f$ -\Delta u = 0 \f$ with the harmonic solution
 *        \f$ u = \exp(\pi \sqrt{d-1}\, x_0) \prod_{i>0} \sin(\pi x_i) \f$ as Dirichlet values on the whole boundary,
 *        assembled with local degrees of freedom (finite element and hybrid control-volume finite element schemes)
 */
#ifndef DUMUX_TEST_LINEAR_PMULTIGRID_LAPLACE_HH
#define DUMUX_TEST_LINEAR_PMULTIGRID_LAPLACE_HH

#include <cmath>
#include <memory>
#include <numbers>
#include <ranges>
#include <vector>

#include <dune/common/fvector.hh>

#include <dumux/common/boundarytypes_.hh>
#include <dumux/common/concepts/ipdata_.hh>
#include <dumux/common/constraintinfo.hh>
#include <dumux/common/indextraits.hh>
#include <dumux/common/numeqvector.hh>
#include <dumux/common/problem.hh>
#include <dumux/common/properties.hh>
#include <dumux/common/properties/model.hh>
#include <dumux/common/typetraits/localdofs_.hh>
#include <dumux/discretization/concepts.hh>
#include <dumux/discretization/cvfe/gridvariablescache_.hh>
#include <dumux/discretization/cvfe/hybrid/gridvariablescache.hh>
#include <dumux/discretization/cvfe/interpolationpointdata.hh>
#include <dumux/discretization/cvfe/quadraturerules.hh>
#include <dumux/discretization/cvfe/variablesadapter.hh>
#include <dumux/discretization/defaultlocaloperator.hh>
#include <dumux/discretization/dirichletconstraints.hh>
#include <dumux/discretization/fem/gridvariablescache.hh>
#include <dumux/discretization/gridvariables.hh>

namespace Dumux::LaplaceTest {

struct Indices
{ static constexpr int solutionIdx = 0; };

struct ModelTraits
{
    using Indices = LaplaceTest::Indices;
    static constexpr int numEq() { return 1; }
};

//! The solution at an interpolation point
template<class PV>
class Variables
{
public:
    using PrimaryVariables = PV;
    using Indices = LaplaceTest::Indices;

    template<class ElementSolution, class Problem, class ElementDiscretization, Concept::LocalDofIpData IpData>
    void update(const ElementSolution& elemSol, const Problem&, const ElementDiscretization&, const IpData& ipData)
    { priVars_ = elemSol[ipData.localDofIndex()]; }

    typename PV::value_type solution() const
    { return priVars_[0]; }

    typename PV::value_type priVar(const int pvIdx) const
    { return priVars_[pvIdx]; }

    const PrimaryVariables& priVars() const
    { return priVars_; }

    typename PV::value_type extrusionFactor() const
    { return 1.0; }

private:
    PrimaryVariables priVars_;
};

//! The flux \f$ -\nabla u \cdot \mathbf{n} \f$ over control volume faces and the Galerkin terms \f$ \nabla u \cdot \nabla w \f$
template<class TypeTag>
class LocalResidual : public DiscretizationDefaultLocalOperator<TypeTag>
{
    using ParentType = DiscretizationDefaultLocalOperator<TypeTag>;
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using Problem = GetPropType<TypeTag, Properties::Problem>;
    using GridDiscretization = GetPropType<TypeTag, Properties::GridGeometry>;
    using ElementDiscretization = typename GridDiscretization::LocalView;
    using GridView = typename GridDiscretization::GridView;
    using Element = typename GridView::template Codim<0>::Entity;
    using ElementVariables = typename GetPropType<TypeTag, Properties::GridVariables>::GridVariablesCache::LocalView;
    using NumEqVector = Dumux::NumEqVector<GetPropType<TypeTag, Properties::PrimaryVariables>>;
    using Vector = Dune::FieldVector<Scalar, GridView::dimensionworld>;

public:
    using ParentType::ParentType;
    using ElementResidualVector = typename ParentType::ElementResidualVector;

    template<class SubControlVolume>
    NumEqVector storageIntegral(const ElementDiscretization&, const ElementVariables&, const SubControlVolume&, bool) const
    { return NumEqVector(0.0); }

    template<class SubControlVolume>
    NumEqVector sourceIntegral(const ElementDiscretization&, const ElementVariables&, const SubControlVolume&) const
    { return NumEqVector(0.0); }

    template<class SubControlVolumeFace>
    NumEqVector fluxIntegral(const ElementDiscretization& elemDisc, const ElementVariables& elemVars,
                             const SubControlVolumeFace& scvf) const
    {
        NumEqVector flux(0.0);
        for (const auto& qpData : CVFE::quadratureRule(elemDisc, scvf))
            flux[0] -= qpData.weight()*(gradient_(elemDisc, elemVars, cache(elemVars, qpData.ipData()))*scvf.unitOuterNormal());
        return flux;
    }

    void addToElementFluxAndSourceResidual(ElementResidualVector& residual,
                                           const Problem&,
                                           const Element& element,
                                           const ElementDiscretization& elemDisc,
                                           const ElementVariables& elemVars) const
    {
        if constexpr (Detail::LocalDofs::hasNonCVLocalDofsInterface<ElementDiscretization>())
        {
            if (std::ranges::empty(nonCVLocalDofs(elemDisc)))
                return;

            for (const auto& qpData : CVFE::quadratureRule(elemDisc, element))
            {
                const auto& ipCache = cache(elemVars, qpData.ipData());
                const auto gradient = gradient_(elemDisc, elemVars, ipCache);
                for (const auto& localDof : nonCVLocalDofs(elemDisc))
                    residual[localDof.index()][0] += qpData.weight()*(gradient*ipCache.gradN(localDof.index()));
            }
        }
    }

private:
    template<class ShapeData>
    Vector gradient_(const ElementDiscretization& elemDisc, const ElementVariables& elemVars, const ShapeData& shapeData) const
    {
        Vector gradient(0.0);
        for (const auto& localDof : localDofs(elemDisc))
            gradient.axpy(elemVars[localDof].solution(), shapeData.gradN(localDof.index()));
        return gradient;
    }
};

//! Dirichlet values of the harmonic solution at all boundary dofs
template<class TypeTag>
class Problem : public Dumux::Experimental::Problem<TypeTag>
{
    using ParentType = Dumux::Experimental::Problem<TypeTag>;
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using GridGeometry = GetPropType<TypeTag, Properties::GridGeometry>;
    using GridView = typename GridGeometry::GridView;
    using GlobalPosition = typename GridView::template Codim<0>::Entity::Geometry::GlobalCoordinate;
    using PrimaryVariables = GetPropType<TypeTag, Properties::PrimaryVariables>;
    using BoundaryTypes = Dumux::Experimental::BoundaryTypes<PrimaryVariables::size()>;
    using ConstraintInfo = Dumux::DirichletConstraintInfo<PrimaryVariables::size()>;
    using DirichletConstraintData = Dumux::DirichletConstraintData<ConstraintInfo, PrimaryVariables,
                                                                   typename IndexTraits<GridView>::GridIndex>;

public:
    explicit Problem(std::shared_ptr<const GridGeometry> gridGeometry)
    : ParentType(gridGeometry)
    {
        appendDirichletConstraints_();
        makeConstraintsConsistent(this->gridDiscretization(), constraints_);
    }

    const std::vector<DirichletConstraintData>& constraints() const
    { return constraints_; }

    //! no flux boundaries: the boundary dofs are constrained
    BoundaryTypes boundaryTypesAtPos(const GlobalPosition&) const
    { return BoundaryTypes{}; }

    static Scalar exactSolution(const GlobalPosition& x)
    {
        using std::exp; using std::sin; using std::sqrt;
        constexpr int dim = GlobalPosition::dimension;
        constexpr Scalar pi = std::numbers::pi_v<Scalar>;
        Scalar u = exp(pi*sqrt(Scalar(dim - 1))*x[0]);
        for (int i = 1; i < dim; ++i)
            u *= sin(pi*x[i]);
        return u;
    }

private:
    void appendDirichletConstraints_()
    {
        auto elemDisc = localView(this->gridDiscretization());
        for (const auto& element : elements(this->gridDiscretization().gridView()))
        {
            elemDisc.bind(element);
            for (const auto& boundaryFace : boundaryFaces(elemDisc))
            {
                for (const auto& localDof : localDofs(elemDisc, boundaryFace))
                {
                    ConstraintInfo info;
                    info.setAll();
                    const PrimaryVariables value(exactSolution(ipData(elemDisc, localDof).global()));
                    constraints_.push_back(DirichletConstraintData{std::move(info), value, localDof.dofIndex()});
                }
            }
        }
    }

    std::vector<DirichletConstraintData> constraints_;
};

} // end namespace Dumux::LaplaceTest

namespace Dumux::Properties {

namespace TTag {
struct LaplaceTest { using InheritsFrom = std::tuple<ModelProperties>; };
} // end namespace TTag

template<class TypeTag>
struct ModelTraits<TypeTag, TTag::LaplaceTest>
{ using type = LaplaceTest::ModelTraits; };

template<class TypeTag>
struct LocalResidual<TypeTag, TTag::LaplaceTest>
{ using type = LaplaceTest::LocalResidual<TypeTag>; };

template<class TypeTag>
struct VolumeVariables<TypeTag, TTag::LaplaceTest>
{ using type = LaplaceTest::Variables<GetPropType<TypeTag, Properties::PrimaryVariables>>; };

template<class TypeTag>
struct Problem<TypeTag, TTag::LaplaceTest>
{ using type = LaplaceTest::Problem<TypeTag>; };

template<class TypeTag>
struct EnableGridGeometryCache<TypeTag, TTag::LaplaceTest> { static constexpr bool value = true; };

template<class TypeTag>
struct GridVariables<TypeTag, TTag::LaplaceTest>
{
private:
    using GG = GetPropType<TypeTag, Properties::GridGeometry>;
    using ElementDiscretization = typename GG::LocalView;
    using Problem = GetPropType<TypeTag, Properties::Problem>;
    using Variables = Dumux::Detail::CVFE::VariablesAdapter<GetPropType<TypeTag, Properties::VolumeVariables>>;
    using IPDataCache = Dumux::CVFE::LocalBasisInterpolationPointData<GG>;
    static constexpr bool enableCache = getPropValue<TypeTag, Properties::EnableGridVolumeVariablesCache>();

    using FECache = Dumux::Experimental::FEGridVariablesCache<
        Dumux::Experimental::FEDefaultGridVariablesCacheTraits<Problem, Variables, IPDataCache>, enableCache>;
    using HybridCache = Dumux::Experimental::CVFE::HybridCVFEGridVariablesCache<
        Dumux::Experimental::CVFE::HybridCVFEDefaultGridVariablesCacheTraits<Problem, Variables, IPDataCache>, enableCache>;
    using CVFECache = Dumux::Experimental::CVFE::CVFEGridVariablesCache<
        Dumux::Experimental::CVFE::CVFEDefaultGridVariablesCacheTraits<Problem, Variables, IPDataCache>, enableCache>;
    using Cache = std::conditional_t<Dumux::Experimental::Concepts::FEElementDiscretization<ElementDiscretization>, FECache,
                  std::conditional_t<Dumux::Experimental::Concepts::HybridElementDiscretization<ElementDiscretization>, HybridCache,
                                     CVFECache>>;
public:
    using type = Dumux::Experimental::GridVariables<GG, Cache>;
};

} // end namespace Dumux::Properties

#endif
