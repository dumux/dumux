// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup AssemblyTests
 * \brief Heat equation \f$ u_t - \Delta u = 0 \f$ assembled with local degrees of freedom
 *        on a grid that is periodic in both directions
 */
#ifndef DUMUX_TEST_ASSEMBLY_MULTISTAGE_PERIODIC_PROPERTIES_HH
#define DUMUX_TEST_ASSEMBLY_MULTISTAGE_PERIODIC_PROPERTIES_HH

#include <ranges>
#include <tuple>
#include <vector>

#include <dune/common/fvector.hh>
#include <dune/grid/spgrid.hh>

#include <dumux/common/boundarytypes_.hh>
#include <dumux/common/concepts/ipdata_.hh>
#include <dumux/common/numeqvector.hh>
#include <dumux/common/problem.hh>
#include <dumux/common/properties.hh>
#include <dumux/common/properties/model.hh>
#include <dumux/common/typetraits/localdofs_.hh>
#include <dumux/discretization/box.hh>
#include <dumux/discretization/cvfe/quadraturerules.hh>
#include <dumux/discretization/defaultlocaloperator.hh>
#include <dumux/discretization/pq2.hh>
#include <dumux/io/grid/gridmanager_sp.hh>

namespace Dumux {

struct HeatModelTraits
{
    struct Indices {};
    static constexpr int numEq() { return 1; }
};

//! The temperature at a local degree of freedom
template<class PV>
class HeatVariables
{
public:
    using PrimaryVariables = PV;

    template<class ElementSolution, class Problem, class ElementDiscretization, Concept::LocalDofIpData IpData>
    void update(const ElementSolution& elemSol, const Problem&, const ElementDiscretization&, const IpData& ipData)
    { priVars_ = elemSol[ipData.localDofIndex()]; }

    double temperature() const { return priVars_[0]; }
    double priVar(int pvIdx) const { return priVars_[pvIdx]; }
    const PrimaryVariables& priVars() const { return priVars_; }
    double extrusionFactor() const { return 1.0; }

private:
    PrimaryVariables priVars_;
};

/*!
 * \brief Control-volume balance of the heat equation for the degrees of freedom with a control volume,
 *        Galerkin form for the others, both with a lumped storage
 */
template<class TypeTag>
class HeatLocalResidual : public DiscretizationDefaultLocalOperator<TypeTag>
{
    using ParentType = DiscretizationDefaultLocalOperator<TypeTag>;
    using Problem = GetPropType<TypeTag, Properties::Problem>;
    using NumEqVector = Dumux::NumEqVector<GetPropType<TypeTag, Properties::PrimaryVariables>>;
    using ElementDiscretization = typename GetPropType<TypeTag, Properties::GridGeometry>::LocalView;
    using Element = typename ElementDiscretization::Element;
    using GlobalPosition = typename Element::Geometry::GlobalCoordinate;

public:
    using ElementResidualVector = typename ParentType::ElementResidualVector;
    using ParentType::ParentType;

    template<class ElementVariables, class Scv>
    NumEqVector storageIntegral(const ElementDiscretization&, const ElementVariables& elemVars, const Scv& scv, bool) const
    { return NumEqVector(elemVars[scv].temperature()*scv.volume()); }

    template<class ElementVariables, class Scvf>
    NumEqVector fluxIntegral(const ElementDiscretization& elemDisc, const ElementVariables& elemVars, const Scvf& scvf) const
    {
        NumEqVector flux(0.0);
        for (const auto& qpData : CVFE::quadratureRule(elemDisc, scvf))
            flux[0] -= qpData.weight()*(gradient_(elemDisc, elemVars, qpData.ipData())*scvf.unitOuterNormal());
        return flux;
    }

    template<class ElementVariables>
    void addToElementStorage(ElementResidualVector& storage, const Problem&, const Element& element,
                             const ElementDiscretization& elemDisc, const ElementVariables& elemVars, bool) const
    {
        if constexpr (Detail::LocalDofs::hasNonCVLocalDofsInterface<ElementDiscretization>())
        {
            if (std::ranges::empty(nonCVLocalDofs(elemDisc)))
                return;

            for (const auto& qpData : CVFE::quadratureRule(elemDisc, element))
            {
                const auto& shapeValues = cache(elemVars, qpData.ipData()).shapeValues();
                for (const auto& localDof : nonCVLocalDofs(elemDisc))
                    storage[localDof.index()][0] += qpData.weight()*shapeValues[localDof.index()][0]*elemVars[localDof].temperature();
            }
        }
    }

    template<class ElementVariables>
    void addToElementFluxAndSourceResidual(ElementResidualVector& residual, const Problem&, const Element& element,
                                           const ElementDiscretization& elemDisc, const ElementVariables& elemVars) const
    {
        if constexpr (Detail::LocalDofs::hasNonCVLocalDofsInterface<ElementDiscretization>())
        {
            if (std::ranges::empty(nonCVLocalDofs(elemDisc)))
                return;

            for (const auto& qpData : CVFE::quadratureRule(elemDisc, element))
            {
                const auto& ipCache = cache(elemVars, qpData.ipData());
                const auto gradU = gradient_(elemDisc, elemVars, qpData.ipData());
                for (const auto& localDof : nonCVLocalDofs(elemDisc))
                    residual[localDof.index()][0] += qpData.weight()*(gradU*ipCache.gradN(localDof.index()));
            }
        }
    }

private:
    template<class ElementVariables, class IpData>
    GlobalPosition gradient_(const ElementDiscretization& elemDisc, const ElementVariables& elemVars, const IpData& ipData) const
    {
        const auto& ipCache = cache(elemVars, ipData);
        GlobalPosition gradU(0.0);
        for (const auto& localDof : localDofs(elemDisc))
            gradU.axpy(elemVars[localDof].temperature(), ipCache.gradN(localDof.index()));
        return gradU;
    }
};

//! No boundaries on a grid that is periodic in all directions
template<class TypeTag>
class HeatProblem : public Experimental::Problem<TypeTag>
{
    using ParentType = Experimental::Problem<TypeTag>;
    using PrimaryVariables = GetPropType<TypeTag, Properties::PrimaryVariables>;
    using GlobalPosition = typename GetPropType<TypeTag, Properties::GridGeometry>::GridView::template Codim<0>::Entity::Geometry::GlobalCoordinate;

public:
    using ParentType::ParentType;

    Dumux::Experimental::BoundaryTypes<1> boundaryTypesAtPos(const GlobalPosition&) const
    { return {}; }
};

namespace Properties {

namespace TTag {
struct Heat { using InheritsFrom = std::tuple<ModelProperties>; };
struct HeatBox { using InheritsFrom = std::tuple<Heat, BoxModel>; };
struct HeatPQ2 { using InheritsFrom = std::tuple<Heat, PQ2HybridModel>; };
} // end namespace TTag

template<class TypeTag>
struct ModelTraits<TypeTag, TTag::Heat> { using type = HeatModelTraits; };

template<class TypeTag>
struct LocalResidual<TypeTag, TTag::Heat> { using type = HeatLocalResidual<TypeTag>; };

template<class TypeTag>
struct Variables<TypeTag, TTag::Heat> { using type = HeatVariables<GetPropType<TypeTag, Properties::PrimaryVariables>>; };

template<class TypeTag>
struct VolumeVariables<TypeTag, TTag::Heat> { using type = GetPropType<TypeTag, Properties::Variables>; };

template<class TypeTag>
struct Grid<TypeTag, TTag::Heat> { using type = Dune::SPGrid<double, 2>; };

template<class TypeTag>
struct Problem<TypeTag, TTag::Heat> { using type = HeatProblem<TypeTag>; };

template<class TypeTag>
struct Transition::UseVariables<TypeTag, TTag::HeatBox> { static constexpr bool value = true; };

} // end namespace Properties
} // end namespace Dumux

#endif
