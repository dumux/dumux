// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup MortarCoupling
 * \brief Predefined projector implementations for mortar-coupling models, and the function
 *        spaces on mortar domains and traces they are built on.
 */
#ifndef DUMUX_MULTIDOMAIN_MORTAR_PROJECTORS_HH
#define DUMUX_MULTIDOMAIN_MORTAR_PROJECTORS_HH

#include <cstddef>
#include <memory>
#include <type_traits>
#include <utility>
#include <vector>

#include <dune/common/exceptions.hh>
#include <dune/common/fmatrix.hh>
#include <dune/common/fvector.hh>
#include <dune/geometry/quadraturerules.hh>
#include <dune/istl/bcrsmatrix.hh>

#include <dumux/common/concepts/entityset_.hh>
#include <dumux/common/concepts/mortartrace_.hh>
#include <dumux/assembly/jacobianpattern.hh>
#include <dumux/geometry/boundingboxtree.hh>
#include <dumux/geometry/geometricentityset.hh>
#include <dumux/geometry/intersectionentityset.hh>

#include <dumux/discretization/method.hh>
#include <dumux/discretization/box/fvgridgeometry.hh>
#include <dumux/discretization/projection/projector.hh>
#include <dumux/discretization/projection/brokenbasis.hh>
#include <dumux/discretization/projection/l2_projection.hh>

#include "projectorinterface.hh"

namespace Dumux::Mortar {

#ifndef DOXYGEN
namespace Detail {

/*!
 * \brief Holds a piecewise constant function space basis over an entity set.
 */
template<Concept::GeometricEntitySet EntitySet>
class BrokenTraceSpace
{
public:
    using Basis = BrokenLagrangeBasis<EntitySet, 0>;

    explicit BrokenTraceSpace(std::shared_ptr<const EntitySet> entitySet)
    : basis_(std::move(entitySet)) {}

    const Basis& basis() const { return basis_; }

private:
    Basis basis_;
};

/*!
 * \brief Holds a function space basis of the requested order over a grid view.
 *
 * At order zero the space is discontinuous, so a basis over a bare entity set suffices. At
 * higher order it is continuous, which needs the adjacency a grid provides; a control-volume
 * finite element geometry over the same grid view supplies it.
 */
template<class GridView, std::size_t order>
class TraceSpace;

template<class GridView>
class TraceSpace<GridView, 0>
: public BrokenTraceSpace<GridViewGeometricEntitySet<GridView>>
{
    using EntitySet = GridViewGeometricEntitySet<GridView>;

public:
    explicit TraceSpace(const GridView& gridView)
    : BrokenTraceSpace<EntitySet>(std::make_shared<EntitySet>(gridView)) {}
};

template<class GridView>
class TraceSpace<GridView, 1>
{
    using GridGeometry = BoxFVGridGeometry<typename GridView::ctype, GridView>;

public:
    using Basis = FEBasisFromCVFEGridDiscretization<GridGeometry>;

    explicit TraceSpace(const GridView& gridView)
    : gridGeometry_(std::make_shared<GridGeometry>(gridView)), basis_(*gridGeometry_) {}

    const Basis& basis() const { return basis_; }

private:
    std::shared_ptr<GridGeometry> gridGeometry_;
    Basis basis_;
};

/*!
 * \brief The function space of the requested order on a mortar trace: piecewise constants
 *        on any trace, continuous spaces only on a trace that is a grid.
 */
template<std::size_t order, Concept::MortarTrace Trace>
auto makeTraceSpace(const Trace& trace)
{
    if constexpr (order == 0)
        return BrokenTraceSpace<typename Trace::EntitySet>{trace.entitySet()};
    else
    {
        static_assert(Concept::GridMortarTrace<Trace>,
                      "A continuous trace space needs a trace whose cells form a grid");
        return TraceSpace<typename Trace::GridView, order>{trace.gridView()};
    }
}

//! The order of the function space on a mortar domain, that of the mortar's own discretization
template<class MortarGridGeometry>
inline constexpr std::size_t mortarSpaceOrder
    = MortarGridGeometry::discMethod == DiscretizationMethods::box ? 1 : 0;

//! The function space on a mortar domain
template<class MortarGridGeometry>
using MortarSpace = TraceSpace<typename MortarGridGeometry::GridView, mortarSpaceOrder<MortarGridGeometry>>;

template<class MortarGridGeometry>
MortarSpace<MortarGridGeometry> makeMortarSpace(const MortarGridGeometry& mortarGridGeometry)
{ return MortarSpace<MortarGridGeometry>{mortarGridGeometry.gridView()}; }

//! Visit the entities a basis is defined on, whether over a grid view or a bare entity set
template<class Basis, class Visitor>
void forEachBasisEntity(const Basis& basis, Visitor&& visit)
{
    if constexpr (requires { basis.gridView(); })
        for (const auto& element : elements(basis.gridView()))
            visit(element);
    else
        for (const auto& entity : basis.entitySet())
            visit(entity);
}

/*!
 * \brief Visit the quadrature points of a scalar function space, with the local view bound
 *        to the entity, the shape function values at the point and the weight including the
 *        integration element. The rule integrates products of two shape functions exactly.
 */
template<class Scalar, class Basis, class Visitor>
void forEachQuadraturePoint(const Basis& basis, Visitor&& visit)
{
    auto localView = basis.localView();
    std::vector<Dune::FieldVector<Scalar, 1>> shapeValues;
    forEachBasisEntity(basis, [&] (const auto& entity) {
        localView.bind(entity);
        const auto& localBasis = localView.tree().finiteElement().localBasis();
        const auto geometry = entity.geometry();
        using Geometry = std::decay_t<decltype(geometry)>;
        const auto& quad = Dune::QuadratureRules<Scalar, Geometry::mydimension>::rule(
            geometry.type(), 2*localBasis.order()
        );
        for (const auto& qp : quad)
        {
            localBasis.evaluateFunction(qp.position(), shapeValues);
            visit(localView, shapeValues, qp.weight()*geometry.integrationElement(qp.position()));
        }
    });
}

//! The mass matrix of the function space on a mortar domain
template<class Scalar, class MortarGridGeometry>
Dune::BCRSMatrix<Dune::FieldMatrix<Scalar, 1, 1>> assembleMortarMassMatrix(const MortarGridGeometry& mortarGridGeometry)
{
    const auto space = makeMortarSpace(mortarGridGeometry);
    Dune::BCRSMatrix<Dune::FieldMatrix<Scalar, 1, 1>> mass;
    getFEJacobianPattern(space.basis()).exportIdx(mass);
    mass = 0.0;
    forEachQuadraturePoint<Scalar>(space.basis(), [&] (const auto& localView, const auto& shapeValues, Scalar weight) {
        for (std::size_t i = 0; i < shapeValues.size(); ++i)
            for (std::size_t j = 0; j < shapeValues.size(); ++j)
                mass[localView.index(i)][localView.index(j)] += weight*shapeValues[i][0]*shapeValues[j][0];
    });
    return mass;
}

//! The integrals of the basis functions of the function space on a mortar domain
template<class Scalar, class MortarGridGeometry>
std::vector<Scalar> mortarBasisIntegrals(const MortarGridGeometry& mortarGridGeometry)
{
    const auto space = makeMortarSpace(mortarGridGeometry);
    std::vector<Scalar> integrals(mortarGridGeometry.numDofs(), 0.0);
    forEachQuadraturePoint<Scalar>(space.basis(), [&] (const auto& localView, const auto& shapeValues, Scalar weight) {
        for (std::size_t i = 0; i < shapeValues.size(); ++i)
            integrals[localView.index(i)] += weight*shapeValues[i][0];
    });
    return integrals;
}

} // end namespace Detail
#endif // DOXYGEN

/*!
 * \ingroup MortarCoupling
 * \brief Default projector for finite-volume schemes and mortars representing primary
 *        variables (e.g. pressure mortars). Projects between a sub-domain trace and the
 *        mortar domain.
 */
template<typename SolutionVector>
class FVDefaultProjector : public Projector<SolutionVector>
{
    using Scalar = typename SolutionVector::field_type;
    using L2Projector = Dumux::Projector<Scalar>;

public:
    /*!
     * \brief Projectors between a mortar domain and the trace of a subdomain on it.
     * \param mortarGridGeometry The grid geometry of the mortar domain
     * \param trace The trace of the subdomain on the mortar domain
     * \param traceDataOrder The order of the data imposed on the trace: zero for one value
     *        per trace cell, one for one value per trace vertex, which needs a trace whose
     *        cells form a grid
     */
    template<typename MortarGridGeometry, Concept::MortarTrace Trace>
    FVDefaultProjector(const MortarGridGeometry& mortarGridGeometry,
                       const Trace& trace,
                       std::size_t traceDataOrder)
    {
        using MortarEntitySet = typename std::remove_cvref_t<decltype(mortarGridGeometry.boundingBoxTree())>::EntitySet;
        IntersectionEntitySet<MortarEntitySet, typename Trace::EntitySet> glue;
        glue.build(mortarGridGeometry.boundingBoxTree(), trace.boundingBoxTree());

        const auto mortarSpace = Detail::makeMortarSpace(mortarGridGeometry);
        const auto residualTraceSpace = Detail::makeTraceSpace<0>(trace);

        // Imposing the mortar on the trace is an L2-projection: the mortar is a function there.
        const auto makeImposition = [&] (const auto& traceSpace) {
            to_ = std::make_unique<L2Projector>(makeProjector(mortarSpace.basis(), traceSpace.basis(), glue));
        };
        if (traceDataOrder == 0)
            makeImposition(residualTraceSpace);
        else if (traceDataOrder == 1)
        {
            if constexpr (Concept::GridMortarTrace<Trace>)
                makeImposition(Detail::makeTraceSpace<1>(trace));
            else
                DUNE_THROW(Dune::NotImplemented, "Trace data per vertex on a trace whose cells form no grid");
        }
        else
            DUNE_THROW(Dune::NotImplemented, "Trace data of order " << traceDataOrder);

        // The residual is not. It is a functional, obtained by testing the trace against the
        // mortar basis, so it is the projection matrix alone and carries no mass inverse.
        // Applying one would leave the interface operator symmetric only in the mortar mass
        // inner product, which the conjugate gradient method driving it does not use.
        from_ = makeProjectionMatricesPair(mortarSpace.basis(), residualTraceSpace.basis(), glue).second.second;
    }

private:
    SolutionVector toTrace_(const SolutionVector& x) const override { return to_->project(x); }

    SolutionVector fromTrace_(const SolutionVector& x) const override
    {
        SolutionVector result(from_.N());
        result = 0.0;
        from_.mv(x, result);
        return result;
    }

    std::unique_ptr<L2Projector> to_;
    typename L2Projector::Matrix from_;
};

} // end namespace Dumux::Mortar

#endif
