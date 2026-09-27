// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup MortarCoupling
 * \brief Linear operator for sequentially solving mortar models.
 */
#ifndef DUMUX_MULTIDOMAIN_MORTAR_INTERFACE_OPERATOR_HH
#define DUMUX_MULTIDOMAIN_MORTAR_INTERFACE_OPERATOR_HH

#include <memory>
#include <utility>
#include <dune/istl/operators.hh>

#include "model.hh"

namespace Dumux::Mortar {

/*!
 * \ingroup MortarCoupling
 * \brief Linear operator for sequentially solving mortar models.
 *
 * Realizes the interface operator of the non-overlapping domain-decomposition algorithm of
 * \cite Boon2023. One application imposes the mortar data on all subdomains, solves them
 * independently, and assembles the jump of the conjugate traces into the residual.
 *
 * \tparam M The mortar model (see Dumux::Mortar::Model)
 */
template<typename M>
class InterfaceOperator
: public Dune::LinearOperator<typename M::SolutionVector, typename M::SolutionVector>
{
    using FieldType =  typename M::SolutionVector::field_type;

public:
    using Model = M;
    using SolutionVector = typename M::SolutionVector;

    //! The interface operator of the given model, which it takes ownership of
    explicit InterfaceOperator(Model&& model)
    : model_{std::make_shared<M>(std::move(model))}
    {}

    //! The interface operator of the given model, shared with other users such as a preconditioner
    explicit InterfaceOperator(std::shared_ptr<Model> model)
    : model_{std::move(model)}
    {}

    //! apply operator to x:  \f$ y = A(x) \f$
    virtual void apply(const SolutionVector& x, SolutionVector& r) const
    {
        r = 0.0;
        model_->setMortar(x);
        model_->solveSubDomains();
        model_->assembleMortarResidual(r);
    }

    //! apply operator to x, scale and add:  \f$ y = y + \alpha A(x) \f$
    virtual void applyscaleadd(FieldType alpha, const SolutionVector& x, SolutionVector& y) const
    {
        SolutionVector yTmp;

        apply(x, yTmp);
        yTmp *= alpha;

        y += yTmp;
    }

    //! Category of the solver (see SolverCategory::Category)
    virtual Dune::SolverCategory::Category category() const
    { return Dune::SolverCategory::sequential; }

private:
    std::shared_ptr<Model> model_;
};

template<typename MortarSolutionVector,
         typename MortarGridGeometry,
         typename... SubDomainGridGeometries>
InterfaceOperator(Model<MortarSolutionVector, MortarGridGeometry, SubDomainGridGeometries...>&&)
-> InterfaceOperator<Model<MortarSolutionVector, MortarGridGeometry, SubDomainGridGeometries...>>;

template<typename MortarSolutionVector,
         typename MortarGridGeometry,
         typename... SubDomainGridGeometries>
InterfaceOperator(std::shared_ptr<Model<MortarSolutionVector, MortarGridGeometry, SubDomainGridGeometries...>>)
-> InterfaceOperator<Model<MortarSolutionVector, MortarGridGeometry, SubDomainGridGeometries...>>;

}  // namespace Dumux::Mortar

#endif
