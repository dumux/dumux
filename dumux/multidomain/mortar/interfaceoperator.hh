// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*****************************************************************************
 *   See the file COPYING for full copying permissions.                      *
 *                                                                           *
 *   This program is free software: you can redistribute it and/or modify    *
 *   it under the terms of the GNU General Public License as published by    *
 *   the Free Software Foundation, either version 3 of the License, or       *
 *   (at your option) any later version.                                     *
 *                                                                           *
 *   This program is distributed in the hope that it will be useful,         *
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of          *
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the            *
 *   GNU General Public License for more details.                            *
 *                                                                           *
 *   You should have received a copy of the GNU General Public License       *
 *   along with this program.  If not, see <http://www.gnu.org/licenses/>.   *
 *****************************************************************************/
/*!
 * \file
 * \ingroup MultiDomain
 * \ingroup MortarCoupling
 * \brief Linear operator for sequentially solving mortar models.
 */
#ifndef DUMUX_MORTAR_INTERFACE_OPERATOR_HH
#define DUMUX_MORTAR_INTERFACE_OPERATOR_HH

#include <memory>
#include <dune/istl/operators.hh>

#include "model.hh"

namespace Dumux::Mortar {

/*!
 * \ingroup MultiDomain
 * \ingroup MortarCoupling
 * \brief Linear operator for sequentially solving mortar models.
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

    explicit InterfaceOperator(Model&& model)
    : model_{std::make_shared<M>(std::move(model))}
    {}

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
