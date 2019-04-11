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
 * \ingroup TODO doc me.
 * \brief TODO doc me.
 */
#ifndef DUMUX_STOKES_DARCY_MORTAR_INTERFACE_OPERATOR_HH
#define DUMUX_STOKES_DARCY_MORTAR_INTERFACE_OPERATOR_HH

#include <dune/istl/operators.hh>

namespace Dumux {

// TODO REPLACE DARCY2 WITH STOKES SOLVER
template<class DarcySolver, class MortarDarcyProjector,
         class DarcySolver2, class MortarDarcyProjector2>
class InterfaceOperator : public Dune::LinearOperator< typename MortarDarcyProjector::MortarSolution,
                                                       typename MortarDarcyProjector::MortarSolution >
{
    using MortarSolution = typename MortarDarcyProjector::MortarSolution;
    using FieldType =  typename MortarSolution::field_type;

    // make sure types of solution vectors of both projections are equal
    static constexpr bool isSameMortarSol = std::is_same<typename MortarDarcyProjector2::MortarSolution,
                                                         MortarSolution>::value;
    static_assert(isSameMortarSol, "Subdomain projectors do not map to same mortar solution type");

public:
    explicit InterfaceOperator(std::shared_ptr<DarcySolver> darcySolver,
                               std::shared_ptr<MortarDarcyProjector> projector,
                               std::shared_ptr<DarcySolver2> darcySolver2,
                               std::shared_ptr<MortarDarcyProjector2> projector2)
    : darcySolver_(darcySolver)
    , darcySolver2_(darcySolver2)
    , projector_(projector)
    , projector2_(projector2)
    {}

    virtual void apply(const MortarSolution& x, MortarSolution& r) const
    {
        darcySolver_->solve();
        darcySolver2_->solve();

        r = projector_->projectInterfacePressures();
        r -= projector2_->projectInterfacePressures();
    }

    //! apply operator to x, scale and add:  \f$ y = y + \alpha A(x) \f$
    virtual void applyscaleadd(FieldType alpha, const MortarSolution& x, MortarSolution& y) const
    {
        MortarSolution yTmp;

        apply(x, yTmp);
        yTmp *= alpha;

        y += yTmp;
    }

    //! Category of the solver (see SolverCategory::Category)
    virtual Dune::SolverCategory::Category category() const
    {
        return Dune::SolverCategory::sequential;
    }

private:
    std::shared_ptr<DarcySolver> darcySolver_;
    std::shared_ptr<DarcySolver2> darcySolver2_;

    std::shared_ptr<MortarDarcyProjector> projector_;
    std::shared_ptr<MortarDarcyProjector2> projector2_;
};

} // end namespace Dumux

#endif
