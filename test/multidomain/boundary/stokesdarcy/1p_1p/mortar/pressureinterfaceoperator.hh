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
#ifndef DUMUX_MORTAR_PRESSURE_INTERFACE_OPERATOR_HH
#define DUMUX_MORTAR_PRESSURE_INTERFACE_OPERATOR_HH

#include <dune/istl/operators.hh>

namespace Dumux {

template<class Solver1, class Projector1,
         class Solver2, class Projector2>
class PressureInterfaceOperator
: public Dune::LinearOperator< typename Projector1::MortarSolution,
                               typename Projector2::MortarSolution >
{
    using MortarSolution = typename Projector1::MortarSolution;
    using FieldType =  typename MortarSolution::field_type;

    // make sure types of solution vectors of both projections are equal
    static constexpr bool isSameMortarSol = std::is_same<typename Projector2::MortarSolution, MortarSolution>::value;
    static_assert(isSameMortarSol, "Subdomain projectors do not map to same mortar solution type");

public:
    explicit PressureInterfaceOperator(std::shared_ptr<Solver1> solver1,
                                       std::shared_ptr<Projector1> projector1,
                                       std::shared_ptr<Solver2> solver2,
                                       std::shared_ptr<Projector2> projector2)
    : solver1_(solver1)
    , solver2_(solver2)
    , projector1_(projector1)
    , projector2_(projector2)
    {}

    //! apply operator to x:  \f$ y = A(x) \f$
    virtual void apply(const MortarSolution& x, MortarSolution& r) const
    {
        // turn into shared_ptr to set it in projectors
        std::shared_ptr<MortarSolution> p = std::make_shared<MortarSolution>(x);
        projector1_->setMortarSolutionPointer(p);
        projector2_->setMortarSolutionPointer(p);

        solver1_->solve();
        solver2_->solve();

        r = projector1_->projectInterfaceFluxes();
        r += projector2_->projectInterfaceFluxes();
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
    { return Dune::SolverCategory::sequential; }

private:
    std::shared_ptr<Solver1> solver1_;
    std::shared_ptr<Solver2> solver2_;

    std::shared_ptr<Projector1> projector1_;
    std::shared_ptr<Projector2> projector2_;
};

} // end namespace Dumux

#endif
