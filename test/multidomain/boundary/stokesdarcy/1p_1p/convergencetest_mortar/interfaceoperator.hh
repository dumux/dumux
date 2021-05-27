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
#ifndef DUMUX_MORTAR_ONEP_INTERFACE_OPERATOR_HH
#define DUMUX_MORTAR_ONEP_INTERFACE_OPERATOR_HH

#include <dune/istl/operators.hh>

#include "mortarvariabletype.hh"
#include "reconstructionhelper.hh"

namespace Dumux {

template<class Solver1, class Reconstructor1,
         class Solver2, class Reconstructor2,
         class MortarSolutionVector>
class OnePMortarInterfaceOperator
: public Dune::LinearOperator< MortarSolutionVector,
                               MortarSolutionVector >
{
    using FieldType =  typename MortarSolutionVector::field_type;
    using Projector = MortarProjectorBase<MortarSolutionVector>;

    using GridGeometry1 = typename Solver1::FVGridGeometry;
    using GridGeometry2 = typename Solver2::FVGridGeometry;
    using CoupledScvfMap1 = typename MortarReconstructionHelper::ElementScvfIndexMap<GridGeometry1>;
    using CoupledScvfMap2 = typename MortarReconstructionHelper::ElementScvfIndexMap<GridGeometry2>;

public:
    template<class MortarGridGeometry>
    explicit OnePMortarInterfaceOperator(std::shared_ptr<Solver1> solver1, std::shared_ptr<Projector> projector1,
                                         std::shared_ptr<Solver2> solver2, std::shared_ptr<Projector> projector2,
                                         const MortarGridGeometry& mortarGG, OnePMortarVariableType mv)
    : solver1_(solver1)
    , solver2_(solver2)
    , projector1_(projector1)
    , projector2_(projector2)
    , variableType_(mv)
    {
        using Helper = MortarReconstructionHelper;
        coupledScvfMap1_ = Helper::findCoupledScvfs(*solver1->gridGeometryPointer(), mortarGG);
        coupledScvfMap2_ = Helper::findCoupledScvfs(*solver2->gridGeometryPointer(), mortarGG);
    }

    //! apply operator to x:  \f$ y = A(x) \f$
    virtual void apply(const MortarSolutionVector& x, MortarSolutionVector& r) const
    {
        // project mortar pressure into sub-domains
        const auto p1 = projector1_->projectMortarToSubDomain(x);
        const auto p2 = projector2_->projectMortarToSubDomain(x);
        solver1_->problemPointer()->setMortarProjection(p1);
        solver2_->problemPointer()->setMortarProjection(p2);

        solver1_->solve();
        solver2_->solve();

        using R1 = Reconstructor1;
        using R2 = Reconstructor2;

        // compute fluxes in sub-domains
        if (variableType_ == OnePMortarVariableType::pressure)
        {
            const auto flux1 = R1::template recoverNormalFlux<MortarSolutionVector>(*solver1_->gridGeometryPointer(),
                                                                                    *solver1_->gridVariablesPointer(),
                                                                                    *solver1_->solutionPointer(),
                                                                                    coupledScvfMap1_);

            const auto flux2 = R2::template recoverNormalFlux<MortarSolutionVector>(*solver2_->gridGeometryPointer(),
                                                                                    *solver2_->gridVariablesPointer(),
                                                                                    *solver2_->solutionPointer(),
                                                                                    coupledScvfMap2_);

            r = projector1_->projectSubDomainToMortar(flux1);
            r += projector2_->projectSubDomainToMortar(flux2);
        }
        else if (variableType_ == OnePMortarVariableType::flux)
        {
            auto pressure1 = R1::template recoverPressure<MortarSolutionVector>(*solver1_->gridGeometryPointer(),
                                                                                *solver1_->gridVariablesPointer(),
                                                                                *solver1_->solutionPointer(),
                                                                                coupledScvfMap1_);

            auto pressure2 = R2::template recoverPressure<MortarSolutionVector>(*solver2_->gridGeometryPointer(),
                                                                                *solver2_->gridVariablesPointer(),
                                                                                *solver2_->solutionPointer(),
                                                                                coupledScvfMap2_);

            if (solver1_->problemPointer()->isOnNegativeMortarSide()) pressure1 *= -1.0;
            if (solver2_->problemPointer()->isOnNegativeMortarSide()) pressure2 *= -1.0;

            r = projector1_->projectSubDomainToMortar(pressure1);
            r += projector2_->projectSubDomainToMortar(pressure2);
        }
        else
            DUNE_THROW(Dune::InvalidStateException, "Unkown mortar variable type");
    }

    //! apply operator to x, scale and add:  \f$ y = y + \alpha A(x) \f$
    virtual void applyscaleadd(FieldType alpha, const MortarSolutionVector& x, MortarSolutionVector& y) const
    {
        MortarSolutionVector yTmp;

        apply(x, yTmp);
        yTmp *= alpha;

        y += yTmp;
    }

    //! Category of the solver (see SolverCategory::Category)
    virtual Dune::SolverCategory::Category category() const
    { return Dune::SolverCategory::sequential; }

    //! Return the coupled scvf map of sub-domain 1
    const CoupledScvfMap1& coupledScvfMap1() const
    { return coupledScvfMap1_; }

    //! Return the coupled scvf map of sub-domain 2
    const CoupledScvfMap2& coupledScvfMap2() const
    { return coupledScvfMap2_; }

    const Projector& projector1() const
    { return *projector1_; }

    const Projector& projector2() const
    { return *projector2_; }

private:
    std::shared_ptr<Solver1> solver1_;
    std::shared_ptr<Solver2> solver2_;

    std::shared_ptr<Projector> projector1_;
    std::shared_ptr<Projector> projector2_;

    OnePMortarVariableType variableType_;

    CoupledScvfMap1 coupledScvfMap1_;
    CoupledScvfMap2 coupledScvfMap2_;
};

} // end namespace Dumux

#endif
