// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*****************************************************************************
 *   See the file COPYING for full copying permissions.                      *
 *                                                                           *
 *   This program is free software: you can redistribute it and/or modify    *
 *   it under the terms of the GNU General Public License as published by    *
 *   the Free Software Foundation, either version 2 of the License, or       *
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
 *
 * \brief Problem file for the single fracture case of Sandve et al. (2012)
 */
#ifndef DUMUX_1D2D_FACET_ANALYTICAL_PROBLEM_HH
#define DUMUX_1D2D_FACET_ANALYTICAL_PROBLEM_HH

#include <dune/geometry/quadraturerules.hh>

#include "analyticfractureproblem.hh"
#include "analyticmatrixproblem.hh"

#include <dumux/mixeddimension/problem.hh>
#include <dumux/mixeddimension/facet/properties.hh>
#include <dumux/mixeddimension/facet/gmshdualfacetgridcreator.hh>
#include <dumux/mixeddimension/facet/mpfa/couplingmanager.hh>

namespace Dumux
{
template <class TypeTag>
class OnePFacetCouplingProblem;

namespace Properties
{
// Type tag of the global Problem and the two sub problems
NEW_TYPE_TAG(OnePFacetCoupling, INHERITS_FROM(MixedDimensionFacetCoupling));

// Set the problem property
SET_TYPE_PROP(OnePFacetCoupling, Problem, Dumux::OnePFacetCouplingProblem<TypeTag>);

// Set the grid creator
SET_TYPE_PROP(OnePFacetCoupling, GridCreator, Dumux::GmshDualFacetGridCreator<TypeTag>);

// Set the two sub-problems of the global problem
SET_TYPE_PROP(OnePFacetCoupling, BulkProblemTypeTag, TTAG(OnePCCMpfaMatrixProblem));
SET_TYPE_PROP(OnePFacetCoupling, LowDimProblemTypeTag, TTAG(OnePCCFractureProblem));

// The coupling manager
SET_TYPE_PROP(OnePFacetCoupling, CouplingManager, CCMpfaFacetCouplingManager<TypeTag>);

// The linear solver to be used
SET_TYPE_PROP(OnePFacetCoupling, LinearSolver, UMFPackBackend<TypeTag>);

// The sub-problems need to know the global problem's type tag
SET_TYPE_PROP(OnePCCMpfaMatrixProblem, GlobalProblemTypeTag, TTAG(OnePFacetCoupling));
SET_TYPE_PROP(OnePCCFractureProblem, GlobalProblemTypeTag, TTAG(OnePFacetCoupling));

// The subproblems inherit the parameter tree from this problem
SET_PROP(OnePCCMpfaMatrixProblem, ParameterTree) : GET_PROP(TTAG(OnePFacetCoupling), ParameterTree) {};
SET_PROP(OnePCCFractureProblem, ParameterTree) : GET_PROP(TTAG(OnePFacetCoupling), ParameterTree) {};

// Set the grids for the two sub-problems
SET_TYPE_PROP(OnePCCMpfaMatrixProblem, Grid, Dune::ALUGrid<2, 2, Dune::cube, Dune::nonconforming>);
SET_TYPE_PROP(OnePCCFractureProblem, Grid, Dune::FoamGrid<1,2>);

}

template <class TypeTag>
class OnePFacetCouplingProblem : public MixedDimensionProblem<TypeTag>
{
    using ParentType = MixedDimensionProblem<TypeTag>;

    using TimeManager = typename GET_PROP_TYPE(TypeTag, TimeManager);

    using BulkProblemTypeTag = typename GET_PROP_TYPE(TypeTag, BulkProblemTypeTag);
    using LowDimProblemTypeTag = typename GET_PROP_TYPE(TypeTag, LowDimProblemTypeTag);

    using BulkGridView = typename GET_PROP_TYPE(BulkProblemTypeTag, GridView);
    using LowDimGridView = typename GET_PROP_TYPE(LowDimProblemTypeTag, GridView);

    using Scalar = typename GET_PROP_TYPE(TypeTag, Scalar);

public:

    OnePFacetCouplingProblem(TimeManager &timeManager, const BulkGridView &bulkGridView, const LowDimGridView &lowDimGridView)
    : ParentType(timeManager, bulkGridView, lowDimGridView)
    {}

    //! use this method to calculate the L2-norms of the errors
    void postTimeStep()
    {
        if (!this->timeManager().willBeFinished())
            return;

        // norm of the errors in the two domains
        const auto m = calculateNorm<BulkProblemTypeTag>(this->bulkProblem(), this->bulkProblem().gridView());
        const auto f = calculateNorm<LowDimProblemTypeTag>(this->lowDimProblem(), this->lowDimProblem().gridView());

        // output to terminal and file
        std::cout << "Matrix h, L2_p, L2_q: " << m[0] << ", " << m[1] << ", " << m[2] << std::endl;
        std::cout << "Fracture h, L2_p, L2_q: " << f[0] << ", " << f[1] << ", " << f[2] << std::endl;

        std::ofstream file;
        file.open(this->name() + ".log", std::ios::app);
        file << m[0] << " " << m[1] << " " << m[2] << " " << f[0] << " " << f[1] << " " << f[2] << "\n";
        file.close();
    }

private:

    template<class T, class SP, class GV>
    std::array<Scalar, 3> calculateNorm(SP& subProblem, const GV& gridView)
    {
        using std::min;
        using std::max;
        using std::sqrt;

        // extract the flux variables type
        using FluxVariables = typename GET_PROP_TYPE(T, FluxVariables);

        // container to store the results
        const auto Li = Scalar(1.0/GET_RUNTIME_PARAM_FROM_GROUP(TypeTag, int, Grid, NumCells));
        std::array<Scalar, 3> result = {Li, 0.0, 0.0};

        // minimum and maximum pressure/flux values (used below to scale the error)
        Scalar uMin = std::numeric_limits<Scalar>::max();
        Scalar uMax = std::numeric_limits<Scalar>::min();
        Scalar qMin = std::numeric_limits<Scalar>::max();
        Scalar qMax = std::numeric_limits<Scalar>::min();
        Scalar LMax = std::numeric_limits<Scalar>::min();
        Scalar sumL = 0.0;
        Scalar sumA = 0.0;

        for (const auto& element : elements(gridView))
        {
            const auto c = element.geometry().center();
            const auto uh = subProblem.model().curSol()[subProblem.elementMapper().index(element)];

            // integrate the pressure error over the element
            const auto eg = element.geometry();
            sumA += eg.volume()*subProblem.extrusionFactorAtPos(c);

            static const bool doIntegration = GET_RUNTIME_PARAM_FROM_GROUP(TypeTag, bool, L2Error, IntegrateElements);
            if (doIntegration)
            {
                static const int order = GET_RUNTIME_PARAM_FROM_GROUP(TypeTag, int, L2Error, QuadratureOrder);
                const auto& rule = Dune::QuadratureRules<Scalar, GV::dimension>::rule(eg.type(), order);
                for (auto&& qp : rule)
                {
                    auto globalPos = eg.global(qp.position());
                    const auto u = subProblem.exact(element, globalPos);

                    uMin = min(u, uMin);
                    uMax = max(u, uMax);

                    result[1] += (uh - u)*(uh - u)*qp.weight()*eg.integrationElement(qp.position())*subProblem.extrusionFactorAtPos(c);
                }
            }
            else
            {
                auto globalPos = c;
                const auto u = subProblem.exact(element, globalPos);
                result[1] += (uh - u)*(uh - u)*eg.volume()*subProblem.extrusionFactorAtPos(c);

                uMin = min(u, uMin);
                uMax = max(u, uMax);
            }

            // evaluate the error in fluxes for this element
            subProblem.couplingManager().setCouplingContext(element);
            auto fvGeometry = localView(subProblem.model().globalFvGeometry());
            auto elemVolVars = localView(subProblem.model().curGlobalVolVars());
            auto elemFluxVarCache = localView(subProblem.model().globalFluxVarsCache());
            fvGeometry.bind(element);
            elemVolVars.bind(element, fvGeometry, subProblem.model().curSol());
            elemFluxVarCache.bind(element, fvGeometry, elemVolVars);

            // sum up the exact and discrete fluxes
            for (const auto& scvf : scvfs(fvGeometry))
            {
                FluxVariables fluxVars;
                fluxVars.init(subProblem, element, fvGeometry, elemVolVars, scvf, elemFluxVarCache);

                auto upwindRule = [] (const auto& volVars) { return 1.0; };
                Scalar qh = fluxVars.advectiveFlux(/*phaseIdx*/0, upwindRule);
                Scalar q = subProblem.exactFlux(element, scvf);

                const Scalar scvfArea = scvf.area()*subProblem.extrusionFactorAtPos(scvf.ipGlobal());
                qh /= scvfArea;
                q /= scvfArea;

                sumL += scvfArea;
                result[2] += scvfArea*(q-qh)*(q-qh);

                qMax = max(qMax, q);
                qMin = min(qMin, q);
            }
        }

        // final error calculation
        result[0] = GET_RUNTIME_PARAM_FROM_GROUP(TypeTag, Scalar, SpatialParams, FractureAperture)/Li;
        result[1] = sqrt(result[1])/(uMax-uMin)/sqrt(sumA);
        result[2] = sqrt(result[2])/(qMax-qMin)/sqrt(sumL);
        return result;
    }
};

} //end namespace

#endif
