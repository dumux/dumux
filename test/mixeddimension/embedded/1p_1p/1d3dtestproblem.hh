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
 * \brief A test problem for the 1d3d coupled problem:
 * a one-dimensional network is embedded in a three-dimensional matrix
 * Comparison with anaytical solutiom (see D'Angelo 2007 PhD thesis)
 */
#ifndef DUMUX_1D3D_TEST_PROBLEM_HH
#define DUMUX_1D3D_TEST_PROBLEM_HH

#include "bloodflowproblem.hh"
#include "tissueproblem.hh"

#include <dumux/mixeddimension/problem.hh>
#include <dumux/mixeddimension/embedded/cellcentered/bboxtreecouplingmanager.hh>
//#include <dumux/mixeddimension/embedded/cellcentered/bboxtreecouplingmanagersimple.hh>
#include <dumux/mixeddimension/integrationpointsource.hh>

namespace Dumux
{
template <class TypeTag>
class TestOneDThreeDProblem;

namespace Properties
{
NEW_TYPE_TAG(TestOneDThreeDProblem, INHERITS_FROM(MixedDimension));
NEW_TYPE_TAG(TestOneDThreeDCCProblem, INHERITS_FROM(TestOneDThreeDProblem));

// Set the problem property
SET_TYPE_PROP(TestOneDThreeDProblem, Problem, Dumux::TestOneDThreeDProblem<TypeTag>);

// Set the coupling manager
//SET_TYPE_PROP(TestOneDThreeDCCProblem, CouplingManager, Dumux::CCBBoxTreeEmbeddedCouplingManagerSimple<TypeTag>);
SET_TYPE_PROP(TestOneDThreeDCCProblem, CouplingManager, Dumux::CCBBoxTreeEmbeddedCouplingManager<TypeTag>);
////////////////////////////////////////////////////////////////////////////
// Set the two sub-problems of the global problem
SET_TYPE_PROP(TestOneDThreeDCCProblem, LowDimProblemTypeTag, TTAG(BloodFlowCCProblem));
SET_TYPE_PROP(TestOneDThreeDCCProblem, BulkProblemTypeTag, TTAG(TissueCCProblem));
////////////////////////////////////////////////////////////////////////////

// publish this problem in the sub problems
SET_TYPE_PROP(BloodFlowCCProblem, GlobalProblemTypeTag, TTAG(TestOneDThreeDCCProblem));
SET_TYPE_PROP(TissueCCProblem, GlobalProblemTypeTag, TTAG(TestOneDThreeDCCProblem));

// The subproblems inherit the parameter tree from this problem
SET_PROP(BloodFlowCCProblem, ParameterTree) : GET_PROP(TTAG(TestOneDThreeDCCProblem), ParameterTree) {};
SET_PROP(TissueCCProblem, ParameterTree) : GET_PROP(TTAG(TestOneDThreeDCCProblem), ParameterTree) {};

// Set the point source type of the subproblems to an id'ed point source
SET_TYPE_PROP(BloodFlowCCProblem, PointSource, Dumux::IntegrationPointSource<TTAG(BloodFlowCCProblem), unsigned int>);
SET_TYPE_PROP(BloodFlowCCProblem, PointSourceHelper, Dumux::IntegrationPointSourceHelper<TTAG(BloodFlowCCProblem)>);
SET_TYPE_PROP(TissueCCProblem, PointSource, Dumux::IntegrationPointSource<TTAG(TissueCCProblem), unsigned int>);
SET_TYPE_PROP(TissueCCProblem, PointSourceHelper, Dumux::IntegrationPointSourceHelper<TTAG(TissueCCProblem)>);

SET_TYPE_PROP(TestOneDThreeDProblem, LinearSolver, ILU0BiCGSTABBackend);

}//end namespace properties

template <class TypeTag>
class TestOneDThreeDProblem : public MixedDimensionProblem<TypeTag>
{
    using ParentType = MixedDimensionProblem<TypeTag>;
    using TimeManager = typename GET_PROP_TYPE(TypeTag, TimeManager);
    using Scalar = typename GET_PROP_TYPE(TypeTag, Scalar);

    // obtain the type tags of the sub problems
    using BulkProblemTypeTag = typename GET_PROP_TYPE(TypeTag, BulkProblemTypeTag);
    using LowDimProblemTypeTag = typename GET_PROP_TYPE(TypeTag, LowDimProblemTypeTag);

    // obtain types from the sub problem type tags
    using BulkProblem = typename GET_PROP_TYPE(BulkProblemTypeTag, Problem);
    using LowDimProblem = typename GET_PROP_TYPE(LowDimProblemTypeTag, Problem);

    using BulkGridView = typename GET_PROP_TYPE(BulkProblemTypeTag, GridView);
    using LowDimGridView = typename GET_PROP_TYPE(LowDimProblemTypeTag, GridView);

    enum { bulkIsBox = GET_PROP_VALUE(BulkProblemTypeTag, ImplicitIsBox) };
    enum { lowDimIsBox = GET_PROP_VALUE(LowDimProblemTypeTag, ImplicitIsBox) };

public:
    TestOneDThreeDProblem(TimeManager &timeManager, const BulkGridView &bulkGridView, const LowDimGridView &lowDimgridView)
    : ParentType(timeManager, bulkGridView, lowDimgridView)
    {
        order_ = GET_RUNTIME_PARAM_FROM_GROUP(TypeTag, int, Problem, NormIntegrationOrder);
        excludeInnerBulk_ = GET_RUNTIME_PARAM_FROM_GROUP(TypeTag, bool, Problem, NormExcludeInnerBulk);
    }

    /*!
     * \brief Called by Dumux::TimeManager in order to do a time
     *        integration on the model. Iterative algorithm for one time step.
     */
    void timeIntegration()
    {
        ParentType::timeIntegration();

        // compute the L2-norm with respect to the analytical solution
        // discrete l2-norm for finite volumes
        auto norm1D = this->calculateNorm(this->lowDimProblem(), this->lowDimGridView(), order_, lowDimIsBox);
        auto norm3D = this->calculateNorm(this->bulkProblem(), this->bulkGridView(), order_, bulkIsBox);
        // normalize the norm
        norm1D.first /= this->calculateNormExactSolution(this->lowDimProblem(), this->lowDimGridView(), order_, lowDimIsBox);
        norm3D.first /= this->calculateNormExactSolution(this->bulkProblem(), this->bulkGridView(), order_, bulkIsBox);

        // ouput result to terminal
        std::cout << "hmax_3d: " << norm3D.second << " "
                  << "hmax_1d: " << norm1D.second << " "
                  << "L2_3d: " << norm3D.first << " "
                  << "L2_1d: " << norm1D.first << std::endl;

        // ... and file.
        file_.open(this->name() + ".log", std::ios::app);
        file_ << norm3D.second << " " << norm1D.second << " " << norm3D.first << " " << norm1D.first << "\n";
        file_.close();
    }

    //! Calculate the L2-Norm of the solution and hmax of the grid
    template<class SubProblem, class GridView>
    std::pair<Scalar, Scalar> calculateNorm(const SubProblem &sp, const GridView &gv, int order, bool isBox)
    {
        typename Dune::PQkLocalFiniteElementCache<Scalar, Scalar, GridView::dimension, 1> feCache;

        // calculate the L2-Norm and hmax
        const auto& solution = sp.model().curSol();
        Scalar LTwoNorm = 0.0;
        Scalar hMax = 0.0;

        // iterate over all elements
        for (const auto& element : elements(gv))
        {
            const unsigned int eIdx = gv.indexSet().index(element);
            const auto geometry = element.geometry();
            const auto center = geometry.center();
            if (int(GridView::dimension) == 3 && excludeInnerBulk_ &&
                std::sqrt(center[0]*center[0] + center[1]*center[1]) < GET_RUNTIME_PARAM_FROM_GROUP(TypeTag, Scalar, SpatialParams, Radius) + 0.01)
                continue;
            hMax = std::max(geometryDiameter(geometry), hMax);

            using Quad = Dune::QuadratureRule<Scalar, GridView::dimension>;
            const Quad &quad = Dune::QuadratureRules<Scalar, GridView::dimension>::rule(geometry.type(), order);
            for(auto&& qp : quad)
            {
                const auto globalPos = geometry.global(qp.position());
                const Scalar pe = sp.exactSolution(globalPos);
                Scalar p(0.0);
                if(!isBox)
                    p = solution[eIdx][0];
                else
                {
                    // do interpolation with ansatz functions
                    std::vector<Dune::FieldVector<Scalar,1> > shapeValues;
                    const auto& localFiniteElement = feCache.get(geometry.type());
                    localFiniteElement.localBasis().evaluateFunction(qp.position(), shapeValues);
                    for (int vIdx = 0; vIdx < shapeValues.size(); ++vIdx)
                        p += shapeValues[vIdx]*solution[sp.model().dofMapper().subIndex(element, vIdx, GridView::dimension)][0];
                }

                LTwoNorm += (p - pe)*(p - pe)*qp.weight()*geometry.integrationElement(qp.position());
            }
        }
        return std::make_pair(std::sqrt(LTwoNorm), hMax);
    }

    template<class SubProblem, class GridView>
    Scalar calculateNormExactSolution(const SubProblem &sp, const GridView &gv, int order, bool isBox)
    {
        // calculate the L2-Norm
        Scalar LTwoNorm = 0.0;

        // iterate over all elements
        for (const auto& element : elements(gv))
        {
            const auto geometry = element.geometry();
            const auto center = geometry.center();
            if (int(GridView::dimension) == 3 && excludeInnerBulk_ &&
                std::sqrt(center[0]*center[0] + center[1]*center[1]) < GET_RUNTIME_PARAM_FROM_GROUP(TypeTag, Scalar, SpatialParams, Radius))
                continue;
            using Quad = Dune::QuadratureRule<Scalar, GridView::dimension>;
            const Quad &quad = Dune::QuadratureRules<Scalar, GridView::dimension>::rule(geometry.type(), order);
            for(auto&& qp : quad)
            {
                const auto globalPos = geometry.global(qp.position());
                const Scalar pe = sp.exactSolution(globalPos);
                LTwoNorm += pe*pe*qp.weight()*geometry.integrationElement(qp.position());
            }
        }
        return std::sqrt(LTwoNorm);
    }

    template<class Geometry>
    Scalar geometryDiameter(const Geometry& geometry)
    {
        auto type = geometry.type();
        if(type.isLine())
        {
            return (geometry.corner(0)-geometry.corner(1)).two_norm();
        }
        else if(type.isHexahedron())
        {
           return (geometry.corner(0)-geometry.corner(1)).two_norm();
        }
        else if(type.isTetrahedron())
        {
            const auto p0 = geometry.corner(0);
            const auto p1 = geometry.corner(1);
            const auto p2 = geometry.corner(2);
            const auto p3 = geometry.corner(3);

            // Compute side length
            const auto a = (p1-p2).two_norm();
            const auto b = (p0-p2).two_norm();
            const auto c = (p0-p1).two_norm();
            const auto aa = (p0-p3).two_norm();
            const auto bb = (p1-p3).two_norm();
            const auto cc = (p2-p3).two_norm();

            // compute some helper variables
            const auto la = a*aa;
            const auto lb = b*bb;
            const auto lc = c*cc;
            const auto s = 0.5*(la+lb+lc);
            const auto area = std::sqrt(s*(s-la)*(s-lb)*(s-lc));
            return area/(3.0*geometry.volume());
        }
        else
        {
            DUNE_THROW(Dune::NotImplemented, "Diameter for geometry " << type << ".");
        }
    }

private:
    std::ofstream file_;
    int order_;
    bool excludeInnerBulk_;
};

} //end namespace Dumux

#endif
