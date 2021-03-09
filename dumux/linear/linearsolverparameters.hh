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
 * \ingroup Linear
 * \brief Generates a parameter tree required for the linear solvers and precondioners of the Dune ISTL
 */

#ifndef DUMUX_LINEAR_SOLVER_PARAMETERS_HH
#define DUMUX_LINEAR_SOLVER_PARAMETERS_HH

#include <string>
#include <array>
#include <vector>

#include <dune/common/parametertree.hh>
#include <dumux/common/parameters.hh>

namespace Dumux {

template<class LinearSolverTraits>
class LinearSolverParameters
{
public:

    //! Translation table for solver parameters
    static const std::vector<std::array<std::string, 2>> dumuxToIstlSolverParams;

    //! Create a tree containing parameters required for the linear solvers and precondioners of the Dune ISTL
    static Dune::ParameterTree createParameterTree(const std::string& paramGroup = "")
    {
        Dune::ParameterTree params;
        setDefaultParameters(params, paramGroup);
        fillValuesForIstlKeys(params, paramGroup);
        return params;
    }

    //! Set some defaults for the solver parameters
    static void setDefaultParameters(Dune::ParameterTree& params, const std::string& paramGroup = "")
    {
        params["restart"] = "10";
        params["maxit"] = "250";
        params["reduction"] = "1e-13";
        params["verbose"] = "0";
        params["preconditioner.iterations"] = "1";
        params["preconditioner.relaxation"] = "1.0";
        params["preconditioner.verbosity"] = "0";
        params["preconditioner.defaultAggregationSizeMode"] = "isotropic";
        params["preconditioner.defaultAggregationDimension"] = std::to_string(LinearSolverTraits::GridView::dimension);
        params["preconditioner.maxLevel"] = "100";
        params["ParameterGroup"] = paramGroup;
        params["preconditioner.ParameterGroup"] = paramGroup;
    }

    //! Iterate over all keys required by the ISTL, translate them to Dumux syntax and add values to tree
    static void fillValuesForIstlKeys(Dune::ParameterTree& params, const std::string& paramGroup = "")
    {
        const auto linearSolverGroups = getParamSubGroups("LinearSolver", paramGroup);
        if (linearSolverGroups.empty()) // no linear solver parameters were specified
            return;

        for (const auto& [dumuxKey, istlKey] : dumuxToIstlSolverParams)
        {
            for (const auto& group : linearSolverGroups)
            {
                const auto fullDumuxKey = group + "." + dumuxKey;
                const auto value = getParam<std::string>(fullDumuxKey, "");
                if (!value.empty())
                {
                    params[istlKey] = value;
                    break; // skip groups with smaller depth in the tree
                }
            }
        }
    }
};

//! Translation table for solver parameters
//! TODO change to constexpr array of std::string_view once we require g++ >= 7.3 (bug in older versions)
template<class LinearSolverTraits>
const std::vector<std::array<std::string, 2>>
LinearSolverParameters<LinearSolverTraits>::dumuxToIstlSolverParams =
{
    // solver params
    {"Verbosity", "verbose"},
    {"MaxIterations", "maxit"},
    {"ResidualReduction", "reduction"},
    {"Type", "type"},
    {"GMResRestart", "restart"}, // cycles before restarting
    {"Restart", "restart"}, // cycles before restarting
    {"MaxOrthogonalizationVectors", "mmax"},

    // preconditioner params
    {"Preconditioner.Verbosity", "preconditioner.verbosity"},
    {"Preconditioner.Type", "preconditioner.type"},
    {"Preconditioner.Iterations", "preconditioner.iterations"},
    {"Preconditioner.Relaxation", "preconditioner.relaxation"},
    {"Preconditioner.ILUOrder", "preconditioner.n"},
    {"Preconditioner.ILUResort", "preconditioner.resort"},
    {"Preconditioner.AmgSmootherRelaxation", "preconditioner.smootherRelaxation"},
    {"Preconditioner.AmgSmootherIterations", "preconditioner.smootherIterations"},
    {"Preconditioner.AmgMaxLevel", "preconditioner.maxLevel"},
    {"Preconditioner.AmgCoarsenTarget", "preconditioner.coarsenTarget"},
    {"Preconditioner.AmgMinCoarseningRate", "preconditioner.minCoarseningRate"},
    {"Preconditioner.AmgAccumulationMode", "preconditioner.accumulationMode"},
    {"Preconditioner.AmgProlongationDampingFactor", "preconditioner.prolongationDampingFactor"},
    {"Preconditioner.AmgAlpha", "preconditioner.alpha"},
    {"Preconditioner.AmgBeta", "preconditioner.beta"},
    {"Preconditioner.AmgAdditive", "preconditioner.additive"},
    {"Preconditioner.AmgGamma", "preconditioner.gamma"},
    {"Preconditioner.AmgPreSmoothingSteps", "preconditioner.preSteps"},
    {"Preconditioner.AmgPostSmoothingSteps", "preconditioner.postSteps"},
    {"Preconditioner.AmgCriterionSymmetric", "preconditioner.criterionSymmetric"},
    {"Preconditioner.AmgStrengthMeasure", "preconditioner.strengthMeasure"},
    {"Preconditioner.AmgDiagonalRowIndex", "preconditioner.diagonalRowIndex"},
    {"Preconditioner.AmgDefaultAggregationSizeMode", "preconditioner.defaultAggregationSizeMode"},
    {"Preconditioner.AmgDefaultAggregationDimension", "preconditioner.defaultAggregationDimension"},
    {"Preconditioner.AmgMaxAggregateDistance", "preconditioner.maxAggregateDistance"},
    {"Preconditioner.AmgMinAggregateSize", "preconditioner.minAggregateSize"},
    {"Preconditioner.AmgMaxAggregateSize", "preconditioner.maxAggregateSize"}
};

} // end namespace Dumux

#endif
