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
 * \ingroup FEMFlux
 * \brief Specialization of the effective stress law for finite element schemes. This computes
 *        the stress tensor resulting from mechanical deformation and the pore pressure.
 */
#ifndef DUMUX_DISCRETIZATION_FEM_EFFECTIVE_STRESS_LAW_HH
#define DUMUX_DISCRETIZATION_FEM_EFFECTIVE_STRESS_LAW_HH

#include <dumux/flux/effectivestresslaw.hh>
#include <dumux/discretization/method.hh>

namespace Dumux {

/*!
 * \ingroup FEMFlux
 * \brief Effective stress law for finite element schemes
 * \tparam StressType type used for the computation of
 *         purely mechanical stresses (i.e. material law)
 * \tparam GridGeometry the finite element grid geometry
 */
template<class StressType, class GridGeometry>
class EffectiveStressLaw<StressType, GridGeometry, DiscretizationMethod::fem>
{
    using FEElementGeometry = typename GridGeometry::LocalView;
    using GridView = typename GridGeometry::GridView;
    using Element = typename GridView::template Codim<0>::Entity;

    static constexpr int dim = GridView::dimension;
    static constexpr int dimWorld = GridView::dimensionworld;
    static_assert(dim == dimWorld, "EffectiveStressLaw not implemented for network/surface grids");
    static_assert(StressType::discMethod == DiscretizationMethod::fem, "The provided stress type must be specialized for fe schemes");

public:
    //! export the type used for scalar values
    using Scalar = typename StressType::Scalar;
    //! export the type used for the stress tensor
    using StressTensor = typename StressType::StressTensor;
    //! export the type used for force vectors
    using ForceVector = typename StressType::ForceVector;
    //! state the discretization method this implementation belongs to
    static constexpr DiscretizationMethod discMethod = DiscretizationMethod::fem;

    //! assembles the stress tensor at a given integration point
    template<class Problem, class ElementSolution, class IpData, class SecondaryVariables>
    static StressTensor stressTensor(const Problem& problem,
                                     const Element& element,
                                     const FEElementGeometry& feGeometry,
                                     const ElementSolution& elemSol,
                                     const IpData& ipData,
                                     const SecondaryVariables& secVars)
    {
        // compute the purely mechanical stress
        auto sigma = StressType::stressTensor(problem, element, feGeometry, elemSol, ipData, secVars);

        // obtain biot coefficient and effective pore pressure
        const auto biotCoeff = problem.spatialParams().biotCoefficient(element, feGeometry, elemSol, ipData, secVars);
        const auto effPress = problem.effectivePorePressure(element, feGeometry, elemSol, ipData, secVars);

        // subtract pore pressure from the diagonal entries
        const auto bcp = biotCoeff*effPress;
        for (int i = 0; i < dim; ++i)
            sigma[i][i] -= bcp;

        return sigma;
    }

    //! assembles the (effective) stress tensor of the solid skeleton at a given integration point
    template<class Problem, class ElementSolution, class IpData, class SecondaryVariables>
    static StressTensor effectiveStressTensor(const Problem& problem,
                                              const Element& element,
                                              const FEElementGeometry& feGeometry,
                                              const ElementSolution& elemSol,
                                              const IpData& ipData,
                                              const SecondaryVariables& secVars)
    { return StressType::stressTensor(problem, element, feGeometry, elemSol, ipData, secVars); }
};

} // end namespace Dumux

#endif
