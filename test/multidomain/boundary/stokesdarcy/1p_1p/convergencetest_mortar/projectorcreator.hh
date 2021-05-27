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
#ifndef DUMUX_MORTAR_PROJECTOR_CREATOR_HH
#define DUMUX_MORTAR_PROJECTOR_CREATOR_HH

#include <dune/common/promotiontraits.hh>
#include <dune/functions/functionspacebases/lagrangebasis.hh>

#include <dumux/discretization/functionspacebasis.hh>
#include <dumux/discretization/method.hh>
#include <dumux/multidomain/glue.hh>

#include "projector.hh"
#include "mortarvariabletype.hh"

namespace Dumux {

//! HACK for staggered. Provide FunctionSpaceBasisTraits.
//! We use piecewise constants here.
template< class GridGeometry >
struct FunctionSpaceBasisTraits<GridGeometry, DiscretizationMethod::staggered>
{ using GlobalBasis = Dune::Functions::LagrangeBasis<typename GridGeometry::GridView, /*order*/0>; };

/*!
 * \ingroup TODO doc me.
 * \brief TODO doc me.
 */
class DefaultMortarProjectorCreator
{
    template<class MortarSolution>
    using BaseProjectorPointer = std::shared_ptr<MortarProjectorBase<MortarSolution>>;

    template<class MortarSolution>
    using ProjectorPointer = std::shared_ptr<MortarProjector<MortarSolution>>;

    template<class MortarSolution>
    using TransposedProjectorPointer = std::shared_ptr<TransposedMortarProjector<MortarSolution>>;
public:
    template<class MortarSolution, class Solver1, class Solver2, class MortarGG>
    static std::pair< BaseProjectorPointer<MortarSolution>,
                      BaseProjectorPointer<MortarSolution> >
    makeProjectors(const Solver1& solver1,
                   const Solver2& solver2,
                   const MortarGG& mortarGG,
                   OnePMortarVariableType mv)
    {
        const auto glue1 = makeGlue(*solver1.gridGeometryPointer(), mortarGG);
        const auto glue2 = makeGlue(*solver2.gridGeometryPointer(), mortarGG);

        using namespace Dune::Functions::BasisFactory;
        // if (mv == OnePMortarVariableType::pressure)
        // {
            // Currently we only support schemes with piecewise constant fluxes
            // const auto fluxBasis1 = makeBasis(solver1.gridGeometryPointer()->gridView(), lagrange<0>());
            // const auto fluxBasis2 = makeBasis(solver2.gridGeometryPointer()->gridView(), lagrange<0>());
            // const auto& mortarBasis = getFunctionSpaceBasis(mortarGG);
            //
            // auto baseProjectors1 = makeProjectorPair(fluxBasis1, mortarBasis, glue1);
            // auto baseProjectors2 = makeProjectorPair(fluxBasis2, mortarBasis, glue2);
            //
            // using P = MortarProjector<MortarSolution>;
            // auto p1 = std::make_shared<P>(std::move(baseProjectors1.second), std::move(baseProjectors1.first));
            // auto p2 = std::make_shared<P>(std::move(baseProjectors2.second), std::move(baseProjectors2.first));
            // return std::make_pair(p1, p2);
        // }
        // else
        // {

            const auto sd1Basis = getFunctionSpaceBasis(*solver1.gridGeometryPointer());
            const auto sd2Basis = getFunctionSpaceBasis(*solver2.gridGeometryPointer());
            const auto& mortarBasis = getFunctionSpaceBasis(mortarGG);

            auto baseProjectorMatrices1 = makeProjectionMatricesPair(sd1Basis, mortarBasis, glue1).second;
            auto baseProjectorMatrices2 = makeProjectionMatricesPair(sd2Basis, mortarBasis, glue2).second;

            using P = TransposedMortarProjector<MortarSolution>;
            auto p1 = std::make_shared<P>(std::move(baseProjectorMatrices1.first), std::move(baseProjectorMatrices1.second));
            auto p2 = std::make_shared<P>(std::move(baseProjectorMatrices2.first), std::move(baseProjectorMatrices2.second));
            return std::make_pair(p1, p2);

        // }
    }
};

} // end namespace Dumux

#endif
