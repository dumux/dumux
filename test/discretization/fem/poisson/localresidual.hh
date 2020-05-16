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
 * \todo TODO: WHICH GROUP?
 * \ingroup TODO: WHICH GROUP?
 * \copydoc FEPoissonLocalResidual
 */
#ifndef DUMUX_FEM_POISSON_LOCALRESIDUAL_HH
#define DUMUX_FEM_POISSON_LOCALRESIDUAL_HH

#include <dumux/common/math.hh>
#include <dumux/assembly/felocalresidual.hh>

namespace Dumux {

/*!
 * \file
 * \todo TODO: WHICH GROUP?
 * \ingroup TODO: WHICH GROUP?
 * \brief The local residual class for a poisson problem
 *        using the finite element method
 * \tparam GV The grid variables type
 */
template<class GV>
class FEPoissonLocalResidual : public FELocalResidual< GV, FEPoissonLocalResidual<GV> >
{
    using ThisType = FEPoissonLocalResidual<GV>;
    using ParentType = FELocalResidual<GV, ThisType>;
    using IpVariables = typename GV::IntegrationPointVariables;

public:
    //! export flux term type
    using typename ParentType::FluxTerm;

    //! pull up base class constructors
    using ParentType::ParentType;

    /*!
     * \brief Calculate the flux term of the equation
     * \param elemSol The element solution vector
     * \param ipData The shape function values/gradients evaluated at the integration point
     * \param ipVars The primary/secondary variables evaluated at the integration point
     */
    template<class ElementSolution, class IpData>
    FluxTerm computeFlux(const ElementSolution& elemSol,
                         const IpData& ipData,
                         const IpVariables& ipVars) const
    {
        if (elemSol.size() != ipData.size())
            DUNE_THROW(Dune::InvalidStateException, "Element solution size mismatch");

        // evaluate gradient in solution
        typename IpData::GlobalPosition gradX(0.0);
        for (unsigned int i = 0; i < ipData.size(); ++i)
        {
            auto tmp = ipData.gradN(i);
            tmp *= elemSol[i];
            gradX += tmp;
        }

        // The flux is tensor*gradX
        FluxTerm result(0.0);
        if (result.size() != 1)
            DUNE_THROW(Dune::InvalidStateException, "Flux term size mismatch");

        result[0] = mv(this->problem_().poissonTensor(), gradX);
        return result;
    }
};

} // end namespace Dumux

#endif
