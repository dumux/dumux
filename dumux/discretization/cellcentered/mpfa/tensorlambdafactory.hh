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
 * \ingroup CCMpfaDiscretization
 * \brief Helper class to be used to obtain lambda functions for the tensors
 *        involved in the laws that describe the different kind of fluxes that
 *        occur in DuMuX models (i.e. advective, diffusive and heat conduction fluxes).
 *        The local systems appearing in Mpfa methods have to be solved subject to
 *        the different tensors. This class returns lambda expressions to be used in the
 *        local systems. The specialization for other discretization methods allows
 *        compatibility with the TPFA scheme, which could be used for one or more of the tensors.
 */
#ifndef DUMUX_DISCRETIZATION_MPFA_TENSOR_LAMBDA_FACTORY_HH
#define DUMUX_DISCRETIZATION_MPFA_TENSOR_LAMBDA_FACTORY_HH

#include <dumux/discretization/methods.hh>
#include <dumux/discretization/cellcentered/mpfa/tensorlambdafactory.hh>

namespace Dumux {

/*!
 * \ingroup CCMpfaDiscretization
 * \brief Helper class to be used to obtain lambda functions for the tensors
 *        involved in the laws that describe the different kind of fluxes that
 *        occur in DuMuX models (i.e. advective, diffusive and heat conduction fluxes).
 *        The local systems appearing in Mpfa methods have to be solved subject to
 *        the different tensors. This class returns lambda expressions to be used in the
 *        local systems. The specialization for discretization methods other than mpfa allows
 *        compatibility with the TPFA scheme, which could be used for one or more of the tensors.
 *        The interfaces of the lambdas are chosen such that all involved tensors can be extracted
 *        with the given arguments.
 */
template<DiscretizationMethod discMethod>
class TensorLambdaFactory;

//! Specialization for mpfa schemes
template<>
class TensorLambdaFactory<DiscretizationMethod::ccmpfa>
{
public:

    //! return void lambda here
    static auto getAdvectionLambda()
    {
        return [] (const auto& problem,
                   const auto& element,
                   const auto& volVars,
                   const auto& fvGeometry,
                   const auto& scv)
               { return volVars.permeability(); };
    }

    //! return void lambda here
    static auto getDiffusionLambda(unsigned int phaseIdx, unsigned int compIdx)
    {
        return [phaseIdx, compIdx] (const auto& problem,
                                    const auto& element,
                                    const auto& volVars,
                                    const auto& fvGeometry,
                                    const auto& scv)
               { return volVars.diffusionCoefficient(phaseIdx, compIdx); };
    }

    //! return void lambda here
    template<class ThermalConductivityModel>
    static auto getHeatConductionLambda()
    {
        return [] (const auto& problem,
                   const auto& element,
                   const auto& volVars,
                   const auto& fvGeometry,
                   const auto& scv)
               { return ThermalConductivityModel::effectiveThermalConductivity(volVars,
                                                                               problem.spatialParams(),
                                                                               element,
                                                                               fvGeometry,
                                                                               scv); };
    }
};

} // end namespace

#endif
