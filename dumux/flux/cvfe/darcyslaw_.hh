// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup CVFEFlux
 * \brief Darcy's law for control-volume finite element schemes, evaluated at an interpolation point
 */
#ifndef DUMUX_FLUX_CVFE_DARCYS_LAW__HH
#define DUMUX_FLUX_CVFE_DARCYS_LAW__HH

#include <type_traits>

#include <dumux/common/math.hh>
#include <dumux/common/parameters.hh>
#include <dumux/common/concepts/ipdata_.hh>
#include <dumux/discretization/cvfe/localdof.hh>
#include <dumux/flux/gravitytermdensity.hh>
#include <dumux/flux/cvfe/coefficients_.hh>

namespace Dumux::Experimental {

/*!
 * \ingroup CVFEFlux
 * \brief Darcy's law for control-volume finite element schemes, evaluated at an interpolation point
 *
 * Each term is interpolated using the shape values of the basis functions.
 *
 * \tparam Scalar the scalar type for scalar physical quantities
 * \tparam GridDiscretization the grid discretization
 * \tparam GravityTermDensity provides the density \f$\rho\f$ in the gravity term,
 *         `GravityTermDensity::density(problem, vars, phaseIdx)`
 */
template<class Scalar, class GridDiscretization, class GravityTermDensity = PhaseDensityInGravityTerm>
class CVFEDarcysLawAtIp
{
    using GridView = typename GridDiscretization::GridView;
    using Element = typename GridView::template Codim<0>::Entity;
    using GlobalPosition = typename Element::Geometry::GlobalCoordinate;
    using ElementDiscretization = typename GridDiscretization::LocalView;

public:

    /*!
     * \brief Returns the flux term of a fluid phase at an interpolation point
     * \note Looks up the cache of the interpolation point, see the overload below.
     */
    template<class Problem, class ElementVariables, Concept::IpData IpData, class Tensor>
    static GlobalPosition fluxTerm(const Problem& problem,
                                   const Element& element,
                                   const ElementDiscretization& elemDisc,
                                   const ElementVariables& elemVars,
                                   const int phaseIdx,
                                   const IpData& ipData,
                                   const Tensor& tensor)
    { return fluxTerm(problem, element, elemDisc, elemVars, phaseIdx, ipData, cache(elemVars, ipData), tensor); }

    /*!
     * \brief Returns the flux term of a fluid phase at an interpolation point
     * \note This term is \f$-\mathbf{T} \left( \nabla p - \rho \mathbf{g} \right)\f$. The tensor
     *       \f$\mathbf{T}\f$ is given per local dof, `tensor(vars)`, and then interpolated alongside the
     *       potential gradient, or per interpolation point, `tensor(elemDisc, elemVars, ipData)`, for
     *       instance an average over the face being integrated. A scheme that upwinds the mobility
     *       passes the permeability only.
     */
    template<class Problem, class ElementVariables, Concept::IpData IpData, class IpCache, class Tensor>
    static GlobalPosition fluxTerm(const Problem& problem,
                                   const Element& element,
                                   const ElementDiscretization& elemDisc,
                                   const ElementVariables& elemVars,
                                   const int phaseIdx,
                                   const IpData& ipData,
                                   const IpCache& ipCache,
                                   const Tensor& tensor)
    {
        static const bool enableGravity = getParamFromGroup<bool>(problem.paramGroup(), "Problem.EnableGravity");

        constexpr bool tensorPerLocalDof = Dumux::Detail::CVFE::LocalDofCoefficient<Tensor, ElementVariables>;
        using InterpolatedTensor = typename std::conditional_t<
            tensorPerLocalDof,
            Dumux::Detail::CVFE::LocalDofValue<Tensor, ElementVariables>,
            std::type_identity<Scalar>
        >::type;

        const auto& shapeValues = ipCache.shapeValues();
        [[maybe_unused]] InterpolatedTensor interpolatedTensor(0.0);
        GlobalPosition gradPotential(0.0);
        Scalar rho(0.0);
        for (const auto& localDof : localDofs(elemDisc))
        {
            const auto& vars = elemVars[localDof];
            const auto shapeValue = shapeValues[localDof.index()][0];

            if constexpr (tensorPerLocalDof)
                interpolatedTensor += shapeValue*tensor(vars);

            if (enableGravity)
                rho += GravityTermDensity::density(problem, vars, phaseIdx)*shapeValue;

            gradPotential.axpy(vars.pressure(phaseIdx), ipCache.gradN(localDof.index()));
        }

        if (enableGravity)
            gradPotential.axpy(-rho, problem.spatialParams().gravity(ipData.global()));

        auto fluxTerm = [&]
        {
            if constexpr (tensorPerLocalDof)
                return mv(interpolatedTensor, gradPotential);
            else
                return mv(tensor(elemDisc, elemVars, ipData), gradPotential);
        }();
        fluxTerm *= -1.0;
        return fluxTerm;
    }
};

} // end namespace Dumux::Experimental

#endif
