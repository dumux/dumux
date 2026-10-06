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

#include <concepts>
#include <type_traits>

#include <dune/common/typetraits.hh>

#include <dumux/common/parameters.hh>
#include <dumux/common/concepts/ipdata_.hh>
#include <dumux/discretization/cvfe/localdof.hh>

namespace Dumux::Experimental {

/*!
 * \ingroup CVFEFlux
 * \brief Darcy's law for control-volume finite element schemes, evaluated at an interpolation point
 *
 * Each term is interpolated using the shape values of the basis functions.
 *
 * \tparam Scalar the scalar type for scalar physical quantities
 * \tparam GridDiscretization the grid discretization
 */
template<class Scalar, class GridDiscretization>
class CVFEDarcysLawAtIp
{
    using GridView = typename GridDiscretization::GridView;
    using Element = typename GridView::template Codim<0>::Entity;
    using GlobalPosition = typename Element::Geometry::GlobalCoordinate;
    using ElementDiscretization = typename GridDiscretization::LocalView;

public:

    /*!
     * \brief Returns the flux term of a fluid phase at an interpolation point
     * \note This term is \f$-\mathbf{T} \left( \nabla p - \rho \mathbf{g} \right)\f$, where the
     *       tensor \f$\mathbf{T}\f$ is given per local dof and interpolated alongside the potential
     *       gradient. A scheme that upwinds the mobility passes the permeability only.
     */
    template<class Problem, class ElementVariables, Concept::IpData IpData, class Tensor>
        requires std::invocable<Tensor, const typename ElementVariables::Variables&>
    static GlobalPosition fluxTerm(const Problem& problem,
                                   const Element& element,
                                   const ElementDiscretization& elemDisc,
                                   const ElementVariables& elemVars,
                                   const int phaseIdx,
                                   const IpData& ipData,
                                   const Tensor& tensor)
    {
        static const bool enableGravity = getParamFromGroup<bool>(problem.paramGroup(), "Problem.EnableGravity");

        const auto& ipCache = cache(elemVars, ipData);
        const auto& shapeValues = ipCache.shapeValues();

        std::decay_t<decltype(tensor(elemVars[0]))> interpolatedTensor(0.0);
        GlobalPosition gradP(0.0);
        Scalar rho(0.0);

        for (const auto& localDof : localDofs(elemDisc))
        {
            const auto& vars = elemVars[localDof];
            const auto shapeValue = shapeValues[localDof.index()][0];

            interpolatedTensor += shapeValue*tensor(vars);

            if (enableGravity)
                rho += vars.density(phaseIdx)*shapeValue;

            // the global shape function gradient
            gradP.axpy(vars.pressure(phaseIdx), ipCache.gradN(localDof.index()));
        }

        if (enableGravity)
            gradP.axpy(-rho, problem.spatialParams().gravity(ipData.global()));

        return fluxTerm_(interpolatedTensor, gradP);
    }

    /*!
     * \brief Returns the flux term of a fluid phase at an interpolation point
     * \note Overload for a tensor that is a property of the interpolation point rather than of
     *       the local dofs, for instance an average over the face being integrated.
     */
    template<class Problem, class ElementVariables, Concept::IpData IpData, class Tensor>
        requires std::invocable<Tensor, const ElementDiscretization&, const ElementVariables&, const IpData&>
    static GlobalPosition fluxTerm(const Problem& problem,
                                   const Element& element,
                                   const ElementDiscretization& elemDisc,
                                   const ElementVariables& elemVars,
                                   const int phaseIdx,
                                   const IpData& ipData,
                                   const Tensor& tensor)
    {
        static const bool enableGravity = getParamFromGroup<bool>(problem.paramGroup(), "Problem.EnableGravity");

        const auto& ipCache = cache(elemVars, ipData);
        const auto& shapeValues = ipCache.shapeValues();

        GlobalPosition gradP(0.0);
        Scalar rho(0.0);

        for (const auto& localDof : localDofs(elemDisc))
        {
            const auto& vars = elemVars[localDof];

            if (enableGravity)
                rho += vars.density(phaseIdx)*shapeValues[localDof.index()][0];

            // the global shape function gradient
            gradP.axpy(vars.pressure(phaseIdx), ipCache.gradN(localDof.index()));
        }

        if (enableGravity)
            gradP.axpy(-rho, problem.spatialParams().gravity(ipData.global()));

        return fluxTerm_(tensor(elemDisc, elemVars, ipData), gradP);
    }

private:
    /*!
     * \brief Returns \f$-\mathbf{T} \nabla \Phi\f$ for the given potential gradient
     * \note An isotropic medium is described by a scalar, which is a valid tensor here but
     *       offers no matrix-vector product, so the two cases are formed differently.
     */
    template<class Tensor>
    static GlobalPosition fluxTerm_(const Tensor& tensor, const GlobalPosition& gradP)
    {
        GlobalPosition fluxTerm(0.0);
        if constexpr (Dune::IsNumber<Tensor>::value)
            fluxTerm.axpy(-tensor, gradP);
        else
        {
            tensor.mv(gradP, fluxTerm);
            fluxTerm *= -1.0;
        }
        return fluxTerm;
    }
};

} // end namespace Dumux::Experimental

#endif
