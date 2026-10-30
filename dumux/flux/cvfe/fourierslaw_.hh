// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup CVFEFlux
 * \brief Fourier's law for control-volume finite element schemes, evaluated at an interpolation point
 */
#ifndef DUMUX_FLUX_CVFE_FOURIERS_LAW__HH
#define DUMUX_FLUX_CVFE_FOURIERS_LAW__HH

#include <dumux/common/math.hh>
#include <dumux/common/concepts/ipdata_.hh>
#include <dumux/flux/cvfe/coefficients_.hh>

namespace Dumux::Experimental {

/*!
 * \ingroup CVFEFlux
 * \brief Fourier's law for control-volume finite element schemes, evaluated at an interpolation point
 *
 * The conductive heat flux is \f$ -\boldsymbol{\lambda} \nabla T \f$ with the temperature
 * gradient of the interpolated temperature.
 *
 * \tparam Scalar the scalar type for scalar physical quantities
 * \tparam GridDiscretization the grid discretization
 */
template<class Scalar, class GridDiscretization>
class CVFEFouriersLawAtIp
{
    using GridView = typename GridDiscretization::GridView;
    using Element = typename GridView::template Codim<0>::Entity;
    using GlobalPosition = typename Element::Geometry::GlobalCoordinate;
    using ElementDiscretization = typename GridDiscretization::LocalView;

public:
    /*!
     * \brief Returns the conductive heat flux term at an interpolation point
     * \note Looks up the cache of the interpolation point, see the overload below.
     */
    template<class Problem, class ElementVariables, Concept::IpData IpData, class Conductivity>
    static GlobalPosition fluxTerm(const Problem& problem,
                                   const Element& element,
                                   const ElementDiscretization& elemDisc,
                                   const ElementVariables& elemVars,
                                   const IpData& ipData,
                                   const Conductivity& conductivity)
    { return fluxTerm(problem, element, elemDisc, elemVars, ipData, cache(elemVars, ipData), conductivity); }

    /*!
     * \brief Returns the conductive heat flux term \f$ -\boldsymbol{\lambda} \nabla T \f$ at an interpolation point
     * \note The thermal conductivity is given per local dof, `conductivity(vars)`, and then interpolated,
     *       or per interpolation point, `conductivity(elemDisc, elemVars, ipData)`, for instance
     *       Experimental::CVFE::faceTensorAverageAtScvf.
     */
    template<class Problem, class ElementVariables, Concept::IpData IpData, class IpCache, class Conductivity>
    static GlobalPosition fluxTerm(const Problem& problem,
                                   const Element& element,
                                   const ElementDiscretization& elemDisc,
                                   const ElementVariables& elemVars,
                                   const IpData& ipData,
                                   const IpCache& ipCache,
                                   const Conductivity& conductivity)
    {
        const auto gradT = Dumux::Detail::CVFE::gradientAtIp(elemDisc, elemVars, ipCache,
                                                             [] (const auto& vars) { return vars.temperature(); });
        auto fluxTerm = mv(Dumux::Detail::CVFE::coefficientAtIp(elemDisc, elemVars, ipData, ipCache, conductivity), gradT);
        fluxTerm *= -1.0;
        return fluxTerm;
    }
};

} // end namespace Dumux::Experimental

#endif
