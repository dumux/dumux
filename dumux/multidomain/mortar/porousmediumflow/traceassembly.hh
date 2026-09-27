// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup MortarCoupling
 * \brief The traces a porous-medium subdomain hands back to the mortar domain.
 */
#ifndef DUMUX_MULTIDOMAIN_MORTAR_POROUSMEDIUMFLOW_TRACE_ASSEMBLY_HH
#define DUMUX_MULTIDOMAIN_MORTAR_POROUSMEDIUMFLOW_TRACE_ASSEMBLY_HH

#include <type_traits>

#include <dune/common/fvector.hh>

#include <dumux/common/math.hh>
#include <dumux/multidomain/mortar/traceassembly.hh>

namespace Dumux::Mortar {

/*!
 * \ingroup MortarCoupling
 * \brief The mass flux of a porous-medium phase through a boundary sub-control volume face,
 *        outward.
 */
template<typename FluxVariables>
auto darcyMassFlux(int phaseIdx = 0)
{
    return traceIntegrand<TraceEntity::subControlVolumeFace>(
        [phaseIdx] (const auto& context, const auto& scvf)
        {
            FluxVariables fluxVars;
            fluxVars.init(
                context.problem(),
                context.element(), context.fvGeometry(),
                context.elemVolVars(), scvf, context.elemFluxVarsCache()
            );
            const auto flux = fluxVars.advectiveFlux(phaseIdx, [phaseIdx] (const auto& volVars) { return volVars.density(phaseIdx)*volVars.mobility(phaseIdx); });
            return Dune::FieldVector<std::decay_t<decltype(flux)>, 1>{flux};
        }
    );
}

/*!
 * \ingroup MortarCoupling
 * \brief The pressure on a boundary sub-control volume face of a cell-centred scheme,
 *        reconstructed by two points from the flux imposed on that face, integrated over
 *        the face.
 *
 * Where the mortar datum is the flux, the conjugate trace is the pressure on the face. The
 * cell value alone differs from it by one cell-to-face distance times the gradient the
 * imposed flux implies, which biases the pressure-continuity residual of the mortar. The
 * gradient follows from the cell's permeability in the direction of the face normal.
 *
 * \note The integrand reads the coupling manager by reference, so it must not outlive it;
 *       it is meant to be passed straight to assembleTrace.
 */
template<typename CouplingManager>
auto reconstructedFacePressure(const CouplingManager& couplingManager, int phaseIdx = 0)
{
    return traceIntegrand<TraceEntity::subControlVolumeFace>(
        [&couplingManager, phaseIdx] (const auto& context, const auto& scvf)
        {
            using Extrusion = typename std::decay_t<decltype(context.fvGeometry().gridGeometry())>::Extrusion;
            const auto& scv = context.fvGeometry().scv(scvf.insideScvIdx());
            const auto& volVars = context.elemVolVars()[scvf.insideScvIdx()];
            const auto distance = (scvf.ipGlobal() - scv.center()).two_norm();
            const auto& normal = scvf.unitOuterNormal();
            const auto normalPermeability = vtmv(normal, volVars.permeability(), normal);
            const auto flux = couplingManager.traceAt(context.element(), scvf)[0];
            const auto mobility = volVars.density(phaseIdx)*volVars.mobility(phaseIdx);
            const auto pressure = volVars.pressure(phaseIdx) - flux*distance/(mobility*normalPermeability);
            return Dune::FieldVector<std::decay_t<decltype(pressure)>, 1>{pressure*Extrusion::area(context.fvGeometry(), scvf)};
        }
    );
}

} // end namespace Dumux::Mortar

#endif
