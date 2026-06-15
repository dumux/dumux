// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup BedloadTransportModel
 * \copydoc Dumux::BedloadFluxVariables
 */
#ifndef DUMUX_BEDLOAD_FLUXVARIABLES_HH
#define DUMUX_BEDLOAD_FLUXVARIABLES_HH

#include <dumux/flux/fluxvariablesbase.hh>
#include <dumux/flux/bedload/bedloadflux.hh>

namespace Dumux {

/*!
 * \ingroup BedloadTransportModel
 * \brief The flux variables class for the bedload transport model.
 *
 */
template<class TypeTag>
class BedloadFluxVariables: public FluxVariablesBase<GetPropType<TypeTag, Properties::Problem>,
                                                     typename GetPropType<TypeTag, Properties::GridGeometry>::LocalView,
                                                     typename GetPropType<TypeTag, Properties::GridVolumeVariables>::LocalView,
                                                     typename GetPropType<TypeTag, Properties::GridFluxVariablesCache>::LocalView>
{
    using Problem = GetPropType<TypeTag, Properties::Problem>;
    using ModelTraits = GetPropType<TypeTag, Properties::ModelTraits>;
    using PrimaryVariables = GetPropType<TypeTag, Properties::PrimaryVariables>;
    using NumEqVector = Dumux::NumEqVector<PrimaryVariables>;

    using GridVariables = GetPropType<TypeTag, Properties::GridVariables>;
    using GridVolumeVariables = typename GridVariables::GridVolumeVariables;
    using GridGeometry = GetPropType<TypeTag, Properties::GridGeometry>;
    using FVElementGeometry = typename GridGeometry::LocalView;
    using GridView = typename GridGeometry::GridView;

    using Element = typename GridView::template Codim<0>::Entity;
    using ElementVolumeVariables = typename GridVolumeVariables::LocalView;
    using SubControlVolumeFace = typename FVElementGeometry::SubControlVolumeFace;

    static constexpr bool enableBedload = ModelTraits::enableBedload();

public:

    /*!
     * \brief Returns the bedload flux computed by the local Lax-Friedrichs approach
     */
    NumEqVector bedloadFlux(const Problem& problem,
                            const Element& element,
                            const FVElementGeometry& fvGeometry,
                            const ElementVolumeVariables& elemVolVars,
                            const SubControlVolumeFace& scvf) const
    {
        if (enableBedload)
            return BedloadFlux<NumEqVector>::flux(problem, element, fvGeometry, elemVolVars, scvf);

        return NumEqVector(0.0);
    }
};

} // end namespace Dumux

#endif
