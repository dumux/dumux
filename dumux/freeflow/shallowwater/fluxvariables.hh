// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup ShallowWaterModels
 * \copydoc Dumux::ShallowWaterFluxVariables
 */
#ifndef DUMUX_FREEFLOW_SHALLOW_WATER_FLUXVARIABLES_HH
#define DUMUX_FREEFLOW_SHALLOW_WATER_FLUXVARIABLES_HH

#include <dumux/common/properties.hh>
#include <dumux/common/numeqvector.hh>
#include <dumux/flux/fluxvariablesbase.hh>

namespace Dumux {

/*!
 * \ingroup ShallowWaterModels
 * \brief The flux variables class for the shallow water model.
 *
 */
template<class TypeTag>
class ShallowWaterFluxVariables
: public FluxVariablesBase<GetPropType<TypeTag, Properties::Problem>,
                           typename GetPropType<TypeTag, Properties::GridGeometry>::LocalView,
                           typename GetPropType<TypeTag, Properties::GridVolumeVariables>::LocalView,
                           typename GetPropType<TypeTag, Properties::GridFluxVariablesCache>::LocalView>
{
    using Problem = GetPropType<TypeTag, Properties::Problem>;
    using ModelTraits = GetPropType<TypeTag, Properties::ModelTraits>;
    using NumEqVector = Dumux::NumEqVector<GetPropType<TypeTag, Properties::PrimaryVariables>>;
    using AdvectionType = GetPropType<TypeTag, Properties::AdvectionType>;
    using ViscousFluxType = GetPropType<TypeTag, Properties::ViscousFluxType>;

    using GridVariables = GetPropType<TypeTag, Properties::GridVariables>;
    using GridVolumeVariables = typename GridVariables::GridVolumeVariables;
    using GridGeometry = GetPropType<TypeTag, Properties::GridGeometry>;
    using FVElementGeometry = typename GridGeometry::LocalView;
    using GridView = typename GridGeometry::GridView;

    using Element = typename GridView::template Codim<0>::Entity;
    using ElementVolumeVariables = typename GridVolumeVariables::LocalView;
    using SubControlVolumeFace = typename FVElementGeometry::SubControlVolumeFace;

    static constexpr bool enableAdvection = ModelTraits::enableAdvection();

public:

    /*!
     * \brief Returns the advective flux computed by the Riemann solver
     *
     */
    NumEqVector advectiveFlux(const Problem& problem,
                              const Element& element,
                              const FVElementGeometry& fvGeometry,
                              const ElementVolumeVariables& elemVolVars,
                              const SubControlVolumeFace& scvf) const
    {
        if (enableAdvection)
            return AdvectionType::flux(problem, element, fvGeometry, elemVolVars, scvf);

        return NumEqVector(0.0);
    }

    /*!
     * \brief Returns the viscous momentum flux
     *
     */
    NumEqVector viscousFlux(const Problem& problem,
                            const Element& element,
                            const FVElementGeometry& fvGeometry,
                            const ElementVolumeVariables& elemVolVars,
                            const SubControlVolumeFace& scvf) const
    {
        // Add viscous momentum flux
        return ViscousFluxType::flux(problem, element, fvGeometry, elemVolVars, scvf);
    }
};

} // end namespace Dumux

#endif
