// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup SpatialParameters
 * \ingroup PNMOnePModel
 * \brief Spatial parameters for single-phase flow in the pore network of a sphere packing
 */
#ifndef DUMUX_PNM_EXTRACTION_SPHERE_PACKING_SPATIAL_PARAMS_HH
#define DUMUX_PNM_EXTRACTION_SPHERE_PACKING_SPATIAL_PARAMS_HH

#include <memory>
#include <vector>

#include <dune/common/exceptions.hh>

#include <dumux/porenetwork/1p/spatialparams.hh>

namespace Dumux::PoreNetwork {

/*!
 * \ingroup SpatialParameters
 * \ingroup PNMOnePModel
 * \brief Single-phase spatial parameters with the throat hydraulic radius from the grid data, as required by
 *        TransmissibilityChareyre, and the throat transmissibility if the grid data has the element parameter
 *        ThroatTransmissibility, as required by TransmissibilityFromSpatialParams
 */
template<class GridGeometry, class Scalar>
class SpherePackingOnePSpatialParams
: public OnePSpatialParams<GridGeometry, Scalar, SpherePackingOnePSpatialParams<GridGeometry, Scalar>>
{
    using ParentType = OnePSpatialParams<GridGeometry, Scalar, SpherePackingOnePSpatialParams<GridGeometry, Scalar>>;
    using Element = typename GridGeometry::GridView::template Codim<0>::Entity;

public:
    template<class GridData>
    SpherePackingOnePSpatialParams(std::shared_ptr<const GridGeometry> gridGeometry, const GridData& gridData)
    : ParentType(gridGeometry)
    , hydraulicRadius_(gridGeometry->gridView().size(0))
    {
        const bool hasTransmissibility = gridData.gridHasElementParameter("ThroatTransmissibility");
        if (hasTransmissibility)
            transmissibility_.resize(gridGeometry->gridView().size(0));
        for (const auto& element : elements(gridGeometry->gridView()))
        {
            const auto eIdx = gridGeometry->elementMapper().index(element);
            hydraulicRadius_[eIdx] = gridData.getParameter(element, "ThroatHydraulicRadius");
            if (hasTransmissibility)
                transmissibility_[eIdx] = gridData.getParameter(element, "ThroatTransmissibility");
        }
    }

    template<class ElementVolumeVariables>
    Scalar throatHydraulicRadius(const Element& element, const ElementVolumeVariables& elemVolVars) const
    { return hydraulicRadius_[this->gridGeometry().elementMapper().index(element)]; }

    template<class ElementVolumeVariables>
    Scalar throatTransmissibility(const Element& element, const ElementVolumeVariables& elemVolVars) const
    {
        if (transmissibility_.empty())
            DUNE_THROW(Dune::InvalidStateException, "The grid data has no element parameter ThroatTransmissibility");
        return transmissibility_[this->gridGeometry().elementMapper().index(element)];
    }

private:
    std::vector<Scalar> hydraulicRadius_;
    std::vector<Scalar> transmissibility_;
};

} // end namespace Dumux::PoreNetwork

#endif
