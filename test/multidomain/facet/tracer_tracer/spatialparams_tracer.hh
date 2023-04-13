// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup FacetTests
 * \brief Definition of the spatial parameters for the tracer problem.
 */

#ifndef DUMUX_TEST_TPFAFACETCOUPLING_TRACER_SPATIALPARAMS_HH
#define DUMUX_TEST_TPFAFACETCOUPLING_TRACER_SPATIALPARAMS_HH

#include <dumux/porousmediumflow/properties.hh>
#include <dumux/porousmediumflow/fvspatialparams1p.hh>

namespace Dumux {

/*!
 * \ingroup FacetTests
 * \brief Definition of the spatial parameters for the tracer problem.
 */
template<class GridGeometry, class Scalar>
class TracerSpatialParams
: public FVPorousMediumFlowSpatialParamsOneP<GridGeometry, Scalar,
                                         TracerSpatialParams<GridGeometry, Scalar>>
{
    using GridView = typename GridGeometry::GridView;
    using FVElementGeometry = typename GridGeometry::LocalView;
    using SubControlVolume = typename FVElementGeometry::SubControlVolume;
    using SubControlVolumeFace = typename FVElementGeometry::SubControlVolumeFace;
    using Element = typename GridView::template Codim<0>::Entity;
    using ParentType = FVPorousMediumFlowSpatialParamsOneP<GridGeometry, Scalar,
                                                       TracerSpatialParams<GridGeometry, Scalar>>;

    static constexpr bool isBox = GridGeometry::discMethod == DiscretizationMethods::box;
    static constexpr int dimWorld = GridView::dimensionworld;
    using GlobalPosition = typename Dune::FieldVector<Scalar, dimWorld>;

public:

    TracerSpatialParams(std::shared_ptr<const GridGeometry> gridGeometry,
                        const std::vector< std::vector<Scalar> >& volumeFluxes,
                        const std::string& paramGroup = "")
    : ParentType(gridGeometry)
    , volumeFlux_(volumeFluxes)
    , porosity_(getParamFromGroup<Scalar>(paramGroup, "SpatialParams.Porosity"))
    , extrusion_(getParamFromGroup<Scalar>(paramGroup, "SpatialParams.Aperture", 1.0))
    {}

    /*!
     * \brief Defines the porosity \f$\mathrm{[-]}\f$.
     *
     * \param globalPos The global position
     */
    Scalar porosityAtPos(const GlobalPosition& globalPos) const
    { return porosity_; }

    //! Returns the extrusion factor
    Scalar extrusionFactorAtPos(const GlobalPosition& globalPos) const
    { return extrusion_; }

    //! Fluid properties that are spatial parameters in the tracer model.
    //! They can possibly vary with space but are usually constants.
    //! Fluid density
    Scalar fluidDensity(const Element &element,
                        const SubControlVolume& scv) const
    { return 1000; }

    //! Fluid molar mass
    Scalar fluidMolarMass(const Element &element,
                          const SubControlVolume& scv) const
    { return 18.0; }

    Scalar fluidMolarMass(const GlobalPosition &globalPos) const
    { return 18.0; }

    //! Velocity field
    template<class ElementVolumeVariables>
    Scalar volumeFlux(const Element &element,
                      const FVElementGeometry& fvGeometry,
                      const ElementVolumeVariables& elemVolVars,
                      const SubControlVolumeFace& scvf) const
    {
        return isBox ? volumeFlux_[this->gridGeometry().elementMapper().index(element)][scvf.index()]
                     : volumeFlux_[scvf.index()][0];
    }

private:
    std::vector< std::vector<Scalar> > volumeFlux_;
    Scalar porosity_;
    Scalar extrusion_;
};

} // end namespace Dumux

#endif
