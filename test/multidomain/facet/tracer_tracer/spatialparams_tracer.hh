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
 * \ingroup TracerTests
 * \brief Definition of the spatial parameters for the tracer problem
 */
#ifndef DUMUX_TEST_TPFAFACETCOUPLING_TRACER_SPATIALPARAMS_HH
#define DUMUX_TEST_TPFAFACETCOUPLING_TRACER_SPATIALPARAMS_HH

#include <dumux/porousmediumflow/properties.hh>
#include <dumux/material/spatialparams/fv1p.hh>

namespace Dumux {

/*!
 * \ingroup TracerTests
 * \brief Definition of the spatial parameters for the tracer problem
 */
template<class FVGridGeometry, class Scalar>
class TracerSpatialParams
: public FVSpatialParamsOneP<FVGridGeometry, Scalar,
                             TracerSpatialParams<FVGridGeometry, Scalar>>
{
    using GridView = typename FVGridGeometry::GridView;
    using FVElementGeometry = typename FVGridGeometry::LocalView;
    using SubControlVolume = typename FVElementGeometry::SubControlVolume;
    using SubControlVolumeFace = typename FVElementGeometry::SubControlVolumeFace;
    using Element = typename GridView::template Codim<0>::Entity;
    using ParentType = FVSpatialParamsOneP<FVGridGeometry, Scalar,
                                           TracerSpatialParams<FVGridGeometry, Scalar>>;

    static constexpr bool isBox = FVGridGeometry::discMethod == DiscretizationMethod::box;
    static constexpr int dimWorld = GridView::dimensionworld;
    using GlobalPosition = typename Dune::FieldVector<Scalar, dimWorld>;

public:

    TracerSpatialParams(std::shared_ptr<const FVGridGeometry> fvGridGeometry,
                        const std::vector< std::vector<Scalar> >& volumeFluxes,
                        const std::string& paramGroup = "")
    : ParentType(fvGridGeometry)
    , volumeFlux_(volumeFluxes)
    , porosity_(getParamFromGroup<Scalar>(paramGroup, "SpatialParams.Porosity"))
    {}

    /*!
     * \brief Define the porosity \f$\mathrm{[-]}\f$.
     *
     * \param globalPos The global position
     */
    Scalar porosityAtPos(const GlobalPosition& globalPos) const
    { return porosity_; }

    //! Fluid properties that are spatial params in the tracer model
    //! They can possible vary with space but are usually constants
    //! fluid density
    Scalar fluidDensity(const Element &element,
                        const SubControlVolume& scv) const
    { return 1000; }

    //! fluid molar mass
    Scalar fluidMolarMass(const Element &element,
                          const SubControlVolume& scv) const
    { return 18.0; }

    Scalar fluidMolarMass(const GlobalPosition &globalPos) const
    { return 18.0; }

    //! velocity field
    template<class ElementVolumeVariables>
    Scalar volumeFlux(const Element &element,
                      const FVElementGeometry& fvGeometry,
                      const ElementVolumeVariables& elemVolVars,
                      const SubControlVolumeFace& scvf) const
    {
        return isBox ? volumeFlux_[this->fvGridGeometry().elementMapper().index(element)][scvf.index()]
                     : volumeFlux_[scvf.index()][0];
    }

private:
    std::vector< std::vector<Scalar> > volumeFlux_;
    Scalar porosity_;
};

} // end namespace Dumux

#endif
