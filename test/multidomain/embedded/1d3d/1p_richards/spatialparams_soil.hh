// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*****************************************************************************
 *   See the file COPYING for full copying permissions.                      *
 *                                                                           *
 *   This program is free software: you can redistribute it and/or modify    *
 *   it under the terms of the GNU General Public License as published by    *
 *   the Free Software Foundation, either version 3 of the License, or       *
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
 * \ingroup EmbeddedTests
 * \brief Definition of the spatial parameters for the tissue problem.
 */

#ifndef DUMUX_SOIL_SPATIAL_PARAMS_HH
#define DUMUX_SOIL_SPATIAL_PARAMS_HH

#include <dumux/common/parameters.hh>
#include <dumux/material/spatialparams/fv.hh>
#include <dumux/material/fluidmatrixinteractions/2p/regularizedvangenuchten.hh>
#include <dumux/material/fluidmatrixinteractions/2p/efftoabslaw.hh>

namespace Dumux {

/*!
 * \ingroup EmbeddedTests
 * \brief Definition of the spatial parameters for the tissue problem.
 */
template<class FVGridGeometry, class Scalar>
class SoilSpatialParams
: public FVSpatialParams<FVGridGeometry,Scalar, SoilSpatialParams<FVGridGeometry, Scalar>>
{
    using ThisType = SoilSpatialParams<FVGridGeometry, Scalar>;
    using ParentType = FVSpatialParams<FVGridGeometry, Scalar, ThisType>;
    using GridView = typename FVGridGeometry::GridView;
    using Element = typename GridView::template Codim<0>::Entity;
    using SubControlVolume = typename FVGridGeometry::SubControlVolume;
    using GlobalPosition = typename Element::Geometry::GlobalCoordinate;

public:
    // export permeability type
    using PermeabilityType = Scalar;
    // export material law type
    using MaterialLaw = EffToAbsLaw<RegularizedVanGenuchten<Scalar>>;
    // export material law params
    using MaterialLawParams = typename MaterialLaw::Params;

    SoilSpatialParams(std::shared_ptr<const FVGridGeometry> fvGridGeometry)
    : ParentType(fvGridGeometry)
    {
        // residual saturations
        materialParams_.setSwr(0.05);
        materialParams_.setSnr(0.0);

        // parameters for the Van Genuchten law
        // alpha and n
        materialParams_.setVgAlpha(2.956e-4);
        materialParams_.setVgn(1.5);

        // perm and poro
        permeability_ = getParam<Scalar>("Soil.SpatialParams.Permeability");
        porosity_ = getParam<Scalar>("Soil.SpatialParams.Porosity");
    }

    /*!
     * \brief Defines the intrinsic permeability \f$\mathrm{[m^2]}\f$.
     *
     * \param element The element
     * \param scv The sub control volume
     * \param elemSol The element solution vector
     */
    template<class ElementSolution>
    PermeabilityType permeability(const Element& element,
                                  const SubControlVolume& scv,
                                  const ElementSolution& elemSol) const
    {
        return permeability_;
    }

    /*!
     * \brief Defines the porosity \f$\mathrm{[-]}\f$.
     *
     * \param element The current finite element
     * \param scv The sub control volume
     * \param elemSol The current element solution vector
     */
    template<class ElementSolution>
    Scalar porosity(const Element& element,
                    const SubControlVolume& scv,
                    const ElementSolution& elemSol) const
    {
        return porosity_;
    }

    /*!
     * \brief Returns the parameters for the material law at a given location
     *
     * This method is not actually required by the Richards model, but provided
     * for the convenience of the RichardsLensProblem
     *
     * \param globalPos A global coordinate vector
     */
    const MaterialLawParams& materialLawParamsAtPos(const GlobalPosition &globalPos) const
    { return materialParams_; }

private:
    MaterialLawParams materialParams_;
    Scalar permeability_;
    Scalar porosity_;
};

} // end namespace Dumux

#endif
