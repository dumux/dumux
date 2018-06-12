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
#ifndef DUMUX_TWOP_TRACER_TEST_SPATIAL_PARAMS_HH
#define DUMUX_TWOP_TRACER_TEST_SPATIAL_PARAMS_HH

#include <dumux/material/spatialparams/fv.hh>
#include <iostream>

namespace Dumux {

/*!
 * \ingroup TracerTests
 * \brief Definition of the spatial parameters for the tracer problem
 */
template<class TypeTag>
class TwoPTracerTestSpatialParams
: public FVSpatialParams<typename GET_PROP_TYPE(TypeTag, FVGridGeometry),
                             typename GET_PROP_TYPE(TypeTag, Scalar),
                             TwoPTracerTestSpatialParams<TypeTag>>
{
    using Scalar = typename GET_PROP_TYPE(TypeTag, Scalar);
    using FVGridGeometry = typename GET_PROP_TYPE(TypeTag, FVGridGeometry);
    using GridView = typename FVGridGeometry::GridView;
    using FVElementGeometry = typename FVGridGeometry::LocalView;
    using SubControlVolume = typename FVElementGeometry::SubControlVolume;
    using SubControlVolumeFace = typename FVElementGeometry::SubControlVolumeFace;
    using Element = typename GridView::template Codim<0>::Entity;
    using ParentType = FVSpatialParams<FVGridGeometry, Scalar, TwoPTracerTestSpatialParams<TypeTag>>;
    using ModelTraits = typename GET_PROP_TYPE(TypeTag, ModelTraits);


    //static const int numPhases = ModelTraits::numPhases();
    //static const int numScvf = FVGridGeometry::numScvf();
    //using FieldMatrix = Dune::FieldMatrix <Scalar, numPhases, numScvf>;
    //using PhaseVector = std::vector<Scalar>;
    //using FieldVector = Dune::FieldVector<PhaseVector, 2>;

    static const int numPhases = ModelTraits::numPhases();
    using ScvfVector  = std::vector<Scalar>;
    using FieldVector = std::vector<ScvfVector, numPhases>;

    static const int dimWorld = GridView::dimensionworld;
    using GlobalPosition = typename Dune::FieldVector<Scalar, dimWorld>;

    //static const int numPhases = 2;
    //using numPhases = typename GET_PROP_TYPE(TypeTag, ModelTraits)::numPhases();
    //using numScvf = typename GET_PROP_TYPE(TypeTag, FVGridGeometry)::numScvf();
    //static const int numScvf = fvGeometry->numScvf();

      // number of rows
      //static const int rows = Jacobian::rows;
      //static_assert((rows == dimRange), "Number of rows and range dimension do not coincide");
      // number of cols
      //static const int cols = Jacobian::cols;
      //static_assert((cols == dimDomain), "Number of columns and domain dimension do not coincide");



public:

    TwoPTracerTestSpatialParams(std::shared_ptr<const FVGridGeometry> fvGridGeometry)
    : ParentType(fvGridGeometry) {}

    /*!
     * \brief Define the porosity \f$\mathrm{[-]}\f$.
     *
     * \param globalPos The global position
     */
    Scalar porosityAtPos(const GlobalPosition& globalPos) const
    { return 0.4; }

    /*!
     * \brief Define the dispersivity.
     *
     * \param element The finite element
     * \param scv The sub-control volume
     * \param elemSol The solution for all dofs of the element
     */
    template<class ElementSolution>
    Scalar dispersivity(const Element &element,
                        const SubControlVolume& scv,
                        const ElementSolution& elemSol) const
    { return 0; }

    //! Fluid properties that are spatial params in the tracer model
    //! They can possible vary with space but are usually constants

    //! fluid density

    Scalar fluidDensity(const int phaseIdx,
                        const Element &element,
                        const SubControlVolume& scv) const
    {
        if (phaseIdx == 0)
            return 1000;
        else
            return 1000;
    }

    // TODO density und molar mass (x2) für verschiedene Phasen differenzieren!

    //! fluid molar mass
    Scalar fluidMolarMass(const int phaseIdx,
                          const Element &element,
                          const SubControlVolume& scv) const
    {
        if (phaseIdx == 0)
            return 18.0;
        else
            return 131.0;
    }

    Scalar fluidMolarMass(const int phaseIdx,
                          const GlobalPosition &globalPos) const
    { if (phaseIdx == 0)
            return 18.0;
      else
            return 18.0;
    }

    //! velocity field
    template<class ElementVolumeVariables>
    Scalar volumeFlux(const Element &element,
                      const FVElementGeometry& fvGeometry,
                      const ElementVolumeVariables& elemVolVars,
                      const SubControlVolumeFace& scvf) const
    {
        return volumeFlux_[scvf.index()];
    }

    void setVolumeFlux(const FieldVector& f)
    { volumeFlux_ = f; }

private:
    FieldVector volumeFlux_;
    //std::vector<Scalar> volumeFlux_;
};

} // end namespace Dumux

#endif
