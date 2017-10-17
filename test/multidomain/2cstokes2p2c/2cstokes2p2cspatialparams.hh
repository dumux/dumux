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
 *
 * \brief Spatial parameters for the
 *        coupling of an isothermal two-component Stokes
 *        and an isothermal two-phase two-component Darcy model.
 */

#ifndef DUMUX_TWOCSTOKES_2P2C_SPATIALPARAMS_HH
#define DUMUX_TWOCSTOKES_2P2C_SPATIALPARAMS_HH

#include <dune/grid/io/file/vtk/common.hh>

#include <dumux/material/spatialparams/implicit.hh>
#include <dumux/material/fluidmatrixinteractions/2p/linearmaterial.hh>
#include <dumux/material/fluidmatrixinteractions/2p/regularizedvangenuchten.hh>
#include <dumux/material/fluidmatrixinteractions/2p/regularizedbrookscorey.hh>
#include <dumux/material/fluidmatrixinteractions/2p/efftoabslaw.hh>

namespace Dumux
{
//forward declaration
template<class TypeTag>
class TwoCStokesTwoPTwoCSpatialParams;

namespace Properties
{
// The spatial parameters TypeTag
NEW_TYPE_TAG(TwoCStokesTwoPTwoCSpatialParams);

// Set the spatial parameters
SET_TYPE_PROP(TwoCStokesTwoPTwoCSpatialParams, SpatialParams,
              TwoCStokesTwoPTwoCSpatialParams<TypeTag>);

// Set the material law parameterized by absolute saturations
SET_TYPE_PROP(TwoCStokesTwoPTwoCSpatialParams,
              MaterialLaw,
              EffToAbsLaw<RegularizedVanGenuchten<typename GET_PROP_TYPE(TypeTag, Scalar)>>);
//               EffToAbsLaw<RegularizedBrooksCorey<typename GET_PROP_TYPE(TypeTag, Scalar)> >);
}


/*!
 * \ingroup ImplicitTestProblems
 * \ingroup TwoPTwoCStokesTwoCModel
 * \brief Definition of the spatial parameters for
 *        the coupling of an isothermal two-component Stokes
 *        and an isothermal two-phase two-component Darcy model.
 */
template<class TypeTag>
class TwoCStokesTwoPTwoCSpatialParams : public ImplicitSpatialParams<TypeTag>
{
    using ParentType = ImplicitSpatialParams<TypeTag>;
    using Scalar = typename GET_PROP_TYPE(TypeTag, Scalar);
    using Problem = typename GET_PROP_TYPE(TypeTag, Problem);
    using ElementSolutionVector = typename GET_PROP_TYPE(TypeTag, ElementSolutionVector);
    using SubControlVolume = typename GET_PROP_TYPE(TypeTag, SubControlVolume);
    using GridView = typename GET_PROP_TYPE(TypeTag, GridView);
    using CoordScalar = typename GridView::ctype;

    enum {
        dim=GridView::dimension,
        dimWorld=GridView::dimensionworld
    };
    using GlobalPosition = Dune::FieldVector<Scalar,dimWorld>;

    using Element = typename GridView::template Codim<0>::Entity;
    using FVElementGeometry = typename GET_PROP_TYPE(TypeTag, FVElementGeometry);

    using MaterialLaw = typename GET_PROP_TYPE(TypeTag, MaterialLaw);
    using MaterialLawParams = typename MaterialLaw::Params;



public:
    // export permeability type
    using PermeabilityType = Scalar;

    /*!
     * \brief Spatial parameters for the
     *        coupling of an isothermal two-component Stokes
     *        and an isothermal two-phase two-component Darcy model.
     *
     * \param gridView The GridView which is used by the problem
     */
    TwoCStokesTwoPTwoCSpatialParams(const Problem& problem, const GridView& gridView)
    : ParentType(problem, gridView)
    {
        porosity_ = GET_RUNTIME_PARAM_FROM_GROUP(TypeTag, Scalar, SpatialParams, Porosity);
        permeability_ = GET_RUNTIME_PARAM_FROM_GROUP(TypeTag, Scalar, SpatialParams, Permeability);
        alphaBJ_ = GET_RUNTIME_PARAM_FROM_GROUP(TypeTag, Scalar, SpatialParams, AlphaBJ);

        // residual saturations
        params_.setSwr(GET_RUNTIME_PARAM_FROM_GROUP(TypeTag, Scalar, SpatialParams, Swr));
        params_.setSnr(GET_RUNTIME_PARAM_FROM_GROUP(TypeTag, Scalar, SpatialParams, Snr));
        // parameters for the vanGenuchten law
        params_.setVgAlpha(GET_RUNTIME_PARAM_FROM_GROUP(TypeTag, Scalar, SpatialParams, VgAlpha));
        params_.setVgn(GET_RUNTIME_PARAM_FROM_GROUP(TypeTag, Scalar, SpatialParams, VgN));
    }

    /*!
     * \brief Returns the intrinsic permeability tensor \f$[m^2]\f$
     *
     * \param element The finite element
     * \param fvGeometry The finite volume geometry of the element
     * \param scvIdx The local index of the sub-control volume
     */
    const PermeabilityType intrinsicPermeability(const Element &element,
                                       const FVElementGeometry &fvGeometry,
                                       const int scvIdx) const
    {
        return permeability_;
    }

    /*!
     * \brief Function for defining the (intrinsic) permeability \f$[m^2]\f$.
     *
     * \param globalPos The global position
     * \return the intrinsic permeability
     */
    PermeabilityType permeabilityAtPos(const GlobalPosition& globalPos) const
    {
        return permeability_;
    }

    /*!
     * \brief Returns the porosity \f$[-]\f$
     *
     * \param element The finite element
     * \param fvGeometry The finite volume geometry of the element
     * \param scvIdx The local index of the sub-control volume
     */
    Scalar porosity(const Element& element,
                    const SubControlVolume& scv,
                    const ElementSolutionVector& elemSol) const
    {
        return porosity_;
    }

    /*!
     * \brief Returns the parameter object for the material law
     *
     * \param element The finite element
     * \param fvGeometry The finite volume geometry of the element
     * \param scvIdx The local index of the sub-control volume
     */
//    const MaterialLawParams& materialLawParams(const Element &element,
//                                               const FVElementGeometry &fvGeometry,
//                                               const int scvIdx) const
    const MaterialLawParams& materialLawParams(const Element &element,
                                               const SubControlVolume &scv,
                                               const ElementSolutionVector &elemSol) const
    {
        return params_;
    }

    /*!
     * \brief Evaluate the Beavers-Joseph coefficient at given position
     *
     * \param globalPos The global position
     *
     * \return Beavers-Joseph coefficient
     */
    Scalar beaversJosephCoeffAtPos(const GlobalPosition &globalPos) const
    {
        return alphaBJ_;
    }

private:
    Scalar permeability_;
    Scalar porosity_;
    Scalar alphaBJ_;
    MaterialLawParams params_;
};

} // end namespace Dumux

#endif // DUMUX_TWOCSTOKES_2P2C_SPATIALPARAMS_HH
