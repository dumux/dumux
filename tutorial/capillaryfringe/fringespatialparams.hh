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
 * \brief spatial parameters for the TwoPTwoCTestProblem
 */
#ifndef DUMUX_FRINGE_SPATIAL_PARAMETERS_HH
#define DUMUX_FRINGE_SPATIAL_PARAMETERS_HH

#include <dumux/material/spatialparams/implicit.hh>
#include <dumux/material/fluidmatrixinteractions/2p/regularizedvangenuchten.hh>
#include <dumux/material/fluidmatrixinteractions/2p/linearmaterial.hh>
#include <dumux/material/fluidmatrixinteractions/2p/efftoabslaw.hh>

namespace Dumux
{
template<class TypeTag>
class FringeSpatialParams;

namespace Properties
{
// The spatial parameters TypeTag
NEW_TYPE_TAG(FringeSpatialParams);

// Set the spatial parameters
SET_TYPE_PROP(FringeSpatialParams, SpatialParams, FringeSpatialParams<TypeTag>);

// Set the material law
SET_PROP(FringeSpatialParams, MaterialLaw)
{
private:
    // define the material law which is parameterized by effective
    // saturations
    using Scalar = typename GET_PROP_TYPE(TypeTag, Scalar);
public:
    // define the material law parameterized by absolute saturations
    using type = EffToAbsLaw<RegularizedVanGenuchten<Scalar>>;
};
} // end namespace Properties

/*!
 * \ingroup TwoPTwoCModel
 * \ingroup ImplicitTestProblems
 * \brief The spatial parameters for the TwoPTwoCTestProblem
 */
template<class TypeTag>
class FringeSpatialParams : public ImplicitSpatialParams<TypeTag>
{
    using ParentType = ImplicitSpatialParams<TypeTag>;
    using Problem = typename GET_PROP_TYPE(TypeTag, Problem);
    using GridView = typename GET_PROP_TYPE(TypeTag, GridView);
    using Scalar = typename GET_PROP_TYPE(TypeTag, Scalar);

    enum {
        dim=GridView::dimension,
        dimWorld=GridView::dimensionworld
    };

    using GlobalPosition = Dune::FieldVector<Scalar, dimWorld>;
    using MaterialLaw = typename GET_PROP_TYPE(TypeTag, MaterialLaw);
    using MaterialLawParams = typename MaterialLaw::Params;

public:
    // export permeability type
    using PermeabilityType = Scalar;

    /*!
     * \brief Constructor
     *
     * \param gridView The DUNE GridView representing the spatial
     *                 domain of the problem.
     */
    FringeSpatialParams(const Problem& problem, const GridView& gridView)
        : ParentType(problem, gridView)
    {

        lensLowerLeft_ = GET_RUNTIME_PARAM(TypeTag, GlobalPosition, LensSpatialParams.LowerLeft);
        lensUpperRight_ = GET_RUNTIME_PARAM(TypeTag, GlobalPosition, LensSpatialParams.UpperRight);

        // residual saturations and parameters for the Van Genuchten law
        materialParams_.setSwr(GET_RUNTIME_PARAM(TypeTag, Scalar, SpatialParams.Swr));
        materialParams_.setSnr(GET_RUNTIME_PARAM(TypeTag, Scalar, SpatialParams.Snr));
        materialParams_.setVgAlpha(GET_RUNTIME_PARAM(TypeTag, Scalar, SpatialParams.VanGenAlpha));
        materialParams_.setVgn(GET_RUNTIME_PARAM(TypeTag, Scalar, SpatialParams.VanGenN));
        materialParams_.setPcLowSw(GET_RUNTIME_PARAM(TypeTag, Scalar, SpatialParams.RegThresholdSw));

        // same for lens: residual saturations and parameters for the Van Genuchten law
        materialParamsLens_.setSwr(GET_RUNTIME_PARAM(TypeTag, Scalar, LensSpatialParams.Swr));
        materialParamsLens_.setSnr(GET_RUNTIME_PARAM(TypeTag, Scalar, LensSpatialParams.Snr));
        materialParamsLens_.setVgAlpha(GET_RUNTIME_PARAM(TypeTag, Scalar, LensSpatialParams.VanGenAlpha));
        materialParamsLens_.setVgn(GET_RUNTIME_PARAM(TypeTag, Scalar, LensSpatialParams.VanGenN));
        materialParamsLens_.setPcLowSw(GET_RUNTIME_PARAM(TypeTag, Scalar, LensSpatialParams.RegThresholdSw));

        permeability_ = GET_RUNTIME_PARAM(TypeTag, Scalar, SpatialParams.Permeability);
        lensPermeability_ = GET_RUNTIME_PARAM(TypeTag, Scalar, LensSpatialParams.Permeability);
        porosity_ = GET_RUNTIME_PARAM(TypeTag, Scalar, SpatialParams.Porosity);
        lensPorosity_ = GET_RUNTIME_PARAM(TypeTag, Scalar, LensSpatialParams.Porosity);
        hasLens_ = GET_RUNTIME_PARAM(TypeTag, bool, LensSpatialParams.HasLens);
    }

    /*!
     * \brief Returns the intrinsic permeability tensor [m^2] at a given location
     *
     * \param globalPos The global position where we evaluate
     */
    PermeabilityType permeabilityAtPos(const GlobalPosition& globalPos) const
    {
        if (isInLens_(globalPos))
            return lensPermeability_;
        else
            return permeability_;

    }

    /*!
     * \brief Returns the porosity [] at a given location
     *
     * \param globalPos The global position where we evaluate
     */
    Scalar porosityAtPos(const GlobalPosition& globalPos) const
    {
        if (isInLens_(globalPos))
            return lensPorosity_;
        else
            return porosity_;
    }

    /*!
     * \brief Returns the parameters for the material law at a given location
     *
     * This method is not actually required by the TwoPTwoC model, but provided
     * for the convenience of the TwoPTwoCLensProblem
     *
     * \param globalPos A global coordinate vector
     */
    const MaterialLawParams& materialLawParamsAtPos(const GlobalPosition &globalPos) const
    {
         if (isInLens_(globalPos))
             return materialParamsLens_;
         else
            return materialParams_;
    }

        /*!
     * \brief Returns the heat capacity \f$[J / (kg K)]\f$ of the rock matrix.
     *
     * This is only required for non-isothermal models.
     *
     * \param globalPos The global position
     */
    Scalar solidHeatCapacityAtPos(const GlobalPosition& globalPos) const
    { return 1480; /*specific heat capacity of granite [J / (kg K)]*/ }

    /*!
     * \brief Returns the mass density \f$[kg / m^3]\f$ of the rock matrix.
     *
     * This is only required for non-isothermal models.
     *
     * \param globalPos The global position
     */
    Scalar solidDensityAtPos(const GlobalPosition& globalPos) const
    { return 1600; /*density of granite [kg/m^3]*/ }

    /*!
     * \brief Returns the thermal conductivity \f$\mathrm{[W/(m K)]}\f$ of the solid
     *
     * This is only required for non-isothermal models.
     *
     * \param globalPos The global position
     */
    Scalar solidThermalConductivityAtPos(const GlobalPosition& globalPos) const
    { return 2; }

private:
    bool isInLens_(const GlobalPosition &globalPos) const
    {
        if(!hasLens_) return false;

        for (int i = 0; i < dimWorld; ++i)
            if (globalPos[i] < lensLowerLeft_[i] - eps_ || globalPos[i] > lensUpperRight_[i] + eps_)
                return false;

        return true;
    }
    static constexpr Scalar eps_ = 1e-7;

    MaterialLawParams materialParams_;
    MaterialLawParams materialParamsLens_;
    Scalar permeability_, lensPermeability_;
    Scalar porosity_, lensPorosity_;

    GlobalPosition lensLowerLeft_;
    GlobalPosition lensUpperRight_;

    bool hasLens_;
};

} // end namespace Dumux

#endif
