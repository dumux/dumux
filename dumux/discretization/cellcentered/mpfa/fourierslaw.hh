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
 * \ingroup CCMpfaDiscretization
 * \brief Fourier's law for cell-centered finite volume schemes with multi-point flux approximation
 */
#ifndef DUMUX_DISCRETIZATION_CC_MPFA_FOURIERS_LAW_HH
#define DUMUX_DISCRETIZATION_CC_MPFA_FOURIERS_LAW_HH

#include <dumux/common/properties.hh>
#include <dumux/discretization/methods.hh>

#include <dumux/discretization/methods.hh>

namespace Dumux
{
//! forward declaration of the method-specific implementation
template<class TypeTag, DiscretizationMethods discMethod>
class FouriersLawImplementation;

/*!
* \ingroup CCMpfaDiscretization
* \brief Fourier's law for cell-centered finite volume schemes with two-point flux approximation
*/
template <class TypeTag>
class FouriersLawImplementation<TypeTag, DiscretizationMethods::CCMpfa>
{
    using Scalar = typename GET_PROP_TYPE(TypeTag, Scalar);
    using Problem = typename GET_PROP_TYPE(TypeTag, Problem);
    using GridView = typename GET_PROP_TYPE(TypeTag, GridView);
    using Element = typename GridView::template Codim<0>::Entity;
    using SubControlVolumeFace = typename GET_PROP_TYPE(TypeTag, SubControlVolumeFace);
    using FVElementGeometry = typename GET_PROP_TYPE(TypeTag, FVElementGeometry);
    using ElementVolumeVariables = typename GET_PROP_TYPE(TypeTag, ElementVolumeVariables);
    using ElementFluxVarsCache = typename GET_PROP_TYPE(TypeTag, ElementFluxVariablesCache);
    using FluxVariablesCache = typename GET_PROP_TYPE(TypeTag, FluxVariablesCache);
    using ThermalConductivityModel = typename GET_PROP_TYPE(TypeTag, ThermalConductivityModel);

    //! Class that fills the cache corresponding to mpfa Darcy's Law
    class MpfaFouriersLawCacheFiller
    {
    public:
        //! Function to fill an MpfaDarcysLawCache of a given scvf
        //! This interface has to be met by any cache filler class for heat conduction quantities
        template<class FluxVariablesCacheFiller>
        static void fill(FluxVariablesCache& scvfFluxVarsCache,
                         const Problem& problem,
                         const Element& element,
                         const FVElementGeometry& fvGeometry,
                         const ElementVolumeVariables& elemVolVars,
                         const SubControlVolumeFace& scvf,
                         const FluxVariablesCacheFiller& fluxVarsCacheFiller)
        {
          // get interaction volume from the flux vars cache filler & upate the cache
          if (fvGeometry.fvGridGeometry().vertexUsesSecondaryInteractionVolume(scvf.vertexIndex()))
              scvfFluxVarsCache.updateHeatConduction(fluxVarsCacheFiller.secondaryInteractionVolume(),
                                                     fluxVarsCacheFiller.secondaryIvLocalFaceData(),
                                                     fluxVarsCacheFiller.secondaryIvDataHandle(),
                                                     scvf);
          else
              scvfFluxVarsCache.updateHeatConduction(fluxVarsCacheFiller.primaryInteractionVolume(),
                                                     fluxVarsCacheFiller.primaryIvLocalFaceData(),
                                                     fluxVarsCacheFiller.primaryIvDataHandle(),
                                                     scvf);
        }
    };

    //! The cache used in conjunction with the mpfa Fourier's Law
    class MpfaFouriersLawCache
    {
        // In the current implementation of the flux variables cache we cannot
        // make a disctinction between dynamic (mpfa-o) and static (mpfa-l)
        // matrix and vector types, as currently the cache class can only be templated
        // by a type tag (and there can only be one). We use a dynamic vector here to
        // make sure it works in case one of the two used interaction volume types uses
        // dynamic types performance is thus lowered for schemes using static types.
        // TODO: this has to be checked thoroughly as soon as a scheme using static types
        //       is implemented. One idea to overcome the performance drop could be only
        //       storing the iv-local index here and obtain tij always from the datahandle
        //       of the fluxVarsCacheContainer
        using GridIndexType = typename GridView::IndexSet::IndexType;
        using Vector = Dune::DynamicVector< Scalar >;
        using Matrix = Dune::DynamicMatrix< Scalar >;
        using Stencil = std::vector< GridIndexType >;

        using PrimaryInteractionVolume = typename GET_PROP_TYPE(TypeTag, PrimaryInteractionVolume);
        using PrimaryIvVector = typename PrimaryInteractionVolume::Traits::Vector;
        using PrimaryIvMatrix = typename PrimaryInteractionVolume::Traits::Matrix;
        using PrimaryStencil = typename PrimaryInteractionVolume::Traits::Stencil;

        static_assert( std::is_convertible<PrimaryIvVector*, Vector*>::value,
                       "The vector type used in primary interaction volumes is not convertible to Dune::DynamicVector!" );
        static_assert( std::is_convertible<PrimaryIvMatrix*, Matrix*>::value,
                       "The matrix type used in primary interaction volumes is not convertible to Dune::DynamicMatrix!" );
        static_assert( std::is_convertible<PrimaryStencil*, Stencil*>::value,
                       "The stencil type used in primary interaction volumes is not convertible to std::vector<GridIndexType>!" );

        using SecondaryInteractionVolume = typename GET_PROP_TYPE(TypeTag, SecondaryInteractionVolume);
        using SecondaryIvVector = typename SecondaryInteractionVolume::Traits::Vector;
        using SecondaryIvMatrix = typename SecondaryInteractionVolume::Traits::Matrix;
        using SecondaryStencil = typename SecondaryInteractionVolume::Traits::Stencil;

        static_assert( std::is_convertible<SecondaryIvVector*, Vector*>::value,
                       "The vector type used in secondary interaction volumes is not convertible to Dune::DynamicVector!" );
        static_assert( std::is_convertible<SecondaryIvMatrix*, Matrix*>::value,
                       "The matrix type used in secondary interaction volumes is not convertible to Dune::DynamicMatrix!" );
        static_assert( std::is_convertible<SecondaryStencil*, Stencil*>::value,
                       "The stencil type used in secondary interaction volumes is not convertible to std::vector<GridIndexType>!" );

        static constexpr int dim = GridView::dimension;
        static constexpr int dimWorld = GridView::dimensionworld;

    public:
        // export filler type
        using Filler = MpfaFouriersLawCacheFiller;

        /*!
         * \brief Update cached objects (transmissibilities)
         *
         * \tparam InteractionVolume The (mpfa scheme-specific) interaction volume
         * \tparam LocalFaceData The object used to store iv-local info on an scvf
         * \tparam DataHandle The object used to store transmissibility matrices etc.
         *
         * \param iv The interaction volume this scvf is embedded in
         * \param localFaceData iv-local info on this scvf
         * \param dataHandle Transmissibility matrix & gravity data of this iv
         * \param scvf The sub-control volume face
         */
        template<class InteractionVolume, class LocalFaceData, class DataHandle>
        void updateHeatConduction(const InteractionVolume& iv,
                                  const LocalFaceData& localFaceData,
                                  const DataHandle& dataHandle,
                                  const SubControlVolumeFace &scvf)
        {
            stencil_ = &iv.stencil();
            switchFluxSign_ = localFaceData.isOutside();

            // store pointer to the temperature vector of this iv
            Tj_ = &dataHandle.temperatures();

            const auto ivLocalIdx = localFaceData.ivLocalScvfIndex();
            if (dim == dimWorld)
                Tij_ = &dataHandle.heatConductionT()[ivLocalIdx];
            else
                Tij_ = localFaceData.isOutside() ? &dataHandle.heatConductionTout()[ivLocalIdx][localFaceData.scvfLocalOutsideScvfIndex()]
                                                 : &dataHandle.heatConductionT()[ivLocalIdx];
        }

        //! Coefficients for the cell (& Dirichlet) unknowns in flux expressions
        const Vector& heatConductionTij() const { return *Tij_; }

        //! The stencil corresponding to the transmissibilities
        const Stencil& heatConductionStencil() const { return *stencil_; }

        //! In the interaction volume-local system of eq we have one unknown per face.
        //! On scvfs on this face, but in "outside" (neighbor) elements of it, we have
        //! to take the negative value of the fluxes due to the flipped normal vector.
        //! This function returns whether or not this scvf is an "outside" face in the iv.
        bool heatConductionSwitchFluxSign() const { return switchFluxSign_; }

        //! The cell (& Dirichlet) temperatures within this interaction volume
        const Vector& temperatures() const { return *Tj_; }

    private:
        bool switchFluxSign_;
        const Stencil* stencil_; //!< The stencil, i.e. the grid indices j
        const Vector* Tij_;      //!< The transmissibilities such that f = Tij*Tj
        const Vector* Tj_;       //!< The interaction-volume wide temperature Tj
    };

public:
    // state the discretization method this implementation belongs to
    static const DiscretizationMethods myDiscretizationMethod = DiscretizationMethods::CCMpfa;

    // state the type for the corresponding cache and its filler
    using Cache = MpfaFouriersLawCache;

    //! Compute the conductive flux across an scvf
    static Scalar flux(const Problem& problem,
                       const Element& element,
                       const FVElementGeometry& fvGeometry,
                       const ElementVolumeVariables& elemVolVars,
                       const SubControlVolumeFace& scvf,
                       const ElementFluxVarsCache& elemFluxVarsCache)
    {
        const auto& fluxVarsCache = elemFluxVarsCache[scvf];
        const auto& tij = fluxVarsCache.heatConductionTij();
        const auto& Tj = fluxVarsCache.temperatures();

        // compute Tij*tj
        Scalar flux = tij*Tj;
        if (fluxVarsCache.heatConductionSwitchFluxSign())
            flux *= -1.0;

        return flux;
    }
};

} // end namespace Dumux

#endif
