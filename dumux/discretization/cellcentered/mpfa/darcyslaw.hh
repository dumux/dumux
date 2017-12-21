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
 * \brief Darcy's Law for cell-centered finite volume schemes
 *        with multi-point flux approximation.
 */
#ifndef DUMUX_DISCRETIZATION_CC_MPFA_DARCYS_LAW_HH
#define DUMUX_DISCRETIZATION_CC_MPFA_DARCYS_LAW_HH

#include <dune/common/densevector.hh>
#include <dune/common/densematrix.hh>

#include <dumux/common/properties.hh>
#include <dumux/common/parameters.hh>
#include <dumux/discretization/methods.hh>

#include <dumux/discretization/methods.hh>

namespace Dumux
{
//! forward declaration of the method-specific implementation
template<class TypeTag, DiscretizationMethods discMethod>
class DarcysLawImplementation;

/*!
 * \ingroup CCMpfaDiscretization
 * \brief Darcy's law for cell-centered finite volume schemes with multi-point flux approximation.
 */
template <class TypeTag>
class DarcysLawImplementation<TypeTag, DiscretizationMethods::CCMpfa>
{
    using Scalar = typename GET_PROP_TYPE(TypeTag, Scalar);
    using Problem = typename GET_PROP_TYPE(TypeTag, Problem);
    using GridView = typename GET_PROP_TYPE(TypeTag, GridView);
    using Element = typename GridView::template Codim<0>::Entity;

    using FVElementGeometry = typename GET_PROP_TYPE(TypeTag, FVElementGeometry);
    using SubControlVolumeFace = typename GET_PROP_TYPE(TypeTag, SubControlVolumeFace);
    using ElementVolumeVariables = typename GET_PROP_TYPE(TypeTag, ElementVolumeVariables);
    using ElementFluxVariablesCache = typename GET_PROP_TYPE(TypeTag, ElementFluxVariablesCache);
    using FluxVariablesCache = typename GET_PROP_TYPE(TypeTag, FluxVariablesCache);

    //! Class that fills the cache corresponding to mpfa Darcy's Law
    class MpfaDarcysLawCacheFiller
    {
    public:
        //! Function to fill an MpfaDarcysLawCache of a given scvf
        //! This interface has to be met by any advection-related cache filler class
        template<class FluxVariablesCacheFiller>
        static void fill(FluxVariablesCache& scvfFluxVarsCache,
                         const Problem& problem,
                         const Element& element,
                         const FVElementGeometry& fvGeometry,
                         const ElementVolumeVariables& elemVolVars,
                         const SubControlVolumeFace& scvf,
                         const FluxVariablesCacheFiller& fluxVarsCacheFiller)
        {
            // get interaction volume related data from the filler class & upate the cache
            if (fvGeometry.fvGridGeometry().vertexUsesSecondaryInteractionVolume(scvf.vertexIndex()))
                scvfFluxVarsCache.updateAdvection(fluxVarsCacheFiller.secondaryInteractionVolume(),
                                                  fluxVarsCacheFiller.secondaryIvLocalFaceData(),
                                                  fluxVarsCacheFiller.secondaryIvDataHandle(),
                                                  scvf);
            else
                scvfFluxVarsCache.updateAdvection(fluxVarsCacheFiller.primaryInteractionVolume(),
                                                  fluxVarsCacheFiller.primaryIvLocalFaceData(),
                                                  fluxVarsCacheFiller.primaryIvDataHandle(),
                                                  scvf);
        }
    };

    //! The cache used in conjunction with the mpfa Darcy's Law
    class MpfaDarcysLawCache
    {
        static constexpr int dim = GridView::dimension;
        static constexpr int dimWorld = GridView::dimensionworld;
        static constexpr int numPhases = GET_PROP_VALUE(TypeTag, NumPhases);

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
        using Vector = Dune::DynamicVector< Scalar >;
        using Matrix = Dune::DynamicMatrix< Scalar >;

        using PrimaryInteractionVolume = typename GET_PROP_TYPE(TypeTag, PrimaryInteractionVolume);
        using PrimaryIvVector = typename PrimaryInteractionVolume::Traits::Vector;
        using PrimaryIvMatrix = typename PrimaryInteractionVolume::Traits::Matrix;

        static_assert( std::is_convertible<PrimaryIvVector*, Vector*>::value,
                       "The vector type used in primary interaction volumes is not convertible to Dune::DynamicVector!" );
        static_assert( std::is_convertible<PrimaryIvMatrix*, Matrix*>::value,
                       "The matrix type used in primary interaction volumes is not convertible to Dune::DynamicMatrix!" );

        using SecondaryInteractionVolume = typename GET_PROP_TYPE(TypeTag, SecondaryInteractionVolume);
        using SecondaryIvVector = typename SecondaryInteractionVolume::Traits::Vector;
        using SecondaryIvMatrix = typename SecondaryInteractionVolume::Traits::Matrix;

        static_assert( std::is_convertible<SecondaryIvVector*, Vector*>::value,
                       "The vector type used in secondary interaction volumes is not convertible to Dune::DynamicVector!" );
        static_assert( std::is_convertible<SecondaryIvMatrix*, Matrix*>::value,
                       "The matrix type used in secondary interaction volumes is not convertible to Dune::DynamicMatrix!" );

    public:
        // export the filler type
        using Filler = MpfaDarcysLawCacheFiller;

        /*!
         * \brief Update cached objects (transmissibilities and gravity)
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
        template< class InteractionVolume, class LocalFaceData, class DataHandle >
        void updateAdvection(const InteractionVolume& iv,
                             const LocalFaceData& localFaceData,
                             const DataHandle& dataHandle,
                             const SubControlVolumeFace &scvf)
        {
            switchFluxSign_ = localFaceData.isOutside();

            for (unsigned int pIdx = 0; pIdx < numPhases; ++pIdx)
                pj_[pIdx] = &dataHandle.pressures(pIdx);

            static const bool enableGravity = getParamFromGroup<bool>(GET_PROP_VALUE(TypeTag, ModelParameterGroup), "Problem.EnableGravity");
            const auto ivLocalIdx = localFaceData.ivLocalScvfIndex();

            // standard grids
            if (dim == dimWorld)
            {
                Tij_ = &dataHandle.advectionT()[ivLocalIdx];

                if (enableGravity)
                    for (unsigned int phaseIdx = 0; phaseIdx < numPhases; ++phaseIdx)
                        g_[phaseIdx] = dataHandle.gravity(phaseIdx)[ivLocalIdx];
            }

            // surface grids
            else
            {
                if (!localFaceData.isOutside())
                {
                    Tij_ = &dataHandle.advectionT()[ivLocalIdx];

                    if (enableGravity)
                        for (unsigned int phaseIdx = 0; phaseIdx < numPhases; ++phaseIdx)
                            g_[phaseIdx] = dataHandle.gravity(phaseIdx)[ivLocalIdx];
                }
                else
                {
                    const auto idxInOutsideFaces = localFaceData.scvfLocalOutsideScvfIndex();
                    Tij_ = &dataHandle.advectionTout()[ivLocalIdx][idxInOutsideFaces];

                    if (enableGravity)
                        for (unsigned int phaseIdx = 0; phaseIdx < numPhases; ++phaseIdx)
                            g_[phaseIdx] = dataHandle.gravityOutside(phaseIdx)[ivLocalIdx][idxInOutsideFaces];
                }
            }
        }

        //! Coefficients for the cell (& Dirichlet) unknowns in flux expressions
        const Vector& advectionTij() const { return *Tij_; }

        //! The cell (& Dirichlet) pressures within this interaction volume
        const Vector& pressures(unsigned int phaseIdx) const { return *pj_[phaseIdx]; }

        //! The gravitational acceleration for a phase on this scvf
        Scalar gravity(unsigned int phaseIdx) const { return g_[phaseIdx]; }

        //! In the interaction volume-local system of eq we have one unknown per face.
        //! On scvfs on this face, but in "outside" (neighbor) elements of it, we have
        //! to take the negative value of the fluxes due to the flipped normal vector.
        //! This function returns whether or not this scvf is an "outside" face in the iv.
        bool advectionSwitchFluxSign() const { return switchFluxSign_; }

    private:
        bool switchFluxSign_;
        const Vector* Tij_;                       //!< The transmissibilities such that f = Tij*pj
        std::array< Scalar, numPhases > g_;       //!< Gravitational flux contribution on this face
        std::array<const Vector*, numPhases> pj_; //!< The interaction-volume wide phase pressures pj
    };

public:
    // state the discretization method this implementation belongs to
    static const DiscretizationMethods myDiscretizationMethod = DiscretizationMethods::CCMpfa;

    // export the type for the corresponding cache
    using Cache = MpfaDarcysLawCache;

    //! Compute the advective flux across an scvf
    static Scalar flux(const Problem& problem,
                       const Element& element,
                       const FVElementGeometry& fvGeometry,
                       const ElementVolumeVariables& elemVolVars,
                       const SubControlVolumeFace& scvf,
                       const unsigned int phaseIdx,
                       const ElementFluxVariablesCache& elemFluxVarsCache)
    {
        const auto& fluxVarsCache = elemFluxVarsCache[scvf];
        const auto& tij = fluxVarsCache.advectionTij();
        const auto& pj = fluxVarsCache.pressures(phaseIdx);

        // compute t_ij*p_j
        Scalar scvfFlux = tij*pj;

        // maybe add gravitational acceleration
        static const bool enableGravity = getParamFromGroup<bool>(GET_PROP_VALUE(TypeTag, ModelParameterGroup), "Problem.EnableGravity");
        if (enableGravity)
            scvfFlux += fluxVarsCache.gravity(phaseIdx);

        // switch the sign if necessary
        if (fluxVarsCache.advectionSwitchFluxSign())
            scvfFlux *= -1.0;

        return scvfFlux;
    }
};

} // end namespace

#endif
