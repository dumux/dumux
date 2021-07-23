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
 * \ingroup CCTpfaFlux
 * \brief Forchheimers's law for cell-centered finite volume schemes with two-point flux approximation
 */
#ifndef DUMUX_DISCRETIZATION_CC_TPFA_FORCHHEIMERS_LAW_HH
#define DUMUX_DISCRETIZATION_CC_TPFA_FORCHHEIMERS_LAW_HH

#include <dune/common/fvector.hh>
#include <dune/common/fmatrix.hh>

#include <dumux/common/math.hh>
#include <dumux/common/parameters.hh>
#include <dumux/common/properties.hh>
#include <dumux/common/typetraits/typetraits.hh>

#include <dumux/discretization/method.hh>
#include <dumux/discretization/extrusion.hh>
#include <dumux/flux/cctpfa/darcyslaw.hh>

namespace Dumux {

// forward declarations
template<class TypeTag, class ForchheimerVelocity, DiscretizationMethod discMethod>
class ForchheimersLawImplementation;

/*!
 * \ingroup CCTpfaFlux
 * \brief Forchheimer's law for cell-centered finite volume schemes with two-point flux approximation
 * \note Forchheimer's law is specialized for network and surface grids (i.e. if grid dim < dimWorld)
 * \tparam Scalar the scalar type for scalar physical quantities
 * \tparam GridGeometry the grid geometry
 * \tparam ForchheimerVelocity class for the calculation of the Forchheimer velocity
 * \tparam isNetwork whether we are computing on a network grid embedded in a higher world dimension
 */
template<class Scalar, class GridGeometry, class ForchheimerVelocity, bool isNetwork>
class CCTpfaForchheimersLaw;

/*!
 * \ingroup CCTpfaFlux
 * \brief Forchheimer's law for cell-centered finite volume schemes with two-point flux approximation
 * \note Forchheimer's law is specialized for network and surface grids (i.e. if grid dim < dimWorld)
 */
template <class TypeTag, class ForchheimerVelocity>
class ForchheimersLawImplementation<TypeTag, ForchheimerVelocity, DiscretizationMethod::cctpfa>
: public CCTpfaForchheimersLaw<GetPropType<TypeTag, Properties::Scalar>,
                               GetPropType<TypeTag, Properties::GridGeometry>,
                               ForchheimerVelocity,
                               (GetPropType<TypeTag, Properties::GridGeometry>::GridView::dimension < GetPropType<TypeTag, Properties::GridGeometry>::GridView::dimensionworld)>
{};

/*!
 * \ingroup CCTpfaFlux
 * \brief Class that fills the cache corresponding to tpfa Forchheimer's Law
 */
template<class GridGeometry>
class TpfaForchheimersLawCacheFiller
{
    using FVElementGeometry = typename GridGeometry::LocalView;
    using SubControlVolumeFace = typename GridGeometry::SubControlVolumeFace;
    using Element = typename GridGeometry::GridView::template Codim<0>::Entity;

public:
    //! Function to fill a TpfaForchheimersLawCache of a given scvf
    //! This interface has to be met by any advection-related cache filler class
    //! TODO: Probably get cache type out of the filler
    template<class FluxVariablesCache, class Problem, class ElementVolumeVariables, class FluxVariablesCacheFiller>
    static void fill(FluxVariablesCache& scvfFluxVarsCache,
                     const Problem& problem,
                     const Element& element,
                     const FVElementGeometry& fvGeometry,
                     const ElementVolumeVariables& elemVolVars,
                     const SubControlVolumeFace& scvf,
                     const FluxVariablesCacheFiller& fluxVarsCacheFiller)
    {
        scvfFluxVarsCache.updateAdvection(problem, element, fvGeometry, elemVolVars, scvf);
    }
};

/*!
 * \ingroup CCTpfaFlux
 * \brief The cache corresponding to tpfa Forchheimer's Law
 */
template<class AdvectionType, class GridGeometry>
class TpfaForchheimersLawCache
{
    using Scalar = typename AdvectionType::Scalar;
    using FVElementGeometry = typename GridGeometry::LocalView;
    using SubControlVolumeFace = typename GridGeometry::SubControlVolumeFace;
    using Element = typename GridGeometry::GridView::template Codim<0>::Entity;
    static constexpr int dimWorld = GridGeometry::GridView::dimensionworld;
    using DimWorldMatrix = Dune::FieldMatrix<Scalar, dimWorld, dimWorld>;

public:
    using Filler = TpfaForchheimersLawCacheFiller<GridGeometry>;

    template<class Problem, class ElementVolumeVariables>
    void updateAdvection(const Problem& problem,
                         const Element& element,
                         const FVElementGeometry& fvGeometry,
                         const ElementVolumeVariables& elemVolVars,
                         const SubControlVolumeFace &scvf)
    {
        tij_ = AdvectionType::calculateTransmissibility(problem, element, fvGeometry, elemVolVars, scvf);
        harmonicMeanSqrtK_ = AdvectionType::calculateHarmonicMeanSqrtPermeability(problem, elemVolVars, scvf);
    }

    const Scalar& advectionTij() const
    { return tij_; }

    const DimWorldMatrix& harmonicMeanSqrtPermeability() const
    { return harmonicMeanSqrtK_; }

private:
    Scalar tij_;
    DimWorldMatrix harmonicMeanSqrtK_;
};

/*!
 * \ingroup CCTpfaFlux
 * \brief Specialization of the CCTpfaForchheimersLaw grids where dim=dimWorld
 */
template<class ScalarType, class GridGeometry, class ForchheimerVelocity>
class CCTpfaForchheimersLaw<ScalarType, GridGeometry, ForchheimerVelocity, /*isNetwork*/ false>
{
    using ThisType = CCTpfaForchheimersLaw<ScalarType, GridGeometry, ForchheimerVelocity, /*isNetwork*/ false>;
    using FVElementGeometry = typename GridGeometry::LocalView;
    using SubControlVolumeFace = typename GridGeometry::SubControlVolumeFace;
    using Extrusion = Extrusion_t<GridGeometry>;
    using GridView = typename GridGeometry::GridView;
    using Element = typename GridView::template Codim<0>::Entity;

    using DimWorldVector = typename ForchheimerVelocity::DimWorldVector;
    using DimWorldMatrix = typename ForchheimerVelocity::DimWorldMatrix;

    using DarcysLaw = CCTpfaDarcysLaw<ScalarType, GridGeometry, /*isNetwork*/ false>;

  public:
    //! state the scalar type of the law
    using Scalar = ScalarType;

    //! state the discretization method this implementation belongs to
    static const DiscretizationMethod discMethod = DiscretizationMethod::cctpfa;

    //! state the type for the corresponding cache
    using Cache = TpfaForchheimersLawCache<ThisType, GridGeometry>;

    /*! \brief Compute the advective flux of a phase across
    *          the given sub-control volume face using the Forchheimer equation.
    *
    *          The flux is given in N*m, and can be converted
    *          into a volume flux (m^3/s) or mass flux (kg/s) by applying an upwind scheme
    *          for the mobility or the product of density and mobility, respectively.
    */
    template<class Problem, class ElementVolumeVariables, class ElementFluxVarsCache>
    static Scalar flux(const Problem& problem,
                       const Element& element,
                       const FVElementGeometry& fvGeometry,
                       const ElementVolumeVariables& elemVolVars,
                       const SubControlVolumeFace& scvf,
                       int phaseIdx,
                       const ElementFluxVarsCache& elemFluxVarsCache)
    {
        // Get the volume flux based on Darcy's law. The value returned by this method needs to be multiplied with the
        // mobility (upwinding).
        Scalar darcyFlux = DarcysLaw::flux(problem, element, fvGeometry, elemVolVars, scvf, phaseIdx, elemFluxVarsCache);
        auto upwindTerm = [phaseIdx](const auto& volVars){ return volVars.mobility(phaseIdx); };
        DimWorldVector darcyVelocity = scvf.unitOuterNormal();
        darcyVelocity *= ForchheimerVelocity::UpwindScheme::apply(elemVolVars, scvf, upwindTerm, darcyFlux, phaseIdx);
        darcyVelocity /= Extrusion::area(scvf);

        const auto velocity = ForchheimerVelocity::velocity(fvGeometry,
                                                            elemVolVars,
                                                            scvf,
                                                            phaseIdx,
                                                            elemFluxVarsCache[scvf].harmonicMeanSqrtPermeability(),
                                                            darcyVelocity);

        Scalar flux = velocity * scvf.unitOuterNormal();
        flux *= Extrusion::area(scvf);
        return flux;
    }

    //! The flux variables cache has to be bound to an element prior to flux calculations
    //! During the binding, the transmissibility will be computed and stored using the method below.
    template<class Problem, class ElementVolumeVariables>
    static Scalar calculateTransmissibility(const Problem& problem,
                                            const Element& element,
                                            const FVElementGeometry& fvGeometry,
                                            const ElementVolumeVariables& elemVolVars,
                                            const SubControlVolumeFace& scvf)
    {
        return DarcysLaw::calculateTransmissibility(problem, element, fvGeometry, elemVolVars, scvf);
    }

    /*! \brief Returns the harmonic mean of \f$\sqrt{K_0}\f$ and \f$\sqrt{K_1}\f$.
     *
     * This is a specialization for scalar-valued permeabilities which returns a tensor with identical diagonal entries.
     */
    template<class Problem, class ElementVolumeVariables>
    static DimWorldMatrix calculateHarmonicMeanSqrtPermeability(const Problem& problem,
                                                                const ElementVolumeVariables& elemVolVars,
                                                                const SubControlVolumeFace& scvf)
    {
        return ForchheimerVelocity::calculateHarmonicMeanSqrtPermeability(problem, elemVolVars, scvf);
    }
};

/*!
 * \ingroup CCTpfaFlux
 * \brief Specialization of the CCTpfaForchheimersLaw grids where dim<dimWorld
 */
template<class ScalarType, class GridGeometry, class ForchheimerVelocity>
class CCTpfaForchheimersLaw<ScalarType, GridGeometry, ForchheimerVelocity, /*isNetwork*/ true>
{
    static_assert(AlwaysFalse<ScalarType>::value, "Forchheimer not implemented for network grids");
};

} // end namespace Dumux

#endif
