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
 * \brief This file contains the data which is required to calculate
 *        volume and mass fluxes of fluid phases over a face of a finite volume by means
 *        of the Darcy approximation. Specializations are provided for the different discretization methods.
 */
#ifndef DUMUX_DISCRETIZATION_CC_TPFA_DARCYS_LAW_HH
#define DUMUX_DISCRETIZATION_CC_TPFA_DARCYS_LAW_HH

#include <memory>

#include <dune/common/float_cmp.hh>

#include <dumux/common/math.hh>
#include <dumux/common/parameters.hh>

#include <dumux/implicit/properties.hh>
#include <dumux/discretization/methods.hh>

namespace Dumux
{

namespace Properties
{
// forward declaration of properties
NEW_PROP_TAG(ProblemEnableGravity);
}

/*!
 * \ingroup DarcysLaw
 * \brief Specialization of Darcy's Law for the CCTpfa method.
 */
template <class TypeTag>
class DarcysLawImplementation<TypeTag, DiscretizationMethods::CCTpfa>
{
    using Implementation = DarcysLawImplementation<TypeTag, DiscretizationMethods::CCTpfa>;
    using Problem = typename GET_PROP_TYPE(TypeTag, Problem);
    using SubControlVolume = typename GET_PROP_TYPE(TypeTag, SubControlVolume);
    using SubControlVolumeFace = typename GET_PROP_TYPE(TypeTag, SubControlVolumeFace);
    using GridView = typename GET_PROP_TYPE(TypeTag, GridView);
    using Scalar = typename GET_PROP_TYPE(TypeTag, Scalar);
    using FVElementGeometry = typename GET_PROP_TYPE(TypeTag, FVElementGeometry);
    using ElementVolumeVariables = typename GET_PROP_TYPE(TypeTag, ElementVolumeVariables);
    using VolumeVariables = typename GET_PROP_TYPE(TypeTag, VolumeVariables);
    using ElementFluxVarsCache = typename GET_PROP_TYPE(TypeTag, ElementFluxVariablesCache);
    using FluxVariablesCache = typename GET_PROP_TYPE(TypeTag, FluxVariablesCache);

    using Element = typename GridView::template Codim<0>::Entity;
    using IndexType = typename GridView::IndexSet::IndexType;

    static const int dim = GridView::dimension;
    static const int dimWorld = GridView::dimensionworld;

    using DimWorldMatrix = Dune::FieldMatrix<Scalar, dimWorld, dimWorld>;
    using GlobalPosition = Dune::FieldVector<Scalar, dimWorld>;

    class TpfaDarcysLawCache
    {
    public:
        void updateAdvection(const Problem& problem,
                             const Element& element,
                             const FVElementGeometry& fvGeometry,
                             const ElementVolumeVariables& elemVolVars,
                             const SubControlVolumeFace &scvf)
        {
            tij_ = Implementation::calculateTransmissibilities(problem, element, fvGeometry, elemVolVars, scvf);
        }

        const Scalar& tij() const
        { return tij_; }

    private:
        Scalar tij_;
    };

    //! Class that fills the cache corresponding to tpfa Darcy's Law
    class TpfaDarcysLawCacheFiller
    {
    public:
        //! Function to fill a TpfaDarcysLawCache of a given scvf
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
            scvfFluxVarsCache.updateAdvection(problem,
                                              element,
                                              fvGeometry,
                                              elemVolVars,
                                              scvf);
        }
    };

public:
    // state the discretization method this implementation belongs to
    static const DiscretizationMethods myDiscretizationMethod = DiscretizationMethods::CCTpfa;

    // state the type for the corresponding cache and its filler
    using Cache = TpfaDarcysLawCache;
    using CacheFiller = TpfaDarcysLawCacheFiller;

    static Scalar flux(const Problem& problem,
                       const Element& element,
                       const FVElementGeometry& fvGeometry,
                       const ElementVolumeVariables& elemVolVars,
                       const SubControlVolumeFace& scvf,
                       int phaseIdx,
                       const ElementFluxVarsCache& elemFluxVarsCache)
    {
        const auto& fluxVarsCache = elemFluxVarsCache[scvf];

        // Get the inside and outside volume variables
        const auto& insideScv = fvGeometry.scv(scvf.insideScvIdx());
        const auto& insideVolVars = elemVolVars[insideScv];
        const auto& outsideVolVars = elemVolVars[scvf.outsideScvIdx()];

        if (GET_PARAM_FROM_GROUP(TypeTag, bool, Problem, EnableGravity))
        {
            // do averaging for the density over all neighboring elements
            const auto rho = [&]()
            {
                // boundaries
                if (scvf.boundary())
                    return insideVolVars.density(phaseIdx);

                // inner faces with two neighboring elements
                else if (scvf.numOutsideScvs() == 1)
                    return (insideVolVars.density(phaseIdx) + outsideVolVars.density(phaseIdx))*0.5;

                // inner faces in networks (general case)
                else
                {
                    Scalar rho(insideVolVars.density(phaseIdx));
                    for (unsigned int i = 0; i < scvf.numOutsideScvs(); ++i)
                    {
                        const auto outsideScvIdx = scvf.outsideScvIdx(i);
                        const auto& outsideVolVars = elemVolVars[outsideScvIdx];
                        rho += outsideVolVars.density(phaseIdx);
                    }
                    return rho/(scvf.numOutsideScvs()+1);
                }
            }();

            // ask for the gravitational acceleration in the inside neighbor
            const auto xInside = insideScv.center();
            const auto gInside = problem.gravityAtPos(xInside);
            const auto hInside = insideVolVars.pressure(phaseIdx) - rho*(gInside*xInside);
            const auto hOutside = [&]()
            {
                // boundaries
                if (scvf.boundary())
                {
                    const auto xOutside = scvf.ipGlobal();
                    const auto gOutside = problem.gravityAtPos(xOutside);
                    return outsideVolVars.pressure(phaseIdx) - rho*(gOutside*xOutside);
                }

                // inner faces with two neighboring elements
                else if (scvf.numOutsideScvs() == 1)
                {
                    const auto xOutside = fvGeometry.scv(scvf.outsideScvIdx()).center();
                    const auto gOutside = problem.gravityAtPos(xOutside);
                    return outsideVolVars.pressure(phaseIdx) - rho*(gOutside*xOutside);
                }

                // inner faces in networks (general case)
                else
                {
                    const auto& insideFluxVarsCache = elemFluxVarsCache[scvf];

                    Scalar sumTi(insideFluxVarsCache.tij());
                    Scalar sumPTi(insideFluxVarsCache.tij()*hInside);
                    for (unsigned int i = 0; i < scvf.numOutsideScvs(); ++i)
                    {
                        const auto outsideScvIdx = scvf.outsideScvIdx(i);
                        const auto& flippedScvf = fvGeometry.flipScvf(scvf.index(), i);
                        const auto& outsideVolVars = elemVolVars[outsideScvIdx];
                        const auto& outsideFluxVarsCache = elemFluxVarsCache[flippedScvf];
                        const auto xOutside = scvf.boundary() ? scvf.ipGlobal() : fvGeometry.scv(outsideScvIdx).center();
                        const auto gOutside = problem.gravityAtPos(xOutside);

                        sumTi += outsideFluxVarsCache.tij();
                        sumPTi += outsideFluxVarsCache.tij()*(outsideVolVars.pressure(phaseIdx) - rho*(gOutside*xOutside));
                    }
                    return sumPTi/sumTi;
                }
            }();

            return fluxVarsCache.tij()*(hInside - hOutside);
        }
        else // no gravity
        {
            const auto pInside = insideVolVars.pressure(phaseIdx);
            const auto pOutside = [&]()
            {
                // Dirichlet boundaries and inner faces with two neighboring elements
                if (scvf.numOutsideScvs() <= 1)
                    return outsideVolVars.pressure(phaseIdx);

                // inner faces in networks (general case)
                else
                {

                    const auto& insideFluxVarsCache = elemFluxVarsCache[scvf];
                    Scalar sumTi(insideFluxVarsCache.tij());
                    Scalar sumPTi(insideFluxVarsCache.tij()*pInside);

                    for (unsigned int i = 0; i < scvf.numOutsideScvs(); ++i)
                    {
                        const auto outsideScvIdx = scvf.outsideScvIdx(i);
                        const auto& flippedScvf = fvGeometry.flipScvf(scvf.index(), i);
                        const auto& outsideVolVars = elemVolVars[outsideScvIdx];
                        const auto& outsideFluxVarsCache = elemFluxVarsCache[flippedScvf];
                        sumTi += outsideFluxVarsCache.tij();
                        sumPTi += outsideFluxVarsCache.tij()*outsideVolVars.pressure(phaseIdx);
                    }
                    return sumPTi/sumTi;
                }
            }();

            return fluxVarsCache.tij()*(pInside - pOutside);
        }
    }

    // The flux variables cache has to be bound to an element prior to flux calculations
    // During the binding, the transmissibilities will be computed and stored using the method below.
    static Scalar calculateTransmissibilities(const Problem& problem,
                                              const Element& element,
                                              const FVElementGeometry& fvGeometry,
                                              const ElementVolumeVariables& elemVolVars,
                                              const SubControlVolumeFace& scvf)
    {
        Scalar tij;

        const auto insideScvIdx = scvf.insideScvIdx();
        const auto& insideScv = fvGeometry.scv(insideScvIdx);
        const auto& insideVolVars = elemVolVars[insideScvIdx];
        // check if we evaluate the permeability in the volume (for discontinuous fields, default)
        // or at the scvf center for analytical permeability fields (e.g. convergence studies)
        auto getPermeability = [&problem](const VolumeVariables& volVars,
                                          const GlobalPosition& scvfIpGlobal)
                               {
                                    if (GET_PROP_VALUE(TypeTag, EvaluatePermeabilityAtScvfIP))
                                        return problem.spatialParams().permeabilityAtPos(scvfIpGlobal);
                                    else
                                        return volVars.permeability();
                               };

        const Scalar ti = calculateOmega_(scvf, getPermeability(insideVolVars, scvf.ipGlobal()),
                                          insideScv, insideVolVars.extrusionFactor());

        // for the boundary (dirichlet) or at branching points we only need ti
        if (scvf.boundary() || scvf.numOutsideScvs() > 1)
        {
            tij = scvf.area()*ti;
        }

        // otherwise we compute a tpfa harmonic mean
        else
        {
            const auto outsideScvIdx = scvf.outsideScvIdx();
            // as we assemble fluxes from the neighbor to our element the outside index
            // refers to the scv of our element, so we use the scv method
            const auto& outsideScv = fvGeometry.scv(outsideScvIdx);
            const auto& outsideVolVars = elemVolVars[outsideScvIdx];

            const Scalar tj = [&]()
            {
                // normal grids
                if (dim == dimWorld)
                    return -1.0*calculateOmega_(scvf, getPermeability(outsideVolVars, scvf.ipGlobal()),
                                            outsideScv, outsideVolVars.extrusionFactor());

                // embedded surface and network grids
                //(the outside normal vector might differ from the inside normal vector)
                else
                    return calculateOmega_(fvGeometry.flipScvf(scvf.index()), getPermeability(outsideVolVars, scvf.ipGlobal()),
                                           outsideScv, outsideVolVars.extrusionFactor());

            }();

            // harmonic mean (check for division by zero!)
            if (ti*tj <= 0.0)
                tij = 0;
            else
                tij = scvf.area()*(ti * tj)/(ti + tj);
        }

        return tij;
    }

private:
    //! compute the transmissibility ti, overload for tensor permeabilites
    static Scalar calculateOmega_(const SubControlVolumeFace& scvf,
                                  const DimWorldMatrix &K,
                                  const SubControlVolume &scv,
                                  Scalar extrusionFactor)
    {
        GlobalPosition Knormal;
        K.mv(scvf.unitOuterNormal(), Knormal);

        auto distanceVector = scvf.ipGlobal();
        distanceVector -= scv.center();
        distanceVector /= distanceVector.two_norm2();

        Scalar omega = Knormal * distanceVector;
        omega *= extrusionFactor;

        return omega;
    }

    //! compute the transmissibility ti, overload for scalar permeabilites
    static Scalar calculateOmega_(const SubControlVolumeFace& scvf,
                                  const Scalar K,
                                  const SubControlVolume &scv,
                                  Scalar extrusionFactor)
    {
        auto distanceVector = scvf.ipGlobal();
        distanceVector -= scv.center();
        distanceVector /= distanceVector.two_norm2();

        Scalar omega = K * (distanceVector * scvf.unitOuterNormal());
        omega *= extrusionFactor;

        return omega;
    }
};

} // end namespace Dumux

#endif
