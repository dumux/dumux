// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup Assembly
 * \ingroup CVFEDiscretization
 * \brief Calculates the element-wise residual for control-volume finite element schemes
 */
#ifndef DUMUX_CVFE_LOCAL_RESIDUAL_HH
#define DUMUX_CVFE_LOCAL_RESIDUAL_HH

#include <dune/common/std/type_traits.hh>
#include <dune/geometry/type.hh>
#include <dune/istl/matrix.hh>

#include <dumux/common/properties.hh>
#include <dumux/common/numeqvector.hh>
#include <dumux/assembly/fvlocalresidual.hh>
#include <dumux/discretization/extrusion.hh>

namespace Dumux::Detail {

template<class Imp, class P, class G, class S, class V>
using TimeInfoInterfaceCVFEDetector = decltype(
    std::declval<Imp>().computeStorage(
        std::declval<P>(), std::declval<G>(), std::declval<S>(), std::declval<V>(), true
    )
);

template<class Imp, class P, class G, class S, class V>
constexpr inline bool hasTimeInfoInterfaceCVFE()
{ return Dune::Std::is_detected<TimeInfoInterfaceCVFEDetector, Imp, P, G, S, V>::value; }

template<class Imp>
using SCVFIsOverlappingDetector = decltype(
    std::declval<Imp>().isOverlapping()
);

template<class Imp>
constexpr inline bool hasScvfIsOverlapping()
{ return Dune::Std::is_detected<SCVFIsOverlappingDetector, Imp>::value; }

} // end namespace Dumux::Detail


namespace Dumux {

/*!
 * \ingroup Assembly
 * \ingroup CVFEDiscretization
 * \brief The element-wise residual for control-volume finite element schemes
 * \tparam TypeTag the TypeTag
 */
template<class TypeTag>
class CVFELocalResidual : public FVLocalResidual<TypeTag>
{
    using ParentType = FVLocalResidual<TypeTag>;
    using Implementation = GetPropType<TypeTag, Properties::LocalResidual>;
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using Problem = GetPropType<TypeTag, Properties::Problem>;
    using GridGeometry = GetPropType<TypeTag, Properties::GridGeometry>;
    using GridView = typename GridGeometry::GridView;
    using Element = typename GridView::template Codim<0>::Entity;
    using ElementBoundaryTypes = GetPropType<TypeTag, Properties::ElementBoundaryTypes>;
    using FVElementGeometry = typename GridGeometry::LocalView;
    using GridVolumeVariables = GetPropType<TypeTag, Properties::GridVolumeVariables>;
    using ElementVolumeVariables = typename GridVolumeVariables::LocalView;
    using VolumeVariables = typename GridVolumeVariables::VolumeVariables;
    using SubControlVolume = typename FVElementGeometry::SubControlVolume;
    using SubControlVolumeFace = typename FVElementGeometry::SubControlVolumeFace;
    using ElementFluxVariablesCache = typename GetPropType<TypeTag, Properties::GridFluxVariablesCache>::LocalView;
    using PrimaryVariables = GetPropType<TypeTag, Properties::PrimaryVariables>;
    using NumEqVector = Dumux::NumEqVector<PrimaryVariables>;
    using Extrusion = Extrusion_t<GridGeometry>;

public:
    using ElementResidualVector = typename ParentType::ElementResidualVector;
    using ParentType::ParentType;

    //! evaluate flux residuals for one sub control volume face and add to residual
    void evalFlux(ElementResidualVector& residual,
                  const Problem& problem,
                  const Element& element,
                  const FVElementGeometry& fvGeometry,
                  const ElementVolumeVariables& elemVolVars,
                  const ElementBoundaryTypes& elemBcTypes,
                  const ElementFluxVariablesCache& elemFluxVarsCache,
                  const SubControlVolumeFace& scvf) const
    {
        const auto flux = evalFlux(problem, element, fvGeometry, elemVolVars, elemBcTypes, elemFluxVarsCache, scvf);
        if (!scvf.boundary())
        {
            const auto& insideScv = fvGeometry.scv(scvf.insideScvIdx());
            const auto& outsideScv = fvGeometry.scv(scvf.outsideScvIdx());
            residual[insideScv.localDofIndex()] += flux;

            // for control-volume finite element schemes with overlapping control volumes
            if constexpr (Detail::hasScvfIsOverlapping<SubControlVolumeFace>())
            {
                if (!scvf.isOverlapping())
                    residual[outsideScv.localDofIndex()] -= flux;
            }
            else
                residual[outsideScv.localDofIndex()] -= flux;
        }
        else
        {
            const auto& insideScv = fvGeometry.scv(scvf.insideScvIdx());
            residual[insideScv.localDofIndex()] += flux;
        }
    }

    //! evaluate flux residuals for one sub control volume face
    NumEqVector evalFlux(const Problem& problem,
                         const Element& element,
                         const FVElementGeometry& fvGeometry,
                         const ElementVolumeVariables& elemVolVars,
                         const ElementBoundaryTypes& elemBcTypes,
                         const ElementFluxVariablesCache& elemFluxVarsCache,
                         const SubControlVolumeFace& scvf) const
    {
        NumEqVector flux(0.0);

        // inner faces
        if (!scvf.boundary())
            flux += this->asImp().computeFlux(problem, element, fvGeometry, elemVolVars, scvf, elemFluxVarsCache);

        // boundary faces
        else
        {
            const auto& scv = fvGeometry.scv(scvf.insideScvIdx());
            const auto& bcTypes = elemBcTypes.get(fvGeometry, scv);

            // Treat Neumann and Robin ("solution dependent Neumann") boundary conditions.
            // For Dirichlet there is no addition to the residual here but they
            // are enforced strongly by replacing the residual entry afterwards.
            if (bcTypes.hasNeumann())
            {
                auto neumannFluxes = problem.neumann(element, fvGeometry, elemVolVars, elemFluxVarsCache, scvf);

                // multiply neumann fluxes with the area and the extrusion factor
                neumannFluxes *= Extrusion::area(fvGeometry, scvf)*elemVolVars[scv].extrusionFactor();

                // only add fluxes to equations for which Neumann is set
                for (int eqIdx = 0; eqIdx < NumEqVector::dimension; ++eqIdx)
                    if (bcTypes.isNeumann(eqIdx))
                        flux[eqIdx] += neumannFluxes[eqIdx];
            }
        }

        return flux;
    }

    using ParentType::evalStorage;
    /*!
     * \brief Compute the storage local residual, i.e. the deviation of the
     *        storage term from zero for instationary problems.
     *
     * \param residual The residual vector to fill
     * \param problem The problem to solve
     * \param element The DUNE Codim<0> entity for which the residual
     *                ought to be calculated
     * \param fvGeometry The finite-volume geometry of the element
     * \param prevElemVolVars The volume averaged variables for all
     *                        sub-control volumes of the element at the previous time level
     * \param curElemVolVars The volume averaged variables for all
     *                       sub-control volumes of the element at the current  time level
     * \param scv The sub control volume the storage term is integrated over
     */
    void evalStorage(ElementResidualVector& residual,
                     const Problem& problem,
                     const Element& element,
                     const FVElementGeometry& fvGeometry,
                     const ElementVolumeVariables& prevElemVolVars,
                     const ElementVolumeVariables& curElemVolVars,
                     const SubControlVolume& scv) const
    {
        const auto& curVolVars = curElemVolVars[scv];
        const auto& prevVolVars = prevElemVolVars[scv];

        // Compute storage with the model specific storage residual
        // This addresses issues #792/#940 in ad-hoc way by additionally providing crude time level information (previous or current)
        // to the low-level interfaces if this is supported by the LocalResidual implementation
        NumEqVector prevStorage = computeStorageImpl_(problem, fvGeometry, scv, prevVolVars, /*previous time level?*/true);
        NumEqVector storage = computeStorageImpl_(problem, fvGeometry, scv, curVolVars, /*previous time level?*/false);

        prevStorage *= prevVolVars.extrusionFactor();
        storage *= curVolVars.extrusionFactor();

        storage -= prevStorage;
        storage *= Extrusion::volume(fvGeometry, scv);
        storage /= this->timeLoop().timeStepSize();

        residual[scv.localDofIndex()] += storage;
    }

private:
    NumEqVector computeStorageImpl_(const Problem& problem,
                                    const FVElementGeometry& fvGeometry,
                                    const SubControlVolume& scv,
                                    const VolumeVariables& volVars,
                                    [[maybe_unused]] bool isPreviousTimeStep) const
    {
        if constexpr (Detail::hasTimeInfoInterfaceCVFE<Implementation, Problem, FVElementGeometry, SubControlVolume, VolumeVariables>())
            return this->asImp().computeStorage(problem, fvGeometry, scv, volVars, isPreviousTimeStep);
        else
            return this->asImp().computeStorage(problem, scv, volVars);
    }
};

} // end namespace Dumux

#endif
