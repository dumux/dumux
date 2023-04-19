// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup PoreNetworkModels
 * \copydoc Dumux::PoreNetwork::BoundaryFlux
 */
#ifndef DUMUX_PNM_BOUNDARYFLUX_HH
#define DUMUX_PNM_BOUNDARYFLUX_HH

#include <algorithm>
#include <numeric>
#include <vector>
#include <type_traits>
#include <unordered_map>
#include <string_view>
#include <iostream>

#include <dune/common/exceptions.hh>
#include <dumux/common/typetraits/problem.hh>
#include <dumux/discretization/cvfe/elementboundarytypes.hh>

namespace Dumux::PoreNetwork {

/*!
 * \ingroup PoreNetworkModels
 * \brief Class for the calculation of fluxes at the boundary
 *        of pore-network models.
 */
template<class GridVariables, class LocalResidual, class SolutionVector>
class BoundaryFlux
{
    using Problem = std::decay_t<decltype(std::declval<LocalResidual>().problem())>;
    using GridGeometry = typename ProblemTraits<Problem>::GridGeometry;
    using GridView = typename GridGeometry::GridView;
    using Element = typename GridView::template Codim<0>::Entity;
    using FVElementGeometry = typename GridGeometry::LocalView;
    using SubControlVolume = typename FVElementGeometry::SubControlVolume;
    using BoundaryTypes = typename ProblemTraits<Problem>::BoundaryTypes;
    using ElementBoundaryTypes = CVFEElementBoundaryTypes<BoundaryTypes>;

    using NumEqVector = typename SolutionVector::block_type;
    static constexpr auto numEq = NumEqVector::dimension;

    //! result struct that holds both the total flux and the flux per pore
    struct Result
    {
        NumEqVector totalFlux;
        std::unordered_map<std::size_t, NumEqVector> fluxPerPore;

        //! make the total flux printable to e.g. std::cout
        friend std::ostream& operator<< (std::ostream& stream, const Result& result)
        {
            stream << result.totalFlux;
            return stream;
        }

        //! allow to get the total flux for a given eqIdx
        const auto& operator[] (int eqIdx) const
        { return totalFlux[eqIdx]; }

        //! make the total flux assignable to NumEqVector through implicit conversion
        operator NumEqVector() const
        { return totalFlux; }
    };

public:
    // export the Scalar type
    using Scalar = typename GridVariables::Scalar;

    BoundaryFlux(const GridVariables& gridVariables,
                 const LocalResidual& localResidual,
                 const SolutionVector& sol)
    : localResidual_(localResidual)
    , gridVariables_(gridVariables)
    , sol_(sol)
    , isStationary_(localResidual.isStationary())
    {
        const auto numDofs = localResidual_.problem().gridGeometry().numDofs();
        isConsidered_.resize(numDofs, false);
        boundaryFluxes_.resize(numDofs);
    }

    /*!
     * \brief Returns the cumulative flux in \f$\mathrm{[\frac{kg}{s}]}\f$, \f$\mathrm{[\frac{mole}{s}]}\f$ or \f$\mathrm{[\frac{J}{s}]}\f$ of several pores for a given list of pore labels to consider
     *
     * \param labels A list of pore labels which will be considered for the flux calculation
     * \param verbose If set true, the fluxes at all individual SCVs are printed
     */
    template<class Label>
    Result getFlux(const std::vector<Label>& labels, const bool verbose = false) const
    {
        // helper lambda to decide which scvs to consider for flux calculation
        auto restriction = [&] (const SubControlVolume& scv)
        {
            const Label poreLabel = localResidual_.problem().gridGeometry().poreLabel(scv.dofIndex());
            return std::any_of(labels.begin(), labels.end(),
                               [&](const Label l){ return l == poreLabel; });
        };

        std::fill(boundaryFluxes_.begin(), boundaryFluxes_.end(), NumEqVector(0.0));
        std::fill(isConsidered_.begin(), isConsidered_.end(), false);

        // sum up the fluxes
        for (const auto& element : elements(localResidual_.problem().gridGeometry().gridView()))
            computeBoundaryFlux_(element, restriction, verbose);

        Result result;
        result.totalFlux = std::accumulate(boundaryFluxes_.begin(), boundaryFluxes_.end(), NumEqVector(0.0));;
        for (int i = 0; i < isConsidered_.size(); ++i)
        {
            if (isConsidered_[i])
                result.fluxPerPore[i] = boundaryFluxes_[i];
        }

        return result;
    }

    /*!
     * \brief Returns the cumulative flux in \f$\mathrm{[\frac{kg}{s}]}\f$, \f$\mathrm{[\frac{mole}{s}]}\f$ or \f$\mathrm{[\frac{J}{s}]}\f$ of several pores at a given location on the boundary
     *
     * \param minMax Consider bBoxMin or bBoxMax by setting "min" or "max"
     * \param coord x, y or z coordinate at which bBoxMin or bBoxMax is evaluated
     * \param verbose If set true, the fluxes at all individual SCVs are printed
     */
    Result getFlux(std::string_view minMax, const int coord, const bool verbose = false) const
    {
        if (!(minMax == "min" || minMax == "max"))
            DUNE_THROW(Dune::InvalidStateException,
                    "second argument must be either 'min' or 'max' (string) !");

        const Scalar eps = 1e-6; //TODO
        auto onMeasuringBoundary = [&] (const Scalar pos)
        {
            return ( (minMax == "min" && pos < localResidual_.problem().gridGeometry().bBoxMin()[coord] + eps) ||
                     (minMax == "max" && pos > localResidual_.problem().gridGeometry().bBoxMax()[coord] - eps) );
        };

        // helper lambda to decide which scvs to consider for flux calculation
        auto restriction = [&] (const SubControlVolume& scv)
        {
            bool considerAllDirections = coord < 0 ? true : false;

            //only consider SCVs on the boundary
            bool considerScv = localResidual_.problem().gridGeometry().dofOnBoundary(scv.dofIndex()) && onMeasuringBoundary(scv.dofPosition()[coord]);

            //check whether a vertex lies on a boundary and also check whether this boundary shall be
            // considered for the flux calculation
            if (considerScv && !considerAllDirections)
            {
                const auto& pos = scv.dofPosition();
                if (!(pos[coord] < localResidual_.problem().gridGeometry().bBoxMin()[coord] + eps || pos[coord] > localResidual_.problem().gridGeometry().bBoxMax()[coord] -eps ))
                considerScv = false;
            }

            return considerScv;
        };

        std::fill(boundaryFluxes_.begin(), boundaryFluxes_.end(), NumEqVector(0.0));
        std::fill(isConsidered_.begin(), isConsidered_.end(), false);

        // sum up the fluxes
        for (const auto& element : elements(localResidual_.problem().gridGeometry().gridView()))
            computeBoundaryFlux_(element, restriction, verbose);

        Result result;
        result.totalFlux = std::accumulate(boundaryFluxes_.begin(), boundaryFluxes_.end(), NumEqVector(0.0));;
        for (int i = 0; i < isConsidered_.size(); ++i)
        {
            if (isConsidered_[i])
                result.fluxPerPore[i] = boundaryFluxes_[i];
        }

        return result;
    }

private:

    /*!
     * \brief Computes the cumulative flux in \f$\mathrm{[\frac{kg}{s}]}\f$, \f$\mathrm{[\frac{mole}{s}]}\f$ or \f$\mathrm{[\frac{J}{s}]}\f$ of an individual pore on the boundary
     *
     * \param element The element
     * \param considerScv A lambda function to decide whether to consider a scv or not
     * \param verbose If set true, the fluxes at all individual SCVs are printed
     */
    template<class RestrictingFunction>
    void computeBoundaryFlux_(const Element& element,
                              RestrictingFunction considerScv,
                              const bool verbose = false) const
    {
        // by default, all coordinate directions are considered for the definition of a boundary

        // make sure FVElementGeometry and volume variables are bound to the element
        const auto fvGeometry = localView(localResidual_.problem().gridGeometry()).bind(element);
        const auto curElemVolVars = localView(gridVariables_.curGridVolVars()).bind(element, fvGeometry, sol_);

        auto prevElemVolVars = localView(gridVariables_.prevGridVolVars());
        if (!isStationary_)
            prevElemVolVars.bindElement(element, fvGeometry, sol_);

        ElementBoundaryTypes elemBcTypes;
        elemBcTypes.update(localResidual_.problem(), element, fvGeometry);
        const auto elemFluxVarsCache = localView(gridVariables_.gridFluxVarsCache()).bindElement(element, fvGeometry, curElemVolVars);
        auto residual = localResidual_.evalFluxAndSource(element, fvGeometry, curElemVolVars, elemFluxVarsCache, elemBcTypes);

        if (!isStationary_)
            residual += localResidual_.evalStorage(element, fvGeometry, prevElemVolVars, curElemVolVars);

        for (auto&& scv : scvs(fvGeometry))
        {
            // compute the boundary flux using the local residual of the element's scv on the boundary
            if (considerScv(scv))
            {
                isConsidered_[scv.dofIndex()] = true;

                // get the type of the boundary condition on the scv
                const auto& bcTypes = elemBcTypes.get(fvGeometry, scv);

                NumEqVector flux(0.0);
                for (int eqIdx = 0; eqIdx < NumEqVector::dimension; ++eqIdx)
                {
                    // Check the type of the boundary condition.
                    // If BC is Dirichlet type, the flux is equal to the local residual of the element's scv on the boundary.
                    // Otherwise the flux is either zero or equal to a source term at the element's scv on the boundary.
                    // Technicaly, the PNM considers source terms instead of Neumann BCs.
                    if (!bcTypes.isDirichlet(eqIdx))
                    {
                        auto source = localResidual_.computeSource(localResidual_.problem(), element, fvGeometry, curElemVolVars, scv);
                        source *= scv.volume() * curElemVolVars[scv].extrusionFactor();
                        flux[eqIdx] = source[eqIdx];
                    }
                    else
                        flux[eqIdx] = residual[scv.indexInElement()][eqIdx];
                }

                // The flux must be subtracted:
                // On an inlet boundary, the flux part of the local residual will be positive, since all fluxes will leave the SCV towards to interior domain.
                // For the domain itself, however, the sign has to be negative, since mass is entering the system.
                boundaryFluxes_[scv.dofIndex()] -= flux;

                if (verbose)
                    std::cout << "SCV of element " << scv.elementIndex()  << " at vertex " << scv.dofIndex() << " has flux: " << flux << std::endl;
            }
        }
    }

    const LocalResidual localResidual_; // store a copy of the local residual
    const GridVariables& gridVariables_;
    const SolutionVector& sol_;
    bool isStationary_;
    mutable std::vector<bool> isConsidered_;
    mutable std::vector<NumEqVector> boundaryFluxes_;
};

} // end Dumux::PoreNetwork

#endif
