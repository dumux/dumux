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
 * \brief This file contains a class for the calculation of fluxes at the boundary
          of pore-network models.
 *
 */
#ifndef DUMUX_PNM_BOUNDARYFLUX_HH
#define DUMUX_PNM_BOUNDARYFLUX_HH

#include <algorithm>
#include <vector>
#include <dumux/discretization/box/elementboundarytypes.hh>

namespace Dumux
{
template<class Assembler>
class PoreNetworkModelBoundaryFlux
{
    using Scalar = typename Assembler::Scalar;
    using GridGeometry = typename Assembler::GridGeometry;
    using GridView = typename GridGeometry::GridView;
    using Element = typename GridView::template Codim<0>::Entity;
    using FVElementGeometry = typename GridGeometry::LocalView;
    using SubControlVolume = typename FVElementGeometry::SubControlVolume;
    using Problem = typename Assembler::Problem;
    using BoundaryTypes = std::decay_t<decltype(std::declval<Problem>().boundaryTypes(std::declval<Element>(), std::declval<SubControlVolume>()))>;
    using ElementBoundaryTypes = BoxElementBoundaryTypes<BoundaryTypes>;

    using SolutionVector = typename Assembler::ResidualType;;
    using NumEqVector = typename SolutionVector::block_type;
    using GridVariables = typename Assembler::GridVariables;

    static constexpr auto dim = GridView::dimension;

public:

    PoreNetworkModelBoundaryFlux(const Assembler& assembler,
                                 const SolutionVector& sol)
    : assembler_(assembler)
    , problem_(assembler.problem())
    , gridVariables_(assembler.gridVariables())
    , sol_(sol)
    , isStationary_(assembler.isStationaryProblem()) {}

    /*!
     * \brief Returns the cumulative flux in \f$\mathrm{[\frac{kg}{s}]}\f$ of several pore throats for a given list of pore labels to consider
     *
     * \param labels A list of pore labels which will be considered for the flux calculation
     * \param verbose If set true, the fluxes at all individual SCVs are printed
     */
    template<class Label>
    NumEqVector getFlux(const std::vector<Label>& labels, const bool verbose = false) const
    {
        // helper lambda to decide which scvs to consider for flux calculation
        auto restriction = [&] (const SubControlVolume& scv)
        {
            const Label poreLabel = problem_.gridGeometry().poreLabel(scv.dofIndex());
            return std::any_of(labels.begin(), labels.end(),
                               [&](const Label l){ return l == poreLabel; });
        };

        NumEqVector flux(0.0);

        // sum up the fluxes
        for (const auto& element : elements(problem_.gridGeometry().gridView()))
            flux += getFlux(element, restriction, verbose);

        return flux;
    }

    /*!
     * \brief Returns the cumulative flux in \f$\mathrm{[\frac{kg}{s}]}\f$ of several pore throats at a given location on the boundary
     *
     * \param minMax Consider bBoxMin or bBoxMax by setting "min" or "max"
     * \param coord x, y or z coordinate at which bBoxMin or bBoxMax is evaluated
     * \param verbose If set true, the fluxes at all individual SCVs are printed
     */
    NumEqVector getFlux(const std::string minMax, const int coord, const bool verbose = false) const
    {
        if(!(minMax == "min" || minMax == "max"))
            DUNE_THROW(Dune::InvalidStateException,
                    "second argument must be either 'min' or 'max' (string) !");

        const Scalar eps = 1e-6; //TODO
        auto onMeasuringBoundary = [&] (const Scalar pos)
        {
            return ( (minMax == "min" && pos < problem_.gridGeometry().bBoxMin()[coord] + eps) ||
                     (minMax == "max" && pos > problem_.gridGeometry().bBoxMax()[coord] - eps) );
        };

        // helper lambda to decide which scvs to consider for flux calculation
        auto restriction = [&] (const SubControlVolume& scv)
        {
            bool considerAllDirections = coord < 0 ? true : false;

            //only consider SCVs on the boundary
            bool considerScv = problem_.gridGeometry().dofOnBoundary(scv.dofIndex()) && onMeasuringBoundary(scv.dofPosition()[coord]);

            //check whether a vertex lies on a boundary and also check whether this boundary shall be
            // considered for the flux calculation
            if(considerScv && !considerAllDirections)
            {
                const auto& pos = scv.dofPosition();
                if (!(pos[coord] < problem_.gridGeometry().bBoxMin()[coord] + eps || pos[coord] > problem_.gridGeometry().bBoxMax()[coord] -eps ))
                considerScv = false;
            }

            return considerScv;
        };

        NumEqVector flux(0.0);

        // sum up the fluxes
        for(const auto& element : elements(problem_.gridGeometry().gridView()))
            flux += getFlux(element, restriction, verbose);

        return flux;
    }

    /*!
     * \brief Returns the cumulative flux in \f$\mathrm{[\frac{kg}{s}]}\f$ of several pore throats at a given location on the boundary
     *
     * \param element The element
     * \param considerScv A lambda function to decide whether to consider a scv or not
     * \param verbose If set true, the fluxes at all individual SCVs are printed
     */
    template<class RestrictingFunction>
    NumEqVector getFlux(const Element& element,
                        RestrictingFunction considerScv,
                        const bool verbose = false) const
    {
        NumEqVector flux(0.0);

        // by default, all coordinate directions are considered for the definition of a boundary

        // make sure FVElementGeometry and volume variables are bound to the element
        auto fvGeometry = localView(problem_.gridGeometry());
        fvGeometry.bind(element);

        auto curElemVolVars = localView(gridVariables_.curGridVolVars());
        curElemVolVars.bind(element, fvGeometry, sol_);

        auto prevElemVolVars = localView(gridVariables_.prevGridVolVars());
        if (!isStationary_)
            prevElemVolVars.bindElement(element, fvGeometry, sol_);

        auto elemFluxVarsCache = localView(gridVariables_.gridFluxVarsCache());
        elemFluxVarsCache.bindElement(element, fvGeometry, curElemVolVars);

        ElementBoundaryTypes elemBcTypes;
        elemBcTypes.update(problem_, element, fvGeometry);

        const auto localResidual = assembler_.localResidual();

        auto residual = localResidual.evalFluxAndSource(element, fvGeometry, curElemVolVars, elemFluxVarsCache, elemBcTypes);

        if (!isStationary_)
            residual += localResidual.evalStorage(element, fvGeometry, prevElemVolVars, curElemVolVars);

        for(auto&& scv : scvs(fvGeometry))
        {
            // compute the boundary flux using the local residual of the element's scv on the boundary
            if(considerScv(scv))
            {
                // The flux must be substracted:
                // On an inlet boundary, the flux part of the local residual will be positive, since all fluxes will leave the SCV towards to interior domain.
                // For the domain itself, however, the sign has to be negative, since mass is entering the system.
                flux -= residual[scv.indexInElement()];

                if(verbose)
                {
                    std::cout << "SCV of element " << scv.elementIndex()  << " at vertex " << scv.dofIndex() << " has flux: " << residual[scv.indexInElement()] << std::endl;
                }
            }
        }
        return flux;
    }

private:
    const Assembler& assembler_;
    const Problem& problem_;
    const GridVariables& gridVariables_;
    const SolutionVector& sol_;

    bool isStationary_;
};


} // end namespace

#endif
