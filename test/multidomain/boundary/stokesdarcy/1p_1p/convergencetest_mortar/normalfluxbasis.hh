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
 * \ingroup Discretization
 * \brief Provides helper aliases and functionality to obtain the types
 *        and instances of Dune::Functions function space bases for the
 *        velocity spaces of different discretization schemes.
 */
#ifndef DUMUX_DISCRETIZATION_VELOCITY_BASIS_HH
#define DUMUX_DISCRETIZATION_VELOCITY_BASIS_HH

#if HAVE_DUNE_FUNCTIONS

#include <dune/istl/bvector.hh>
#include <dune/common/fvector.hh>
#include <dune/common/exceptions.hh>
#include <dune/functions/functionspacebases/lagrangedgbasis.hh>
#include <dune/functions/functionspacebases/raviartthomasbasis.hh>
#include <dune/functions/gridfunctions/discreteglobalbasisfunction.hh>

#include <dumux/discretization/method.hh>

namespace Dumux {

/*!
 * \ingroup Discretization
 * \brief Traits class that contains information on the function
 *        space basis used by a discretization scheme (underlying
 *        a grid geometry class) for the velocities. Specializations
 *        of the traits for different schemes are provided below.
 */
template< class GridGeometry,
          DiscretizationMethod dm = GridGeometry::discMethod >
struct VelocityFunctionSpaceBasisTraits;

/*!
 * \ingroup Discretization
 * \brief Creates a Dune::Functions object of the function space
 *        basis for the velocity.
 */
template<class GridGeometry>
typename VelocityFunctionSpaceBasisTraits<GridGeometry>::GlobalBasis
getVelocityFunctionSpaceBasis(const GridGeometry& gridGeometry)
{ return {gridGeometry.gridView()}; }

/*!
 * \ingroup Discretization
 * \brief Creates the coefficient vector. TODO doc properly. Spec for tpfa
 */
template<class FluxVariables, class GridGeometry, class GridVariables, class SolutionVector,
         std::enable_if_t<GridGeometry::discMethod == DiscretizationMethod::cctpfa, int> = 0>
std::array< Dune::BlockVector< Dune::FieldVector<typename SolutionVector::field_type, 1> >,
            FluxVariables::numPhases >
getVelocityCoefficientVector(const typename VelocityFunctionSpaceBasisTraits<GridGeometry>::GlobalBasis& velocityBasis,
                             const GridGeometry& gridGeometry,
                             const GridVariables& gridVariables,
                             const SolutionVector& x)
{
    static constexpr int numPhases = FluxVariables::numPhases;

    using FluxRange = Dune::FieldVector<typename SolutionVector::field_type, 1>;
    using CoefficientVector = Dune::BlockVector<FluxRange>;
    using ResultType = std::array< CoefficientVector, numPhases >;

    ResultType coefficients;
    std::for_each(coefficients.begin(),
                  coefficients.end(),
                  [&velocityBasis] (auto& c) { c.resize(velocityBasis.size()); });

    std::vector<bool> visited(velocityBasis.size(), false);
    for (const auto& element : elements(gridGeometry.gridView()))
    {
        auto basisLocalView = velocityBasis.localView();
        auto fvGeometry = localView(gridGeometry);
        auto elemVolVars = localView(gridVariables.curGridVolVars());
        auto elemFluxVarsCache = localView(gridVariables.gridFluxVarsCache());

        basisLocalView.bind(element);
        fvGeometry.bind(element);
        elemVolVars.bind(element, fvGeometry, x);
        elemFluxVarsCache.bind(element, fvGeometry, elemVolVars);

        unsigned int localDofCounter = 0;
        assert(basisLocalView.tree().finiteElement().localBasis().size() == fvGeometry.numScvf());
        for (const auto& scvf : scvfs(fvGeometry))
        {
            const auto idx = basisLocalView.index(localDofCounter++);
            if (!visited[idx])
            {
                visited[idx] = true;
                const auto& insideVolVars = elemVolVars[ fvGeometry.scv(scvf.insideScvIdx()) ];
                const auto& problem = elemVolVars.gridVolVars().problem();

                FluxVariables fluxVars;
                fluxVars.init(problem, element, fvGeometry, elemVolVars, scvf, elemFluxVarsCache);

                // lambda to compute flux across a non-neumann face
                const auto scvfArea = insideVolVars.extrusionFactor()*scvf.area();
                auto computeFlux = [&] (int pIdx)
                {
                    auto up = [pIdx] (const auto& vv) { return vv.mobility(pIdx); };
                    return fluxVars.advectiveFlux(pIdx, up)/scvfArea;
                };

                if (scvf.boundary())
                {
                    const auto bcTypes = problem.boundaryTypes(element, scvf);
                    if (bcTypes.hasOnlyNeumann())
                    {
                        const auto neumannFluxes = problem.neumann(element, fvGeometry, elemVolVars, scvf);

                        // TODO check for zeroes and throw if not zero
                        if (numPhases > 1)
                            DUNE_THROW(Dune::InvalidStateException, "Don't know how to interpret Neumann fluxes");

                        // std::cout << "NEUMFLUX VS EXACT " << neumannFluxes[0] << " - " << problem.exactFlux(scvf.ipGlobal()) << std::endl;
                        for (unsigned int pIdx = 0; pIdx < numPhases; ++pIdx)
                            coefficients[pIdx][idx] = neumannFluxes[0]/insideVolVars.density(pIdx);
                    }
                    else if (bcTypes.hasOnlyDirichlet())
                    {
                        // std::cout << "diri flux VS EXACT " << computeFlux(0) << " - " << problem.exactFlux(scvf.ipGlobal()) << std::endl;

                        for (unsigned int pIdx = 0; pIdx < numPhases; ++pIdx)
                            coefficients[pIdx][idx] = computeFlux(pIdx);
                    }
                    else
                        DUNE_THROW(Dune::NotImplemented, "Mixed boundary types for cc schemes");
                }
                else
                {
                    // std::cout << "flux VS EXACT " << computeFlux(0) << " - " << problem.exactFlux(scvf.ipGlobal()) << std::endl;

                    const int sign = scvf.insideScvIdx() > scvf.outsideScvIdx() ? 1.0 : -1.0;
                    for (unsigned int pIdx = 0; pIdx < numPhases; ++pIdx)
                        coefficients[pIdx][idx] = sign*computeFlux(pIdx);
                }
            }
        }
    }

    return coefficients;
}

/*!
 * \ingroup Discretization
 * \brief Creates the coefficient vector. TODO doc properly. Spec for mpfa
 */
template<class FluxVariables, class GridGeometry, class GridVariables, class SolutionVector,
         std::enable_if_t<GridGeometry::discMethod != DiscretizationMethod::cctpfa, int> = 0>
std::array< Dune::BlockVector< Dune::FieldVector<typename SolutionVector::field_type, GridGeometry::GridView::dimension> >,
            FluxVariables::numPhases >
getVelocityCoefficientVector(const typename VelocityFunctionSpaceBasisTraits<GridGeometry>::GlobalBasis& velocityBasis,
                             const GridGeometry& gridGeometry,
                             const GridVariables& gridVariables,
                             const SolutionVector& x)
{ DUNE_THROW(Dune::NotImplemented, "Velocity coefficient vector computation for provided scheme"); }


///////////////////////////////////////////////////////////////////
// Specializations of the VelocityFunctionSpaceBasisTraits class //
///////////////////////////////////////////////////////////////////

//! Traits specialization: cc schemes use Raviart-Thomas basis of order 0
template< class GridGeometry >
struct VelocityFunctionSpaceBasisTraits<GridGeometry, DiscretizationMethod::cctpfa>
{ using GlobalBasis = Dune::Functions::RaviartThomasBasis<typename GridGeometry::GridView, /*order*/0>; };

//! Traits specialization: cc schemes use Raviart-Thomas basis of order 0
template< class GridGeometry >
struct VelocityFunctionSpaceBasisTraits<GridGeometry, DiscretizationMethod::ccmpfa>
{ using GlobalBasis = Dune::Functions::RaviartThomasBasis<typename GridGeometry::GridView, /*order*/0>; };

//! Traits specialization: staggered uses Raviart-Thomas basis of order 0
template< class GridGeometry >
struct VelocityFunctionSpaceBasisTraits<GridGeometry, DiscretizationMethod::staggered>
{ using GlobalBasis = Dune::Functions::RaviartThomasBasis<typename GridGeometry::GridView, /*order*/0>; };

//! Traits specialization: box scheme uses a discontinous first order lagrange basis
template< class GridGeometry >
struct VelocityFunctionSpaceBasisTraits<GridGeometry, DiscretizationMethod::box>
{ using GlobalBasis = Dune::Functions::LagrangeDGBasis<typename GridGeometry::GridView, /*order*/1>; };

//! Traits specialization: TODO what to put here? Or do not specialize at all?
// template< class GridGeometry >
// struct VelocitySpaceBasisTraits<GridGeometry, DiscretizationMethod::fem>
// { using GlobalBasis = typename GridGeometry::FEBasis; };

} // end namespace Dumux

#endif // HAVE_DUNE_FUNCTIONS
#endif // DUMUX_DISCRETIZATION_VELOCITY_SPACE_BASIS_HH
