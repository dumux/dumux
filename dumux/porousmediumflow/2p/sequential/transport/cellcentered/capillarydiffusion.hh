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
 * \ingroup SequentialTwoPModel
 * \brief  Class for defining the diffusive capillary pressure term of a 2p saturation equation.
 */
#ifndef DUMUX_CAPILLARYDIFFUSION_HH
#define DUMUX_CAPILLARYDIFFUSION_HH

#include <dumux/porousmediumflow/2p/sequential/transport/cellcentered/diffusivepart.hh>
#include "properties.hh"

#include <dumux/common/deprecated.hh>

namespace Dumux {
/*!
 * \ingroup SequentialTwoPModel
 * \brief  Class for defining the diffusive capillary pressure term of a saturation equation.
 *
 * Defines the diffusive capillary pressure term of the form
 *
 * \f[
 * \bar \lambda \boldsymbol K \textbf{grad} \, p_c,
 * \f]
 *
 * where \f$ \bar \lambda = \lambda_w f_n = \lambda_n f_w \f$ and \f$ \lambda \f$ is a phase mobility
 * and \f$ f \f$ a phase fractional flow function,
 * \f$ \boldsymbol K \f$ is the intrinsic permeability and \f$ p_c = p_c(S_w) \f$ the capillary pressure.
 *
 * \tparam TypeTag The Type Tag
 */
template<class TypeTag>
class CapillaryDiffusion: public DiffusivePart<TypeTag>
{
private:
    using GridView = typename GetPropType<TypeTag, Properties::GridGeometry>::GridView;
      using Scalar = GetPropType<TypeTag, Properties::Scalar>;
      using Problem = GetPropType<TypeTag, Properties::Problem>;
      using Indices = typename GetPropType<TypeTag, Properties::ModelTraits>::Indices;

      using SpatialParams = GetPropType<TypeTag, Properties::SpatialParams>;

      using FluidSystem = GetPropType<TypeTag, Properties::FluidSystem>;
      using FluidState = GetPropType<TypeTag, Properties::FluidState>;

      using CellData = GetPropType<TypeTag, Properties::CellData>;

      using BoundaryTypes = GetPropType<TypeTag, Properties::SequentialBoundaryTypes>;
      using SolutionTypes = GetProp<TypeTag, Properties::SolutionTypes>;
      using PrimaryVariables = typename SolutionTypes::PrimaryVariables;

    enum
    {
        dim = GridView::dimension, dimWorld = GridView::dimensionworld
    };
    enum
    {
        wPhaseIdx = Indices::wPhaseIdx, nPhaseIdx = Indices::nPhaseIdx,
        pressEqIdx = Indices::pressureEqIdx
    };

    using Intersection = typename GridView::Intersection;
    using DimVector = Dune::FieldVector<Scalar, dim>;
    using GlobalPosition = Dune::FieldVector<Scalar, dimWorld>;
    using DimMatrix = Dune::FieldMatrix<Scalar, dim, dim>;

public:

    /*!
     * \brief Returns capillary diffusion term for current element face
     *
     *  \param flux        Flux vector (gets the flux from the function)
     *  \param intersection  Intersection of two grid elements/global boundary
     *  \param satI           saturation of current element
     *  \param satJ           saturation of neighbor element
     *  \param pcGradient     gradient of capillary pressure between element I and J
     */
    void getFlux (DimVector& flux, const Intersection& intersection, Scalar satI, Scalar satJ,
                  const DimVector& pcGradient) const
    {
        auto element = intersection.inside();
        // get global coordinate of cell center
        const GlobalPosition& globalPos = element.geometry().center();

        int globalIdxI = problem_.variables().index(element);
        CellData& CellDataI = problem_.variables().cellData(globalIdxI);

        // get geometry type of face
        //Dune::GeometryType faceGT = isIt->geometryInInside().type();

        Scalar temperature = problem_.temperature(element);
        Scalar referencePressure = problem_.referencePressure(element);

        //get lambda_bar = lambda_n*f_w
        Scalar mobBar = 0;
        Scalar mobilityWI = 0;
        Scalar mobilityNwI = 0;

        // old material law interface is deprecated: Replace this by
        // const auto& fluidMatrixInteraction = spatialParams.fluidMatrixInteractionAtPos(element.geometry().center());
        // after the release of 3.3, when the deprecated interface is no longer supported
        const auto fluidMatrixInteraction = Deprecated::makePcKrSw(Scalar{}, problem_.spatialParams(), element);

        if (preComput_)
        {
            mobilityWI = CellDataI.mobility(wPhaseIdx);
            mobilityNwI = CellDataI.mobility(nPhaseIdx);
        }
        else
        {
            FluidState fluidState;
            fluidState.setPressure(wPhaseIdx, referencePressure);
            fluidState.setPressure(nPhaseIdx, referencePressure);
            fluidState.setTemperature(temperature);
            mobilityWI = fluidMatrixInteraction.krw(satI);
            mobilityWI /= FluidSystem::viscosity(fluidState, wPhaseIdx);
            mobilityNwI = fluidMatrixInteraction.krn(satI);
            mobilityNwI /= FluidSystem::viscosity(fluidState, nPhaseIdx);
        }

        DimMatrix meanPermeability(0);

        if (intersection.neighbor())
        {
            // access neighbor
            auto neighbor = intersection.outside();

            int globalIdxJ = problem_.variables().index(neighbor);
            CellData& cellDataJ = problem_.variables().cellData(globalIdxJ);

            // neighbor cell center in global coordinates
            const GlobalPosition& globalPosNeighbor = neighbor.geometry().center();

            // distance vector between barycenters
            DimVector distVec = globalPosNeighbor - globalPos;

            // compute distance between cell centers
            Scalar dist = distVec.two_norm();

            DimVector unitDistVec(distVec);
            unitDistVec /= dist;

            // get permeability
            problem_.spatialParams().meanK(meanPermeability,
                    problem_.spatialParams().intrinsicPermeability(element),
                    problem_.spatialParams().intrinsicPermeability(neighbor));


            Scalar mobilityWJ = 0;
            Scalar mobilityNwJ = 0;
            //get lambda_bar = lambda_n*f_w
            if(preComput_)
            {
                mobilityWJ = cellDataJ.mobility(wPhaseIdx);
                mobilityNwJ = cellDataJ.mobility(nPhaseIdx);
            }
            else
            {
                FluidState fluidState;
                fluidState.setPressure(wPhaseIdx, referencePressure);
                fluidState.setPressure(nPhaseIdx, referencePressure);
                fluidState.setTemperature(temperature);

                // old material law interface is deprecated: Replace this by
                // const auto& fluidMatrixInteractionNeighbor = spatialParams.fluidMatrixInteractionAtPos(neighbor.geometry().center());
                // after the release of 3.3, when the deprecated interface is no longer supported
                const auto fluidMatrixInteractionNeighbor = Deprecated::makePcKrSw(Scalar{}, problem_.spatialParams(), neighbor);

                mobilityWJ = fluidMatrixInteractionNeighbor.krw(satJ);
                mobilityWJ /= FluidSystem::viscosity(fluidState, wPhaseIdx);
                mobilityNwJ = fluidMatrixInteractionNeighbor.krn(satJ);
                mobilityNwJ /= FluidSystem::viscosity(fluidState, nPhaseIdx);
            }
            Scalar mobilityWMean = 0.5*(mobilityWI + mobilityWJ);
            Scalar mobilityNwMean = 0.5*(mobilityNwI + mobilityNwJ);
            mobBar = mobilityWMean*mobilityNwMean/(mobilityWMean+mobilityNwMean);
         }//end intersection with neighbor
        else
        {
            BoundaryTypes bcTypes;
            problem_.boundaryTypes(bcTypes, intersection);
            if (bcTypes.isNeumann(pressEqIdx))
            {
                PrimaryVariables priVars;
                problem_.neumann(priVars, intersection);
                if (priVars[wPhaseIdx] == 0)
                {
                    flux = 0;
                    return;
                }
            }
            // get permeability
            problem_.spatialParams().meanK(meanPermeability,
                    problem_.spatialParams().intrinsicPermeability(element));

            Scalar mobilityWJ = 0;
            Scalar mobilityNwJ = 0;

            //calculate lambda_n*f_w at the boundary
            FluidState fluidState;
            fluidState.setPressure(wPhaseIdx, referencePressure);
            fluidState.setPressure(nPhaseIdx, referencePressure);
            fluidState.setTemperature(temperature);
            mobilityWJ = fluidMatrixInteraction.krw(satJ);
            mobilityWJ /= FluidSystem::viscosity(fluidState, wPhaseIdx);
            mobilityNwJ = fluidMatrixInteraction.krn(satJ);
            mobilityNwJ /= FluidSystem::viscosity(fluidState, nPhaseIdx);

            Scalar mobWMean = 0.5 * (mobilityWI + mobilityWJ);
            Scalar mobNwMean = 0.5 * (mobilityNwI + mobilityNwJ);

            mobBar = mobWMean * mobNwMean / (mobWMean + mobNwMean);
        }

        // set result to K*grad(pc)
        meanPermeability.mv(pcGradient, flux);

        // set result to f_w*lambda_n*K*grad(pc)
        flux *= mobBar;
    }

    /*!
     * \brief Constructs a CapillaryDiffusion object
     *
     *  \param problem A problem class object
     */
    CapillaryDiffusion (Problem& problem)
    : DiffusivePart<TypeTag>(problem), problem_(problem), preComput_(getPropValue<TypeTag, Properties::PrecomputedConstRels>())
    {}

private:
    Problem& problem_;//problem data
    const bool preComput_;//if preCompute = true the mobilities are taken from the variable object,
                          //if preCompute = false new mobilities will be taken (for implicit Scheme)
};
} // end namespace Dumux

#endif
