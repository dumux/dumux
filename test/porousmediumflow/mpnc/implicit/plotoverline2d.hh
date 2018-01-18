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
 * \brief Plot variables over a line specified by two arguments.
 *        These output files are meant for visualization with another
 *        program (matlab, gnuplot...)
 *
 */
#ifndef DUMUX_PLOTOVERLINE_2D_HH
#define DUMUX_PLOTOVERLINE_2D_HH

#include <dumux/common/valgrind.hh>
#include <dumux/common/properties.hh>

#include <dune/common/fvector.hh>
#include <dune/common/fmatrix.hh>
#include <dune/grid/common/rangegenerators.hh>

#include <iostream>
#include <fstream>

namespace Dumux
{
namespace Properties
{
    NEW_PROP_TAG(Scalar);
    NEW_PROP_TAG(Problem);
    NEW_PROP_TAG(GridView);
    NEW_PROP_TAG(DofMapper);
    NEW_PROP_TAG(FluidSystem);
    NEW_PROP_TAG(ElementSolutionVector);
    NEW_PROP_TAG(SolutionVector);
    NEW_PROP_TAG(FVElementGeometry);
    NEW_PROP_TAG(TwoPIAIndices);
    NEW_PROP_TAG(NumEq);
    NEW_PROP_TAG(MaterialLaw);
    NEW_PROP_TAG(ElementVolumeVariables);
    NEW_PROP_TAG(AwnSurface);
    NEW_PROP_TAG(AwnSurfaceParams);
}

template<class TypeTag>
class PlotOverLine2D
{
    using Scalar = typename GET_PROP_TYPE(TypeTag, Scalar);
    using Problem = typename GET_PROP_TYPE(TypeTag, Problem);
    using GridView = typename GET_PROP_TYPE(TypeTag, GridView);

    using FVElementGeometry = typename GET_PROP_TYPE(TypeTag, FVGridGeometry)::LocalView;
    using DofMapper = typename GET_PROP_TYPE(TypeTag, DofMapper);
    using ElementSolutionVector = typename GET_PROP_TYPE(TypeTag, ElementSolutionVector);
    using SolutionVector = typename GET_PROP_TYPE(TypeTag, SolutionVector);
    using MaterialLaw = typename GET_PROP_TYPE(TypeTag, MaterialLaw);
    using aterialLawParams = typename MaterialLaw::Params;

    using ElementVolumeVariables = typename GET_PROP_TYPE(TypeTag, ElementVolumeVariables);
    using FluidSystem = typename GET_PROP_TYPE(TypeTag, FluidSystem);

    enum {
        wPhaseIdx = FluidSystem::wPhaseIdx,
        nPhaseIdx = FluidSystem::nPhaseIdx,
        sPhaseIdx = FluidSystem::sPhaseIdx,
        wCompIdx  = FluidSystem::wCompIdx,
        nCompIdx  = FluidSystem::nCompIdx,

        // Grid and world dimension
        dim         = GridView::dimension,
        dimWorld    = GridView::dimensionworld,
    };

    using GlobalPosition = Dune::FieldVector<Scalar, dimWorld>;

public:
    /*!
     * \brief A function that writes results over a line (like paraview's plotOverline into a text file.)
     *
     *        The writer needs to be called in postTimeStep().
     *
     *        This function puts output variables (TemperaturePhase, Saturation, t, tIndex, ...)
     *        over space (1D, over a line) into a text file,
     *        so they can be read in by another program like matlab.
     *        The file can be found by the extension: dat
     */
    void write(const Problem & problem,
               const GlobalPosition & pointOne,
               const GlobalPosition & pointTwo,
               const bool appendData,
               const std::string appendOutputName = "")
    {
        static_assert(dim==2, "this implements plot over Line: so far this works only for 2D");

        FVElementGeometry fvGeometry;
        ElementVolumeVariables elemVolVars;

        // so many vertices in the domain
        const unsigned int numGlobalVerts = problem.gridView().size(dim);

        // check whether a vertex was already visited by true / false vector
        std::vector<bool> isVisited(numGlobalVerts);
        std::fill(isVisited.begin(), isVisited.end(), false);


        // filename of the output file
        std::string fileName = problem.name();
        fileName += appendOutputName;
        fileName += ".dat";
        std::ofstream dataFile;
        const unsigned int timeStepIndex = problem.timeManager().timeStepIndex();

        // Writing a header into the output
        if (timeStepIndex == 0 || !appendData)
        {
            dataFile.open(fileName.c_str());
            dataFile << "# This is a DuMuX output file for further processing with the preferred graphics program of your choice. \n";

            dataFile << "# This output file was written from "<< __FILE__ << ", line " <<__LINE__ << "\n";
            dataFile << "# This output file was generated from code compiled at " << __TIME__ <<", "<< __DATE__<< "\n";
            dataFile << "\n";
            dataFile << "# Header\n";
            dataFile << "#timestep time x y Sw Tw Tn Ts xH2On xH2OnEquil xN2w xN2wEquil\n";
            dataFile.close();
        }

        // Looping over all elements of the domain
        for (const auto& element : elements(problem.gridView()))
        {
            // updating the volume variables
            fvGeometry.update(problem.gridView(), element);
            elemVolVars.update(problem, element, fvGeometry, false);

            // number of scv
            const unsigned int numScv = fvGeometry.numScv;

            for (unsigned int scvIdx=0; scvIdx < numScv; ++scvIdx)
            {
                // find some global identification
                const unsigned int vIdxGlobal = problem.vertexMapper().subIndex(element, scvIdx, dim);

                // only write out if the vertex was not already visited
                if (isVisited[vIdxGlobal])
                    continue;

                isVisited[vIdxGlobal] = true ;

                // Getting the spatial coordinate
                const GlobalPosition & globalPosCurrent
                    = fvGeometry.subContVol[scvIdx].global;

                // write output if the current location is between two specified points
                if (isBetween(globalPosCurrent, pointOne, pointTwo))
                {
                    const Scalar time         = problem.timeManager().time();
                    const Scalar saturationW  = elemVolVars[scvIdx].saturation(wPhaseIdx);
                    const Scalar Tw           = elemVolVars[scvIdx].temperature(wPhaseIdx);
                    const Scalar Tn           = elemVolVars[scvIdx].temperature(nPhaseIdx);
                    const Scalar Ts           = elemVolVars[scvIdx].temperature(sPhaseIdx);
                    const Scalar xH2On        = elemVolVars[scvIdx].moleFraction(nPhaseIdx, wCompIdx);
                    const Scalar xH2OnEquil   = elemVolVars[scvIdx].xEquil(nPhaseIdx, wCompIdx);
                    const Scalar xN2w         = elemVolVars[scvIdx].moleFraction(wPhaseIdx, nCompIdx);
                    const Scalar xN2wEquil    = elemVolVars[scvIdx].xEquil(wPhaseIdx, nCompIdx);

                    // actual output into the text file
                    // This could be done much more efficiently
                    // if writing to the text file would be done all at once.
                    dataFile.open(fileName.c_str(), std::ios::app);
                    dataFile << timeStepIndex
                            <<" "
                            << time
                            << " "
                            << globalPosCurrent[0]
                            << " "
                            << globalPosCurrent[1];

                    dataFile <<" "
                                << saturationW
                                <<" " << Tw
                                <<" " << Tn
                                <<" " << Ts
                                <<" " << xH2On
                                <<" " << xH2OnEquil
                                <<" " << xN2w
                                <<" " << xN2wEquil
                                ;
                    dataFile <<"\n";
                    dataFile.close();
                }
            }
        }
        return ;
    }

    /*!
     * \brief   Check whether the current point is on a line between two points
     */
    bool isBetween(const GlobalPosition & globalPosCurrent,
                         const GlobalPosition & pointOne,
                         const GlobalPosition & pointTwo) const
    {
        return (    pointOne[0] - globalPosCurrent[0] <= 1e-5
            and pointOne[1] - globalPosCurrent[1] <= 1e-5
            and globalPosCurrent[0] - pointTwo[0] <= 1e-5
            and globalPosCurrent[1] - pointTwo[1] <= 1e-5 );
    }
};

}
#endif
