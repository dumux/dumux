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
#include <dumux/common/propertysystem.hh>

#include <dune/common/version.hh>
#include <dune/common/fvector.hh>
#include <dune/common/fmatrix.hh>

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
    typedef typename GET_PROP_TYPE(TypeTag, Scalar)               Scalar;
    typedef typename GET_PROP_TYPE(TypeTag, Problem)              Problem;
    typedef typename GET_PROP_TYPE(TypeTag, GridView)             GridView;
    typedef typename GridView::template Codim<0>::Iterator        ElementIterator;

    typedef typename GET_PROP_TYPE(TypeTag, FVElementGeometry)    FVElementGeometry;
    typedef typename GET_PROP_TYPE(TypeTag, DofMapper)            DofMapper;
    typedef typename GET_PROP_TYPE(TypeTag, ElementSolutionVector) ElementSolutionVector;
    typedef typename GET_PROP_TYPE(TypeTag, SolutionVector)       SolutionVector;
    typedef typename GET_PROP_TYPE(TypeTag, MaterialLaw)          MaterialLaw;
    typedef typename MaterialLaw::Params                          aterialLawParams;

    typedef typename GET_PROP_TYPE(TypeTag, ElementVolumeVariables) ElementVolumeVariables;
    typedef typename GET_PROP_TYPE(TypeTag, FluidSystem)          FluidSystem;

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

    typedef Dune::FieldVector<Scalar, dimWorld> GlobalPosition;

public:
    /*
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

        // Looping over all elements of the domain
        ElementIterator eEndIt = problem.gridView().template end<0>();
        for (ElementIterator eIt = problem.gridView().template begin<0>() ; eIt != eEndIt; ++eIt)
        {
            // updating the volume variables
            fvGeometry.update(problem.gridView(), *eIt);
            elemVolVars.update(problem, *eIt, fvGeometry, false);

            // number of scv
            const unsigned int numScv = fvGeometry.numScv;

            for (unsigned int scvIdx=0; scvIdx < numScv; ++scvIdx)
            {
                // find some global identification
#if DUNE_VERSION_NEWER(DUNE_COMMON, 2, 4)
                const unsigned int vIdxGlobal = problem.vertexMapper().subIndex(*eIt, scvIdx, dim);
#else
                const unsigned int vIdxGlobal = problem.vertexMapper().map(*eIt, scvIdx, dim);
#endif
                // only write out if the vertex was not already visited
                if (isVisited[vIdxGlobal])
                    continue;

                isVisited[vIdxGlobal] = true ;

                // Getting the spatial coordinate
                const GlobalPosition & globalPosCurrent
                    = fvGeometry.subContVol[scvIdx].global;

                std::ofstream dataFile;
                const unsigned int timeStepIndex = problem.timeManager().timeStepIndex() ;

                // filename of the output file
                std::string fileName = problem.name();

                // filename consists of problem name + some function argument + .dat
                // this way several plot overlines can be written from one simulation
                fileName += appendOutputName;
                fileName +=".dat";

                // Writing a header into the output
                if (timeStepIndex == 0)
                {
                    dataFile.open(fileName.c_str());
                    dataFile << "# This is a DuMuX output file for further processing with the preferred graphics program of your choice. \n";

                    dataFile << "# This output file was written from "<< __FILE__ << ", line " <<__LINE__ << "\n";
                    dataFile << "# This output file was generated from code compiled at " << __TIME__ <<", "<< __DATE__<< "\n";
                    dataFile << "\n";
                    dataFile << "# Header\n";
                    dataFile << "#timestep\t time\t\t \t\t x \t\t y  \t\tSw \t\t\t Tw\t\t Tn\t Ts \t xH2On \t xH2OnEquil \t xN2w \txN2wEquil\n";
                    dataFile.close();
                }

                // write output if the current location is between two specified points
                if (isBetween(globalPosCurrent, pointOne, pointTwo))
                {
                    const Scalar time         = problem.timeManager().time();
                    const Scalar saturationW  = elemVolVars[scvIdx].fluidState().saturation(wPhaseIdx);
                    const Scalar Tw           = elemVolVars[scvIdx].fluidState().temperature(wPhaseIdx);
                    const Scalar Tn           = elemVolVars[scvIdx].fluidState().temperature(nPhaseIdx);
                    const Scalar Ts           = elemVolVars[scvIdx].fluidState().temperature(sPhaseIdx);
                    const Scalar xH2On        = elemVolVars[scvIdx].fluidState().moleFraction(nPhaseIdx, wCompIdx);
                    const Scalar xH2OnEquil   = elemVolVars[scvIdx].xEquil(nPhaseIdx, wCompIdx);
                    const Scalar xN2w         = elemVolVars[scvIdx].fluidState().moleFraction(wPhaseIdx, nCompIdx);
                    const Scalar xN2wEquil    = elemVolVars[scvIdx].xEquil(wPhaseIdx, nCompIdx);

                    // actual output into the text file
                    // This could be done much more efficiently
                    // if writing to the text file would be done all at once.
                    dataFile.open(fileName.c_str() , std::ios::app);
                    dataFile << timeStepIndex
                            <<"\t\t\t"
                            << time
                            << "\t\t\t"
                            << globalPosCurrent[0]
                            << "\t\t\t"
                            << globalPosCurrent[1];

                    dataFile <<"\t\t\t"
                                << saturationW
                                <<"\t\t" << Tw
                                <<"\t\t" << Tn
                                <<"\t\t" << Ts
                                <<"\t\t" << xH2On
                                <<"\t\t" << xH2OnEquil
                                <<"\t\t" << xN2w
                                <<"\t\t" << xN2wEquil
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
    const bool isBetween(const GlobalPosition & globalPosCurrent,
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
