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
#ifndef DUMUX_FVMPFAL2PFABOUND3DVELOCITIES2P_HH
#define DUMUX_FVMPFAL2PFABOUND3DVELOCITIES2P_HH

#include "fvmpfal3dpressure2p.hh"
#include "fvmpfal3dvelocity2p.hh"

/**
 * @file
 * @brief  Base class for defining an instance of a numerical diffusion model
 */

namespace Dumux
{
template<class TypeTag> class FvMpfaL3dPressureVelocity2p: public FvMpfaL3dPressure2p<TypeTag>
{
    typedef FvMpfaL3dPressure2p<TypeTag> ParentType;

    typedef typename GET_PROP_TYPE(TypeTag, GridView) GridView;

    enum
        {
            dim = GridView::dimension, dimWorld = GridView::dimensionworld
        };

    typedef typename GET_PROP_TYPE(TypeTag, Problem) Problem;
    typedef typename GET_PROP_TYPE(TypeTag, CellData) CellData;

    typedef typename GridView::template Codim<dim>::Iterator VertexIterator;
    typedef typename GET_PROP_TYPE(TypeTag, MPFAInteractionVolume) InteractionVolume;

public:
    FvMpfaL3dPressureVelocity2p(Problem& problem) :
        ParentType(problem), problem_(problem), velocity_(problem)
    {}

    void calculateVelocity();

public:

    void updateVelocity()
    {
        this->updateMaterialLaws();

        //reset velocities
        int size = problem_.gridView().size(0);
        for (int i = 0; i < size; i++)
        {
            CellData& cellData = problem_.variables().cellData(i);
            cellData.fluxData().resetVelocity();
        }

        calculateVelocity();
    }

    void initialize(bool solveTwice = true)
    {
        ParentType::initialize();
        velocity_.initialize();
        calculateVelocity();

        return;
    }

    void update()
    {
        ParentType::update();
        calculateVelocity();
    }

    //! \brief Write data files
    /*  \param name file name */
    template<class MultiWriter>
    void addOutputVtkFields(MultiWriter &writer)
    {
        ParentType::addOutputVtkFields(writer);
        velocity_.addOutputVtkFields(writer);
    }

private:
    Problem& problem_;
    FvMpfaL3dVelocity2p<TypeTag> velocity_;
};
// end of template

// only for 3-D general hexahedron
template<class TypeTag>
void FvMpfaL3dPressureVelocity2p<TypeTag>::calculateVelocity()
{
    // run through all vertices
    VertexIterator vItEnd = problem_.gridView().template end<dim>();
    for (VertexIterator vIt = problem_.gridView().template begin<dim>(); vIt != vItEnd; ++vIt)
    {
        int globalVertIdx = problem_.variables().index(*vIt);

        InteractionVolume& interactionVolume = this->interactionVolumes_.interactionVolume(globalVertIdx);

        // inner interactionvolume
        if (interactionVolume.isInnerVolume())
        {
            velocity_.calculateInnerInteractionVolumeVelocity(interactionVolume, this->interactionVolumes_, this->transmissibilityCalculator_);
        }
        // at least one face on boundary! (boundary interactionvolume)
        else
        {
            velocity_.calculateBoundaryInteractionVolumeVelocity(interactionVolume, this->interactionVolumes_, this->transmissibilityCalculator_);
        } // end boundaries

    } // end vertex iterator

    return;
}

}
// end of Dune namespace
#endif
