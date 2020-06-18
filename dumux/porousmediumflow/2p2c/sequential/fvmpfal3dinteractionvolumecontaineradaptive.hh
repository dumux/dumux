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
 * \ingroup SequentialTwoPTwoCModel
 * \brief Interaction volume container for compositional adaptive 3-d (using MPFA L-method).
 */
#ifndef DUMUX_FVMPFAL3D_2P2CINTERACTIONVOLUMECONTAINER_ADAPTIVE_HH
#define DUMUX_FVMPFAL3D_2P2CINTERACTIONVOLUMECONTAINER_ADAPTIVE_HH

// dumux environment
#include <dumux/porousmediumflow/2p/sequential/diffusion/mpfa/lmethod/3dinteractionvolumecontaineradaptive.hh>

namespace Dumux {

/*!
 * \ingroup SequentialTwoPTwoCModel
 * \brief Interaction volume container for compositional adaptive 3-d (using MPFA L-method) Model
 *
 * Container class which stores MPFA-interaction-volume information for each vertex of a DUNE grid.
 *
 * For the compositional case, we predominantly apply a TPFA, which is the reason why the model loops
 * over all cells and their interfaces. This means that the in the case of a hanging node, the mpfa
 * has to be build from an interface, so the mpfa implementation needs to be access, which is facilitated
 * by this class.
 * Provides methods that determines which MPFA case needs to be considered.
 */
template<class TypeTag>
class FvMpfaL3d2P2CInteractionVolumeContainerAdaptive : public FvMpfaL3dInteractionVolumeContainerAdaptive<TypeTag>
{
    using ParentType = FvMpfaL3dInteractionVolumeContainerAdaptive<TypeTag>;
    using GridView = typename GetPropType<TypeTag, Properties::GridGeometry>::GridView;

    enum
    {
        dim = GridView::dimension, dimWorld = GridView::dimensionworld
    };

    using Scalar = GetPropType<TypeTag, Properties::Scalar>;using Problem = GetPropType<TypeTag, Properties::Problem>;

    using Indices = typename GetPropType<TypeTag, Properties::ModelTraits>::Indices;

    using BoundaryTypes = GetPropType<TypeTag, Properties::SequentialBoundaryTypes>;
    using PrimaryVariables = typename GetProp<TypeTag, Properties::SolutionTypes>::PrimaryVariables;

    using GridTypeIndices = GetPropType<TypeTag, Properties::GridTypeIndices>;

    using Element = typename GridView::Traits::template Codim<0>::Entity;
    using ElementGeometry = typename Element::Geometry;
    using Vertex = typename GridView::Traits::template Codim<dim>::Entity;

    using IntersectionIterator = typename GridView::IntersectionIterator;
    using Intersection = typename GridView::Intersection;
    using IntersectionGeometry = typename Intersection::Geometry;

    using GlobalPosition = Dune::FieldVector<Scalar, dimWorld>;
    using DimMatrix = Dune::FieldMatrix<Scalar, dim, dim>;

    using DimVector = Dune::FieldVector<Scalar, dim>;

    enum
    {
        pressureEqIdx = Indices::pressureEqIdx,
    };
    enum
    {
        innerEdgeFace = 2, innerSideFace = 1
    };
    enum
    {
        realFaceArea = 0, fluxFaceArea = 1
    };

public:
    //! Type for storing an MPFA-interaction-volume. (Usually of type FvMpfaL3dInteractionVolume or FvMpfaL3dInteractionVolumeAdaptive)
    using InteractionVolume = GetPropType<TypeTag, Properties::MPFAInteractionVolume>;

    using GlobalInteractionVolumeVector = std::vector<InteractionVolume>;
    using FaceAreaVector = std::vector<Dune::FieldVector<Dune::FieldVector<Scalar, 2>, 2*dim> >;

    void storeBoundaryInteractionVolume(InteractionVolume& interactionVolume, const Vertex& vertex);

    inline int getMpfaCase8cells(const IntersectionIterator& isIt,
                                    const int localidxLarge,
                                    InteractionVolume& interactionVolume,
                                    bool& properFluxDirection);
    inline int getMpfaCase6cells(const IntersectionIterator& isIt,
                                    InteractionVolume& interactionVolume,
                                    bool& properFluxDirection);
    inline int getMpfaCase2or4cells(const IntersectionIterator& isIt,
                                    InteractionVolume& interactionVolume,
                                    bool& properFluxDirection);

    FvMpfaL3d2P2CInteractionVolumeContainerAdaptive(Problem& problem)
        : ParentType(problem), problem_(problem)
    {
        if (dim != 3)
        {
            DUNE_THROW(Dune::NotImplemented, "Dimension not supported!");
        }
    }
private:
    Problem& problem_;
};

/*!
 * \brief Overwrites the method from the base class FvMpfaL3dInteractionVolumeContainerAdaptive
 * On each boundary, a TPFA is used in compositional models. Therefore we do not need to store interaction
 * volume containers on the boundary cells.
 */
template<class TypeTag>
void FvMpfaL3d2P2CInteractionVolumeContainerAdaptive<TypeTag>::storeBoundaryInteractionVolume(InteractionVolume& interactionVolume, const Vertex& vertex)
{
    return;
}

/*!
 * \brief Determine the subVolumeFaceIdx for a given intersection and interactionVolume with 8 cells
 *
 * The MPFA is about to be calculated through an intersection, and to do so, its place in the
 * local indexing scheme, i.e. its subVolumeFaceIdx, has to be found. This method
 * investigates the case (see FvMpfaL3dInteractionVolumeAdaptive.HangingNodeTypes ) if
 * 8 cells are present in the current interaction region: A Interaction region where the non-adaptive
 * MPFA-model is applied.
 * This requires a local Index of the "large" cell (where the hanging node rests) to get the right
 * subVolumeFaceIdx.
 *  \param isIt The iterator of the intersection the mpfa should be calculated for
 *  \param localIdxLarge The (local) Index of the large cell (on which the hanging node lives)
 *  \param interactionVolume The interaction volume (FvMpfaL3dInteractionVolumeAdaptive) of interest
 *  \param properFluxDirection Indicates whether the flux through the intersection aligns with its normal
 *  \return The Subvolume Face Idx required by the methods in FvMpfaL3dTransmissibilityCalculator
 */
template<class TypeTag>
inline int FvMpfaL3d2P2CInteractionVolumeContainerAdaptive<TypeTag>::getMpfaCase8cells(const IntersectionIterator& isIt,
                                const int localIdxLarge,
                                InteractionVolume& interactionVolume,
                                bool& properFluxDirection)
{
    int subVolumeFaceIdx = -1;

    /************** determine subVolumeFace of interest */
    int localFaceIdxInside = isIt->indexInInside();
    switch(localIdxLarge)
    {
    case 0:{    //cell I = 7
        if(localFaceIdxInside == 0) //left
        {
            subVolumeFaceIdx = 6;
            properFluxDirection = false;
        }
        else if(localFaceIdxInside == 2) //front
            subVolumeFaceIdx = 5;
        else if(localFaceIdxInside == 4) //bottom
            subVolumeFaceIdx = 10;
        break;
        }
    case 1:{    //cell I = 6
        if(localFaceIdxInside == 1) //right
            subVolumeFaceIdx = 6;
        else if(localFaceIdxInside == 2)    //front
        {
            subVolumeFaceIdx = 7;
            properFluxDirection = false;
        }
        else if(localFaceIdxInside == 4)    //bottom
        {
            subVolumeFaceIdx = 11;
            properFluxDirection = false;
        }
        break;
        }
    case 2:{    // cellI = 5
        if(localFaceIdxInside == 0) // left
            subVolumeFaceIdx = 4;
        else if(localFaceIdxInside == 3)    //rear
        {
            subVolumeFaceIdx = 5;
            properFluxDirection = false;
        }
        else if(localFaceIdxInside == 4)    //bottom
        {
            subVolumeFaceIdx = 9;
            properFluxDirection = false;
        }
        break;
        }
    case 3:{    //cellI = 4
        if(localFaceIdxInside == 1) //right
        {
            subVolumeFaceIdx = 4;
            properFluxDirection = false;
        }
        else if(localFaceIdxInside == 3) //rear
            subVolumeFaceIdx = 7;
        else if(localFaceIdxInside == 4)    //botm
            subVolumeFaceIdx = 8;
        break;
        }
    case 4:{    //cellI = 3
        if(localFaceIdxInside == 0) //left
            subVolumeFaceIdx = 2;
        else if(localFaceIdxInside == 2) //front
        {
            subVolumeFaceIdx = 1;
            properFluxDirection = false;
        }
        else if(localFaceIdxInside == 5) //top
        {
            subVolumeFaceIdx = 10;
            properFluxDirection = false;
        }
        break;
        }
    case 5:{    //cellI = 2
        if(localFaceIdxInside == 1) //right
        {
            subVolumeFaceIdx = 2;
            properFluxDirection = false;
        }
        else if(localFaceIdxInside == 2) //front
            subVolumeFaceIdx = 3;
        else if(localFaceIdxInside == 5) //top
            subVolumeFaceIdx = 11;
        break;
        }
    case 6:{    //cellI = 1
        if(localFaceIdxInside == 0) // left
        {
            subVolumeFaceIdx = 0;
            properFluxDirection = false;
        }
        else if(localFaceIdxInside == 3) //rear
            subVolumeFaceIdx = 1;
        else if(localFaceIdxInside == 5) //top
            subVolumeFaceIdx = 9;
        break;
        }
    case 7:{    //cellI= 0
        if(localFaceIdxInside == 1) //right
            subVolumeFaceIdx = 0;
        else if(localFaceIdxInside == 3)//rear
        {
            subVolumeFaceIdx = 3;
            properFluxDirection = false;
        }
        else if(localFaceIdxInside == 5)//top
        {
            subVolumeFaceIdx = 8;
            properFluxDirection = false;
        }
        break;
        }
    }
    return subVolumeFaceIdx;
}

/*!
 * \brief Determine the subVolumeFaceIdx for a given intersection and interactionVolume with 6 cells
 *
 * The MPFA is about to be calculated through an intersection, and to do so, its place in the
 * local indexing scheme, i.e. its subVolumeFaceIdx, has to be found. This method
 * investigates the case (see FvMpfaL3dInteractionVolumeAdaptive.HangingNodeTypes ) if
 * 6 cells are present in the current interaction region.
 * \param isIt The iterator of the intersection the mpfa should be calculated for
 * \param interactionVolume The interaction volume (FvMpfaL3dInteractionVolumeAdaptive) of interest
 * \param properFluxDirection Indicates whether the flux through the intersection aligns with its normal
 * \return The Subvolume Face Idx required by the methods in FvMpfaL3dTransmissibilityCalculator
 */
template<class TypeTag>
inline int FvMpfaL3d2P2CInteractionVolumeContainerAdaptive<TypeTag>::getMpfaCase6cells(const IntersectionIterator& isIt,
                                InteractionVolume& interactionVolume,
                                bool& properFluxDirection)
{
    int mapI(-5), mapJ(-5);
    int mapI2(-5), mapJ2(-5);
    // search for cell I
    for (int i= 0; i<interactionVolume.getElementNumber(); i++)
    {
        if(isIt->inside() == interactionVolume.getSubVolumeElement(i))
        {
            if(mapI == -5)
                mapI = i;
            else
                mapI2 = i;
        }

        if(isIt->outside() == interactionVolume.getSubVolumeElement(i))
        {
            if(mapJ == -5)
                mapJ = i;
            else
                mapJ2 =i;
        }
    }

    for(int passionfruit = 0; passionfruit <=1; passionfruit++) // loop at most twice
    {
        if(mapI== 0 || mapJ == 0)
        {
            if(mapI== 2)
                return 3;
            else if (mapJ== 2){
                properFluxDirection = false;
                return 3;
            }
            else if(mapJ== 1)
                return 0;
            else if (mapI==1){
                properFluxDirection = false;
                return 0;
            }
        }
        // no "else if" because there is also 4 - 0 treated under I or J == 4
        if(mapI== 1 || mapJ == 1)
        {
            if (mapI== 3){
                properFluxDirection = false;
                return 1;
            }
            else if(mapJ== 3)
                return 1;
            else if(mapI == 5){
                properFluxDirection = false;
                return 9;
            }
            else if(mapJ == 5)
                return 9;
        }
        else if(mapI== 2 || mapJ == 2)
        {
            if(mapI == 6){
                properFluxDirection = false;
                return 11;
            }
            else if(mapJ ==6)
                return 11;
        }
        else if(mapI== 3 || mapJ == 3)
        {
            if (mapJ== 7){
                properFluxDirection = false;
                return 10;
            }
            else if(mapI== 7)
                return 10;
        }
        else if(mapI== 4 || mapJ == 4)
        {
            // this has to be subVolFaceIdx 8 because 7 would mean
            // a case that should be modelled by tpfa
            if (mapI== 0){
                properFluxDirection = false;
                return 8;
            }
            else if(mapJ== 0)
                return 8;
            else if(mapI == 6){
                properFluxDirection = false;
                return 7;
            }
            else if(mapJ == 6)
                return 7;
        }
        else if(mapI== 5 || mapJ == 5)
        {
            if (mapJ== 7){
                properFluxDirection = false;
                return 5;
            }
            else if(mapI== 7)
                return 5;
        }
        else if(mapI== 6 || mapJ == 6)
        {
            if (mapI== 7){
                properFluxDirection = false;
                return 6;
            }
            else if(mapJ== 7)
                return 6;
        }
        // configuration not found: I or J are twice in the interaction volume: investigate other localIdx
        if(mapI2 != -5)
            mapI = mapI2;
        else if(mapJ2 != -5)
            mapJ = mapJ2;
    }

    Dune::dgrave << " Could not find "<< interactionVolume.getHangingNodeType() <<" case  configuration for I = "
            << problem_.variables().index(isIt->inside()) << " localIdx " << mapI << " , "
            << problem_.variables().index(isIt->outside()) << " localIdx " << mapJ << std::endl;

    return -1;
}

/*!
 * \brief Determine the subVolumeFaceIdx for a given intersection and interactionVolume with 2/4 cells
 *
 * The MPFA is about to be calculated through an intersection, and to do so, its place in the
 * local indexing scheme, i.e. its subVolumeFaceIdx, has to be found. This method
 * investigates the case (see FvMpfaL3dInteractionVolumeAdaptive.HangingNodeTypes ) if
 * 2 or 4 cells are present in the current interaction region.
 *  \param isIt The iterator of the intersection the mpfa should be calculated for
 *  \param interactionVolume The interaction volume (FvMpfaL3dInteractionVolumeAdaptive) of interest
 *  \param properFluxDirection Indicates whether the flux through the intersection aligns with its normal
 *  \return The Subvolume Face Idx required by the methods in FvMpfaL3dTransmissibilityCalculator
 */
template<class TypeTag>
inline int FvMpfaL3d2P2CInteractionVolumeContainerAdaptive<TypeTag>::getMpfaCase2or4cells(const IntersectionIterator& isIt,
                                InteractionVolume& interactionVolume,
                                bool& properFluxDirection)
{
    int mapI(-5), mapJ(-5);
    int mapI2(-5), mapJ2(-5);
    // search for cell I
    for (int i= 0; i<interactionVolume.getElementNumber(); i++)
    {
        if(isIt->inside() == interactionVolume.getSubVolumeElement(i))
        {
            if(mapI == -5)
                mapI = i;
            else
                mapI2 = i;
        }

        if(isIt->outside() == interactionVolume.getSubVolumeElement(i))
        {
            if(mapJ == -5)
                mapJ = i;
            else
                mapJ2 =i;
        }
    }

    for(int passionfruit = 0; passionfruit <=1; passionfruit++) // loop at most twice
    {

        if(mapI== 0 || mapJ == 0)
        {
            if (mapJ== 2){
                properFluxDirection = false;
                return 3;
            }
            else if(mapI== 2)
                return 3;
            else if (mapJ== 4){
                properFluxDirection = false;
                return 8;
            }
            else if(mapI== 4)
                return 8;
            //it has to be case 0
            else if (mapJ == 0)
                properFluxDirection = false;

            return 0;
        }
        else if(mapI== 1 || mapJ == 1)
        {
            if (mapI== 3){
                properFluxDirection = false;
                return 1;
            }
            else if(mapJ== 3)
                return 1;
            else if(mapI == 5){
                properFluxDirection = false;
                return 9;
            }
            else if(mapJ == 5)
                return 9;
        }
        else if(mapI== 2 || mapJ == 2)
        {
            if (mapJ== 3){
                properFluxDirection = false;
                return 2;
            }
            else if(mapI== 3)
                return 2;
            else if(mapI == 6){
                properFluxDirection = false;
                return 11;
            }
            else if(mapJ ==6)
                return 11;
        }
        else if(mapI== 3 || mapJ == 3)
        {
            if (mapJ== 7){
                properFluxDirection = false;
                return 10;
            }
            else if(mapI== 7)
                return 10;
        }
        else if(mapI== 4 || mapJ == 4)
        {
            if (mapJ== 5){
                properFluxDirection = false;
                return 4;
            }
            else if(mapI== 5)
                return 4;
            else if(mapI == 6){
                properFluxDirection = false;
                return 7;
            }
            else if(mapJ == 6)
                return 7;
        }
        else if(mapI== 5 || mapJ == 5)
        {
            if (mapJ== 7){
                properFluxDirection = false;
                return 5;
            }
            else if(mapI== 7)
                return 5;
        }
        else if(mapI== 6 || mapJ == 6)
        {
            if (mapI== 7){
                properFluxDirection = false;
                return 6;
            }
            else if(mapJ== 7)
                return 6;
        }
        // configuration not found: I or J are twice in the interaction volume: investigate other localIdx
        if(mapI2 != -5)
            mapI = mapI2;
        else if(mapJ2 != -5)
            mapJ = mapJ2;
    }

    Dune::dgrave << " Could not find "<< interactionVolume.getHangingNodeType() <<" case  configuration for I = "
            << problem_.variables().index(isIt->inside()) << " localIdx " << mapI << " and "
            << problem_.variables().index(isIt->outside()) << " localIdx " << mapJ << std::endl;
    return -1;
}


} // end namespace Dumux
#endif
