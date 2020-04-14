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

#ifndef DUMUX_FVMPFAL3DINTERACTIONVOLUME_ADAPTIVE_HH
#define DUMUX_FVMPFAL3DINTERACTIONVOLUME_ADAPTIVE_HH

/**
 * @file
 * @brief  Class including the information of an interaction volume of a MPFA 3D method that does not change with time
 * @author Markus Wolff
 */

#include "linteractionvolume3d.hh"

namespace Dumux
{
//! \cond \private
// Mapper for local interaction volume indices (see doc/docextra/3dmpfa).
class IndexTranslatorAdaptive: public IndexTranslator
{
public:
    enum
    {
        subVolumeTotalNum = 8,
        fluxFacesTotalNum = 12,
        fluxFacesNumOnSubVolume = 3,
        fluxEdgesTotalNum = 6,
        edgesNumOnFluxFace = 2
    };

    static int getOldElemIdxFromNewFaceIdxto0(int zeroFaceIdx, int elementIdx)
    {
        return oldElemIdxFromNewFaceIdxto0_[zeroFaceIdx][elementIdx];
    }

    static int getNewElemIdxFromOldFaceIdxto0(int zeroFaceIdx, int elementIdx)
    {
        return newElemIdxFromOldFaceIdxto0_[zeroFaceIdx][elementIdx];
    }

    static int getOldFaceIdxFromNewIdxto0(int zeroFaceIdx, int fIdx)
    {
        return oldFaceIdxFromNewIdxto0_[zeroFaceIdx][fIdx];
    }

    static int getNewFaceIdxFromOldIdxto0(int zeroFaceIdx, int fIdx)
    {
        return newFaceIdxFromOldIdxto0_[zeroFaceIdx][fIdx];
    }

    static int getOldEdgeIdxFromNewFaceIdxto0(int zeroFaceIdx, int edgeIdx)
    {
        return oldEdgeIdxFromNewFaceIdxto0_[zeroFaceIdx][edgeIdx];
    }

    static int getNewEdgeIdxFromOldFaceIdxto0(int zeroFaceIdx, int edgeIdx)
    {
        return newEdgeIdxFromOldFaceIdxto0_[zeroFaceIdx][edgeIdx];
    }


private:
    static const int oldElemIdxFromNewFaceIdxto0_[fluxFacesTotalNum][subVolumeTotalNum];
    static const int newElemIdxFromOldFaceIdxto0_[fluxFacesTotalNum][subVolumeTotalNum];
    static const int oldFaceIdxFromNewIdxto0_[fluxFacesTotalNum][fluxFacesTotalNum];
    static const int newFaceIdxFromOldIdxto0_[fluxFacesTotalNum][fluxFacesTotalNum];
    static const int oldEdgeIdxFromNewFaceIdxto0_[fluxFacesTotalNum][fluxEdgesTotalNum];
    static const int newEdgeIdxFromOldFaceIdxto0_[fluxFacesTotalNum][fluxEdgesTotalNum];
};

const int IndexTranslatorAdaptive::oldElemIdxFromNewFaceIdxto0_[fluxFacesTotalNum][subVolumeTotalNum] =
{
        {0, 1, 2, 3, 4, 5, 6, 7},
        {1, 3, 0, 2, 5, 7, 4, 6},
        {3, 2, 1, 0, 7, 6, 5, 4},
        {2, 0, 3, 1, 6, 4, 7, 5},
        {5, 4, 7, 6, 1, 0, 3, 2},
        {7, 5, 6, 4, 3, 1, 2, 0},
        {6, 7, 4, 5, 2, 3, 0, 1},
        {4, 6, 5, 7, 0, 2, 1, 3},
        {0, 4, 1, 5, 2, 6, 3, 7},
        {1, 5, 3, 7, 0, 4, 2, 6},
        {3, 7, 2, 6, 1, 5, 0, 4},
        {2, 6, 0, 4, 3, 7, 1, 5}
};

const int IndexTranslatorAdaptive::newElemIdxFromOldFaceIdxto0_[fluxFacesTotalNum][subVolumeTotalNum] =
{
    {0, 1, 2, 3, 4, 5, 6, 7},
    {2, 0, 3, 1, 6, 4, 7, 5},
    {3, 2, 1, 0, 7, 6, 5, 4},
    {1, 3, 0, 2, 5, 7, 4, 6},
    {5, 4, 7, 6, 1, 0, 3, 2},
    {7, 5, 6, 4, 3, 1, 2, 0},
    {6, 7, 4, 5, 2, 3, 0, 1},
    {4, 6, 5, 7, 0, 2, 1, 3},
    {0, 2, 4, 6, 1, 3, 5, 7},
    {4, 0, 6, 2, 5, 1, 7, 3},
    {6, 4, 2, 0, 7, 5, 3, 1},
    {2, 6, 0, 4, 3, 7, 1, 5}
};

const int IndexTranslatorAdaptive::oldFaceIdxFromNewIdxto0_[fluxFacesTotalNum][fluxFacesTotalNum] =
{
        {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11},
        {1, 2, 3, 0, 5, 6, 7, 4, 9, 10, 11, 8},
        {2, 3, 0, 1, 6, 7, 4, 5, 10, 11, 8, 9},
        {3, 0, 1, 2, 7, 4, 5, 6, 11, 8, 9, 10},
        {4, 7, 6, 5, 0, 3, 2, 1, 9, 8, 11, 10},
        {5, 4, 7, 6, 1, 0, 3, 2, 10, 9, 8, 11},
        {6, 5, 4, 7, 2, 1, 0, 3, 11, 10, 9, 8},
        {7, 6, 5, 4, 3, 2, 1, 0, 8, 11, 10, 9},
        {8, 4, 9, 0, 11, 6, 10, 2, 3, 7, 5, 1},
        {9, 5, 10, 1, 8, 7, 11, 3, 0, 4, 6, 2},
        {10, 6, 11, 2, 9, 4, 8, 0, 1, 5, 7, 3},
        {11, 7, 8, 3, 10, 5, 9, 1, 2, 6, 4, 0}
};

const int IndexTranslatorAdaptive::newFaceIdxFromOldIdxto0_[fluxFacesTotalNum][fluxFacesTotalNum] =
{
    {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11},
    {3, 0, 1, 2, 7, 4, 5, 6, 11, 8, 9, 10},
    {2, 3, 0, 1, 6, 7, 4, 5, 10, 11, 8, 9},
    {1, 2, 3, 0, 5, 6, 7, 4, 9, 10, 11, 8},
    {4, 7, 6, 5, 0, 3, 2, 1, 9, 8, 11, 10},
    {5, 4, 7, 6, 1, 0, 3, 2, 10, 9, 8, 11},
    {6, 5, 4, 7, 2, 1, 0, 3, 11, 10, 9, 8},
    {7, 6, 5, 4, 3, 2, 1, 0, 8, 11, 10, 9},
    {3, 11, 7, 8, 1, 10, 5, 9, 0, 2, 6, 4},
    {8, 3, 11, 7, 9, 1, 10, 5, 4, 0, 2, 6},
    {7, 8, 3, 11, 5, 9, 1, 10, 6, 4, 0, 2},
    {11, 7, 8, 3, 10, 5, 9, 1, 2, 6, 4, 0}
};

const int IndexTranslatorAdaptive::oldEdgeIdxFromNewFaceIdxto0_[fluxFacesTotalNum][fluxEdgesTotalNum] =
{
        {0, 1, 2, 3, 4, 5},
        {0, 1, 3, 4, 5, 2},
        {0, 1, 4, 5, 2, 3},
        {0, 1, 5, 2, 3, 4},
        {1, 0, 2, 5, 4, 3},
        {1, 0, 3, 2, 5, 4},
        {1, 0, 4, 3, 2, 5},
        {1, 0, 5, 4, 3, 2},
        {2, 4, 5, 1, 3, 0},
        {3, 5, 2, 1, 4, 0},
        {4, 2, 3, 1, 5, 0},
        {5, 3, 4, 1, 2, 0}
};

const int IndexTranslatorAdaptive::newEdgeIdxFromOldFaceIdxto0_[fluxFacesTotalNum][fluxEdgesTotalNum] =
{
    {0, 1, 2, 3, 4, 5},
    {0, 1, 5, 2, 3, 4},
    {0, 1, 4, 5, 2, 3},
    {0, 1, 3, 4, 5, 2},
    {1, 0, 2, 5, 4, 3},
    {1, 0, 3, 2, 5, 4},
    {1, 0, 4, 3, 2, 5},
    {1, 0, 5, 4, 3, 2},
    {5, 3, 0, 4, 1, 2},
    {5, 3, 2, 0, 4, 1},
    {5, 3, 1, 2, 0, 4},
    {5, 3, 4, 1, 2, 0}
};
//! \endcond

//! \ingroup IMPET mpfa
/*! \brief Class including the information of a 3d interaction volume of an adaptive MPFA L-method that does not change with time.
 *
 * Includes information needed to calculate the transmissibility matrices of an L-interaction-volume.
 *
 */
template<class TypeTag>
class FvMpfaL3dInteractionVolumeAdaptive:public FvMpfaL3dInteractionVolume<TypeTag>
{
private:
    using ParentType = FvMpfaL3dInteractionVolume<TypeTag>;
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using Grid = GetPropType<TypeTag, Properties::Grid>;
    using GridView = typename GetPropType<TypeTag, Properties::GridGeometry>::GridView;
    using Problem = GetPropType<TypeTag, Properties::Problem>;

    enum
    {
        dim = GridView::dimension,
        dimWorld = GridView::dimensionworld,
    };

    using Element = typename GridView::template Codim<0>::Entity;
    using ElementSeed = typename Grid::template Codim<0>::EntitySeed;

    using BoundaryTypes = GetPropType<TypeTag, Properties::SequentialBoundaryTypes>;
    using SolutionTypes = GetProp<TypeTag, Properties::SolutionTypes>;
    using PrimaryVariables = typename SolutionTypes::PrimaryVariables;

    using DimVector = Dune::FieldVector<Scalar, dim>;
    using FieldVectorVector = Dune::FieldVector<DimVector, dim>;
    using FieldVectorVector2 = Dune::FieldVector<DimVector, 2>;
    using FieldVectorVectorVector = Dune::FieldVector<FieldVectorVector2, dim>;
    using IndexVector = Dune::FieldVector<int, dim>;
    using BCTypeVector = std::vector<BoundaryTypes>;
    using BCVector = std::vector<PrimaryVariables>;

public:

    //! \cond \private
    enum FaceTypes
    {
        inside = 1,
        boundary = 0,
        outside = -1,
    };

    //!\copydoc
    enum
    {
        subVolumeTotalNum = IndexTranslatorAdaptive::subVolumeTotalNum,
        fluxFacesTotalNum = IndexTranslatorAdaptive::fluxFacesTotalNum,
        fluxEdgesTotalNum = IndexTranslatorAdaptive::fluxEdgesTotalNum
    };
    //! \endcond

    //! The different hanging node interaction volume types (see dissertation M. Wolff, http://elib.uni-stuttgart.de/opus/volltexte/2013/8661/)
    enum HangingNodeTypes
    {
        noHangingNode = -1,//!< regular interaction volume
        twoSmallCells = 0,//!< hanging-node interaction volume of type 5 or 7
        fourSmallCellsFace = 1,//!<  hanging-node interaction volume of type 1
        fourSmallCellsEdge = 2,//!<  hanging-node interaction volume of type 3
        fourSmallCellsDiag = 3,//!<  hanging-node interaction volume of type 4
        sixSmallCells = 4//!<  hanging-node interaction volume of type 2 or 6
    };

    //! Constructs a FvMpfaL3dInteractionVolumeAdaptive object
    /**
     */
    FvMpfaL3dInteractionVolumeAdaptive(const Grid& grid)
    : ParentType(grid)
    , hangingNodeType_(noHangingNode)
    {}

    //!\copydoc FvMpfaL3dInteractionVolume::reset()
    void reset()
    {
        hangingNodeType_ = noHangingNode;
        existingLevel_.clear();
    }

    //!\copydoc FvMpfaL3dInteractionVolume::setSubVolumeElement()
    void setSubVolumeElement(const Element& element, int subVolumeIdx)
    {
        ParentType::setSubVolumeElement(element, subVolumeIdx);
        existingLevel_.insert(element.level());
    }

    //! Store the type of hanging-node-interaction volume
    /*!
     * \param hNType the type of hanging-node-interaction volume of type FvMpfaL3dInteractionVolumeAdaptive<TypeTag>::HangingNodeTypes
     */
    void setHangingNodeType(int hNType)
    {
        hangingNodeType_ = hNType;
    }

    //! Check if elements in the interaction volume are of the same grid level
    /*!
     * \return <tt>true</tt> if all elements are on the same grid level
     */
    bool sameLevel()
    {
        if (!isHangingNodeVolume())
            return existingLevel_.size() < 2;
        else
        {
            return existingLevel_.size() < 3;
        }
    }

    //! Check if an element of a certain grid level is stored
    /*!
     * \param level the dune grid level
     *
     * \return <tt>true</tt> if an element of a certain grid level is stored
     */
    bool hasLevel(int level)
    {
        return  existingLevel_.find(level) != existingLevel_.end();
    }

    //! Check whether the interaction volume is a hanging-node volume
    /*!
     * \return <tt>true</tt> if the interaction volume is a hanging-node volume
     */
    bool isHangingNodeVolume()
    {
        return hangingNodeType_ != noHangingNode;
    }

    //! The type of the interaction volume as type of FvMpfaL3dInteractionVolumeAdaptive<TypeTag>::HangingNodeTypes
    int getHangingNodeType()
    {
        return hangingNodeType_;
    }

    //!\copydoc FvMpfaL3dInteractionVolume::printInteractionVolumeInfo()
    void printInteractionVolumeInfo()
    {
        ParentType::printInteractionVolumeInfo();

            if (isHangingNodeVolume())
                std::cout<<"hanging node type: "<<hangingNodeType_<<"\n";
    }

//    void printInteractionVolumeInfoToFile(std::ofstream& dataFile)
//    {
//        ParentType::printInteractionVolumeInfoToFile(dataFile);
//
//            if (isHangingNodeVolume())
//                dataFile<<"hanging node type: "<<hangingNodeType_<<"\n";
//    }

private:
    int hangingNodeType_;
    std::set<int> existingLevel_;
};
}
#endif
