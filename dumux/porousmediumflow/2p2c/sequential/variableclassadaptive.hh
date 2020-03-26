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
 * \brief Base class holding the variables for sequential models.
 */
#ifndef DUMUX_VARIABLECLASS2P2C_ADAPTIVE_HH
#define DUMUX_VARIABLECLASS2P2C_ADAPTIVE_HH

#include <dune/common/power.hh>

// for  parallelization
#include <dumux/porousmediumflow/sequential/variableclassadaptive.hh>

namespace Dumux {
/*!
 * \ingroup SequentialTwoPTwoCModel
 * \brief Class holding additionally mpfa data of adaptive compositional models.
 *
 * This class provides the possibility to store and load the transmissibilities (and associated infos)
 * of the mpfa method per irregular face. This class provides the storage container and access methods
 * for both 2D and 3D implementation.
 * Note that according to the number of half-edges (sub-faces) regarded, one ore more transmissibility
 * can be stored and loaded.
 *
 * \tparam TypeTag The Type Tag
 */
template<class TypeTag>
class VariableClass2P2CAdaptive: public VariableClassAdaptive<TypeTag>
{
private:
    using ParentType = VariableClassAdaptive<TypeTag>;
    using BaseType = VariableClass<TypeTag>;

    using Problem = GetPropType<TypeTag, Properties::Problem>;
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using GridView = typename GetPropType<TypeTag, Properties::GridGeometry>::GridView;
    using CellData = GetPropType<TypeTag, Properties::CellData>;
    using AdaptedValues = typename CellData::AdaptedValues;

    using Grid = typename GridView::Grid;
    using IdType = typename Grid::LocalIdSet::IdType;
    using IntersectionIterator = typename GridView::IntersectionIterator;
    using Intersection = typename GridView::Intersection;
    enum
    {
        dim = GridView::dimension
    };
    enum    //!< for first and second half edge (2D) or subvolume face (3D)
    {
        first = 0,
        second = 1,
        diagonal1 = 2,
        diagonal2 = 3
    };
    // convenience shortcuts for Vectors/Matrices
    using GlobalPosition = Dune::FieldVector<Scalar, GridView::dimensionworld>;
    using TransmissivityMatrix = Dune::FieldVector<Scalar,dim+1>;

protected:
    /** in the 2D case, we need to store 1 additional cell. In 3D, we store
    * dim-1=2 cells for one interaction volume, and 8 cells if four interaction volumes
    * are regarded.
    */
    const static int storageRequirement = Dune::StaticPower<dim-1, dim>::power;

    /*!
     * \brief Storage object for data related to the MPFA method
     *
     * This Object stores the transmissibility coefficients
     * for all interaction regions/volumes of an irregular interface
     * in an h-adaptive simulation. It is valid for 2D/3D
     */
    struct mpfaData
    {
        TransmissivityMatrix T1_[(dim-1)*2];
        int globalIdx3_[storageRequirement];
        GlobalPosition globalPos3_[storageRequirement];
        std::vector<IntersectionIterator> secondHalfEdgeIntersection_;
        int interactionRegionsStored;

        //! Constructor for the local storage object of mpfa data
        mpfaData()
        {
            interactionRegionsStored = 0;
        }

        /*!
         * \brief Stores an intersection for the 2D implementation
         *
         * This method also provides the information that both half-edges (2D) are
         * regarded and information was stored: Two interaction regions are applied
         * \param is23 Intersection pointing to 3rd cell of additional interaction region
         */
        void setIntersection(IntersectionIterator& is23)
        {
            secondHalfEdgeIntersection_.push_back(is23);
            interactionRegionsStored = 2;
        }
        //! Acess method to the stored intersection (only 2D)
        const IntersectionIterator& getIntersection()
        {
            return secondHalfEdgeIntersection_[0];
        }
    };
    std::map<IdType, mpfaData> irregularInterfaceMap_; //!< Storage container for mpfa information
    const Grid& grid_; //!< The grid

public:
    /*!
     * \brief Constructs a VariableClass object
     *
     *  @param gridView a DUNE gridview object corresponding to diffusion and transport equation
     */
    VariableClass2P2CAdaptive(const GridView& gridView) :
        ParentType(gridView), grid_(gridView.grid())
    {}

    /*!
     * \brief Resizes sequential variable vectors
     *
     * Method that change the size of the vectors for h-adaptive simulations, and clears recently
     * stored transmissibility matrices for the newly adapted grid.
     *
     * \param size Size of the current (refined and coarsened) grid
     */
    void adaptVariableSize(int size)
    {
        BaseType::adaptVariableSize(size);
        // clear mapper
        irregularInterfaceMap_.clear();
    }
    /*!
     * \brief Reconstruct missing primary variables (where elements are created/deleted)
     *
     * To reconstruct the solution in father elements, problem properties might
     * need to be accessed.
     *
     * \param problem The current problem
     */
    void reconstructPrimVars(Problem& problem)
    {
        ParentType::reconstructPrimVars(problem);

        problem.pressureModel().adaptPressure();
    }

    /*!
     * \brief Stores Mpfa Data of one interaction region on an intersection
     *
     * The method stores information of ONE interaction region (Transmissitivity
     * as well as details of the 3rd cell in the region) into a storage container.
     * The key to each element is the index of the intersection, seen from the smaller
     * cell (only this is unique).
     * If we arrive from the "wrong" (i.e. non-unique) direction, we invert fluxes.
     *
     * \param irregularIs The current irregular intersection
     * \param T1 Transmissitivity matrix for flux calculations
     * \param globalPos3 The position of the 3rd cell of the interaction region
     * \param globalIdx3 The index of the 3rd cell of the interaction region
     */
    void storeMpfaData(const typename GridView::Intersection & irregularIs,
                       const TransmissivityMatrix& T1,
                       const GlobalPosition& globalPos3,
                       const int& globalIdx3)
    {
        IdType intersectionID = grid_.localIdSet().subId(
            irregularIs.inside(), irregularIs.indexInInside(), 1);
        // mapping is only unique from smaller cell (if *inside and not *outside)
        if (irregularIs.inside().level() < irregularIs.outside().level())
        {
            // IS is regarded from larger cell: get the unique number as seen from smaller
            intersectionID = grid_.localIdSet().subId(
                irregularIs.outside(), irregularIs.indexInOutside(), 1);

            // store as if it was seen from smaller: change i & j
            irregularInterfaceMap_[intersectionID].T1_[first][2] = - T1[0];
            irregularInterfaceMap_[intersectionID].T1_[first][1] = - T1[1];
            irregularInterfaceMap_[intersectionID].T1_[first][0] = - T1[2];
        }
        else // we look from smaller cell = unique interface
            // proceed with numbering according to Aavatsmark, seen from cell i
            irregularInterfaceMap_[intersectionID].T1_[first] = T1;

        irregularInterfaceMap_[intersectionID].globalPos3_[0] = globalPos3;
        irregularInterfaceMap_[intersectionID].globalIdx3_[0] = globalIdx3;
    }

    /*!
     * \brief Stores Mpfa Data on an intersection for both half-edges
     *
     * The method stores information of both interaction regions (Transmissitivity
     * as well as details of the 3rd cells of both regions) into a storage container.
     * The key to each element is the index of the intersection, seen from the smaller
     * cell (only this is unique).
     * If we arrive from the "wrong" (i.e. non-unique) direction, we invert fluxes.
     *
     * \param irregularIs The current irregular intersection
     * \param secondHalfEdgeIntersectionIt Iterator to the intersection connecting the second interaction region
     * \param T1 Transmissitivity matrix for flux calculations: unique interaction region
     * \param T1_secondHalfEdge Second transmissitivity matrix for flux calculations for non-unique interaction region
     * \param globalPos3 The position of the 3rd cell of the first interaction region
     * \param globalIdx3 The index of the 3rd cell of the first interaction region
     */
    void storeMpfaData(const typename GridView::Intersection & irregularIs,
                        IntersectionIterator& secondHalfEdgeIntersectionIt,
                        const TransmissivityMatrix& T1,
                        const TransmissivityMatrix& T1_secondHalfEdge,
                        const GlobalPosition& globalPos3,
                        const int& globalIdx3)
    {
        IdType intersectionID
                = grid_.localIdSet().subId(irregularIs.inside(),
                                           irregularIs.indexInInside(), 1);

        // mapping is only unique from smaller cell (if *inside and not *outside)
        if (irregularIs.inside().level() < irregularIs.outside().level())
        {
            // IS is regarded from larger cell: get the unique number as seen from smaller
            intersectionID
                = grid_.localIdSet().subId(irregularIs.outside(),
                                           irregularIs.indexInOutside(), 1);

            // store as if it was seen from smaller: change i & j
            irregularInterfaceMap_[intersectionID].T1_[first][2] = - T1[0];
            irregularInterfaceMap_[intersectionID].T1_[first][1] = - T1[1];
            irregularInterfaceMap_[intersectionID].T1_[first][0] = - T1[2];

            irregularInterfaceMap_[intersectionID].T1_[second][2] = - T1_secondHalfEdge[0];
            irregularInterfaceMap_[intersectionID].T1_[second][1] = - T1_secondHalfEdge[1];
            irregularInterfaceMap_[intersectionID].T1_[second][0] = - T1_secondHalfEdge[2];

        }
        else // we look from smaller cell = unique interface
        {
            // proceed with numbering according to Aavatsmark, seen from cell i
            irregularInterfaceMap_[intersectionID].T1_[first] = T1;
            irregularInterfaceMap_[intersectionID].T1_[second] = T1_secondHalfEdge;
        }

        irregularInterfaceMap_[intersectionID].globalPos3_[0] = globalPos3;
        irregularInterfaceMap_[intersectionID].globalIdx3_[0] = globalIdx3;
        // second half edge
        irregularInterfaceMap_[intersectionID].setIntersection(secondHalfEdgeIntersectionIt);
    }

    /*!
     * \brief Stores 3D Mpfa Data on an intersection
     *
     * The method stores information of the interaction region (Transmissitivity
     * as well as details of the 3rd & 4th cell in the region) into a storage container.
     * The key to each element is the index of the intersection, seen from the smaller
     * cell (only this index is unique).
     * If we arrive from the "wrong" (i.e. non-unique) direction, we invert fluxes.
     *
     * \param irregularIs The current irregular intersection
     * \param T1 Transmissitivity matrix for flux calculations
     * \param globalPos3 The position of the 3rd cell of the interaction region
     * \param globalIdx3 The index of the 3rd cell of the interaction region
     * \param globalPos4 The position of the 4th cell of the interaction region
     * \param globalIdx4 The index of the 4th cell of the interaction region
     * \param subFaceIdx The index of the subface (up to 4 subfaces possible in 3D)
     */
    void storeMpfaData3D(const typename GridView::Intersection & irregularIs,
                        const TransmissivityMatrix& T1,
                        const GlobalPosition& globalPos3,
                        const int& globalIdx3,
                        const GlobalPosition& globalPos4,
                        const int& globalIdx4,
                        int subFaceIdx = 0)
    {
        // if the second interaction Region (subfaceIdx=1) needs to be stored,
        // we need an offset for the second globalPos and globalIdxes
        const int offset = subFaceIdx * 2;

        IdType intersectionID
                = grid_.localIdSet().subId(irregularIs.inside(),
                                           irregularIs.indexInInside(), 1);

        // mapping is only unique from smaller cell (if *inside and not *outside)
        if (irregularIs.inside().level() < irregularIs.outside().level())
        {
            // IS is regarded from larger cell: get the unique ID as seen from smaller
            intersectionID
                = grid_.localIdSet().subId(irregularIs.outside(),
                                           irregularIs.indexInOutside(), 1);

            // store as if it was seen from smaller: change i & j
            irregularInterfaceMap_[intersectionID].T1_[subFaceIdx][0] = -T1[1];
            irregularInterfaceMap_[intersectionID].T1_[subFaceIdx][1] = -T1[0];
            irregularInterfaceMap_[intersectionID].T1_[subFaceIdx][2] = -T1[2];
            irregularInterfaceMap_[intersectionID].T1_[subFaceIdx][3] = -T1[3];
        }
        else
        {
            // proceed with numbering according to case2
            irregularInterfaceMap_[intersectionID].T1_[subFaceIdx][0] = T1[0];
            irregularInterfaceMap_[intersectionID].T1_[subFaceIdx][1] = T1[1];
            irregularInterfaceMap_[intersectionID].T1_[subFaceIdx][2] = T1[2];
            irregularInterfaceMap_[intersectionID].T1_[subFaceIdx][3] = T1[3];
        }
        irregularInterfaceMap_[intersectionID].globalPos3_[offset+0] = globalPos3;
        irregularInterfaceMap_[intersectionID].globalIdx3_[offset+0] = globalIdx3;
        irregularInterfaceMap_[intersectionID].globalPos3_[offset+1] = globalPos4;
        irregularInterfaceMap_[intersectionID].globalIdx3_[offset+1] = globalIdx4;

        using std::max;
        irregularInterfaceMap_[intersectionID].interactionRegionsStored
            = max(irregularInterfaceMap_[intersectionID].interactionRegionsStored, subFaceIdx+1);
    }

    /*!
     * \brief Weigths the transmissivity coefficient by the flux area (3D)
     *
     * Each cell face could contain up to 4 interaction regions in 3D. If fewer interaction
     * regions are considered (e.g. at a boundary), the flux goes through a higher area than
     * was regarded in the calculation of the interaction region. Therefore, the transmissibility
     * coefficients have to be weighted (scaled).
     * If no specific face Index is given, all present Transmissibility matrices are weighted
     * homogeneously by the given weight.
     * \param irregularIs The current irregular intersection
     * \param weight The weighting factor
     * \param subFaceIdx The index of the subface (up to 4 subfaces possible in 3D)
     */
    void performTransmissitivityWeighting(const typename GridView::Intersection & irregularIs,
                       Scalar weight, int subFaceIdx = -1)
    {
        IdType intersectionID
                = grid_.localIdSet().subId(irregularIs.inside(),
                                           irregularIs.indexInInside(), 1);

        // mapping is only unique from smaller cell (if *inside and not *outside)
        if (irregularIs.inside().level() < irregularIs.outside().level())
        {
            // IS is regarded from larger cell: get the unique ID as seen from smaller
            intersectionID
                = grid_.localIdSet().subId(irregularIs.outside(),
                                           irregularIs.indexInOutside(), 1);
        }

        // for subFaceIdx == -1, we weight all subfaces equally
        if(subFaceIdx == -1)
        {
            for(int i = 0; i < irregularInterfaceMap_[intersectionID].interactionRegionsStored; i++)
                irregularInterfaceMap_[intersectionID].T1_[i] *= weight;
        }
        else    //weight specifically
            irregularInterfaceMap_[intersectionID].T1_[subFaceIdx] *= weight;
    }

    /*!
     * \brief Provides access to stored Mpfa Data on an intersection for both half-edges
     *
     * The method gets information of both interaction regions (Transmissitivity
     * as well as details of the 3rd cells of both regions) from a storage container.
     * The key to each element is the index of the intersection, seen from the smaller
     * cell (only this is unique).
     * If we arrive from the "wrong" (i.e. non-unique) direction, we invert fluxes.
     *
     * \param irregularIs The current irregular intersection
     * \param secondHalfEdgeIntersectionIt Iterator to the intersection connecting the second interaction region
     * \param T1 Transmissitivity matrix for flux calculations: unique interaction region
     * \param T1_secondHalfEdge Second transmissitivity matrix for flux calculations for non-unique interaction region
     * \param globalPos3 The position of the 3rd cell of the first interaction region
     * \param globalIdx3 The index of the 3rd cell of the first interaction region
     */
    int getMpfaData(const Intersection& irregularIs,
                        IntersectionIterator& secondHalfEdgeIntersectionIt,
                        TransmissivityMatrix& T1,
                        TransmissivityMatrix& T1_secondHalfEdge,
                        GlobalPosition& globalPos3,
                        int& globalIdx3)
    {
        IdType intersectionID
                = grid_.localIdSet().subId(irregularIs.inside(),
                                           irregularIs.indexInInside(), 1);
        // mapping is only unique from smaller cell (if *inside and not *outside)
        if (irregularIs.inside().level() < irregularIs.outside().level())
        {
            // IS is regarded from larger cell: get the unique number as seen from smaller
            intersectionID
                = grid_.localIdSet().subId(irregularIs.outside(),
                                           irregularIs.indexInOutside(), 1);

            // check if T1ransmissibility matrix was stored for that IF
            if (irregularInterfaceMap_.find(intersectionID) == irregularInterfaceMap_.end())
                return 0; // no stored data!

            // If data is stored, it is so as if the IF is regarded from smaller cell.
            // since we are looking from larger cell, cells i & j have to be changed
            // Additionally, flux points in opposite direction: - sign
            T1[0] = -irregularInterfaceMap_[intersectionID].T1_[first][2];
            T1[1] = -irregularInterfaceMap_[intersectionID].T1_[first][1];
            T1[2] = -irregularInterfaceMap_[intersectionID].T1_[first][0];
            globalPos3 = irregularInterfaceMap_[intersectionID].globalPos3_[0];
            globalIdx3 = irregularInterfaceMap_[intersectionID].globalIdx3_[0];
            //second half edge
            if(irregularInterfaceMap_[intersectionID].interactionRegionsStored == 2)
            {
                secondHalfEdgeIntersectionIt = irregularInterfaceMap_[intersectionID].getIntersection();
                T1_secondHalfEdge[0] = -irregularInterfaceMap_[intersectionID].T1_[second][2];
                T1_secondHalfEdge[1] = -irregularInterfaceMap_[intersectionID].T1_[second][1];
                T1_secondHalfEdge[2] = -irregularInterfaceMap_[intersectionID].T1_[second][0];
                return 2;
            }
            return 1;
        }
        // check if T1ransmissibility matrix was stored for that IF
        if (irregularInterfaceMap_.find(intersectionID) == irregularInterfaceMap_.end())
            return 0; // no stored data!

        T1[0] = irregularInterfaceMap_[intersectionID].T1_[first][0];
        T1[1] = irregularInterfaceMap_[intersectionID].T1_[first][1];
        T1[2] = irregularInterfaceMap_[intersectionID].T1_[first][2];
        globalPos3 = irregularInterfaceMap_[intersectionID].globalPos3_[0];
        globalIdx3 = irregularInterfaceMap_[intersectionID].globalIdx3_[0];
        //second half edge
        if(irregularInterfaceMap_[intersectionID].interactionRegionsStored == 2)
        {
            secondHalfEdgeIntersectionIt = irregularInterfaceMap_[intersectionID].getIntersection();
            T1_secondHalfEdge = irregularInterfaceMap_[intersectionID].T1_[second];
            return 2;
        }
        return 1;
    }

    /*!
     * \brief Provides access to stored 3D Mpfa Data on an intersection for up to 4 subfaces
     *
     * The method gets information of up to 4 interaction regions (Transmissitivity
     * as well as details of the 3rd & 4th cells of the regions) from a storage container.
     * The key to each element is the index of the intersection, seen from the smaller
     * cell (only this is unique).
     * If we arrive from the "wrong" (i.e. non-unique) direction, we invert fluxes.
     *
     * \param irregularIs The current irregular intersection
     * \param T1 Transmissitivity matrix for flux calculations: unique interaction region
     * \param globalPos3 The position of the 3rd cell of the first interaction region
     * \param globalIdx3 The index of the 3rd cell of the first interaction region
     * \param globalPos4 The position of the 4th cell of the interaction region
     * \param globalIdx4 The index of the 4th cell of the interaction region
     * \param subFaceIdx The index of the subface (up to 4 subfaces possible in 3D)
     */
    int getMpfaData3D(const Intersection& irregularIs,
                        TransmissivityMatrix& T1,
                        GlobalPosition& globalPos3,
                        int& globalIdx3,
                        GlobalPosition& globalPos4,
                        int& globalIdx4,
                        int subFaceIdx = 0)
    {
        // if the second interaction Region (subfaceIdx=1) needs to be stored,
        // we need an offset for the second globalPos and globalIdxes
        const int offset = subFaceIdx * 2;

        IdType intersectionID
                = grid_.localIdSet().subId(irregularIs.inside(),
                                           irregularIs.indexInInside(), 1);
        // mapping is only unique from smaller cell (if *inside and not *outside)
        if (irregularIs.inside().level() < irregularIs.outside().level())
        {
            // IS is regarded from larger cell: get the unique number as seen from smaller
            intersectionID
                = grid_.localIdSet().subId(irregularIs.outside(),
                                           irregularIs.indexInOutside(), 1);

            // check if T1ransmissibility matrix was stored for that IF
            if (irregularInterfaceMap_.find(intersectionID) == irregularInterfaceMap_.end())
                return 0; // no stored data!

            // If data is stored, it is so as if the IF is regarded from smaller cell.
            // since we are looking from larger cell, cells i & j have to be changed
            // Additionally, flux points in opposite direction: - sign
            T1[0] = -irregularInterfaceMap_[intersectionID].T1_[subFaceIdx][1];
            T1[1] = -irregularInterfaceMap_[intersectionID].T1_[subFaceIdx][0];
            T1[2] = -irregularInterfaceMap_[intersectionID].T1_[subFaceIdx][2];
            T1[3] = -irregularInterfaceMap_[intersectionID].T1_[subFaceIdx][3];

        }
        else
        {
            // check if T1ransmissibility matrix was stored for that IF
            if (irregularInterfaceMap_.find(intersectionID) == irregularInterfaceMap_.end())
                return 0; // no stored data!

            T1[0] = irregularInterfaceMap_[intersectionID].T1_[subFaceIdx][0];
            T1[1] = irregularInterfaceMap_[intersectionID].T1_[subFaceIdx][1];
            T1[2] = irregularInterfaceMap_[intersectionID].T1_[subFaceIdx][2];
            T1[3] = irregularInterfaceMap_[intersectionID].T1_[subFaceIdx][3];
        }

        // return what does not depend on direction: additional cells
        globalPos3 = irregularInterfaceMap_[intersectionID].globalPos3_[offset+0];
        globalIdx3 = irregularInterfaceMap_[intersectionID].globalIdx3_[offset+0];
        globalPos4 = irregularInterfaceMap_[intersectionID].globalPos3_[offset+1];
        globalIdx4 = irregularInterfaceMap_[intersectionID].globalIdx3_[offset+1];

        return irregularInterfaceMap_[intersectionID].interactionRegionsStored;
    }
    //@}

};
} // end namespace Dumux
#endif
