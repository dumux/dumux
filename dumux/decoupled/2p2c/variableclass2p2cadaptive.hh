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
#ifndef DUMUX_VARIABLECLASS2P2C_ADAPTIVE_HH
#define DUMUX_VARIABLECLASS2P2C_ADAPTIVE_HH


// for  parallelization
#include <dumux/decoupled/common/variableclassadaptive.hh>

/**
 * @file
 * @brief  Base class holding the variables for sequential models.
 * @author Markus Wolff
 */

namespace Dumux
{
/*!
 * \ingroup Adaptive2p2c
 */
//! Class holding additionally mpfa data of adaptive compositional models.
/*!
 * This class provides the possibility to store and load the transmissibilities (and associated infos)
 * of the mpfa method per irregular face. While this class provides the storage container for both 2D
 * and 3D implementation, only the 2D storage methods are included here.
 * Note that according to the number of half-edges (sub-faces) regarded, one ore two transmissibility
 * can be stored and loaded.
 *
 * @tparam TypeTag The Type Tag
 */
template<class TypeTag>
class VariableClass2P2CAdaptive: public VariableClassAdaptive<TypeTag>
{
private:
    typedef VariableClassAdaptive<TypeTag> ParentType;
    typedef VariableClass<TypeTag> BaseType;

    typedef typename GET_PROP_TYPE(TypeTag, Problem) Problem;
    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
    typedef typename GET_PROP_TYPE(TypeTag, GridView) GridView;
    typedef typename GET_PROP_TYPE(TypeTag, CellData) CellData;
    typedef typename CellData::AdaptedValues AdaptedValues;

    typedef typename GridView::Grid Grid;
    typedef typename Grid::LocalIdSet::IdType IdType;
    typedef typename GridView::IntersectionIterator IntersectionIterator;
    typedef typename GridView::Intersection Intersection;
    enum
    {
        dim = GridView::dimension, dimWorld = GridView::dimensionworld
    };
    enum    //!< for first and second half edge (2D) or subvolume face (3D)
    {
        first = 0, second = 1
    };
    // convenience shortcuts for Vectors/Matrices
    typedef Dune::FieldVector<Scalar, GridView::dimensionworld> GlobalPosition;
    typedef Dune::FieldVector<Scalar,dim+1> TransmissivityMatrix;

protected:
    /** in the 2D case, we need to store 1 additional cell. In 3D, we store
    * dim-1=2 cells for one interaction volume, and 4 if two interaction volumes
    * are regarded.
    */
    const static int storageRequirement = (dim-1)*(dim-1);
    //! Storage object for data related to the MPFA method
    /**
     * This Struct stores the transmissibility coefficients
     * for the two half-eges of an irregular interface (one
     * near a hanging node) in an h-adaptive simulation.
     */
    struct mpfaData
    {
        TransmissivityMatrix T1_[2];
        int globalIdx3_[storageRequirement];
        GlobalPosition globalPos3_[storageRequirement];
        std::vector<IntersectionIterator> secondHalfEdgeIntersection_;
        bool hasSecondHalfEdge;

        //! Constructor for the local storage object of mpfa data
        mpfaData()
        {
            hasSecondHalfEdge = false;
        }
        //! stores an intersection
        /** This also provides the information that both half-edges are
         * regarded and information was stored.
         * \param is23 Intersection pointing to 3rd cell of additional interaction region
         */
        void setIntersection(IntersectionIterator& is23)
        {
            secondHalfEdgeIntersection_.push_back(is23);
            hasSecondHalfEdge = true;
        };
        //! Acess method to the stored intersection
        const IntersectionIterator& getIntersection()
        {
            return secondHalfEdgeIntersection_[0];
        };
    };
    std::map<IdType, mpfaData> irregularInterfaceMap_; //!< Storage container for mpfa information
    const Grid& grid_; //!< The grid

public:
    //! Constructs a VariableClass object
    /**
     *  @param gridView a DUNE gridview object corresponding to diffusion and transport equation
     */
    VariableClass2P2CAdaptive(const GridView& gridView) :
        ParentType(gridView), grid_(gridView.grid())
    {}

    //! Resizes decoupled variable vectors
    /*! Method that change the size of the vectors for h-adaptive simulations, and clears recently
     * stored transmissibility matrices for the newly adapted grid.
     *
     *\param size Size of the current (refined and coarsened) grid
     */
    void adaptVariableSize(int size)
    {
        BaseType::adaptVariableSize(size);
        // clear mapper
        irregularInterfaceMap_.clear();
    }
    /*!
     * Reconstruct missing primary variables (where elements are created/deleted)
     *
     * To reconstruct the solution in father elements, problem properties might
     * need to be accessed.
     *
     * @param problem The current problem
     */
    void reconstructPrimVars(Problem& problem)
    {
        ParentType::reconstructPrimVars(problem);

        problem.pressureModel().adaptPressure();
    }

    //! Stores Mpfa Data on an intersection
    /** The method stores information to the interaction region (Transmissitivity
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
    	IdType intersectionID
				= grid_.localIdSet().subId(*irregularIs.inside(),
											irregularIs.indexInInside(), 1);
    	// mapping is only unique from smaller cell (if *inside and not *outside)
    	if (irregularIs.inside().level() < irregularIs.outside().level())
		{
    		// IS is regarded from larger cell: get the unique number as seen from smaller
    		intersectionID
				= grid_.localIdSet().subId(*irregularIs.outside(),
											irregularIs.indexInOutside(), 1);

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
        return;
    }
    //! Stores Mpfa Data on an intersection for both half-edges
    /** The method stores information of both interaction regions (Transmissitivity
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
                = grid_.localIdSet().subId(*irregularIs.inside(),
                                            irregularIs.indexInInside(), 1);

        // mapping is only unique from smaller cell (if *inside and not *outside)
        if (irregularIs.inside().level() < irregularIs.outside().level())
        {
            // IS is regarded from larger cell: get the unique number as seen from smaller
            intersectionID
                = grid_.localIdSet().subId(*irregularIs.outside(),
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


        return;
    }

    //! Provides acess to stored Mpfa Data on an intersection for both half-edges
    /** The method gets information of both interaction regions (Transmissitivity
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
    int getMpfaData(const Intersection& irregularIs,
                        IntersectionIterator& secondHalfEdgeIntersectionIt,
                        TransmissivityMatrix& T1,
                        TransmissivityMatrix& T1_secondHalfEdge,
                        GlobalPosition& globalPos3,
                        int& globalIdx3)
    {
        IdType intersectionID
                = grid_.localIdSet().subId(*irregularIs.inside(),
                                            irregularIs.indexInInside(), 1);
        // mapping is only unique from smaller cell (if *inside and not *outside)
        if (irregularIs.inside().level() < irregularIs.outside().level())
        {
            // IS is regarded from larger cell: get the unique number as seen from smaller
            intersectionID
                = grid_.localIdSet().subId(*irregularIs.outside(),
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
            if(irregularInterfaceMap_[intersectionID].hasSecondHalfEdge)
            {
                secondHalfEdgeIntersectionIt = irregularInterfaceMap_[intersectionID].getIntersection();
                T1_secondHalfEdge[0] = -irregularInterfaceMap_[intersectionID].T1_[second][2];
                T1_secondHalfEdge[1] = -irregularInterfaceMap_[intersectionID].T1_[second][1];
                T1_secondHalfEdge[2] = -irregularInterfaceMap_[intersectionID].T1_[second][0];
    //          Dune::dinfo << "mpfa Info retrieved for isID " << intersectionID
    //                  << "at coordinate " << irregularIs.geometry().center() << " from GlobalIdx " << this->index(*irregularIs.inside())<<std::endl;
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
        if(irregularInterfaceMap_[intersectionID].hasSecondHalfEdge)
        {
            secondHalfEdgeIntersectionIt = irregularInterfaceMap_[intersectionID].getIntersection();
            T1_secondHalfEdge = irregularInterfaceMap_[intersectionID].T1_[second];

    //      Dune::dinfo << "mpfa Info retrieved for isID " << intersectionID
    //              << "at coordinate " << irregularIs.geometry().center() << " from GlobalIdx " << this->index(*irregularIs.inside())<<std::endl;
            return 2;
        }
        return 1;
    }
	//@}

};
}
#endif
