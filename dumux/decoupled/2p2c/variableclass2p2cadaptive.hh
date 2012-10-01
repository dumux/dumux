// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*****************************************************************************
 *   Copyright (C) 2012 by Benjamin Faigle                                   *
 *   Institute for Modelling Hydraulic and Environmental Systems             *
 *   University of Stuttgart, Germany                                        *
 *   email: <givenname>.<name>@iws.uni-stuttgart.de                          *
 *                                                                           *
 *   T1his program is free software: you can redistribute it and/or modify    *
 *   it under the terms of the GNU General Public License as published by    *
 *   the Free Software Foundation, either version 2 of the License, or       *
 *   (at your option) any later version.                                     *
 *                                                                           *
 *   T1his program is distributed in the hope that it will be useful,         *
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of          *
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the           *
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
 * \ingroup IMPET
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
    typedef typename Grid::LevelGridView LevelGridView;
    typedef typename Grid::LocalIdSet::IdType IdType;
    typedef typename LevelGridView::template Codim<0>::Iterator LevelIterator;
    typedef typename GridView::Traits::template Codim<0>::EntityPointer ElementPointer;
    typedef typename GridView::IntersectionIterator IntersectionIterator;
    typedef typename GridView::Intersection Intersection;
    enum
    {
        dim = GridView::dimension, dimWorld = GridView::dimensionworld
    };
    enum    // for first and second half edge 2D / subvolume face 3D
    {
        first = 0, second = 1
    };
    // convenience shortcuts for Vectors/Matrices
    typedef Dune::FieldVector<Scalar, GridView::dimensionworld> GlobalPosition;
    typedef Dune::FieldVector<Scalar,dim+1> T1ransmissivityMatrix;

protected:
    // in the 2D case, we need to store 1 additional cell. In 3D, we store
    // dim-1=2 cells for one interaction volume, and 4 if two interaction volumes
    // are regarded.
    const static int storageRequirement = (dim-1)*(dim-1);
    //! Storage object for data related to the MPFA method
    /*
     * This Struct stores the transmissibility coefficients
     * for the two half-eges of an irregular interface (one
     * near a hanging node) in an h-adaptive simulation.
     */
    struct mpfaData
    {
        T1ransmissivityMatrix T1_[2];
        int globalIdx3_[storageRequirement];
        GlobalPosition globalPos3_[storageRequirement];
        std::vector<IntersectionIterator> secondHalfEdgeIntersection_;
        bool hasSecondHalfEdge;

        mpfaData()
        {
            hasSecondHalfEdge = false;
        }

        void setIntersection(IntersectionIterator& is23)
        {
            secondHalfEdgeIntersection_.push_back(is23);
            hasSecondHalfEdge = true;
        };
        const IntersectionIterator& getIntersection()
        {
            return secondHalfEdgeIntersection_[0];
        };
    };
    std::map<IdType, mpfaData> irregularInterfaceMap_;

    const Grid& grid_;

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

    void storeMpfaData(const typename GridView::Intersection & irregularIs,
    					const T1ransmissivityMatrix& T1,
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

    void storeMpfaData(const typename GridView::Intersection & irregularIs,
                        IntersectionIterator& secondHalfEdgeIntersectionIt,
                        const T1ransmissivityMatrix& T1,
                        const T1ransmissivityMatrix& T11_secondHalfEdge,
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

            irregularInterfaceMap_[intersectionID].T1_[second][2] = - T11_secondHalfEdge[0];
            irregularInterfaceMap_[intersectionID].T1_[second][1] = - T11_secondHalfEdge[1];
            irregularInterfaceMap_[intersectionID].T1_[second][0] = - T11_secondHalfEdge[2];

        }
        else // we look from smaller cell = unique interface
        {
            // proceed with numbering according to Aavatsmark, seen from cell i
            irregularInterfaceMap_[intersectionID].T1_[first] = T1;
            irregularInterfaceMap_[intersectionID].T1_[second] = T11_secondHalfEdge;
        }

        irregularInterfaceMap_[intersectionID].globalPos3_[0] = globalPos3;
        irregularInterfaceMap_[intersectionID].globalIdx3_[0] = globalIdx3;
        // second half edge
        irregularInterfaceMap_[intersectionID].setIntersection(secondHalfEdgeIntersectionIt);


        return;
    }

    int getMpfaData(const Intersection& irregularIs,
                        IntersectionIterator& secondHalfEdgeIntersectionIt,
                        T1ransmissivityMatrix& T1,
                        T1ransmissivityMatrix& T1_secondHalfEdge,
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
