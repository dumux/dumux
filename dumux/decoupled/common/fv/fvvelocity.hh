// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*****************************************************************************
 *   Copyright (C) 2011 by Markus Wolff                                      *
 *   Institute of Hydraulic Engineering                                      *
 *   University of Stuttgart, Germany                                        *
 *   email: <givenname>.<name>@iws.uni-stuttgart.de                          *
 *                                                                           *
 *   This program is free software: you can redistribute it and/or modify    *
 *   it under the terms of the GNU General Public License as published by    *
 *   the Free Software Foundation, either version 2 of the License, or       *
 *   (at your option) any later version.                                     *
 *                                                                           *
 *   This program is distributed in the hope that it will be useful,         *
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of          *
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the           *
 *   GNU General Public License for more details.                            *
 *                                                                           *
 *   You should have received a copy of the GNU General Public License       *
 *   along with this program.  If not, see <http://www.gnu.org/licenses/>.   *
 *****************************************************************************/
#ifndef DUMUX_FVVELOCITY_HH
#define DUMUX_FVVELOCITY_HH

// dumux environment
#include "dumux/common/math.hh"
#include <dumux/decoupled/common/impetproperties.hh>

/**
 * @file
 * @brief  Finite volume velocity reconstruction
 * @author Markus Wolff
 */

namespace Dumux
{

//! The finite volume base class for the solution of a pressure equation
/*! \ingroup multiphase
 *  TODO:
 * \tparam TypeTag The Type Tag
 */
template<class TypeTag> class FVVelocity
{
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(GridView)) GridView;
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(Scalar)) Scalar;
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(Problem)) Problem;

    typedef typename GET_PROP_TYPE(TypeTag, PTAG(CellData)) CellData;

    typedef typename GET_PROP_TYPE(TypeTag, Velocity) Velocity;

    enum
    {
        dim = GridView::dimension, dimWorld = GridView::dimensionworld
    };

    // typedefs to abbreviate several dune classes...
    typedef typename GridView::Traits::template Codim<0>::Entity Element;
    typedef typename GridView::template Codim<0>::Iterator ElementIterator;
    typedef typename GridView::IntersectionIterator IntersectionIterator;
    typedef typename GridView::Intersection Intersection;


public:

    //function which iterates through the grid and calculates the global velocity field
    void calculateVelocity();

    //! \brief Write data files
    /*  \param name file name */
    template<class MultiWriter>
    void addOutputVtkFields(MultiWriter &writer)
    {
        velocity_.addOutputVtkFields(writer);
    }

    //@}

    //! Constructs a FVVelocity object
    /**
     * \param problem a problem class object
     */
    FVVelocity(Problem& problem) :
        problem_(problem), velocity_(problem)
    {}

private:
    Problem& problem_;
    Velocity velocity_;
};


//! function which reconstructs a global velocity field
template<class TypeTag>
void FVVelocity<TypeTag>::calculateVelocity()
{
    Dune::FieldVector<Scalar, dim> vel(0.);

    ElementIterator eItEnd = problem_.gridView().template end<0> ();
    for (ElementIterator eIt = problem_.gridView().template begin<0> (); eIt != eItEnd; ++eIt)
    {
        // cell information
        int globalIdxI = problem_.variables().index(*eIt);
        CellData& cellDataI = problem_.variables().cellData(globalIdxI);


        /*****  flux term ***********/
        // iterate over all faces of the cell
        IntersectionIterator isItEnd = problem_.gridView().template iend(*eIt);
        for (IntersectionIterator isIt = problem_.gridView().template ibegin(*eIt); isIt != isItEnd; ++isIt)
        {
            /************* handle interior face *****************/
            if (isIt->neighbor())
            {
                int globalIdxJ = problem_.variables().index(*(isIt->outside()));

                //calculate only from one side
                if (globalIdxI > globalIdxJ)
                    continue;

                velocity_.calculateVelocity(*isIt, cellDataI);
            }   // end neighbor


            /************* boundary face ************************/
            else
            {
                velocity_.calculateVelocityOnBoundary(*isIt, cellDataI);
            }
        } //end interfaces loop
    } // end grid traversal

    return;
}

}//end namespace Dumux
#endif
