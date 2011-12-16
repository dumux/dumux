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
#ifndef DUMUX_VariableClass2PGridAdapt_GRIDADAPT_HH
#define DUMUX_VariableClass2PGridAdapt_GRIDADAPT_HH

#include <dumux/decoupled/2p/variableclass2p.hh>
#include <dune/grid/utility/persistentcontainer.hh>

/**
 * @file
 * @brief  Class including the variables and data of discretized data of the constitutive relations
 * @author Markus Wolff, Benjamin Faigle
 */

namespace Dumux
{
/*!
 * \ingroup IMPES
 */
//! Class including the variables and data of discretized data of the constitutive relations.
/*! The variables of two-phase flow, which are one pressure and one saturation are stored in this class.
 *
 * \tparam TypeTag The Type Tag
 */
template<class TypeTag>
class VariableClass2PGridAdapt: public VariableClass2P<TypeTag>
{
private:
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(Scalar)) Scalar;
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(GridView)) GridView;
    typedef typename GridView::Grid                         Grid;
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(Problem)) Problem;

    typedef typename GET_PROP(TypeTag, PTAG(SolutionTypes)) SolutionTypes;
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(TransportSolutionType)) TransportSolutionType;

    typedef typename GET_PROP_TYPE(TypeTag, PTAG(TwoPIndices)) Indices;

    typedef VariableClass2P<TypeTag> ParentClass;

    enum
    {
        dim = GridView::dimension, dimWorld = GridView::dimensionworld
    };
    enum
    {
        pw = Indices::pressureW, pn = Indices::pressureNW, pglobal = Indices::pressureGlobal,
    };
    enum
    {
        wPhaseIdx = Indices::wPhaseIdx, nPhaseIdx = Indices::nPhaseIdx
    };

    static const int pressureType_ = GET_PROP_VALUE(TypeTag, PTAG(PressureFormulation));

    typedef typename GridView::Traits::template Codim<0>::Entity Element;
    typedef typename Grid::template Codim<0>::EntityPointer ElementPointer;

    typedef Dune::FieldVector<Scalar, dim> LocalPosition;
    typedef Dune::FieldVector<Scalar, dimWorld> GlobalPosition;
    typedef typename GridView::IndexSet IndexSet;
    typedef typename Grid::LocalIdSet IdSet;
    typedef typename IdSet::IdType IdType;

    typedef typename Grid::LevelGridView LevelGridView;
    typedef typename LevelGridView::template Codim<0>::Iterator LevelIterator;


    typedef typename GridView::template Codim<0>::Iterator ElementIterator;
    typedef typename GridView::IntersectionIterator IntersectionIterator;

    //*******************************
    // TODO: Doc me! Besseren Namen finden!
    struct RestrictedValue
    {
        Scalar saturation;
        Scalar press;
        Scalar volCorr;
        int count;
        RestrictedValue()
        {
            saturation = 0.;
            press = 0.;
            count = 0;
            volCorr = 0;
        }
    };

    const Grid& grid_;
    Dune::PersistentContainer<Grid, RestrictedValue> restrictionmap_;

public:

    //! Constructs a VariableClass object
    /**
     *  @param gridView a DUNE gridview object corresponding to diffusion and transport equation
     */

    VariableClass2PGridAdapt(const GridView& gridView) :
                VariableClass2P<TypeTag> (gridView), grid_(gridView.grid()), restrictionmap_(grid_,0) //codim 0 map
    {}


    /*!
     * Store primary variables
     *
     * To reconstruct the solution in father elements, problem properties might
     * need to be accessed.
     *
     * @param problem The current problem
     */
    void storePrimVars(const Problem& problem)
    {
        // loop over all levels of the grid
        for (int level=grid_.maxLevel(); level>=0; level--)
        {
            //get grid view on level grid
            LevelGridView levelView = grid_.levelView(level);
            for (LevelIterator it = levelView.template begin<0>();
                    it!=levelView.template end<0>(); ++it)
            {
                //get your map entry
                RestrictedValue &rv = restrictionmap_[*it];

                // put your value in the map
                if (it->isLeaf())
                {
                    // get index
                    int indexI = this->index(*it);

                    rv.saturation = this->saturation()[indexI][0];
                    rv.press = this->pressure()[indexI][0];
                    rv.volCorr = this->volumecorrection(indexI);
                    rv.count = 1;
                }
                //Average in father
                if (it.level()>0)
                {
                    ElementPointer epFather = it->father();
                    RestrictedValue& rvf = restrictionmap_[*epFather];
                    rvf.saturation += rv.saturation/rv.count;
                    rvf.press += rv.press/rv.count;
                    rvf.volCorr += rv.volCorr/rv.count;
                    rvf.count += 1;
                }
            }
        }
    }

    /*!
     * Reconstruct missing primary variables (where elements are created/deleted)
     *
     * To reconstruct the solution in father elements, problem properties might
     * need to be accessed.
     *
     * @param problem The current problem
     */
    void reconstructPrimVars(const Problem& problem)
    {
        restrictionmap_.reserve();

        for (int level=0; level<=grid_.maxLevel(); level++)
        {
            LevelGridView levelView = grid_.levelView(level);
            for (LevelIterator it = levelView.template begin <0>();
                    it!=levelView.template end <0>(); ++it)
            {
                if (!it->isNew())
                {
                    //entry is in map, write in leaf
                    if (it->isLeaf())
                    {
                        RestrictedValue &rv = restrictionmap_[*it];
                        int newIdxI = this->index(*it);
                        this->saturation()[newIdxI][0] = rv.saturation/rv.count;
                        this->pressure()[newIdxI][0] = rv.press/rv.count;
                        this->volumecorrection(newIdxI) = rv.volCorr/rv.count;
                    }
                }
                else
                {
                    // value is not in map, interpolate from father element
                    if (it.level()>0)
                    {
                        ElementPointer ep = it->father();
                        RestrictedValue& rvf = restrictionmap_[*ep];
                        if (it->isLeaf())
                        {
                            int newIdxI = this->index(*it);
                            this->saturation()[newIdxI][0] = rvf.saturation/rvf.count;
                            this->pressure()[newIdxI][0] = rvf.press/rvf.count;
                            this->volumecorrection(newIdxI) = rvf.volCorr/rvf.count;
                        }
                        else
                        {
                            //create new entry
                            RestrictedValue& rv = restrictionmap_[*it];
                            rv.saturation =rvf.saturation/rvf.count;
                            rv.press =rvf.press/rvf.count;
                            rv.volCorr =rvf.volCorr/rvf.count;
                            rv.count = 1;
                        }
                    }
                }
            }
        }
        // reset entries in restrictionmap
        restrictionmap_.clear();
    }
};
}
#endif
