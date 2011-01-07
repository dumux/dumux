// $Id$
/*****************************************************************************
 *   Copyright (C) 2010 by Markus Wolff                                      *
 *   Copyright (C) 2009 by Andreas Lauser                                    *
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
/*!
 * \file
 * \brief Base class for all single phase diffusion problem
 * @author Markus Wolff
 */
#ifndef DUMUX_DIFFUSIONPROBLEM_1P_HH
#define DUMUX_DIFFUSIONPROBLEM_1P_HH

#include <dumux/decoupled/common/onemodelproblem.hh>
#include <dumux/decoupled/common/variableclass.hh>
#include <dumux/decoupled/1p/1pproperties.hh>

namespace Dumux
{
/*!
 * \ingroup OnePhase
 *
 * \brief  Base class for all single phase diffusion problem
 *
 * @tparam TypeTag The Type Tag
 * @tparam Implementation The Problem implementation
 */
template<class TypeTag, class Implementation>
class DiffusionProblem1P: public OneModelProblem<TypeTag, Implementation>
{
    typedef OneModelProblem<TypeTag, Implementation> ParentType;

    typedef typename GET_PROP_TYPE(TypeTag, PTAG(GridView)) GridView;
    typedef typename GridView::Grid Grid;typedef typename GET_PROP_TYPE(TypeTag, PTAG(Scalar)) Scalar;

    // material properties
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(Fluid)) Fluid;
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(SpatialParameters)) SpatialParameters;

    enum
    {
        dim = Grid::dimension, dimWorld = Grid::dimensionworld
    };

    typedef Dune::FieldVector<Scalar, dimWorld> GlobalPosition;

public:
    /*!
     * \brief The constructor
     *
     * \param gridView The grid view
     * \param verbose Output flag for the time manager.
     */
    DiffusionProblem1P(const GridView &gridView, bool verbose = true) :
        ParentType(gridView, verbose), gravity_(0)
    {
        spatialParameters_ = new SpatialParameters(gridView);
        newSpatialParams_ = true;
        gravity_ = 0;
        if (GET_PROP_VALUE(TypeTag, PTAG(EnableGravity)))
            gravity_[dim - 1] = -9.81;
    }
    /*!
     * \brief The constructor
     *
     * \param gridView The grid view
     * \param spatialParameters SpatialParameters instantiation
     * \param verbose Output flag for the time manager.
     */
    DiffusionProblem1P(const GridView &gridView, SpatialParameters &spatialParameters, bool verbose = true) :
        ParentType(gridView, verbose), gravity_(0), spatialParameters_(&spatialParameters)
    {
        newSpatialParams_ = false;
        gravity_ = 0;
        if (GET_PROP_VALUE(TypeTag, PTAG(EnableGravity)))
            gravity_[dim - 1] = -9.81;
    }

    virtual ~DiffusionProblem1P()
    {
        if (newSpatialParams_)
        {
            delete spatialParameters_;
        }
    }

    /*!
     * \name Problem parameters
     */
    // \{

    void timeIntegration()
    {
        // set the initial condition of the model
        ParentType::init();

        //end simulation -> no time dependent problem!
        this->timeManager().setFinished();

        return;
    }

    void serialize()
    {
        return;
    }
    void deserialize(double t)
    {
        return;
    }

    /*!
     * \brief Returns the temperature within the domain.
     *
     * This method MUST be overwritten by the actual problem.
     */
    Scalar temperature() const
    {
        return this->asImp_()->temperature();
    }
    ;

    /*!
     * \brief Returns the acceleration due to gravity.
     *
     * If the <tt>EnableGravity</tt> property is true, this means
     * \f$\boldsymbol{g} = ( 0,\dots,\ -9.81)^T \f$, else \f$\boldsymbol{g} = ( 0,\dots, 0)^T \f$
     */
    const GlobalPosition &gravity() const
    {
        return gravity_;
    }

    /*!
     * \brief Returns the spatial parameters object.
     */
    SpatialParameters &spatialParameters()
    {
        return *spatialParameters_;
    }

    /*!
     * \brief Returns the spatial parameters object.
     */
    const SpatialParameters &spatialParameters() const
    {
        return *spatialParameters_;
    }

    // \}

private:
    GlobalPosition gravity_;

    // fluids and material properties
    SpatialParameters* spatialParameters_;
    bool newSpatialParams_;
};

}

#endif
