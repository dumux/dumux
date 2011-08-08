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
 * \brief Base class for stationary sequential 2-phase problems
 * @author Markus Wolff
 */
#ifndef DUMUX_DIFFUSIONPROBLEM_2P_HH
#define DUMUX_DIFFUSIONPROBLEM_2P_HH

#include <dumux/decoupled/common/onemodelproblem.hh>
#include <dumux/decoupled/2p/variableclass2p.hh>
#include <dumux/material/fluidsystems/2p_system.hh>
#include <dumux/decoupled/2p/2pproperties.hh>

namespace Dumux
{
/*!
 * \ingroup IMPESproblems
 * \ingroup Pressure2p
 * \brief  Base class for stationary sequential 2-phase problems
 *
 * @tparam TypeTag The Type Tag
 * @tparam Implementation The Problem implementation
 */
template<class TypeTag>
class DiffusionProblem2P: public OneModelProblem<TypeTag>
{
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(Problem)) Implementation;
    typedef OneModelProblem<TypeTag> ParentType;

    typedef typename GET_PROP_TYPE(TypeTag, PTAG(GridView)) GridView;
    typedef typename GridView::Grid Grid;typedef typename GET_PROP_TYPE(TypeTag, PTAG(Scalar)) Scalar;

    typedef typename GET_PROP_TYPE(TypeTag, PTAG(TimeManager)) TimeManager;

    // material properties
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(FluidSystem)) FluidSystem;
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(SpatialParameters)) SpatialParameters;

    typedef typename GridView::Traits::template Codim<0>::Entity Element;

    enum
    {
        dim = Grid::dimension, dimWorld = Grid::dimensionworld
    };

    typedef Dune::FieldVector<Scalar, dimWorld> GlobalPosition;

    // private!! copy constructor
    DiffusionProblem2P(const DiffusionProblem2P&)
    {}

public:
    /*!
     * \brief The constructor
     *
     * \param gridView The grid view
     */
    DiffusionProblem2P(const GridView &gridView)
    : ParentType(gridView), gravity_(0)
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
     */
    DiffusionProblem2P(const GridView &gridView, SpatialParameters &spatialParameters)
    : ParentType(gridView), gravity_(0), spatialParameters_(&spatialParameters)
    {
        newSpatialParams_ = false;
        gravity_ = 0;
        if (GET_PROP_VALUE(TypeTag, PTAG(EnableGravity)))
            gravity_[dim - 1] = -9.81;
    }

    virtual ~DiffusionProblem2P()
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
     * \param element The element
     *
     */
    Scalar temperature(const Element& element) const
    {
        return this->asImp_().temperatureAtPos(element.geometry().center());
    }

    /*!
     * \brief Returns the temperature within the domain.
     *
     * \param globalPos The position of the center of an element
     *
     */
    Scalar temperatureAtPos(const GlobalPosition& globalPos) const
    {
        // Throw an exception (there is no initial condition)
        DUNE_THROW(Dune::InvalidStateException,
                   "The problem does not provide "
                   "a temperatureAtPos() method.");
    }

    /*!
     * \brief Returns the reference pressure for evaluation of constitutive relations.
     *
     * \param element The element
     *
     */
    Scalar referencePressure(const Element& element) const
    {
        return this->asImp_().referencePressureAtPos(element.geometry().center());
    }

    /*!
     * \brief Returns the reference pressure for evaluation of constitutive relations.
     *
     * \param globalPos The position of the center of an element
     *
     */
    Scalar referencePressureAtPos(const GlobalPosition& globalPos) const
    {
        // Throw an exception (there is no initial condition)
        DUNE_THROW(Dune::InvalidStateException,
                   "The problem does not provide "
                   "a referencePressureAtPos() method.");
    }

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
    //! Returns the implementation of the problem (i.e. static polymorphism)
    Implementation &asImp_()
    { return *static_cast<Implementation *>(this); }

    //! \copydoc Dumux::IMPETProblem::asImp_()
    const Implementation &asImp_() const
    { return *static_cast<const Implementation *>(this); }

    TimeManager timeManager_;

    GlobalPosition gravity_;

    // fluids and material properties
    SpatialParameters* spatialParameters_;
    bool newSpatialParams_;
};

}

#endif
