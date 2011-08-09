/*****************************************************************************
 *   Copyright (C) 2010 by Benjamin Faigle, Markus Wolff                     *
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
 * \brief Base class for sequential 2p2c compositional problems
 */
#ifndef DUMUX_IMPETPROBLEM_2P2C_HH
#define DUMUX_IMPETPROBLEM_2P2C_HH

#include "boundaryconditions2p2c.hh"
#include <dumux/decoupled/common/impet.hh>
#include <dumux/decoupled/common/impetproblem.hh>
#include <dumux/decoupled/2p2c/variableclass2p2c.hh>
#include <dumux/decoupled/2p2c/2p2cproperties.hh>


namespace Dumux
{
/*!
 * \ingroup IMPEC
 * \ingroup IMPECproblems
 * \brief  Base class for all compositional 2-phase problems which use an impet algorithm
 *
 * Differs from .../2p/impes/impesproblem2p.hh in the includes
 * (usage of the compositional properties and variableclass) and
 * empty functions for the boundary formulation. Because of the latter,
 * only the functions of the currently used formulation has to be specified
 * in the specific problem.
 */
template<class TypeTag>
class IMPETProblem2P2C : public IMPETProblem<TypeTag>
{
    typedef IMPETProblem<TypeTag> ParentType;
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(Problem)) Implementation;

    typedef typename GET_PROP_TYPE(TypeTag, PTAG(TimeManager)) TimeManager;

    typedef typename GET_PROP_TYPE(TypeTag, PTAG(GridView)) GridView;
    typedef typename GridView::Grid                         Grid;
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(Scalar))   Scalar;

    typedef typename GridView::Traits::template Codim<0>::Entity Element;

    // material properties
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(FluidSystem))            FluidSystem;
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(SpatialParameters))    SpatialParameters;


    enum {
        dim = Grid::dimension,
        dimWorld = Grid::dimensionworld
    };

    typedef Dune::FieldVector<Scalar, dimWorld>      GlobalPosition;

public:
    /*!
     * \brief The constructor
     *
     * \param gridView The grid view
     * \param verbose Output flag for the time manager.
     */
    IMPETProblem2P2C(TimeManager &timeManager, const GridView &gridView)
        : ParentType(timeManager, gridView),
        gravity_(0)
    {
        newSpatialParams_ = true;
        spatialParameters_ = new SpatialParameters(gridView);

        gravity_ = 0;
        if (GET_PARAM(TypeTag, bool, EnableGravity))
            gravity_[dim - 1] = - 9.81;
    }
    /*!
     * \brief The constructor
     *
     * \param gridView The grid view
     * \param spatialParameters SpatialParameters instantiation
     * \param verbose Output flag for the time manager.
     */
    IMPETProblem2P2C(TimeManager &timeManager, const GridView &gridView, SpatialParameters &spatialParameters)
        : ParentType(timeManager, gridView),
        gravity_(0),spatialParameters_(&spatialParameters)
    {
        newSpatialParams_ = false;
        gravity_ = 0;
        if (GET_PARAM(TypeTag, bool, EnableGravity))
            gravity_[dim - 1] = - 9.81;
    }

    /*!
     * \brief The constructor
     *
     * \param gridView The grid view
     * \param verbose Output flag for the time manager.
     */
    IMPETProblem2P2C(const GridView &gridView, bool verbose = true)
    DUNE_DEPRECATED // use IMPETProblem2P2C(TimeManager&, const GridView &)
        : ParentType(gridView, verbose),
        gravity_(0)
    {
        newSpatialParams_ = true;
        spatialParameters_ = new SpatialParameters(gridView);

        gravity_ = 0;
        if (GET_PARAM(TypeTag, bool, EnableGravity))
            gravity_[dim - 1] = - 9.81;
    }
    /*!
     * \brief The constructor
     *
     * \param gridView The grid view
     * \param spatialParameters SpatialParameters instantiation
     * \param verbose Output flag for the time manager.
     */
    IMPETProblem2P2C(const GridView &gridView, SpatialParameters &spatialParameters, bool verbose = true)
    DUNE_DEPRECATED // use IMPETProblem2P2C(TimeManager&, const GridView &)
        : ParentType(gridView, verbose),
        gravity_(0),spatialParameters_(&spatialParameters)
    {
        newSpatialParams_ = false;
        gravity_ = 0;
        if (GET_PARAM(TypeTag, bool, EnableGravity))
            gravity_[dim - 1] = - 9.81;
    }

    virtual ~IMPETProblem2P2C()
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

    /*!
     * \copydoc Dumux::IMPESProblem2P::temperature()
     */
    Scalar temperature() const
    { return this->asImp_()->temperature(); };

    /*!
     * \copydoc Dumux::IMPESProblem2P::gravity()
     */
    const GlobalPosition &gravity() const
    { return gravity_; }

    //! Saturation initial condition (dimensionless)
    /*! The problem is initialized with the following saturation. Both
     * phases are assumed to contain an equilibrium concentration of the
     * correspondingly other component.
     */
    Scalar initSat(const GlobalPosition& globalPos, const Element& element) const
    {
        DUNE_THROW(Dune::NotImplemented, "please specify initial saturation in the problem!");
        return NAN;
    }
    //! Concentration initial condition (dimensionless)
    /*! The problem is initialized with the following concentration.
     */
    Scalar initConcentration(const GlobalPosition& globalPos, const Element& element) const
    {
        DUNE_THROW(Dune::NotImplemented, "please specify initial Concentration in the problem!");
        return NAN;
    }
    // \}

    /*!
     * \name Access functions
     */
    // \{
    /*!
     * \copydoc Dumux::IMPESProblem2P::spatialParameters()
     */
    SpatialParameters &spatialParameters()
    { return *spatialParameters_; }

    /*!
     * \copydoc Dumux::IMPESProblem2P::spatialParameters()
     */
    const SpatialParameters &spatialParameters() const
    { return *spatialParameters_; }

    // \}

private:
    //! Returns the implementation of the problem (i.e. static polymorphism)
    Implementation &asImp_()
    { return *static_cast<Implementation *>(this); }

    //! \copydoc Dumux::IMPETProblem::asImp_()
    const Implementation &asImp_() const
    { return *static_cast<const Implementation *>(this); }

    GlobalPosition gravity_;

    // fluids and material properties
    SpatialParameters*  spatialParameters_;
    bool newSpatialParams_;
};

}

#endif
