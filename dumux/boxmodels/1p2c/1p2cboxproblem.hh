/*****************************************************************************
 *   Copyright (C) 2009 by Andreas Lauser                                    *
 *   Institute of Hydraulic Engineering                                      *
 *   University of Stuttgart, Germany                                        *
 *   email: <givenname>.<name>@iws.uni-stuttgart.de                          *
 *                                                                           *
 *   This program is free software; you can redistribute it and/or modify    *
 *   it under the terms of the GNU General Public License as published by    *
 *   the Free Software Foundation; either version 2 of the License, or       *
 *   (at your option) any later version, as long as this copyright notice    *
 *   is included in its original form.                                       *
 *                                                                           *
 *   This program is distributed WITHOUT ANY WARRANTY.                       *
 *****************************************************************************/
/*!
 * \file
 *
 * \brief Base class for all problems which use the single-phase,
 *        two-component box model
 */
#ifndef DUMUX_1P2C_BOX_PROBLEM_HH
#define DUMUX_1P2C_BOX_PROBLEM_HH

#include <dumux/boxmodels/tags.hh>
#include <dumux/boxmodels/boxscheme/boxproblem.hh>
#include <dumux/material/twophaserelations.hh>


namespace Dune
{
/*!
 * \ingroup OnePTwoCProblems
 * \brief  Base class for all problems which use the single-phase, two-component box model
 *
 * \todo Please doc me more!
 */
template<class TypeTag, class Implementation>
class OnePTwoCBoxProblem : public BoxProblem<TypeTag, Implementation>
{
    typedef BoxProblem<TypeTag, Implementation> ParentType;

    typedef typename GET_PROP_TYPE(TypeTag, PTAG(GridView)) GridView;
    typedef typename GridView::Grid                         Grid;
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(Scalar))   Scalar;

    // material properties
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(Fluid))    Fluid;
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(Soil))     Soil;

    enum {
        dim = Grid::dimension,
        dimWorld = Grid::dimensionworld
    };

    typedef Dune::FieldVector<Scalar, dimWorld>      GlobalPosition;

public:
    OnePTwoCBoxProblem(const GridView &gridView)
        : ParentType(gridView),
          gravity_(0)
    {
        if (GET_PROP_VALUE(TypeTag, PTAG(EnableGravity)))
            gravity_[dim - 1] = - 9.81;
    }

    /*!
     * \name Problem parameters
     */
    // \{

    /*!
     * \brief Returns the temperature within the domain.
     *
     * This method MUST be overwritten by the actual problem.
     */
    Scalar temperature() const
    { return asImp_()->temperature(); };

    /*!
     * \brief Returns the acceleration due to gravity.
     *
     * If the <tt>EnableGravity</tt> property is true, this means
     * \f$\boldsymbol{g} = ( 0,\dots,\ -9.81)^T \f$, else \f$\boldsymbol{g} = ( 0,\dots, 0)^T \f$
     */
    const GlobalPosition &gravity() const
    { return gravity_; }

    /*! 
     * \brief Fluid properties of the wetting phase.
     */
    const Fluid &fluid() const
    { return fluid_; }

    /*! 
     * \brief Returns the soil properties object.
     */
    Soil &soil()
    { return soil_; }

    /*! 
     * \copydoc soil()
     */
    const Soil &soil() const
    { return soil_; }
    
    // \}

private:
    //! Returns the implementation of the problem (i.e. static polymorphism)
    Implementation *asImp_()
    { return static_cast<Implementation *>(this); }
    
    //! \copydoc asImp_()
    const Implementation *asImp_() const 
    { return static_cast<const Implementation *>(this); }

    // fluids and material properties
    Fluid           fluid_;
    Soil            soil_;

    GlobalPosition  gravity_;
};

}

#endif
