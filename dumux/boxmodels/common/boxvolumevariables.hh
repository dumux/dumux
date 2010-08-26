/*****************************************************************************
 *   Copyright (C) 2010 by Andreas Lauser                                    *
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
 * \brief Base class for the model specific classes which provide
 *        access to all volume averaged quantities.
 */
#ifndef DUMUX_BOX_VOLUME_VARIABLES_HH
#define DUMUX_BOX_VOLUME_VARIABLES_HH

#include "boxproperties.hh"

namespace Dumux
{

/*!
 * \brief Base class for the model specific classes which provide
 *        access to all volume averaged quantities.
 */
template <class TypeTag>
class BoxVolumeVariables
{
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(VolumeVariables)) Implementation;

    typedef typename GET_PROP_TYPE(TypeTag, PTAG(Scalar)) Scalar;
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(Problem)) Problem;
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(GridView)) GridView;
    typedef typename GridView::template Codim<0>::Entity Element;
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(FVElementGeometry)) FVElementGeometry;

    typedef typename GET_PROP_TYPE(TypeTag, PTAG(PrimaryVariables)) PrimaryVariables;

public:
    // default constructor
    BoxVolumeVariables()
    { evalPoint_ = 0; };

    // copy constructor
    BoxVolumeVariables(const BoxVolumeVariables &v)
    { 
        primaryVars_ = v.primaryVars_;
        evalPoint_ = 0;
    };

    // assignment operator
    BoxVolumeVariables &operator=(const BoxVolumeVariables &v)
    {
        evalPoint_ = 0;
        return *this;
    };

    /*!
     * \brief Sets the evaluation point used in the by the local jacobian.
     */
    void setEvalPoint(const Implementation *ep)
    { evalPoint_ = ep; }

    /*!
     * \brief Returns the evaluation point used in the by the local jacobian.
     */
    const Implementation &evalPoint() const
    { return (evalPoint_ == 0)?asImp_():evalPoint_; }

    /*!
     * \brief Update all quantities for a given control volume.
     */
    void update(const PrimaryVariables &priVars,
                const Problem &problem,
                const Element &element,
                const FVElementGeometry &elemGeom,
                int scvIdx,
                bool isOldSol)
    {
        primaryVars_ = priVars;
    }

    /*!
     * \brief Return the vector of primary variables
     */
    const PrimaryVariables &primaryVars() const
    { return primaryVars_; }

    /*!
     * \brief Return a component of primary variable vector
     */
    Scalar primaryVars(int pvIdx) const
    { return primaryVars_[pvIdx]; }

protected:
    const Implementation &asImp_() const
    { return *static_cast<Implementation*>(this); }
    Implementation &asImp_()
    { return *static_cast<Implementation*>(this); }

    PrimaryVariables primaryVars_;

    // the evaluation point of the local jacobian
    const Implementation *evalPoint_;
};

} // end namepace

#endif
