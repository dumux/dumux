// $Id$
/*****************************************************************************
 *   Copyright (C) 2009 by Karin Erbertseder                                 *
 *   Copyright (C) 2009 by Andreas Lauser                                    *
 *   Copyright (C) 2008 by Bernd Flemisch                                    *
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
 * \brief Quantities required by the single-phase, two-component box
 *        model defined on a vertex.
 */
#ifndef DUMUX_1P2C_VOLUME_VARIABLES_HH
#define DUMUX_1P2C_VOLUME_VARIABLES_HH

#include "1p2cfluidstate.hh"

namespace Dumux
{

/*!
 * \brief Contains the quantities which are are constant within a
 *        finite volume in the single-phase, two-component model.
 */
template <class TypeTag>
class OnePTwoCVolumeVariables
{
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(Scalar)) Scalar;
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(GridView)) GridView;
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(Problem)) Problem;
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(OnePTwoCIndices)) Indices;
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(VolumeVariables)) Implementation;
    typedef OnePTwoCFluidState<TypeTag> FluidState;

    typedef typename GridView::template Codim<0>::Entity Element;

    enum {
        numEq = GET_PROP_VALUE(TypeTag, PTAG(NumEq)),
        numPhases = GET_PROP_VALUE(TypeTag, PTAG(NumPhases)),
        numComponents = GET_PROP_VALUE(TypeTag, PTAG(NumComponents)),

        dim = GridView::dimension,
        dimWorld = GridView::dimensionworld,

        konti = Indices::konti,
        transport = Indices::transport
    };


    typedef typename GET_PROP(TypeTag, PTAG(ReferenceElements)) RefElemProp;
    typedef typename RefElemProp::Container ReferenceElements;

    typedef typename GET_PROP_TYPE(TypeTag, PTAG(FVElementGeometry)) FVElementGeometry;

    typedef typename GET_PROP_TYPE(TypeTag, PTAG(PrimaryVariables)) PrimaryVariables;
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(FluidSystem)) FluidSystem;

    typedef Dune::FieldVector<Scalar, dimWorld> GlobalPosition;
    typedef Dune::FieldVector<Scalar, dim> LocalPosition;

public:
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

        porosity = problem.spatialParameters().porosity(element, elemGeom, scvIdx);
        tortuosity = problem.spatialParameters().tortuosity(element, elemGeom, scvIdx);
        dispersivity = problem.spatialParameters().dispersivity(element, elemGeom, scvIdx);

        Scalar temperature = problem.temperature(element, elemGeom, scvIdx);
        fluidState_.update(priVars, temperature);
        pressure = fluidState_.phasePressure(konti);
        molefraction = fluidState_.moleFrac(konti, transport);

        density = FluidSystem::phaseDensity(konti, temperature, pressure, fluidState_);
        molarDensity = FluidSystem::molarDensity(konti, temperature, pressure, fluidState_);
        viscosity = FluidSystem::phaseViscosity(konti, temperature, pressure, fluidState_);
        diffCoeff = FluidSystem::diffCoeff(konti, transport, konti, temperature, pressure, fluidState_);
    }

    /*!
     * \brief Sets the evaluation point used in the by the local jacobian.
     */
    void setEvalPoint(const Implementation *ep)
    { }

    /*!
     * \brief Return the vector of primary variables
     */
    const PrimaryVariables &primaryVars() const
    { return primaryVars_; }

    Scalar porosity;
    Scalar density;
    Scalar viscosity;
    Scalar tortuosity;
    Dune::FieldVector<Scalar,dim> dispersivity;
    Scalar diffCoeff;
    Scalar molefraction;
    Scalar pressure;
    Scalar molarDensity;
    FluidState fluidState_;

protected:
    PrimaryVariables primaryVars_;
};

}

#endif
