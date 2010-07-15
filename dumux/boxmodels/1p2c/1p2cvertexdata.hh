// $Id: 1p2cvertexdata.hh 3838 2010-07-15 08:31:53Z bernd $
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
#ifndef DUMUX_1P2C_VERTEX_DATA_HH
#define DUMUX_1P2C_VERTEX_DATA_HH

#include "1p2cfluidstate.hh"

namespace Dumux
{

/*!
 * \brief Contains the quantities which are are constant within a
 *        finite volume in the single-phase, two-component model.
 */
template <class TypeTag>
class OnePTwoCVertexData
{
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(Scalar))   Scalar;
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(GridView)) GridView;
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(Problem))  Problem;
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(OnePTwoCIndices)) Indices;
    typedef OnePTwoCFluidState<TypeTag> FluidState;

    typedef typename GridView::template Codim<0>::Entity Element;

    enum {
        numEq         = GET_PROP_VALUE(TypeTag, PTAG(NumEq)),
        numPhases     = GET_PROP_VALUE(TypeTag, PTAG(NumPhases)),
        numComponents = GET_PROP_VALUE(TypeTag, PTAG(NumComponents)),

        dim           = GridView::dimension,
        dimWorld      = GridView::dimensionworld,

        konti = Indices::konti,
        transport = Indices::transport
    };

    typedef typename GET_PROP(TypeTag, PTAG(SolutionTypes))     SolutionTypes;
    typedef typename GET_PROP(TypeTag, PTAG(ReferenceElements)) RefElemProp;
    typedef typename RefElemProp::Container                     ReferenceElements;

    typedef typename GET_PROP_TYPE(TypeTag, PTAG(FVElementGeometry)) FVElementGeometry;

    typedef typename SolutionTypes::PrimaryVarVector  PrimaryVarVector;
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(FluidSystem)) FluidSystem;

    typedef Dune::FieldVector<Scalar, dimWorld>  GlobalPosition;
    typedef Dune::FieldVector<Scalar, dim>       LocalPosition;

public:
    /*!
     * \brief Update all quantities for a given control volume.
     */
    void update(const PrimaryVarVector  &sol,
                const Element           &element,
                const FVElementGeometry &elemGeom,
                int                      vertIdx,
                const Problem           &problem,
                bool                     isOldSol)
    {
        porosity = problem.spatialParameters().porosity(element, elemGeom, vertIdx);
        tortuosity = problem.spatialParameters().tortuosity(element, elemGeom, vertIdx);
        dispersivity = problem.spatialParameters().dispersivity(element, elemGeom, vertIdx);

        Scalar temperature = problem.temperature(element, elemGeom, vertIdx);
        fluidState_.update(sol, temperature);
        pressure = fluidState_.phasePressure(konti);
        molefraction = fluidState_.moleFrac(konti, transport);

        density = FluidSystem::phaseDensity(konti, temperature, pressure, fluidState_);
        molarDensity = FluidSystem::molarDensity(konti, temperature, pressure, fluidState_);
        viscosity = FluidSystem::phaseViscosity(konti, temperature, pressure, fluidState_);
        diffCoeff = FluidSystem::diffCoeff(konti, transport, konti, temperature, pressure, fluidState_);
    }

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

};

}

#endif
