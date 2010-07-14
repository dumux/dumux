// $Id: richardsvertexdata.hh 3784 2010-06-24 13:43:57Z bernd $
/*****************************************************************************
 *   Copyright (C) 2009 by Onur Dogan                                        *
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
 * \brief Quantities required by the richards box model defined on a vertex.
 */
#ifndef DUMUX_RICHARDS_VERTEX_DATA_HH
#define DUMUX_RICHARDS_VERTEX_DATA_HH

#include "richardsproperties.hh"

namespace Dumux
{

/*!
 * \ingroup RichardsBoxModel
 * \brief Contains the quantities which are are constant within a
 *        finite volume in the Richards model.
 */
template <class TypeTag>
class RichardsVertexData
{
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(Scalar))   Scalar;
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(GridView)) GridView;

    typedef typename GridView::template Codim<0>::Entity Element;

    enum {
        dim           = GridView::dimension,
        dimWorld      = GridView::dimensionworld,
    };

    typedef typename GET_PROP(TypeTag, PTAG(ReferenceElements)) RefElemProp;
    typedef typename RefElemProp::Container                     ReferenceElements;

    typedef typename GET_PROP_TYPE(TypeTag, PTAG(Problem)) Problem;
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(FVElementGeometry)) FVElementGeometry;
    typedef typename GET_PROP(TypeTag, PTAG(SolutionTypes))     SolutionTypes;
    typedef typename SolutionTypes::PrimaryVarVector            PrimaryVarVector;
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(RichardsIndices))  Indices;

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
        typedef Indices I;

        // coordinates of the vertex
        const GlobalPosition &global = element.geometry().corner(vertIdx);
        const LocalPosition   &local =
            ReferenceElements::general(element.type()).position(vertIdx,
                                                                GridView::dimension);

        /* pc = pNreference - pw || pc = 0 for computing Sw */
        pNreference = problem.pNreference();
        pW = sol[I::pW];
        if (pW >= pNreference)
            pC = 0.0;
        else
            pC = pNreference-pW;

        dSwdpC = problem.materialLaw().dSdP(pC,
                                            global,
                                            element,
                                            local);
        Sw = problem.materialLaw().saturationW(pC,
                                               global,
                                               element,
                                               local);

        temperature = problem.temperature(element, elemGeom, vertIdx);
        mobilityW = problem.materialLaw().mobW(Sw,
                                               global,
                                               element,
                                               local,
                                               temperature,
                                               pW);
        densityW = problem.wettingPhase().density(temperature,
                                                  pW);
        porosity = problem.soil().porosity(global,
                                           element,
                                           local);
    }

    Scalar pNreference;
    Scalar pW;
    Scalar pC;
    Scalar Sw;
    Scalar dSwdpC;

    Scalar densityW;
    Scalar mobilityW;
    Scalar porosity;
    Scalar temperature;
};

}

#endif
