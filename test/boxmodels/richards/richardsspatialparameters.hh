// $Id$
/*****************************************************************************
 *   Copyright (C) 2010 by Bernd Flemisch                                    *
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
#ifndef DUMUX_RICHARDSSPATIALPARAMETERS_HH
#define DUMUX_RICHARDSSPATIALPARAMETERS_HH

#include <dumux/material/spatialparameters/boxspatialparameters.hh>
#include <dumux/material/fluidmatrixinteractions/2p/regularizedlinearmaterial.hh>
#include <dumux/material/fluidmatrixinteractions/2p/efftoabslaw.hh>

/**
 * @file
 * @brief  Class for defining spatial parameters
 * @author Bernd Flemisch, Klaus Mosthaf, Markus Wolff
 */

namespace Dumux
{

/** \todo Please doc me! */

template<class TypeTag>
class RichardsSpatialParameters : public BoxSpatialParameters<TypeTag>
{
    typedef BoxSpatialParameters<TypeTag> ParentType;
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(Grid)) Grid;
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(GridView)) GridView;
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(Scalar)) Scalar;
    typedef typename Grid::ctype CoordScalar;

    enum {
        dim=GridView::dimension,
        dimWorld=GridView::dimensionworld,
    };

    typedef Dune::FieldVector<CoordScalar,dimWorld> GlobalPosition;

    typedef typename GridView::template Codim<0>::Entity Element;
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(FVElementGeometry)) FVElementGeometry;

    // define the material law which is parameterized by effective saturations
    typedef RegularizedLinearMaterial<Scalar> EffectiveLaw;

public:
    // define the material law parameterized by absolute saturations
    typedef EffToAbsLaw<EffectiveLaw> MaterialLaw;
    typedef typename MaterialLaw::Params MaterialLawParams;

    RichardsSpatialParameters(const GridView& gridView)
        : ParentType(gridView)
    {
        // residual saturations
        materialParams_.setSwr(0.05);
        materialParams_.setSnr(0.0);

        // parameters for the linear law
        materialParams_.setEntryPC(0);
        materialParams_.setMaxPC(0);

        perm_ = 1e-11;
        porosity_ = 0.4;
    }

    Scalar intrinsicPermeability(const Element           &element,
                                 const FVElementGeometry &fvElemGeom,
                                 int                      scvIdx) const
    {
        return perm_;
    }

    Scalar porosity(const Element           &element,
                    const FVElementGeometry &fvElemGeom,
                    int                      scvIdx) const
    {
        return porosity_;
    }

    // return the brooks-corey context depending on the position
    const MaterialLawParams& materialLawParams(const Element           &element,
                                                const FVElementGeometry &fvElemGeom,
                                                int                      scvIdx) const
    {
        return materialParams_;
    }

private:
    Scalar perm_;
    Scalar porosity_;
    MaterialLawParams materialParams_;
};

} // end namespace
#endif

