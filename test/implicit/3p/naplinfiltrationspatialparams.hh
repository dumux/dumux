// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*****************************************************************************
 *   Copyright (C) 2011      by Holger Class                                 *
 *   Copyright (C) 2008-2010 by Andreas Lauser                               *
 *   Copyright (C) 2008-2009 by Klaus Mosthaf                                *
 *   Institute for Modelling Hydraulic and Environmental Systems             *
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
 *
 * \brief Definition of the spatial parameters for the kuevette problem.
 */
#ifndef DUMUX_NAPLINFILTRATION_SPATIAL_PARAMS_HH
#define DUMUX_NAPLINFILTRATION_SPATIAL_PARAMS_HH

#include <dumux/implicit/3p3c/3p3cindices.hh>
#include <dumux/material/spatialparams/implicitspatialparams.hh>
#include <dumux/material/fluidmatrixinteractions/3p/parkerVanGen3p.hh>
#include <dumux/material/fluidmatrixinteractions/3p/parkerVanGen3pparams.hh>

namespace Dumux
{

//forward declaration
template<class TypeTag>
class InfiltrationSpatialParams;

namespace Properties
{
// The spatial parameters TypeTag
NEW_TYPE_TAG(InfiltrationSpatialParams);

// Set the spatial parameters
SET_TYPE_PROP(InfiltrationSpatialParams, SpatialParams, Dumux::InfiltrationSpatialParams<TypeTag>);

// Set the material Law
SET_PROP(InfiltrationSpatialParams, MaterialLaw)
{
 private:
    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
 public:
    // define the material law
    typedef ParkerVanGen3P<Scalar> type;
};
}

/*!
 * \ingroup ThreePThreeCModel
 *
 * \brief Definition of the spatial parameters for the infiltration problem
 */
template<class TypeTag>
class InfiltrationSpatialParams : public ImplicitSpatialParams<TypeTag>
{
    typedef ImplicitSpatialParams<TypeTag> ParentType;

    typedef typename GET_PROP_TYPE(TypeTag, Grid) Grid;
    typedef typename GET_PROP_TYPE(TypeTag, GridView) GridView;
    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
    typedef typename Grid::ctype CoordScalar;
    enum {
        dimWorld=GridView::dimensionworld
    };

    typedef Dune::FieldVector<CoordScalar,dimWorld> GlobalPosition;

    typedef typename GET_PROP_TYPE(TypeTag, FVElementGeometry) FVElementGeometry;
    typedef typename GridView::template Codim<0>::Entity Element;


public:
    //get the material law from the property system
    typedef typename GET_PROP_TYPE(TypeTag, MaterialLaw) MaterialLaw;
    typedef typename MaterialLaw::Params MaterialLawParams;

    /*!
     * \brief The constructor
     *
     * \param gv The grid view
     */
    InfiltrationSpatialParams(const GridView &gv)
        : ParentType(gv)
    {
        // intrinsic permeabilities
        fineK_ = GET_RUNTIME_PARAM(TypeTag, Scalar, permeability);
        coarseK_ = GET_RUNTIME_PARAM(TypeTag, Scalar, permeability);

        // porosities
        porosity_ = GET_RUNTIME_PARAM(TypeTag, Scalar, porosity);

        // residual saturations
        materialParams_.setSwr(0.12);
        materialParams_.setSwrx(0.12);
        materialParams_.setSnr(0.07);
        materialParams_.setSgr(0.03);

        // parameters for the 3phase van Genuchten law
        materialParams_.setVgAlpha(GET_RUNTIME_PARAM(TypeTag, Scalar, vanGenuchtenAlpha));
        materialParams_.setVgn(GET_RUNTIME_PARAM(TypeTag, Scalar, vanGenuchtenN));
        materialParams_.setKrRegardsSnr(false);

        // parameters for adsorption
        materialParams_.setKdNAPL(0.);
        materialParams_.setRhoBulk(1500.);
    }

    ~InfiltrationSpatialParams()
    {}

    /*!
     * \brief Apply the intrinsic permeability tensor to a pressure
     *        potential gradient.
     *
     * \param element The current finite element
     * \param fvElemGeom The current finite volume geometry of the element
     * \param scvIdx The index of the sub-control volume
     */
    const Scalar intrinsicPermeability(const Element &element,
                                       const FVElementGeometry &fvElemGeom,
                                       int scvIdx) const
    {
        const GlobalPosition &pos = fvElemGeom.subContVol[scvIdx].global;
        if (isFineMaterial_(pos))
            return fineK_;
        return coarseK_;
    }

    /*!
     * \brief Define the porosity \f$[-]\f$ of the spatial parameters
     *
     * \param element The finite element
     * \param fvElemGeom The finite volume geometry
     * \param scvIdx The local index of the sub-control volume where
     *                    the porosity needs to be defined
     */
    double porosity(const Element &element,
                    const FVElementGeometry &fvElemGeom,
                    int scvIdx) const
    {
        return porosity_;
    }


    /*!
     * \brief return the parameter object for the material law which depends on the position
     *
     * \param element The current finite element
     * \param fvElemGeom The current finite volume geometry of the element
     * \param scvIdx The index of the sub-control volume
     */
    const MaterialLawParams& materialLawParams(const Element &element,
                                               const FVElementGeometry &fvElemGeom,
                                               int scvIdx) const
    {
        return materialParams_;
    }

    const MaterialLawParams& materialLawParams() const
    {
        return materialParams_;
    }
private:
    bool isFineMaterial_(const GlobalPosition &pos) const
    { return
            70. <= pos[0] && pos[0] <= 85. &&
            7.0 <= pos[1] && pos[1] <= 7.50;
    };

    Scalar fineK_;
    Scalar coarseK_;

    Scalar porosity_;

    MaterialLawParams materialParams_;
};

}

#endif
