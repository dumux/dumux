/*****************************************************************************
 *   Copyright (C) 2008-2009 by Markus Wolff                                 *
 *   Copyright (C) 2010 by Klaus Mosthaf                                     *
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
#ifndef BUCKLEYLEVERETT_SPATIALPARAMETERS_HH
#define BUCKLEYLEVERETT_SPATIALPARAMETERS_HH

#include <dumux/material/spatialparameters/fvspatialparameters.hh>
#include <dumux/material/fluidmatrixinteractions/2p/linearmaterial.hh>
#include <dumux/material/fluidmatrixinteractions/2p/regularizedbrookscorey.hh>
#include <dumux/material/fluidmatrixinteractions/2p/efftoabslaw.hh>

namespace Dumux
{

//forward declaration
template<class TypeTag>
class BuckleyLeverettSpatialParams;

namespace Properties
{
// The spatial parameters TypeTag
NEW_TYPE_TAG(BuckleyLeverettSpatialParams);

// Set the spatial parameters
SET_TYPE_PROP(BuckleyLeverettSpatialParams, SpatialParameters, Dumux::BuckleyLeverettSpatialParams<TypeTag>);

// Set the material law
SET_PROP(BuckleyLeverettSpatialParams, MaterialLaw)
{
private:
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(Scalar)) Scalar;
    typedef RegularizedBrooksCorey<Scalar> RawMaterialLaw;
public:
    typedef EffToAbsLaw<RawMaterialLaw> type;
};
}

/** \todo Please doc me! */

template<class TypeTag>
class BuckleyLeverettSpatialParams: public FVSpatialParameters<TypeTag>
{
    typedef FVSpatialParameters<TypeTag> ParentType;
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(Grid)) Grid;
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(GridView)) GridView;
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(Scalar)) Scalar;
    typedef typename Grid::ctype CoordScalar;

    enum
        {dim=Grid::dimension, dimWorld=Grid::dimensionworld, numEq=1};
    typedef typename Grid::Traits::template Codim<0>::Entity Element;

    typedef Dune::FieldVector<CoordScalar, dimWorld> GlobalPosition;
    typedef Dune::FieldMatrix<Scalar,dim,dim> FieldMatrix;

public:
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(MaterialLaw)) MaterialLaw;
    typedef typename MaterialLaw::Params MaterialLawParams;

    Scalar intrinsicPermeability(const Element &element) const
    {
        return constPermeability_;
    }

    Scalar porosity(const Element &element) const
    {
        return porosity_;
    }


    // return the parameter object for the Brooks-Corey material law which depends on the position
    const MaterialLawParams& materialLawParams(const Element &element) const
    {
            return materialLawParams_;
    }


    BuckleyLeverettSpatialParams(const GridView& gridView)
    :ParentType(gridView)
    {
        Dumux::InterfaceSoilProperties interfaceSoilProps("interface_BL.xml");

        // residual saturations
        materialLawParams_.setSwr(interfaceSoilProps.ISP_ResidualSaturationWetting);
        materialLawParams_.setSnr(interfaceSoilProps.ISP_ResidualSaturationNonWetting);

        porosity_ = interfaceSoilProps.ISP_Porosity;

        // parameters for the Brooks-Corey Law
        // entry pressures
        materialLawParams_.setPe(interfaceSoilProps.ISP_BrooksCoreyEntryPressure);
        // Brooks-Corey shape parameters
        materialLawParams_.setLambda(interfaceSoilProps.ISP_BrooksCoreyLambda);

        constPermeability_ = interfaceSoilProps.ISP_Permeability;

        // parameters for the linear
        // entry pressures function
//        materialLawParams_.setEntryPC(0);
//        materialLawParams_.setMaxPC(0);
//        materialLawParams_.setEntryPC(2);
//        materialLawParams_.setMaxPC(3);


    }

private:
    MaterialLawParams materialLawParams_;
    Scalar constPermeability_;
    Scalar porosity_;

};

} // end namespace
#endif
