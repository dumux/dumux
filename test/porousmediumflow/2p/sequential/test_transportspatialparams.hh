// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*****************************************************************************
 *   See the file COPYING for full copying permissions.                      *
 *                                                                           *
 *   This program is free software: you can redistribute it and/or modify    *
 *   it under the terms of the GNU General Public License as published by    *
 *   the Free Software Foundation, either version 2 of the License, or       *
 *   (at your option) any later version.                                     *
 *                                                                           *
 *   This program is distributed in the hope that it will be useful,         *
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of          *
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the            *
 *   GNU General Public License for more details.                            *
 *                                                                           *
 *   You should have received a copy of the GNU General Public License       *
 *   along with this program.  If not, see <http://www.gnu.org/licenses/>.   *
 *****************************************************************************/
/*!
 * \file
 *
 * \brief spatial parameters for the explicit transport test
 */
#ifndef TEST_TRANSPORT_SPATIALPARAMS_HH
#define TEST_TRANSPORT_SPATIALPARAMS_HH

#include <dumux/material/spatialparams/fvspatialparams.hh>
#include <dumux/material/fluidmatrixinteractions/2p/linearmaterial.hh>
#include <dumux/material/fluidmatrixinteractions/2p/regularizedbrookscorey.hh>
#include <dumux/material/fluidmatrixinteractions/2p/efftoabslaw.hh>

namespace Dumux
{

//forward declaration
template<class TypeTag>
class TestTransportSpatialParams;

namespace Properties
{
// The spatial parameters TypeTag
NEW_TYPE_TAG(TestTransportSpatialParams);

// Set the spatial parameters
SET_TYPE_PROP(TestTransportSpatialParams, SpatialParams, Dumux::TestTransportSpatialParams<TypeTag>);

// Set the material law
SET_PROP(TestTransportSpatialParams, MaterialLaw)
{
private:
    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
    typedef LinearMaterial<Scalar> RawMaterialLaw;
public:
    typedef EffToAbsLaw<RawMaterialLaw> type;
};
}

/*!
 * \ingroup IMPETtests
 * \brief spatial parameters for the explicit transport test
 */
template<class TypeTag>
class TestTransportSpatialParams: public FVSpatialParams<TypeTag>
{
    typedef FVSpatialParams<TypeTag> ParentType;
    typedef typename GET_PROP_TYPE(TypeTag, Grid) Grid;
    typedef typename GET_PROP_TYPE(TypeTag, GridView) GridView;
    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;

    typedef typename Grid::Traits::template Codim<0>::Entity Element;

public:
    typedef typename GET_PROP_TYPE(TypeTag, MaterialLaw) MaterialLaw;
    typedef typename MaterialLaw::Params MaterialLawParams;

    Scalar intrinsicPermeability (const Element& element) const
    {
        return 1e-5;
    }

    double porosity(const Element& element) const
    {
        return 0.2;
    }


    // return the parameter object for the Brooks-Corey material law which depends on the position
    const MaterialLawParams& materialLawParams(const Element &element) const
    {
            return materialLawParams_;
    }


    TestTransportSpatialParams(const GridView& gridView)
    : ParentType(gridView)
    {
        // residual saturations
        materialLawParams_.setSwr(0.0);
        materialLawParams_.setSnr(0.0);

        // parameters for the linear entry pressures function
        materialLawParams_.setEntryPc(0);
        materialLawParams_.setMaxPc(0);
    }

private:
    MaterialLawParams materialLawParams_;
};

} // end namespace
#endif
