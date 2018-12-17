// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*****************************************************************************
 *   See the file COPYING for full copying permissions.                      *
 *                                                                           *
 *   This program is free software: you can redistribute it and/or modify    *
 *   it under the terms of the GNU General Public License as published by    *
 *   the Free Software Foundation, either version 3 of the License, or       *
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
 * \ingroup SequentialTwoPTests
 * \brief spatial parameters for the explicit transport test
 */
#ifndef TEST_TRANSPORT_SPATIALPARAMS_HH
#define TEST_TRANSPORT_SPATIALPARAMS_HH

#include <dumux/material/spatialparams/sequentialfv.hh>
#include <dumux/material/fluidmatrixinteractions/2p/linearmaterial.hh>
#include <dumux/material/fluidmatrixinteractions/2p/regularizedbrookscorey.hh>
#include <dumux/material/fluidmatrixinteractions/2p/efftoabslaw.hh>

namespace Dumux
{
/*!
 * \ingroup SequentialTwoPTests
 * \brief spatial parameters for the explicit transport test
 */
//forward declaration
template<class TypeTag>
class TestTransportSpatialParams;

namespace Properties
{
// The spatial parameters TypeTag
NEW_TYPE_TAG(TestTransportSpatialParams);

// Set the spatial parameters
SET_TYPE_PROP(TestTransportSpatialParams, SpatialParams, TestTransportSpatialParams<TypeTag>);

// Set the material law
template<class TypeTag>
struct MaterialLaw<TypeTag, TTag::TestTransportSpatialParams>
{
private:
    using Scalar = typename GET_PROP_TYPE(TypeTag, Scalar);
    using RawMaterialLaw = LinearMaterial<Scalar>;
public:
    using type = EffToAbsLaw<RawMaterialLaw>;
};
}

/*!
 * \ingroup IMPETtests
 * \brief spatial parameters for the explicit transport test
 */
template<class TypeTag>
class TestTransportSpatialParams: public SequentialFVSpatialParams<TypeTag>
{
    using ParentType = SequentialFVSpatialParams<TypeTag>;
    using Problem = typename GET_PROP_TYPE(TypeTag, Problem);
    using Grid = typename GET_PROP_TYPE(TypeTag, Grid);
    using GridView = typename GET_PROP_TYPE(TypeTag, GridView);
    using Scalar = typename GET_PROP_TYPE(TypeTag, Scalar);

    using Element = typename Grid::Traits::template Codim<0>::Entity;

public:
    using MaterialLaw = typename GET_PROP_TYPE(TypeTag, MaterialLaw);
    using MaterialLawParams = typename MaterialLaw::Params;

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


    TestTransportSpatialParams(const Problem& problem)
    : ParentType(problem)
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
