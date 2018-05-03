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
 * \ingroup SequentialTwoPTests
 * \brief spatial parameters for the sequential 2p test
 */
#ifndef TEST_IMPES_ADAPTIVE_SPATIALPARAMS_HH
#define TEST_IMPES_ADAPTIVE_SPATIALPARAMS_HH

#include <dumux/material/spatialparams/sequentialfv.hh>
#include <dumux/material/fluidmatrixinteractions/2p/linearmaterial.hh>
#include <dumux/material/fluidmatrixinteractions/2p/regularizedbrookscorey.hh>
#include <dumux/material/fluidmatrixinteractions/2p/efftoabslaw.hh>

namespace Dumux
{

//forward declaration
template<class TypeTag>
class TestIMPESAdaptiveSpatialParams;

namespace Properties
{
// The spatial parameters TypeTag
NEW_TYPE_TAG(TestIMPESAdaptiveSpatialParams);

// Set the spatial parameters
SET_TYPE_PROP(TestIMPESAdaptiveSpatialParams, SpatialParams, TestIMPESAdaptiveSpatialParams<TypeTag>);

// Set the material law
SET_PROP(TestIMPESAdaptiveSpatialParams, MaterialLaw)
{
private:
    using Scalar = typename GET_PROP_TYPE(TypeTag, Scalar);
    using RawMaterialLaw = RegularizedBrooksCorey<Scalar>;
public:
    using type = EffToAbsLaw<RawMaterialLaw>;
};
}

/*!
 * \ingroup SequentialTwoPTests
 * \brief spatial parameters for the sequential 2p test
 */
template<class TypeTag>
class TestIMPESAdaptiveSpatialParams: public SequentialFVSpatialParams<TypeTag>
{
    using Problem = typename GET_PROP_TYPE(TypeTag, Problem);
    using Grid = typename GET_PROP_TYPE(TypeTag, Grid);
    using GridView = typename GET_PROP_TYPE(TypeTag, GridView);
    using Scalar = typename GET_PROP_TYPE(TypeTag, Scalar);
    using ParentType = SequentialFVSpatialParams<TypeTag>;
    using CoordScalar = typename Grid::ctype;

    enum
        {dimWorld=Grid::dimensionworld};
    using Element = typename Grid::Traits::template Codim<0>::Entity;

    using GlobalPosition = typename Element::Geometry::GlobalCoordinate;

public:
    using MaterialLaw = typename GET_PROP_TYPE(TypeTag, MaterialLaw);
    using MaterialLawParams = typename MaterialLaw::Params;


    Scalar intrinsicPermeability (const Element& element) const
    {
        return 1.0e-7;
    }

    double porosity(const Element& element) const
    {
        return 0.2;
    }


    // return the parameter object for the Brooks-Corey material law which depends on the position
//    const MaterialLawParams& materialLawParamsAtPos(const GlobalPosition& globalPos) const
    const MaterialLawParams& materialLawParams(const Element& element) const
    {
            return materialLawParams_;
    }


    TestIMPESAdaptiveSpatialParams(const Problem& problem)
    : ParentType(problem)
    {
        // residual saturations
        materialLawParams_.setSwr(0.2);
        materialLawParams_.setSnr(0.2);

//        // parameters for the Brooks-Corey Law
//        // entry pressures
        materialLawParams_.setPe(0);
//        // Brooks-Corey shape parameters
        materialLawParams_.setLambda(2);

        // parameters for the linear
        // entry pressures function
//        materialLawParams_.setEntryPc(0);
//        materialLawParams_.setMaxPc(0);
    }

private:
    MaterialLawParams materialLawParams_;
};

} // end namespace
#endif
