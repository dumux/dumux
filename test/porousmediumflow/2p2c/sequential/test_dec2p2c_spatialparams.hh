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
 * \ingroup SequentialTwoPTwoCTests
 * \brief spatial parameters for the sequential 2p2c test
 */
#ifndef TEST_2P2C_SPATIALPARAMS_HH
#define TEST_2P2C_SPATIALPARAMS_HH

#include <dumux/porousmediumflow/2p2c/sequential/properties.hh>
#include <dumux/material/spatialparams/sequentialfv.hh>
#include <dumux/material/fluidmatrixinteractions/2p/linearmaterial.hh>
#include <dumux/material/fluidmatrixinteractions/2p/efftoabslaw.hh>

namespace Dumux
{
//forward declaration
template<class TypeTag>
class Test2P2CSpatialParams;

namespace Properties
{
// The spatial parameters TypeTag
NEW_TYPE_TAG(Test2P2CSpatialParams);

// Set the spatial parameters
SET_TYPE_PROP(Test2P2CSpatialParams, SpatialParams, Test2P2CSpatialParams<TypeTag>);

// Set the material law
SET_PROP(Test2P2CSpatialParams, MaterialLaw)
{
private:
    using Scalar = typename GET_PROP_TYPE(TypeTag, Scalar);
    using RawMaterialLaw = LinearMaterial<Scalar>;
public:
    using type = EffToAbsLaw<RawMaterialLaw>;
};
}

/*!
 * \ingroup SequentialTwoPTwoCTests
 * \brief spatial parameters for the sequential 2p2c test
 */
template<class TypeTag>
class Test2P2CSpatialParams : public SequentialFVSpatialParams<TypeTag>
{
    using GridView = typename GET_PROP_TYPE(TypeTag, GridView);
    using Scalar = typename GET_PROP_TYPE(TypeTag, Scalar);

    enum { dim = GridView::dimension };
    using Element = typename GridView::Traits::template Codim<0>::Entity;

    using FieldMatrix = Dune::FieldMatrix<Scalar, dim, dim>;

public:
    using MaterialLaw = typename GET_PROP_TYPE(TypeTag, MaterialLaw);
    using MaterialLawParams = typename MaterialLaw::Params;

    const FieldMatrix& intrinsicPermeability (const Element& element) const
    {
        return constPermeability_;
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


    Test2P2CSpatialParams(const GridView& gridView) : SequentialFVSpatialParams<TypeTag>(gridView),
            constPermeability_(0)
    {
        // residual saturations
        materialLawParams_.setSwr(0);
        materialLawParams_.setSnr(0);

        materialLawParams_.setEntryPc(0);
        materialLawParams_.setMaxPc(10000);

        for(int i = 0; i < dim; i++)
        {
            constPermeability_[i][i] = 1e-12;
        }
    }

private:
    MaterialLawParams materialLawParams_;
    FieldMatrix constPermeability_;

};

} // end namespace
#endif
