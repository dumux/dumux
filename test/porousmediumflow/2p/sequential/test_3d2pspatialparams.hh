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
 * \brief spatial parameters for the 2p 3d test
 */
#ifndef TEST_3D2P_SPATIALPARAMETERS_HH
#define TEST_3D2P_SPATIALPARAMETERS_HH

#include <dumux/material/spatialparams/sequentialfv.hh>
#include <dumux/material/fluidmatrixinteractions/2p/regularizedbrookscorey.hh>
#include <dumux/material/fluidmatrixinteractions/2p/linearmaterial.hh>
#include <dumux/material/fluidmatrixinteractions/2p/efftoabslaw.hh>


namespace Dumux
{
/*!
 * \ingroup SequentialTwoPTests
 * \brief spatial parameters for the 2p 3d test
 */
//forward declaration
template<class TypeTag>
class Test3d2pSpatialParams;

namespace Properties
{
// The spatial parameters TypeTag
NEW_TYPE_TAG(Test3d2pSpatialParams);

// Set the spatial parameters
SET_TYPE_PROP(Test3d2pSpatialParams, SpatialParams, Test3d2pSpatialParams<TypeTag>);

// Set the material law
template<class TypeTag>
struct MaterialLaw<TypeTag, TTag::Test3d2pSpatialParams>
{
private:
    using Scalar = typename GET_PROP_TYPE(TypeTag, Scalar);
    using RawMaterialLaw = RegularizedBrooksCorey<Scalar>;
//    using RawMaterialLaw = LinearMaterial<Scalar>;
public:
    using type = EffToAbsLaw<RawMaterialLaw>;
};
}

/*!
 *
 * \ingroup IMPESTests
 * \brief spatial parameters for the 2p test using MPFAL 3D method
 */
template<class TypeTag>
class Test3d2pSpatialParams: public SequentialFVSpatialParams<TypeTag>
{
    using ParentType = SequentialFVSpatialParams<TypeTag>;
    using Problem = typename GET_PROP_TYPE(TypeTag, Problem);
    using Grid = typename GET_PROP_TYPE(TypeTag, Grid);
    using GridView = typename GET_PROP_TYPE(TypeTag, GridView);
    using Scalar = typename GET_PROP_TYPE(TypeTag, Scalar);
    using CoordScalar = typename Grid::ctype;

    enum
        {dim=Grid::dimension, dimWorld=Grid::dimensionworld, numEq=1};
    using Element = typename Grid::Traits::template Codim<0>::Entity;

    using GlobalPosition = typename Element::Geometry::GlobalCoordinate;
    using LocalPosition = Dune::FieldVector<CoordScalar, dim>;
    using FieldMatrix = Dune::FieldMatrix<Scalar, dim, dim>;

public:
    using MaterialLaw = typename GET_PROP_TYPE(TypeTag, MaterialLaw);
    using MaterialLawParams = typename MaterialLaw::Params;

    void update (Scalar saturationW, const Element& element)
    {
    }

    const FieldMatrix& intrinsicPermeabilityAtPos (const GlobalPosition& globalPos) const
    {
        return constPermeability_;
    }

    double porosity(const Element& element) const
    {
#if PROBLEM == 1
        return 0.3;
#else
        return 0.2;
#endif
    }


    // return the parameter object for the Brooks-Corey material law which depends on the position
    const MaterialLawParams& materialLawParams(const Element &element) const
    {
            return materialLawParams_;
    }


    Test3d2pSpatialParams(const Problem& problem)
    : ParentType(problem), constPermeability_(0)
    {



//        // parameters for the Brooks-Corey Law
#if PROBLEM == 1
        // residual saturations
        materialLawParams_.setSwr(0.);
        materialLawParams_.setSnr(0.);

        // entry pressures
        materialLawParams_.setPe(5000);
#else
        // residual saturations
        materialLawParams_.setSwr(0.2);
        materialLawParams_.setSnr(0.2);

        // entry pressures
        materialLawParams_.setPe(0);
#endif

        // Brooks-Corey shape parameters
        materialLawParams_.setLambda(2);

#if PROBLEM == 1
        for(int i = 0; i < dim; i++)
        {
            constPermeability_[i][i] = 1e-10;
        }
#else
        for(int i = 0; i < dim; i++)
        {
            constPermeability_[i][i] = 1e-7;
        }
#endif
    }

private:
    MaterialLawParams materialLawParams_;
    FieldMatrix constPermeability_;

};

} // end namespace
#endif
