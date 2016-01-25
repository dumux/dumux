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
 * \brief spatial parameters for the 2p 3d test
 */
#ifndef TEST_3D2P_SPATIALPARAMETERS_HH
#define TEST_3D2P_SPATIALPARAMETERS_HH

#include <dumux/material/spatialparams/fv.hh>
#include <dumux/material/fluidmatrixinteractions/2p/regularizedbrookscorey.hh>
#include <dumux/material/fluidmatrixinteractions/2p/linearmaterial.hh>
#include <dumux/material/fluidmatrixinteractions/2p/efftoabslaw.hh>


namespace Dumux
{

//forward declaration
template<class TypeTag>
class Test3d2pSpatialParams;

namespace Properties
{
// The spatial parameters TypeTag
NEW_TYPE_TAG(Test3d2pSpatialParams);

// Set the spatial parameters
SET_TYPE_PROP(Test3d2pSpatialParams, SpatialParams, Dumux::Test3d2pSpatialParams<TypeTag>);

// Set the material law
SET_PROP(Test3d2pSpatialParams, MaterialLaw)
{
private:
    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
    typedef RegularizedBrooksCorey<Scalar> RawMaterialLaw;
//    typedef LinearMaterial<Scalar> RawMaterialLaw;
public:
    typedef EffToAbsLaw<RawMaterialLaw> type;
};
}

/*!
 *
 * \ingroup IMPESTests
 * \brief spatial parameters for the 2p test using MPFAL 3D method
 */
template<class TypeTag>
class Test3d2pSpatialParams: public FVSpatialParams<TypeTag>
{
    typedef FVSpatialParams<TypeTag> ParentType;
    typedef typename GET_PROP_TYPE(TypeTag, Grid) Grid;
    typedef typename GET_PROP_TYPE(TypeTag, GridView) GridView;
    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
    typedef typename Grid::ctype CoordScalar;

    enum
        {dim=Grid::dimension, dimWorld=Grid::dimensionworld, numEq=1};
    typedef typename Grid::Traits::template Codim<0>::Entity Element;

    typedef Dune::FieldVector<CoordScalar, dimWorld> GlobalPosition;
    typedef Dune::FieldVector<CoordScalar, dim> LocalPosition;
    typedef Dune::FieldMatrix<Scalar,dim,dim> FieldMatrix;

public:
    typedef typename GET_PROP_TYPE(TypeTag, MaterialLaw) MaterialLaw;
    typedef typename MaterialLaw::Params MaterialLawParams;

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


    Test3d2pSpatialParams(const GridView& gridView)
    : ParentType(gridView), constPermeability_(0)
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
