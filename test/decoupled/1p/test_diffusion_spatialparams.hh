// $Id: test_2p_spatialparams.hh 3456 2010-04-09 12:11:51Z mwolff $
/*****************************************************************************
 *   Copyright (C) 2008-2009 by Markus Wolff                                 *
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
/*!
 * \file
 *
 * \brief spatial parameters for the test problem for diffusion models.
 */
#ifndef TEST_DIFFUSION_SPATIALPARAMETERS_HH
#define TEST_DIFFUSION_SPATIALPARAMETERS_HH


#include <dumux/material/fluidmatrixinteractions/2p/linearmaterial.hh>
#include <dumux/material/fluidmatrixinteractions/2p/efftoabslaw.hh>

namespace Dumux
{

/*!
 * \ingroup IMPETtests
 * \brief spatial parameters for the test problem for diffusion models.
 */
template<class TypeTag>
class TestDiffusionSpatialParams
{
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(Grid)) Grid;
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(GridView)) GridView;
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(Scalar)) Scalar;
    typedef typename Grid::ctype CoordScalar;

    enum
        {dim=Grid::dimension, dimWorld=Grid::dimensionworld, numEq=1};
    typedef typename Grid::Traits::template Codim<0>::Entity Element;

    typedef Dune::FieldVector<CoordScalar, dimWorld> GlobalPosition;
    typedef Dune::FieldVector<CoordScalar, dim> LocalPosition;
    typedef Dune::FieldMatrix<Scalar,dim,dim> FieldMatrix;

    typedef LinearMaterial<Scalar>                        RawMaterialLaw;
public:
    typedef EffToAbsLaw<RawMaterialLaw>               MaterialLaw;
    typedef typename MaterialLaw::Params MaterialLawParams;

    void update (Scalar saturationW, const Element& element)
    {

    }

    const FieldMatrix& intrinsicPermeability (const GlobalPosition& globalPos, const Element& element) const
    {
        double rt = globalPos[0]*globalPos[0]+globalPos[1]*globalPos[1];
        permeability_[0][0] = (delta_*globalPos[0]*globalPos[0] + globalPos[1]*globalPos[1])/rt;
        permeability_[0][1] = permeability_[1][0] = -(1.0 - delta_)*globalPos[0]*globalPos[1]/rt;
        permeability_[1][1] = (globalPos[0]*globalPos[0] + delta_*globalPos[1]*globalPos[1])/rt;

        return permeability_;
    }

    double porosity(const GlobalPosition& globalPos, const Element& element) const
    {
        return 0.2;
    }


    // return the parameter object for the Brooks-Corey material law which depends on the position
    const MaterialLawParams& materialLawParams(const GlobalPosition& globalPos, const Element &element) const
    {
            return materialLawParams_;
    }

    void setDelta(const double delta)
    {
        delta_ = delta;
    }

    TestDiffusionSpatialParams(const GridView& gridView)
    : permeability_(0)
    {
        // residual saturations
        materialLawParams_.setSwr(0.0);
        materialLawParams_.setSnr(0.0);

        // parameters for the linear entry pressure function
        materialLawParams_.setEntryPC(0);
        materialLawParams_.setMaxPC(0);
    }

private:
    MaterialLawParams materialLawParams_;
    mutable FieldMatrix permeability_;
    double delta_;
};

} // end namespace
#endif
