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
 * \brief spatial parameters for the test problem for diffusion models.
 */
#ifndef TEST1_FVCA6_SPATIALPARAMETERS_HH
#define TEST1_FVCA6_SPATIALPARAMETERS_HH

#include <dumux/material/spatialparams/sequentialfv.hh>
#include <dumux/material/fluidmatrixinteractions/2p/linearmaterial.hh>
#include <dumux/material/fluidmatrixinteractions/2p/efftoabslaw.hh>

namespace Dumux
{

//forward declaration
template<class TypeTag>
class TestDiffusionSpatialParams3d;

namespace Properties
{
// The spatial parameters TypeTag
NEW_TYPE_TAG(TestDiffusionSpatialParams3d);

// Set the spatial parameters
SET_TYPE_PROP(TestDiffusionSpatialParams3d, SpatialParams, TestDiffusionSpatialParams3d<TypeTag>);

// Set the material law
SET_PROP(TestDiffusionSpatialParams3d, MaterialLaw)
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
 * \brief spatial parameters for the test problem for diffusion models.
 */
template<class TypeTag>
class TestDiffusionSpatialParams3d: public SequentialFVSpatialParams<TypeTag>
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

    using GlobalPosition = Dune::FieldVector<CoordScalar, dimWorld>;
    using FieldMatrix = Dune::FieldMatrix<Scalar, dim, dim>;

public:
    using MaterialLaw = typename GET_PROP_TYPE(TypeTag, MaterialLaw);
    using MaterialLawParams = typename MaterialLaw::Params;

    const FieldMatrix& intrinsicPermeabilityAtPos (const GlobalPosition& globalPos) const
    {
        return permeability_;
    }

    double porosity(const Element& element) const
    {
        return 1.0;
    }


    // return the parameter object for the Brooks-Corey material law which depends on the position
    const MaterialLawParams& materialLawParams(const Element &element) const
    {
            return materialLawParams_;
    }

    TestDiffusionSpatialParams3d(const Problem& problem)
    : ParentType(problem), permeability_(0)
    {
        // residual saturations
        materialLawParams_.setSwr(0.0);
        materialLawParams_.setSnr(0.0);

        // parameters for the linear entry pressure function
        materialLawParams_.setEntryPc(0);
        materialLawParams_.setMaxPc(0);

        // permeability values
        permeability_[0][0] = permeability_[1][1] = permeability_[2][2] = 1.0;
        permeability_[0][1] = permeability_[1][0] = permeability_[1][2] = permeability_[2][1] = 0.5;
    }

private:
    MaterialLawParams materialLawParams_;
    FieldMatrix permeability_;
};

} // end namespace
#endif
