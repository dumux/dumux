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
/**
 * \file
 *
 * \brief Definition of the spatial params properties for the Forchheimer problem
 *
 */

#ifndef DUMUX_FORCHHEIMER_SPATIAL_PARAMS_HH
#define DUMUX_FORCHHEIMER_SPATIAL_PARAMS_HH

#include <dumux/material/spatialparams/implicit.hh>
#include <dumux/material/fluidmatrixinteractions/2p/regularizedlinearmaterial.hh>
#include <dumux/material/fluidmatrixinteractions/2p/efftoabslaw.hh>

#include <dumux/porousmediumflow/mpnc/implicit/model.hh>
#include <dumux/material/fluidmatrixinteractions/mp/mplinearmaterial.hh>
#include <dumux/material/fluidmatrixinteractions/mp/2padapter.hh>

namespace Dumux
{

//forward declaration
template<class TypeTag>
class ForchheimerSpatialParams;

namespace Properties
{
// The spatial parameters TypeTag
NEW_TYPE_TAG(ForchheimerSpatialParams);

// Set the spatial parameters
SET_TYPE_PROP(ForchheimerSpatialParams, SpatialParams, ForchheimerSpatialParams<TypeTag>);

// Set the material Law
SET_PROP(ForchheimerSpatialParams, MaterialLaw)
{
private:
    typedef typename GET_PROP_TYPE(TypeTag, FluidSystem) FluidSystem;
    enum {wPhaseIdx = FluidSystem::wPhaseIdx};
    // define the material law
    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
    typedef RegularizedLinearMaterial<Scalar> EffMaterialLaw;
    typedef EffToAbsLaw<EffMaterialLaw> TwoPMaterialLaw;
public:
    typedef TwoPAdapter<wPhaseIdx, TwoPMaterialLaw> type;
};
}

/**
 * \ingroup MPNCModel
 * \ingroup ImplicitTestProblems
 * \brief Definition of the spatial params properties for the Forchheimer problem
 *
 */
template<class TypeTag>
class ForchheimerSpatialParams : public ImplicitSpatialParams<TypeTag>
{
    typedef ImplicitSpatialParams<TypeTag> ParentType;
    typedef typename GET_PROP_TYPE(TypeTag, Grid) Grid;
    typedef typename GET_PROP_TYPE(TypeTag, GridView) GridView;
    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
    typedef typename GET_PROP_TYPE(TypeTag, SolutionVector) SolutionVector;
    typedef typename GET_PROP_TYPE(TypeTag, FVElementGeometry) FVElementGeometry;
    typedef typename GET_PROP_TYPE(TypeTag, FluidSystem) FluidSystem;
    typedef typename Grid::ctype CoordScalar;

    enum {dimWorld=GridView::dimensionworld};
    enum {wPhaseIdx = FluidSystem::wPhaseIdx};

    typedef typename GridView::template Codim<0>::Entity Element;
    typedef Dune::FieldVector<CoordScalar,dimWorld> GlobalPosition;

public:
    typedef typename GET_PROP_TYPE(TypeTag, MaterialLaw) MaterialLaw;
    typedef typename MaterialLaw::Params MaterialLawParams;

    ForchheimerSpatialParams(const GridView &gridView)
        : ParentType(gridView)
    {
        // intrinsic permeabilities
        K_ = 1e-12;

        // the porosity
        porosity_ = 0.3;

        // residual saturations
        materialParams_.setSwr(0.0);
        materialParams_.setSnr(0.0);

        // parameters for the linear law, i.e. minimum and maximum
        // pressures
        materialParams_.setEntryPc(0.0);
        materialParams_.setMaxPc(0.0);
    }

    ~ForchheimerSpatialParams()
    {}

    /*!
     * \brief Update the spatial parameters with the flow solution
     *        after a timestep.
     *
     * \param globalSol The current solution vector
     */
    void update(const SolutionVector &globalSol)
    {
    }

    /*!
     * \brief Returns the intrinsic permeability tensor.
     *
     * \param element       The current finite element
     * \param fvGeometry    The current finite volume geometry of the element
     * \param scvIdx        The index sub-control volume where the
     *                      intrinsic permeability is given.
     */
    Scalar intrinsicPermeability(const Element &element,
                                 const FVElementGeometry &fvGeometry,
                                 const unsigned int scvIdx) const
    {
        return K_;
    }

    /*!
     * \brief Define the porosity \f$[-]\f$ of the soil
     *
     * \param element     The finite element
     * \param fvGeometry  The finite volume geometry
     * \param scvIdx      The local index of the sub-control volume where
     *                    the porosity needs to be defined
     */
    double porosity(const Element &element,
                    const FVElementGeometry &fvGeometry,
                    const unsigned int scvIdx) const
    {
        return porosity_;
    }

    /*!
     * \brief Function for defining the parameters needed by constitutive relationships (kr-sw, pc-sw, etc.).
     *
     * \param pos The global position of the sub-control volume.
     * \return the material parameters object
     */
    const MaterialLawParams& materialLawParamsAtPos(const GlobalPosition &pos) const
    {
        return materialParams_;
    }


    /*!
     * \brief Apply the Forchheimer coefficient for inertial forces
     *        calculation.
     *
     *        Source: Ward, J.C. 1964 Turbulent flow in porous media. ASCE J. Hydraul. Div 90.
     *        Actually the Forchheimer coefficient is also a function of the dimensions of the
     *        porous medium. Taking it as a constant is only a first approximation
     *        (Nield, Bejan, Convection in porous media, 2006, p. 10)
     *
     * \param element The current finite element
     * \param fvGeometry The current finite volume geometry of the element
     * \param scvIdx The index sub-control volume face where the
     *                      intrinsic velocity ought to be calculated.
     *
     */
    Scalar forchCoeff(const Element &element,
                    const FVElementGeometry &fvGeometry,
                    const unsigned int scvIdx) const
    {
        // If there are better measures / estimates / values available than this default number:
        // here is the place to implement it
        return ParentType::forchCoeff(element,
                                      fvGeometry,
                                      scvIdx);
    }

    Scalar K_;
    Scalar porosity_;
    MaterialLawParams materialParams_;
};

}

#endif
