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
 * \brief Definition of the spatial parameters for the non-isothermal Richards problems.
 */
#ifndef DUMUX_RICHARDSNI_SPATIAL_PARAMS_HH
#define DUMUX_RICHARDSNI_SPATIAL_PARAMS_HH

#include <dumux/porousmediumflow/richards/implicit/model.hh>
#include <dumux/material/fluidmatrixinteractions/2p/efftoabslaw.hh>
#include <dumux/material/fluidmatrixinteractions/2p/linearmaterial.hh>
#include <dumux/material/fluidmatrixinteractions/2p/regularizedvangenuchten.hh>
#include <dumux/material/spatialparams/implicit.hh>

namespace Dumux
{

/*!
 * \ingroup RichardsModel
 * \ingroup ImplicitTestProblems
 *
 * \brief Definition of the spatial parameters for the RichardsNI problems.
 */

//forward declaration
template<class TypeTag>
class RichardsNISpatialParams;

namespace Properties
{
// The spatial parameters TypeTag
NEW_TYPE_TAG(RichardsNISpatialParams);

// Set the spatial parameters
SET_TYPE_PROP(RichardsNISpatialParams, SpatialParams, RichardsNISpatialParams<TypeTag>);

// Set the material law
SET_PROP(RichardsNISpatialParams, MaterialLaw)
{
private:
    // define the material law which is parameterized by effective
    // saturations
    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
    typedef RegularizedVanGenuchten<Scalar> EffectiveLaw;
public:
    // define the material law parameterized by absolute saturations
    typedef EffToAbsLaw<EffectiveLaw> type;
};
}


template<class TypeTag>
class RichardsNISpatialParams : public ImplicitSpatialParams<TypeTag>
{
    typedef ImplicitSpatialParams<TypeTag> ParentType;
    typedef typename GET_PROP_TYPE(TypeTag, Grid) Grid;
    typedef typename GET_PROP_TYPE(TypeTag, GridView) GridView;
    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;

    typedef typename GET_PROP_TYPE(TypeTag, SolutionVector) SolutionVector;


    typedef typename GET_PROP_TYPE(TypeTag, FVElementGeometry) FVElementGeometry;
    typedef typename GridView::template Codim<0>::Entity Element;

    //typedef LinearMaterial<Scalar> EffMaterialLaw;
public:

    typedef typename GET_PROP_TYPE(TypeTag, MaterialLaw) MaterialLaw;
    typedef typename MaterialLaw::Params MaterialLawParams;

    RichardsNISpatialParams(const GridView &gridView)
        : ParentType(gridView)
    {
        permeability_ = 1e-10;
        porosity_ = 0.4;

        // heat conductivity of granite
        lambdaSolid_ = 2.8;

        // residual saturations

        // residual saturations
        materialParams_.setSwr(0.05);
        materialParams_.setSnr(0.0);

        // parameters for the Van Genuchten law
        // alpha and n

        materialParams_.setVgAlpha(0.0037);
        materialParams_.setVgn(4.7);
    }

    ~RichardsNISpatialParams()
    {}


    /*!
     * \brief Update the spatial parameters with the flow solution
     *        after a timestep.
     *
     * \param globalSolution the global solution vector
     */
    void update(const SolutionVector &globalSolution)
    {
    }

    /*!
     * \brief Define the intrinsic permeability \f$\mathrm{[m^2]}\f$.
     *
     * \param element The current finite element
     * \param fvGeometry The current finite volume geometry of the element
     * \param scvIdx The index of the sub-control volume
     */
    const Scalar intrinsicPermeability(const Element &element,
            const FVElementGeometry &fvGeometry,
            const int scvIdx) const
    {
        return permeability_;
    }

    /*!
     * \brief Define the porosity \f$\mathrm{[-]}\f$.
     *
     * \param element The finite element
     * \param fvGeometry The finite volume geometry
     * \param scvIdx The local index of the sub-control volume where
     */
    double porosity(const Element &element,
            const FVElementGeometry &fvGeometry,
            const int scvIdx) const
    {
        return porosity_;
    }

        /*!
     * \brief return the parameter object for the Brooks-Corey material law which depends on the position
     *
     * \param element The current finite element
     * \param fvGeometry The current finite volume geometry of the element
     * \param scvIdx The index of the sub-control volume
     */
    const MaterialLawParams& materialLawParams(const Element &element,
            const FVElementGeometry &fvGeometry,
            const int scvIdx) const
    {
        return materialParams_;
    }

    /*!
     * \brief Returns the heat capacity \f$[J / (kg K)]\f$ of the rock matrix.
     *
     * This is only required for non-isothermal models.
     *
     * \param element The finite element
     * \param fvGeometry The finite volume geometry
     * \param scvIdx The local index of the sub-control volume
     */
    Scalar solidHeatCapacity(const Element &element,
                             const FVElementGeometry &fvGeometry,
                             const int scvIdx) const
    {
        return 790; // specific heat capacity of granite [J / (kg K)]
    }

    /*!
     * \brief Returns the mass density \f$[kg / m^3]\f$ of the rock matrix.
     *
     * This is only required for non-isothermal models.
     *
     * \param element The finite element
     * \param fvGeometry The finite volume geometry
     * \param scvIdx The local index of the sub-control volume
     */
    Scalar solidDensity(const Element &element,
                        const FVElementGeometry &fvGeometry,
                        const int scvIdx) const
    {
        return 2700; // density of granite [kg/m^3]
    }

    /*!
     * \brief Returns the thermal conductivity \f$\mathrm{[W/(m K)]}\f$ of the porous material.
     *
     * \param element The finite element
     * \param fvGeometry The finite volume geometry
     * \param scvIdx The local index of the sub-control volume where
     *                    the heat capacity needs to be defined
     */
    Scalar solidThermalConductivity(const Element &element,
                                    const FVElementGeometry &fvGeometry,
                                    const int scvIdx) const
    {
        return lambdaSolid_;
    }


private:

    MaterialLawParams materialParams_;
    Scalar permeability_;
    Scalar porosity_;
    Scalar lambdaSolid_;
};

}

#endif
