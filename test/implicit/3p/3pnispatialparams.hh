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
 * \brief Definition of the spatial parameters for the 3pni problems.
 */
#ifndef DUMUX_THREEPNI_SPATIAL_PARAMS_HH
#define DUMUX_THREEPNI_SPATIAL_PARAMS_HH

#include <dumux/implicit/3p3c/3p3cindices.hh>
#include <dumux/material/spatialparams/implicitspatialparams.hh>
#include <dumux/material/fluidmatrixinteractions/3p/parkervangen3p.hh>
#include <dumux/material/fluidmatrixinteractions/3p/parkervangen3pparams.hh>

namespace Dumux
{

/*!
 * \ingroup ThreePNIModel
 * \ingroup ImplicitTestProblems
 *
 * \brief Definition of the spatial parameters for the 3pni problems.
 */

//forward declaration
template<class TypeTag>
class ThreePNISpatialParams;

namespace Properties
{
// The spatial parameters TypeTag
NEW_TYPE_TAG(ThreePNISpatialParams);

// Set the spatial parameters
SET_TYPE_PROP(ThreePNISpatialParams, SpatialParams, Dumux::ThreePNISpatialParams<TypeTag>);

// Set the material Law
SET_TYPE_PROP(ThreePNISpatialParams, MaterialLaw, ParkerVanGen3P<typename GET_PROP_TYPE(TypeTag, Scalar)>);
}


template<class TypeTag>
class ThreePNISpatialParams : public ImplicitSpatialParams<TypeTag>
{
    typedef ImplicitSpatialParams<TypeTag> ParentType;
    typedef typename GET_PROP_TYPE(TypeTag, Grid) Grid;
    typedef typename GET_PROP_TYPE(TypeTag, GridView) GridView;
    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;

    typedef typename GET_PROP_TYPE(TypeTag, SolutionVector) SolutionVector;


    typedef typename GET_PROP_TYPE(TypeTag, FVElementGeometry) FVElementGeometry;
    typedef typename GridView::template Codim<0>::Entity Element;

public:
  
    typedef typename GET_PROP_TYPE(TypeTag, MaterialLaw) MaterialLaw;
    typedef typename MaterialLaw::Params MaterialLawParams;
    
    ThreePNISpatialParams(const GridView &gridView)
        : ParentType(gridView)
    {
        permeability_ = 1e-10;
        porosity_ = 0.4;

        // heat conductivity of granite
        lambdaSolid_ = 2.8;
	
	
        // residual saturations
        materialParams_.setSwr(0.12);
        materialParams_.setSwrx(0.12);
        materialParams_.setSnr(0.10);
        materialParams_.setSgr(0.01);

        // parameters for the 3phase van Genuchten law
        materialParams_.setVgAlpha(0.5);
        materialParams_.setVgn(4.0);
        materialParams_.setKrRegardsSnr(true);

        // parameters for adsorption
        materialParams_.setKdNAPL(0.);
        materialParams_.setRhoBulk(1500.);

    }

    ~ThreePNISpatialParams()
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


    bool useTwoPointGradient(const Element &element,
                             const int vertexI,
                             const int vertexJ) const
    {
        return false;
    }

    /*!
     * \brief Returns the heat capacity \f$[J/m^3 K]\f$ of the rock matrix.
     *
     * This is only required for non-isothermal models.
     *
     * \param element The finite element
     * \param fvGeometry The finite volume geometry
     * \param scvIdx The local index of the sub-control volume where
     *                    the heat capacity needs to be defined
     */
    Scalar heatCapacity(const Element &element,
                        const FVElementGeometry &fvGeometry,
                        const int scvIdx) const
    {
        return
            790 // specific heat capacity of granite [J / (kg K)]
            * 2700 // density of granite [kg/m^3]
            * (1 - porosity(element, fvGeometry, scvIdx));
    }

    /*!
     * \brief Returns the thermal conductivity \f$[W/m^2]\f$ of the porous material.
     *
     * \param element The finite element
     * \param fvGeometry The finite volume geometry
     * \param scvIdx The local index of the sub-control volume where
     *                    the heat capacity needs to be defined
     */
    Scalar thermalConductivitySolid(const Element &element,
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
