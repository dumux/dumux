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
 * \brief Definition of the spatial parameters for the heterogeneous
 *        problem which uses the non-isothermal or isothermal CO2
 *        fully implicit model.
 */

#ifndef DUMUX_HETEROGENEOUS_SPATIAL_PARAMS_HH
#define DUMUX_HETEROGENEOUS_SPATIAL_PARAMS_HH

#include <dumux/material/spatialparams/implicitspatialparams.hh>
#include <dumux/material/fluidmatrixinteractions/2p/regularizedbrookscorey.hh>
#include <dumux/material/fluidmatrixinteractions/2p/efftoabslaw.hh>

#include <dumux/implicit/co2/co2model.hh>

namespace Dumux
{

//forward declaration
template<class TypeTag>
class HeterogeneousSpatialParams;

namespace Properties
{
// The spatial parameters TypeTag
NEW_TYPE_TAG(HeterogeneousSpatialParams);

// Set the spatial parameters
SET_TYPE_PROP(HeterogeneousSpatialParams, SpatialParams, Dumux::HeterogeneousSpatialParams<TypeTag>);

// Set the material Law
SET_PROP(HeterogeneousSpatialParams, MaterialLaw)
{
private:
    // define the material law which is parameterized by effective
    // saturations
    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
    typedef RegularizedBrooksCorey<Scalar> EffMaterialLaw;
public:
    // define the material law parameterized by absolute saturations
    typedef EffToAbsLaw<EffMaterialLaw> type;
};
}

/*!
 * \ingroup CO2Model
 * \ingroup ImplicitTestProblems
 * \brief Definition of the spatial parameters for the heterogeneous
 *        problem which uses the non-isothermal or isothermal CO2
 *        fully implicit model.
 */
template<class TypeTag>
class HeterogeneousSpatialParams : public ImplicitSpatialParams<TypeTag>
{
    typedef ImplicitSpatialParams<TypeTag> ParentType;
    typedef typename GET_PROP_TYPE(TypeTag, Grid) Grid;
    typedef typename GET_PROP_TYPE(TypeTag, GridCreator) GridCreator;
    typedef typename GET_PROP_TYPE(TypeTag, GridView) GridView;
    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
    typedef typename GET_PROP_TYPE(TypeTag, FluidSystem) FluidSystem;

    enum { dim=GridView::dimension };

    typedef typename GET_PROP_TYPE(TypeTag, FVElementGeometry) FVElementGeometry;
    typedef typename GridView::template Codim<0>::Entity Element;
    typedef typename GridView::template Codim<0>::Iterator ElementIterator;

public:
    typedef typename GET_PROP_TYPE(TypeTag, MaterialLaw) MaterialLaw;
    typedef typename MaterialLaw::Params MaterialLawParams;

    /*!
     * \brief The constructor
     *
     * \param gridView The grid view
     */
    HeterogeneousSpatialParams(const GridView &gridView)
        : ParentType(gridView), gridView_(gridView)
    {
        /*
         * Layer Index Setup:
         */
        barrierTop_ = 1;
        barrierMiddle_ = 2;
        reservoir_ = 3;

        // heat conductivity of granite
        lambdaSolid_ = 2.8;

        //Set the permeability for the layers
        barrierTopK_ = 1e-17; //sqm
        barrierMiddleK_ = 1e-15; //sqm
        reservoirK_ = 1e-14; //sqm

        //Set the effective porosity of the layers
        barrierTopPorosity_ = 0.001;
        barrierMiddlePorosity_ = 0.05;
        reservoirPorosity_ = 0.2;

        // Same material parameters for every layer
        materialParams_.setSwr(0.2);
        materialParams_.setSwr(0.05);
        materialParams_.setLambda(2.0);
        materialParams_.setPe(1e4);
    }

    ~HeterogeneousSpatialParams()
    {}

    /*!
     * \brief Reads layer information from the grid
     *
     */
    void setParams()
    {
        int numElements = gridView_.size(0);
        paramIdx_.resize(numElements);

        ElementIterator eIt = gridView_.template begin<0>();
        const ElementIterator eEndIt = gridView_.template end<0>();
        for (; eIt != eEndIt; ++eIt)
        {
            int eIdx = gridView_.indexSet().index(*eIt);
            int param = GridCreator::parameters(*eIt)[0];
            paramIdx_[eIdx] = param;
        }
    }

    /*!
     * \brief Returns the scalar intrinsic permeability \f$[m^2]\f$
     *
     * \param element The finite element
     * \param fvGeometry The finite volume geometry of the element
     * \param scvIdx The local index of the sub-control volume
     */
    const Scalar intrinsicPermeability(const Element &element,
                                       const FVElementGeometry &fvGeometry,
                                       int scvIdx) const
    {
        //Get the global index of the element
        int eIdx = gridView_.indexSet().index(element);

        if (paramIdx_[eIdx] == barrierTop_)
            return barrierTopK_;
        else if (paramIdx_[eIdx] == barrierMiddle_)
            return barrierMiddleK_;
        else
            return reservoirK_;
    }

    /*!
     * \brief Returns the porosity \f$[-]\f$
     *
     * \param element The finite element
     * \param fvGeometry The finite volume geometry of the element
     * \param scvIdx The local index of the sub-control volume
     */
    Scalar porosity(const Element &element,
                    const FVElementGeometry &fvGeometry,
                    int scvIdx) const
    {
        //Get the global index of the element
        int eIdx = gridView_.indexSet().index(element);

        if (paramIdx_[eIdx] == barrierTop_)
            return barrierTopPorosity_;
        else if (paramIdx_[eIdx] == barrierMiddle_)
            return barrierMiddlePorosity_;
        else
            return reservoirPorosity_;
    }


    /*!
     * \brief Returns the parameter object for the Brooks-Corey material law
     *
     * \param element The finite element
     * \param fvGeometry The finite volume geometry of the element
     * \param scvIdx The local index of the sub-control volume
     */
    const MaterialLawParams& materialLawParams(const Element &element,
                                                const FVElementGeometry &fvGeometry,
                                                int scvIdx) const
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
     * \brief Returns the thermal conductivity \f$\mathrm{[W/(m K)]}\f$ of the solid
     *
     * This is only required for non-isothermal models.
     *
     * \param element The finite element
     * \param fvGeometry The finite volume geometry of the element
     * \param scvIdx The local index of the sub-control volume
     */
    Scalar solidThermalConductivity(const Element &element,
                                    const FVElementGeometry &fvGeometry,
                                    const int scvIdx) const
    {
        return lambdaSolid_;
    }



private:

    int barrierTop_;
    int barrierMiddle_;
    int reservoir_;


    Scalar barrierTopPorosity_;
    Scalar barrierMiddlePorosity_;
    Scalar reservoirPorosity_;

    Scalar barrierTopK_;
    Scalar barrierMiddleK_;
    Scalar reservoirK_;
    Scalar lambdaSolid_;

    MaterialLawParams materialParams_;

    const GridView gridView_;
    std::vector<int> paramIdx_;
};

}

#endif
