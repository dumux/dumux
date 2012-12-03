// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*****************************************************************************
 *   Copyright (C) 2008-2009 by Klaus Mosthaf                                *
 *   Copyright (C) 2008-2009 by Andreas Lauser                               *
 *   Institute for Modelling Hydraulic and Environmental Systems             *
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
 * \brief Definition of the spatial parameters for the injection
 *        problem which uses the isothermal CO2 box model
 */

#ifndef DUMUX_HETEROGENEOUS_SPATIAL_PARAMS_HH
#define DUMUX_HETEROGENEOUS_SPATIAL_PARAMS_HH

#include <dumux/material/spatialparams/implicitspatialparams.hh>
#include <dumux/material/fluidmatrixinteractions/2p/linearmaterial.hh>
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
 * \ingroup BoxTestProblems
 * \brief Definition of the spatial parameters for the injection
 *        problem which uses the isothermal CO2 box model
 */
template<class TypeTag>
class HeterogeneousSpatialParams : public ImplicitSpatialParams<TypeTag>
{
    typedef ImplicitSpatialParams<TypeTag> ParentType;
    typedef typename GET_PROP_TYPE(TypeTag, Grid) Grid;
    typedef typename GET_PROP_TYPE(TypeTag, GridView) GridView;
    typedef Dune::GridPtr<Grid> GridPointer;
    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
    typedef typename Grid::ctype CoordScalar;
    typedef typename GET_PROP_TYPE(TypeTag, FluidSystem) FluidSystem;

    enum {
        dim=GridView::dimension,
        dimWorld=GridView::dimensionworld,

        lPhaseIdx = FluidSystem::lPhaseIdx
    };

    typedef Dune::FieldVector<CoordScalar,dimWorld> GlobalPosition;
    typedef Dune::FieldVector<CoordScalar,dimWorld> Vector;

    typedef typename GET_PROP_TYPE(TypeTag, FluxVariables) FluxVariables;
    typedef typename GET_PROP_TYPE(TypeTag, ElementVolumeVariables) ElementVolumeVariables;

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
        : ParentType(gridView)
    {
        /*
         * Layer Index Setup:
         * barrierTop = 4
         * barrierMiddle = 5
         * reservoir = 6
         */
        barrierTop_ = 1;
        barrierMiddle_ = 2;
        reservoir_ = 3;

        //Set the permeability tensor for the layers
        //Scalar anisotropy = 0.1;
        //for (int i = 0; i < dim; i++)
            barrierTopK_ = 1e-17; //sqm
        //barrierTopK_[dim-1][dim-1] = barrierTopK_[0][0]*anisotropy;

        //for (int i = 0; i < dim; i++)
            barrierMiddleK_ = 1e-15; //sqm
        //barrierMiddleK_[dim-1][dim-1] = barrierMiddleK_[0][0]*anisotropy;

     //   for (int i = 0; i < dim; i++)
            reservoirK_ = 1e-14; //sqm
        //reservoirK_[dim-1][dim-1] = reservoirK_[0][0]*anisotropy;

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

    void setParams(GridPointer *gridPtr)
    {
        gridPtr_ = gridPtr;
        int numElems = (*gridPtr_)->leafView().size(0);
        paramIdx_.resize(numElems);

        ElementIterator elemIt = (*gridPtr_)->leafView().template begin<0>();
        const ElementIterator elemEndIt = (*gridPtr_)->leafView().template end<0>();
        for (; elemIt != elemEndIt; ++elemIt)
        {
            int elemIdx = (*gridPtr_)->leafView().indexSet().index(*elemIt);
            int param = (*gridPtr_).parameters(*elemIt)[0];
            paramIdx_[elemIdx] = param;
            //std::cout<<"param: "<<paramIdx_[elemIdx]<<std::endl;
        }

    }

    /*!
     * \brief Apply the intrinsic permeability tensor to a pressure
     *        potential gradient.
     *
     * \param element The current finite element
     * \param fvElemGeom The current finite volume geometry of the element
     * \param scvIdx The index of the sub-control volume
     */
    const Scalar intrinsicPermeability(const Element &element,
                                       const FVElementGeometry &fvElemGeom,
                                       int scvIdx) const
    {
        int elemIdx = (*gridPtr_)->leafView().indexSet().index(element); //Get the global index of the element

        if (paramIdx_[elemIdx] == barrierTop_)
            return barrierTopK_;
        else if (paramIdx_[elemIdx] == barrierMiddle_)
            return barrierMiddleK_;
        else
            return reservoirK_;
    }

    /*!
     * \brief Define the porosity \f$[-]\f$ of the spatial parameters
     *
     * \param element The finite element
     * \param fvElemGeom The finite volume geometry
     * \param scvIdx The local index of the sub-control volume where
     *                    the porosity needs to be defined
     */
    Scalar porosity(const Element &element,
                    const FVElementGeometry &fvElemGeom,
                    int scvIdx) const
    {
        int elemIdx = (*gridPtr_)->leafView().indexSet().index(element); //Get the global index of the element

        if (paramIdx_[elemIdx] == barrierTop_)
            return barrierTopPorosity_;
        else if (paramIdx_[elemIdx] == barrierMiddle_)
            return barrierMiddlePorosity_;
        else
            return reservoirPorosity_;
    }


    /*!
     * \brief return the parameter object for the Brooks-Corey material law which depends on the position
     *
    * \param element The current finite element
    * \param fvElemGeom The current finite volume geometry of the element
    * \param scvIdx The index of the sub-control volume
    */
    const MaterialLawParams& materialLawParams(const Element &element,
                                                const FVElementGeometry &fvElemGeom,
                                                int scvIdx) const
    {

        return materialParams_;
    }

    /*!
     * \brief Returns the heat capacity \f$[J/m^3 K]\f$ of the rock matrix.
     *
     * This is only required for non-isothermal models.
     *
     * \param element The finite element
     * \param fvElemGeom The finite volume geometry
     * \param scvIdx The local index of the sub-control volume where
     *                    the heat capacity needs to be defined
     */
    double heatCapacity(const Element &element,
                        const FVElementGeometry &fvElemGeom,
                        int scvIdx) const
    {
        return
            790 // specific heat capacity of granite [J / (kg K)]
            * 2700 // density of granite [kg/m^3]
            * (1 - porosity(element, fvElemGeom, scvIdx));
    }

    /*!
     * \brief Calculate the heat flux \f$[W/m^2]\f$ through the
     *        rock matrix based on the temperature gradient \f$[K / m]\f$
     *
     * This is only required for non-isothermal models.
     *
     * \param heatFlux The resulting heat flux vector
     * \param fluxDat The flux variables
     * \param vDat The volume variables
     * \param tempGrad The temperature gradient
     * \param element The current finite element
     * \param fvElemGeom The finite volume geometry of the current element
     * \param scvfIdx The local index of the sub-control volume face where
     *                    the matrix heat flux should be calculated
     */
    void matrixHeatFlux(Vector &heatFlux,
                        const FluxVariables &fluxDat,
                        const ElementVolumeVariables &vDat,
                        const Vector &tempGrad,
                        const Element &element,
                        const FVElementGeometry &fvElemGeom,
                        int scvfIdx) const
    {
        static const Scalar lWater = 0.6;
        static const Scalar lGranite = 2.8;

        // arithmetic mean of the liquid saturation and the porosity
        const int i = fvElemGeom.subContVolFace[scvfIdx].i;
        const int j = fvElemGeom.subContVolFace[scvfIdx].j;
        Scalar Sl = std::max<Scalar>(0.0, (vDat[i].saturation(lPhaseIdx) +
                                           vDat[j].saturation(lPhaseIdx)) / 2);
        Scalar poro = (porosity(element, fvElemGeom, i) +
                       porosity(element, fvElemGeom, j)) / 2;

        Scalar lsat = pow(lGranite, (1-poro)) * pow(lWater, poro);
        Scalar ldry = pow(lGranite, (1-poro));

        // the heat conductivity of the matrix. in general this is a
        // tensorial value, but we assume isotropic heat conductivity.
        Scalar heatCond = ldry + sqrt(Sl) * (ldry - lsat);

        // the matrix heat flux is the negative temperature gradient
        // times the heat conductivity.
        heatFlux = tempGrad;
        heatFlux *= -heatCond;
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

    MaterialLawParams materialParams_;

    GridPointer *gridPtr_;
    std::vector<int> paramIdx_;
};

}

#endif
