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
 * \brief Definition of the spatial parameters for the SAGD problem.
 */

#ifndef DUMUX_SAGD_SPATIAL_PARAMS_HH
#define DUMUX_SAGD_SPATIAL_PARAMS_HH

#include <dumux/material/spatialparams/implicit.hh>

#include <dumux/porousmediumflow/3pwateroil/indices.hh>

#include <dumux/material/fluidmatrixinteractions/3p/regularizedparkervangen3p.hh>
#include <dumux/material/fluidmatrixinteractions/3p/regularizedparkervangen3pparams.hh>
#include <dumux/material/fluidmatrixinteractions/3p/efftoabslaw.hh>

namespace Dumux
{

//forward declaration
template<class TypeTag>
class SagdSpatialParams;

namespace Properties
{
// The spatial parameters TypeTag
NEW_TYPE_TAG(SagdSpatialParams);

// Set the spatial parameters
SET_TYPE_PROP(SagdSpatialParams, SpatialParams, Dumux::SagdSpatialParams<TypeTag>);

// Set the material Law
SET_PROP(SagdSpatialParams, MaterialLaw)
{
 private:
    // define the material law which is parameterized by effective
    // saturations
    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
    typedef RegularizedParkerVanGen3P<Scalar> EffectiveLaw;
 public:
    // define the material law parameterized by absolute saturations
    typedef EffToAbsLaw<EffectiveLaw> type;
};
}

/*!
 * \ingroup ThreePThreeCNIModel
 *
 * \brief Definition of the spatial parameters for the SAGD problem
 */
template<class TypeTag>
class SagdSpatialParams : public ImplicitSpatialParams<TypeTag>
{
    typedef ImplicitSpatialParams<TypeTag> ParentType;

    typedef typename GET_PROP_TYPE(TypeTag, Grid) Grid;
    typedef typename GET_PROP_TYPE(TypeTag, GridView) GridView;
    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
    typedef typename Grid::ctype CoordScalar;
    enum {
        dim=GridView::dimension,
        dimWorld=GridView::dimensionworld
    };

    typedef typename GET_PROP_TYPE(TypeTag, Indices) Indices;
    enum {
        wPhaseIdx = Indices::wPhaseIdx,
        nPhaseIdx = Indices::nPhaseIdx
    };

    typedef Dune::FieldVector<CoordScalar,dimWorld> GlobalPosition;
    typedef Dune::FieldVector<CoordScalar,dimWorld> DimVector;


    typedef typename GET_PROP_TYPE(TypeTag, SolutionVector) SolutionVector;

    typedef typename GET_PROP_TYPE(TypeTag, FluxVariables) FluxVariables;
    typedef typename GET_PROP_TYPE(TypeTag, ElementVolumeVariables) ElementVolumeVariables;

    typedef typename GET_PROP_TYPE(TypeTag, FVElementGeometry) FVElementGeometry;
    typedef typename GridView::template Codim<0>::Entity Element;


public:
    typedef typename GET_PROP_TYPE(TypeTag, MaterialLaw) MaterialLaw;
    typedef typename MaterialLaw::Params MaterialLawParams;

    /*!
     * \brief The constructor
     *
     * \param gridView The grid view
     */
    SagdSpatialParams(const GridView &gridView)
        : ParentType(gridView)
    {
        layerBottom_ = 35.0;

        // intrinsic permeabilities
        fineK_ = 1e-16;
        coarseK_ = 4e-14;

        // porosities
        finePorosity_ = 0.10;
        coarsePorosity_ = 0.1;

        // heat conductivity of granite
        //lambdaSolid_ = 2.8;

        // specific heat capacities
        fineHeatCap_ = 850.;
        coarseHeatCap_ = 850;

        // residual saturations
        fineMaterialParams_.setSwr(0.1);
        fineMaterialParams_.setSwrx(0.12);  //Total liquid Residual Saturation
        fineMaterialParams_.setSnr(0.09);   //Residual of NAPL if there is no water
        fineMaterialParams_.setSgr(0.01);
        coarseMaterialParams_.setSwr(0.1);
        coarseMaterialParams_.setSwrx(0.12);
        coarseMaterialParams_.setSnr(0.09);
        coarseMaterialParams_.setSgr(0.01);

        // parameters for the 3phase van Genuchten law
        fineMaterialParams_.setVgn(4.0);
        coarseMaterialParams_.setVgn(4.0);
        fineMaterialParams_.setVgAlpha(1.);
        coarseMaterialParams_.setVgAlpha(1.);

        coarseMaterialParams_.setKrRegardsSnr(false);
        fineMaterialParams_.setKrRegardsSnr(false);
    }

    ~SagdSpatialParams()
    {}

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
        const GlobalPosition &pos = fvElemGeom.subContVol[scvIdx].global;
        if (isFineMaterial_(pos))
            return fineK_;
        return coarseK_;
    }

    /*!
     * \brief Define the porosity \f$[-]\f$ of the spatial parameters
     *
     * \param element The finite element
     * \param fvGeometry The finite volume geometry
     * \param scvIdx The local index of the sub-control volume where
     *                    the porosity needs to be defined
     */
    Scalar porosity(const Element &element,
                    const FVElementGeometry &fvGeometry,
                    const int scvIdx) const
    {
        const GlobalPosition &pos = fvGeometry.subContVol[scvIdx].global;
        if (isFineMaterial_(pos))
            return finePorosity_;
        else
            return coarsePorosity_;
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
        const GlobalPosition &pos = fvGeometry.subContVol[scvIdx].global;
        if (isFineMaterial_(pos))
            return fineMaterialParams_;
        else
            return coarseMaterialParams_;
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
        const GlobalPosition &pos = fvGeometry.subContVol[scvIdx].global;
        if (isFineMaterial_(pos))
            return fineHeatCap_ * 2650 // density of sand [kg/m^3]
                * (1 - porosity(element, fvGeometry, scvIdx));
        else
            return coarseHeatCap_ * 2650 // density of sand [kg/m^3]
                * (1 - porosity(element, fvGeometry, scvIdx));
    }

    /*!
     * \brief Calculate the heat flux \f$[W/m^2]\f$ through the
     *        rock matrix based on the temperature gradient \f$[K / m]\f$
     *
     * This is only required for non-isothermal models.
     *
     * \param heatFlux The resulting heat flux vector
     * \param fluxDat The flux variables
     * \param elemVolVars The volume variables
     * \param tempGrad The temperature gradient
     * \param element The current finite element
     * \param fvGeometry The finite volume geometry of the current element
     * \param faceIdx The local index of the sub-control volume face where
     *                    the matrix heat flux should be calculated
     */
    void matrixHeatFlux(DimVector &heatFlux,
                        const FluxVariables &fluxDat,
                        const ElementVolumeVariables &elemVolVars,
                        const DimVector &tempGrad,
                        const Element &element,
                        const FVElementGeometry &fvGeometry,
                        int faceIdx) const
    {
        static const Scalar ldry = 0.35;
        static const Scalar lSw1 = 1.8;
        static const Scalar lSn1 = 0.65;

        // arithmetic mean of the liquid saturation and the porosity
        const int i = fvGeometry.subContVolFace[faceIdx].i;
        const int j = fvGeometry.subContVolFace[faceIdx].j;
        const Scalar Sw = std::max(0.0, (elemVolVars[i].saturation(wPhaseIdx) +
                                         elemVolVars[j].saturation(wPhaseIdx)) / 2);
        const Scalar Sn = std::max(0.0, (elemVolVars[i].saturation(nPhaseIdx) +
                                         elemVolVars[j].saturation(nPhaseIdx)) / 2);

        // the heat conductivity of the matrix. in general this is a
        // tensorial value, but we assume isotropic heat conductivity.
        const Scalar heatCond = ldry + std::sqrt(Sw) * (lSw1-ldry) + std::sqrt(Sn) * (lSn1-ldry);

        // the matrix heat flux is the negative temperature gradient
        // times the heat conductivity.
        heatFlux = tempGrad;
        heatFlux *= -heatCond;
    }

private:
    bool isFineMaterial_(const GlobalPosition &pos) const
    {
        return pos[dim-1] > layerBottom_;
    };

    Scalar layerBottom_;
    Scalar lambdaSolid_;

    Scalar fineK_;
    Scalar coarseK_;

    Scalar finePorosity_;
    Scalar coarsePorosity_;

    Scalar fineHeatCap_;
    Scalar coarseHeatCap_;

    MaterialLawParams fineMaterialParams_;
    MaterialLawParams coarseMaterialParams_;
};

}

#endif
