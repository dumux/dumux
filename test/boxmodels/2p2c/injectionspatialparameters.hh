// $Id$
/*****************************************************************************
 *   Copyright (C) 2008-2009 by Klaus Mosthaf                                *
 *   Copyright (C) 2008-2009 by Andreas Lauser                               *
 *   Institute of Hydraulic Engineering                                      *
 *   University of Stuttgart, Germany                                        *
 *   email: <givenname>.<name>@iws.uni-stuttgart.de                          *
 *                                                                           *
 *   This program is free software; you can redistribute it and/or modify    *
 *   it under the terms of the GNU General Public License as published by    *
 *   the Free Software Foundation; either version 2 of the License, or       *
 *   (at your option) any later version, as long as this copyright notice    *
 *   is included in its original form.                                       *
 *                                                                           *
 *   This program is distributed WITHOUT ANY WARRANTY.                       *
 *****************************************************************************/
#ifndef DUMUX_INJECTION_SPATIAL_PARAMETERS_HH
#define DUMUX_INJECTION_SPATIAL_PARAMETERS_HH

#include <dumux/material/spatialparameters/boxspatialparameters.hh>
#include <dumux/material/fluidmatrixinteractions/2p/linearmaterial.hh>
#include <dumux/material/fluidmatrixinteractions/2p/regularizedbrookscorey.hh>
#include <dumux/material/fluidmatrixinteractions/2p/efftoabslaw.hh>

namespace Dumux
{

/**
 * \brief Definition of the soil properties for the injection problem
 *
 */
template<class TypeTag>
class InjectionSpatialParameters : public BoxSpatialParameters<TypeTag>
{
    typedef BoxSpatialParameters<TypeTag> ParentType;
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(Grid)) Grid;
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(GridView)) GridView;
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(Scalar)) Scalar;
    typedef typename Grid::ctype CoordScalar;
    enum {
        dim=GridView::dimension,
        dimWorld=GridView::dimensionworld,
    };

    typedef typename GET_PROP_TYPE(TypeTag, PTAG(FluidSystem)) FluidSystem;
    enum {
        lPhaseIdx = FluidSystem::wPhaseIdx,
        gPhaseIdx = FluidSystem::nPhaseIdx,
    };

    typedef Dune::FieldVector<CoordScalar,dim> LocalPosition;
    typedef Dune::FieldVector<CoordScalar,dimWorld> GlobalPosition;
    typedef Dune::FieldVector<CoordScalar,dimWorld> Vector;


    typedef typename GET_PROP_TYPE(TypeTag, PTAG(SolutionVector)) SolutionVector;

    typedef typename GET_PROP_TYPE(TypeTag, PTAG(SecondaryVars)) SecondaryVars;
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(FluxVars)) FluxVars;
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(ElementSecondaryVars)) ElementSecondaryVars;

    typedef typename GET_PROP_TYPE(TypeTag, PTAG(FVElementGeometry)) FVElementGeometry;
    typedef typename GridView::template Codim<0>::Entity Element;

    typedef RegularizedBrooksCorey<Scalar> EffMaterialLaw;
    //typedef LinearMaterial<Scalar> EffMaterialLaw;
public:
    typedef EffToAbsLaw<EffMaterialLaw> MaterialLaw;
    typedef typename MaterialLaw::Params MaterialLawParams;

    InjectionSpatialParameters(const GridView &gv)
        : ParentType(gv)
    {
        layerBottom_ = 22.0;

        // intrinsic permeabilities
        fineK_ = 1e-13;
        coarseK_ = 1e-12;

        // porosities
        finePorosity_ = 0.3;
        coarsePorosity_ = 0.3;

        // residual saturations
        fineMaterialParams_.setSwr(0.2);
        fineMaterialParams_.setSnr(0.0);
        coarseMaterialParams_.setSwr(0.2);
        coarseMaterialParams_.setSnr(0.0);

        // parameters for the Brooks-Corey law
        fineMaterialParams_.setPe(1e4);
        coarseMaterialParams_.setPe(1e4);
        fineMaterialParams_.setAlpha(2.0);
        coarseMaterialParams_.setAlpha(2.0);
    }

    ~InjectionSpatialParameters()
    {}


    /*!
     * \brief Update the spatial parameters with the flow solution
     *        after a timestep.
     *
     * \TODO call interface
     */
    void update(const SolutionVector &globalSolution)
    {
    };

    /*!
     * \brief Apply the intrinsic permeability tensor to a pressure
     *        potential gradient.
     *
     * \param element The current finite element
     * \param fvElemGeom The current finite volume geometry of the element
     * \param scvfIdx The index sub-control volume face where the
     *                      intrinsic velocity ought to be calculated.
     */
    const Scalar intrinsicPermeability(const Element           &element,
                                       const FVElementGeometry &fvElemGeom,
                                       int scvIdx) const
    {
        const GlobalPosition &pos = fvElemGeom.subContVol[scvIdx].global;
        if (isFineMaterial_(pos))
            return fineK_;
        return coarseK_;
    }

    /*!
     * \brief Define the porosity \f$[-]\f$ of the soil
     *
     * \param vDat The data defined on the sub-control volume
     * \param element The finite element
     * \param fvElemGeom The finite volume geometry
     * \param scvIdx The local index of the sub-control volume where
     *                    the porosity needs to be defined
     */
    double porosity(const Element           &element,
                    const FVElementGeometry &fvElemGeom,
                    int scvIdx) const
    {
        const GlobalPosition &pos = fvElemGeom.subContVol[scvIdx].global;
        if (isFineMaterial_(pos))
            return finePorosity_;
        return coarsePorosity_;
    }


    // return the brooks-corey context depending on the position
    const MaterialLawParams& materialLawParams(const Element           &element,
                                                const FVElementGeometry &fvElemGeom,
                                                int scvIdx) const
    {
        const GlobalPosition &pos = fvElemGeom.subContVol[scvIdx].global;
        if (isFineMaterial_(pos))
            return fineMaterialParams_;
        return coarseMaterialParams_;
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
    double heatCapacity(const Element           &element,
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
     * \param heatFlux The result vector
     * \param tempGrad The temperature gradient
     * \param element The current finite element
     * \param fvElemGeom The finite volume geometry of the current element
     * \param scvfIdx The local index of the sub-control volume face where
     *                    the matrix heat flux should be calculated
     */
    void matrixHeatFlux(Vector                  &heatFlux,
                        const FluxVars          &fluxDat,
                        const ElementSecondaryVars   &vDat,
                        const Vector            &tempGrad,
                        const Element           &element,
                        const FVElementGeometry &fvElemGeom,
                        int scvfIdx) const
    {
        static const Scalar lWater = 0.6;
        static const Scalar lGranite = 2.8;

        // arithmetic mean of the liquid saturation and the porosity
        const int i = fvElemGeom.subContVolFace[scvfIdx].i;
        const int j = fvElemGeom.subContVolFace[scvfIdx].j;
        Scalar Sl   = std::max(0.0, (vDat[i].saturation(lPhaseIdx) +
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
        heatFlux  = tempGrad;
        heatFlux *= -heatCond;
    }

private:
    bool isFineMaterial_(const GlobalPosition &pos) const
    { return pos[dim-1] > layerBottom_; };

    Scalar fineK_;
    Scalar coarseK_;
    Scalar layerBottom_;

    Scalar finePorosity_;
    Scalar coarsePorosity_;

    MaterialLawParams fineMaterialParams_;
    MaterialLawParams coarseMaterialParams_;
};

}

#endif
