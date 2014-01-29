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
#ifndef DUMUX_TWOCNISTOKES2P2CNISPATIALPARAMS_HH
#define DUMUX_TWOCNISTOKES2P2CNISPATIALPARAMS_HH

#include <dune/grid/io/file/vtk/common.hh>

#include <dumux/material/spatialparams/implicitspatialparams.hh>
#include <dumux/material/fluidmatrixinteractions/2p/regularizedvangenuchten.hh>
#include <dumux/material/fluidmatrixinteractions/2p/regularizedbrookscorey.hh>
#include <dumux/material/fluidmatrixinteractions/2p/efftoabslaw.hh>
//#include <dumux/io/plotfluidmatrixlaw.hh>

namespace Dumux
{

//forward declaration
template<class TypeTag>
class TwoCNIStokesTwoPTwoCNISpatialParams;

namespace Properties
{
// The spatial parameters TypeTag
NEW_TYPE_TAG(TwoCNIStokesTwoPTwoCNISpatialParams);

// Set the spatial parameters
SET_TYPE_PROP(TwoCNIStokesTwoPTwoCNISpatialParams, SpatialParams,
        TwoCNIStokesTwoPTwoCNISpatialParams<TypeTag>);

// Set the material Law
SET_PROP(TwoCNIStokesTwoPTwoCNISpatialParams, MaterialLaw)
{
private:
    // define the material law which is parameterized by effective
    // saturations
    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
    typedef RegularizedVanGenuchten<Scalar> EffMaterialLaw;
public:
    // define the material law parameterized by absolute saturations
    typedef EffToAbsLaw<EffMaterialLaw> type;
};
}


/*!
 * \brief docme
 */
template<class TypeTag>
class TwoCNIStokesTwoPTwoCNISpatialParams : public ImplicitSpatialParams<TypeTag>
{
    typedef ImplicitSpatialParams<TypeTag> ParentType;
    typedef typename GET_PROP_TYPE(TypeTag, GridView) GridView;
    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
    typedef typename GridView::ctype CoordScalar;

    enum {
        dim=GridView::dimension,
        dimWorld=GridView::dimensionworld
    };

    typedef Dune::FieldVector<CoordScalar,dim> LocalPosition;
    typedef Dune::FieldVector<CoordScalar,dimWorld> GlobalPosition;
    typedef Dune::FieldVector<CoordScalar,dimWorld> DimVector;

    typedef typename GridView::IndexSet IndexSet;
    typedef typename GridView::template Codim<0>::Iterator ElementIterator;
    typedef typename GridView::template Codim<0>::Entity Element;

    typedef typename GET_PROP_TYPE(TypeTag, FVElementGeometry) FVElementGeometry;
    typedef typename GET_PROP_TYPE(TypeTag, FluxVariables) FluxVariables;
    typedef typename GET_PROP_TYPE(TypeTag, ElementVolumeVariables) ElementVolumeVariables;

    typedef std::vector<Scalar> PermeabilityType;
    typedef typename GET_PROP(TypeTag, ParameterTree) ParameterTree;
    typedef typename GridView::template Codim<dim>::Iterator VertexIterator;

public:
    typedef typename GET_PROP_TYPE(TypeTag, MaterialLaw) MaterialLaw;
    typedef typename MaterialLaw::Params MaterialLawParams;
    typedef std::vector<MaterialLawParams> MaterialLawParamsVector;

    /*!
     * \brief docme
     *
     * \param gridView The GridView which is used by the problem.
     *
     */
    TwoCNIStokesTwoPTwoCNISpatialParams(const GridView& gridView)
        : ParentType(gridView),
          permeability_(gridView.size(dim), 0.0),
          vanGenuchtenAlpha_(gridView.size(dim), 0.0),
          indexSet_(gridView.indexSet())
    {
        try
        {
            soilType_ = GET_RUNTIME_PARAM_FROM_GROUP(TypeTag, int, SpatialParams, SoilType);
            xMaterialInterface_ = GET_RUNTIME_PARAM_FROM_GROUP(TypeTag, Scalar, SpatialParams, MaterialInterfaceX);

            // porosities
            coarsePorosity_ = GET_RUNTIME_PARAM_FROM_GROUP(TypeTag, Scalar, SpatialParams.Coarse, Porosity1);
            mediumPorosity_ = GET_RUNTIME_PARAM_FROM_GROUP(TypeTag, Scalar, SpatialParams.Medium, Porosity2);
            finePorosity_ = GET_RUNTIME_PARAM_FROM_GROUP(TypeTag, Scalar, SpatialParams.Fine, Porosity3);

            // intrinsic permeabilities
            coarsePermeability_ = GET_RUNTIME_PARAM_FROM_GROUP(TypeTag, Scalar, SpatialParams.Coarse, Permeability1);
            mediumPermeability_ = GET_RUNTIME_PARAM_FROM_GROUP(TypeTag, Scalar, SpatialParams.Medium, Permeability2);
            finePermeability_ = GET_RUNTIME_PARAM_FROM_GROUP(TypeTag, Scalar, SpatialParams.Fine, Permeability3);

            // thermal conductivity of the solid material
            coarseLambdaSolid_ = GET_RUNTIME_PARAM_FROM_GROUP(TypeTag, Scalar, SpatialParams.Coarse, LambdaSolid1);
            mediumLambdaSolid_ = GET_RUNTIME_PARAM_FROM_GROUP(TypeTag, Scalar, SpatialParams.Medium, LambdaSolid2);
            fineLambdaSolid_ = GET_RUNTIME_PARAM_FROM_GROUP(TypeTag, Scalar, SpatialParams.Fine, LambdaSolid3);

            if (soilType_ != 0)
            {
                // residual saturations
                coarseParams_.setSwr(GET_RUNTIME_PARAM_FROM_GROUP(TypeTag, Scalar, SpatialParams.Coarse, Swr1));
                coarseParams_.setSnr(GET_RUNTIME_PARAM_FROM_GROUP(TypeTag, Scalar, SpatialParams.Coarse, Snr1));
                mediumParams_.setSwr(GET_RUNTIME_PARAM_FROM_GROUP(TypeTag, Scalar, SpatialParams.Medium, Swr2));
                mediumParams_.setSnr(GET_RUNTIME_PARAM_FROM_GROUP(TypeTag, Scalar, SpatialParams.Medium, Snr2));
                fineParams_.setSwr(GET_RUNTIME_PARAM_FROM_GROUP(TypeTag, Scalar, SpatialParams.Fine, Swr3));
                fineParams_.setSnr(GET_RUNTIME_PARAM_FROM_GROUP(TypeTag, Scalar, SpatialParams.Fine, Snr3));

                // parameters for the vanGenuchten law
                coarseParams_.setVgAlpha(GET_RUNTIME_PARAM_FROM_GROUP(TypeTag, Scalar, SpatialParams.Coarse, VgAlpha1));
                coarseParams_.setVgn(GET_RUNTIME_PARAM_FROM_GROUP(TypeTag, Scalar, SpatialParams.Coarse, VgN1));
                mediumParams_.setVgAlpha(GET_RUNTIME_PARAM_FROM_GROUP(TypeTag, Scalar, SpatialParams.Medium, VgAlpha2));
                mediumParams_.setVgn(GET_RUNTIME_PARAM_FROM_GROUP(TypeTag, Scalar, SpatialParams.Medium, VgN2));
                fineParams_.setVgAlpha(GET_RUNTIME_PARAM_FROM_GROUP(TypeTag, Scalar, SpatialParams.Fine, VgAlpha3));
                fineParams_.setVgn(GET_RUNTIME_PARAM_FROM_GROUP(TypeTag, Scalar, SpatialParams.Fine, VgN3));
            }
        }
        catch (Dumux::ParameterException &e) {
            std::cerr << e << ". Abort!\n";
            exit(1) ;
        }
        catch (...) {
            std::cerr << "Unknown exception thrown!\n";
            exit(1);
        }
    }

    ~TwoCNIStokesTwoPTwoCNISpatialParams()
    { }

    /*!
     * \brief Apply the intrinsic permeability tensor to a pressure
     *        potential gradient.
     *
     * \param element       The current finite element
     * \param fvGeometry    The current finite volume geometry of the element
     * \param scvIdx       The index sub-control volume face where the
     *                      intrinsic velocity ought to be calculated.
     */
    const Scalar intrinsicPermeability(const Element &element,
                                 	   const FVElementGeometry &fvGeometry,
                                 	   const int scvIdx) const
    {
		const GlobalPosition &globalPos = element.geometry().corner(scvIdx);

		if (checkSoilType(globalPos) == 1)
			return coarsePermeability_;
		if (checkSoilType(globalPos) == 2)
			return mediumPermeability_;
		if (checkSoilType(globalPos) == 3)
			return finePermeability_;
		if (checkSoilType(globalPos) == 4)
			return finePermeability_;
		else
			return mediumPermeability_;
    }

    //! \copydoc Dumux::ImplicitSpatialParamsOneP::porosity()
    Scalar porosity(const Element &element,
                    const FVElementGeometry &fvGeometry,
                    const int scvIdx) const
    {
        const GlobalPosition &globalPos = element.geometry().corner(scvIdx);

        if (checkSoilType(globalPos) == 1)
            return coarsePorosity_;
        if (checkSoilType(globalPos) == 2)
            return mediumPorosity_;
        if (checkSoilType(globalPos) == 3)
            return finePorosity_;
        else
            return mediumPorosity_;
    }


    //! \copydoc Dumux::ImplicitSpatialParams::materialLawParams()
    const MaterialLawParams& materialLawParams(const Element &element,
                                               const FVElementGeometry &fvGeometry,
                                               const int scvIdx) const
    {
		const GlobalPosition &globalPos = element.geometry().corner(scvIdx);
		if (checkSoilType(globalPos)==1)
			return coarseParams_;
		if (checkSoilType(globalPos)==2)
			return mediumParams_;
		if (checkSoilType(globalPos)==3)
			return fineParams_;
//            if (checkSoilType(globalPos)==4)
//                return leverettJParams_;
		else
			return mediumParams_;
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
            // 830 for sand
            790 // specific heat capacity of granite [J / (kg K)]
            * 2700 // density of granite [kg/m^3]
            * (1 - porosity(element, fvGeometry, scvIdx));
    }

    /*!
     * \brief docme
     *
     * \param element The finite element
     * \param fvGeometry The finite-volume geometry
     * \param scvIdx The local subcontrolvolume index
     *
     */
    Scalar thermalConductivitySolid(const Element &element,
                                    const FVElementGeometry &fvGeometry,
                                    const int scvIdx) const
    {
        const GlobalPosition &pos = element.geometry().corner(scvIdx);

        if (checkSoilType(pos) == 1)
            return coarseLambdaSolid_;
        if (checkSoilType(pos) == 2)
            return mediumLambdaSolid_;
        if (checkSoilType(pos) == 3)
            return fineLambdaSolid_;
        else
            return mediumLambdaSolid_;
    }

    /*!
     * \brief docme
     *
     * \param pos The position of the boundary face's integration point in global coordinates
     *
     */
    const unsigned checkSoilType(const GlobalPosition &pos) const
    {

// one soil, which can be chosen as runtime parameter:
// 1: coarse
// 2: medium
// 3: fine
// 4: LeverettJ (x < xMaterialInterface)
		return soilType_;
    }

    // this is called from the coupled problem
    // and creates a gnuplot output of the Pc-Sw curve
    /*!
     * \brief docme
     */
    void plotMaterialLaw()
    {
//        if (soilType_ == 0)
//        {
//            std::cout << "Material law plot not possible for heterogeneous media!\n";
//            return;
//        }
//        if (GET_RUNTIME_PARAM_FROM_GROUP(TypeTag, bool, SpatialParams.Coarse, PlotMaterialLaw1))
//        {
//            PlotFluidMatrixLaw<MaterialLaw> coarseFluidMatrixLaw_;
//            coarseFluidMatrixLaw_.plotpC(coarseParams_,
//                    GET_RUNTIME_PARAM_FROM_GROUP(TypeTag, Scalar, SpatialParams.Coarse, Swr1),
//                    1.0 - GET_RUNTIME_PARAM_FROM_GROUP(TypeTag, Scalar, SpatialParams.Coarse, Snr1),
//                    "pcSw_coarse");
//        }
//        if (GET_RUNTIME_PARAM_FROM_GROUP(TypeTag, bool, SpatialParams.Medium, PlotMaterialLaw2))
//        {
//            PlotFluidMatrixLaw<MaterialLaw> mediumFluidMatrixLaw_;
//            mediumFluidMatrixLaw_.plotpC(mediumParams_,
//                    GET_RUNTIME_PARAM_FROM_GROUP(TypeTag, Scalar, SpatialParams.Medium, Swr2),
//                    1.0 - GET_RUNTIME_PARAM_FROM_GROUP(TypeTag, Scalar, SpatialParams.Medium, Snr2),
//                    "pcSw_medium");
//            PlotFluidMatrixLaw<MaterialLaw> mediumFluidMatrixLaw2_;
//            mediumFluidMatrixLaw_.plotkrw(mediumParams_,
//                    GET_RUNTIME_PARAM_FROM_GROUP(TypeTag, Scalar, SpatialParams.Medium, Swr2),
//                    1.0 - GET_RUNTIME_PARAM_FROM_GROUP(TypeTag, Scalar, SpatialParams.Medium, Snr2));
//            PlotFluidMatrixLaw<MaterialLaw> mediumFluidMatrixLaw3_;
//            mediumFluidMatrixLaw_.plotkrn(mediumParams_,
//                    GET_RUNTIME_PARAM_FROM_GROUP(TypeTag, Scalar, SpatialParams.Medium, Swr2),
//                    1.0 - GET_RUNTIME_PARAM_FROM_GROUP(TypeTag, Scalar, SpatialParams.Medium, Snr2));
//        }
//        if (GET_RUNTIME_PARAM_FROM_GROUP(TypeTag, bool, SpatialParams.Fine, PlotMaterialLaw3))
//        {
//            PlotFluidMatrixLaw<MaterialLaw> fineFluidMatrixLaw_;
//            fineFluidMatrixLaw_.plotpC(fineParams_,
//                    GET_RUNTIME_PARAM_FROM_GROUP(TypeTag, Scalar, SpatialParams.Fine, Swr3),
//                    1.0 - GET_RUNTIME_PARAM_FROM_GROUP(TypeTag, Scalar, SpatialParams.Fine, Snr3),
//                    "pcSw_fine");
//        }
//        if (GET_RUNTIME_PARAM_FROM_GROUP(TypeTag, bool, SpatialParams.LeverettJ, PlotMaterialLawJ))
//        {
//            PlotFluidMatrixLaw<MaterialLaw> leverettJFluidMatrixLaw_;
//            leverettJFluidMatrixLaw_.plotpC(leverettJParams_,
//                    GET_RUNTIME_PARAM_FROM_GROUP(TypeTag, Scalar, SpatialParams.Medium, Swr2),
//                    1.0 - GET_RUNTIME_PARAM_FROM_GROUP(TypeTag, Scalar, SpatialParams.Medium, Snr2),
//                    "pcSw_leverett");
//        }
    }

private:

    unsigned soilType_;

    Scalar coarsePermeability_;
    Scalar mediumPermeability_;
    Scalar finePermeability_;

    Scalar coarsePorosity_;
    Scalar mediumPorosity_;
    Scalar finePorosity_;

    Scalar coarseLambdaSolid_;
    Scalar mediumLambdaSolid_;
    Scalar fineLambdaSolid_;

    Scalar xMaterialInterface_;
    MaterialLawParamsVector materialLawParams_;
    PermeabilityType permeability_;
    PermeabilityType vanGenuchtenAlpha_;
    const IndexSet& indexSet_;

    MaterialLawParams coarseParams_;
    MaterialLawParams mediumParams_;
    MaterialLawParams fineParams_;
};

} // end namespace
#endif

