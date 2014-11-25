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
 * \brief Spatial parameters for the
 *        coupling of an isothermal two-component Stokes
 *        and an isothermal two-phase two-component Darcy model.
 */

#ifndef DUMUX_TWOCSTOKES_2P2C_SPATIALPARAMS_HH
#define DUMUX_TWOCSTOKES_2P2C_SPATIALPARAMS_HH

#include <dune/grid/io/file/vtk/common.hh>

//#include <dumux/material/spatialparameters/gstatrandompermeability.hh>
#include <dumux/material/spatialparams/implicitspatialparams.hh>
#include <dumux/material/fluidmatrixinteractions/2p/linearmaterial.hh>
#include <dumux/material/fluidmatrixinteractions/2p/regularizedvangenuchten.hh>
//#include <dumux/material/fluidmatrixinteractions/2p/linearizedregvangenuchtenparams.hh>
//#include <dumux/material/fluidmatrixinteractions/2p/regularizedbrookscorey.hh>
#include <dumux/material/fluidmatrixinteractions/2p/efftoabslaw.hh>
//#include <dumux/io/plotfluidmatrixlaw.hh>

namespace Dumux
{

//forward declaration
template<class TypeTag>
class TwoCStokesTwoPTwoCSpatialParams;

namespace Properties
{
// The spatial parameters TypeTag
NEW_TYPE_TAG(TwoCStokesTwoPTwoCSpatialParams);

// Set the spatial parameters
SET_TYPE_PROP(TwoCStokesTwoPTwoCSpatialParams, SpatialParams,
        TwoCStokesTwoPTwoCSpatialParams<TypeTag>);

// Set the material Law
SET_PROP(TwoCStokesTwoPTwoCSpatialParams, MaterialLaw)
{
private:
    // define the material law which is parameterized by effective
    // saturations
    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
    // typedef RegularizedBrooksCorey<Scalar> EffMaterialLaw;
    // typedef RegularizedVanGenuchten<Scalar, LinearizedRegVanGenuchtenParams<Scalar, TypeTag> > EffMaterialLaw;
    typedef RegularizedVanGenuchten<Scalar> EffMaterialLaw;
public:
    // define the material law parameterized by absolute saturations
    typedef EffToAbsLaw<EffMaterialLaw> type;
};
}


/*!
 * \ingroup ImplicitTestProblems
 * \ingroup MultidomainProblems
 * \brief Definition of the spatial parameters for
 *        the coupling of an isothermal two-component Stokes
 *        and an isothermal two-phase two-component Darcy model.
 */
template<class TypeTag>
class TwoCStokesTwoPTwoCSpatialParams : public ImplicitSpatialParams<TypeTag>
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
     * \brief Spatial parameters for the
     *        coupling of an isothermal two-component Stokes
     *        and an isothermal two-phase two-component Darcy model.
     *
     * \param gridView The GridView which is used by the problem
     */
    TwoCStokesTwoPTwoCSpatialParams(const GridView& gridView)
        : ParentType(gridView),
//          permeability_(gridView.size(dim), 0.0),
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

    /*!
     * \brief The destructor
     */
    ~TwoCStokesTwoPTwoCSpatialParams()
    {}

    /*!
     * \brief Returns the intrinsic permeability tensor \f$[m^2]\f$
     *
     * \param element The finite element
     * \param fvGeometry The finite volume geometry of the element
     * \param scvIdx The local index of the sub-control volume
     */
    const Scalar intrinsicPermeability(const Element &element,
                                       const FVElementGeometry &fvGeometry,
                                       const int scvIdx) const
    {
//        // heterogeneous parameter field computed with GSTAT
//        if (soilType_ == 0)
//            return permeability_[indexSet_.index(
//                *(element.template subEntity<dim> (scvIdx)))];
//        // soil can be chosen by setting the parameter soiltype accordingly
//        else
//        {
            const GlobalPosition &globalPos = element.geometry().corner(scvIdx);

            if (checkSoilType(globalPos) == 1)
                return coarsePermeability_;
            if (checkSoilType(globalPos) == 3)
                return finePermeability_;
            else
                return mediumPermeability_;
//        }
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
                    const int scvIdx) const
    {
        const GlobalPosition &globalPos = element.geometry().corner(scvIdx);

        if (checkSoilType(globalPos) == 1)
            return coarsePorosity_;
        if (checkSoilType(globalPos) == 3)
            return finePorosity_;
        else
            return mediumPorosity_;
    }


    /*!
     * \brief Returns the parameter object for the material law
     *
     * \param element The finite element
     * \param fvGeometry The finite volume geometry of the element
     * \param scvIdx The local index of the sub-control volume
     */
    const MaterialLawParams& materialLawParams(const Element &element,
                                               const FVElementGeometry &fvGeometry,
                                               const int scvIdx) const
    {
//        if (soilType_ == 0)
//            return materialLawParams_[indexSet_.index(
//                *(element.template subEntity<dim> (scvIdx)))];
//        else
//        {
            const GlobalPosition &globalPos = element.geometry().corner(scvIdx);
            if (checkSoilType(globalPos)==1)
                return coarseParams_;
            if (checkSoilType(globalPos)==3)
                return fineParams_;
            else
                return mediumParams_;
//        }
    }


    /*!
     * \brief Returns the effective heat capacity \f$[J/m^3 K]\f$
     *
     * This is only required for non-isothermal models. This function does not
     * return the specific heat capacity, but an effective heat capacity, which is
     * \f$c_\textrm{p,eff,s} = c_\textrm{p,s} \varrho_\textrm{s} \left(1 - \phi\right)\f$
     *
     * \param element The finite element
     * \param fvGeometry The finite volume geometry of the element
     * \param scvIdx The local index of the sub-control volume
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
     * \brief Returns the thermal conductivity \f$[W/m^2]\f$ of the solid
     *
     * This is only required for non-isothermal models.
     *
     * \param element The finite element
     * \param fvGeometry The finite volume geometry of the element
     * \param scvIdx The local index of the sub-control volume
     */
    Scalar thermalConductivitySolid(const Element &element,
                                    const FVElementGeometry &fvGeometry,
                                    const int scvIdx) const
    {
        const GlobalPosition &pos = element.geometry().corner(scvIdx);

        if (checkSoilType(pos) == 1)
            return coarseLambdaSolid_;
        if (checkSoilType(pos) == 3)
            return fineLambdaSolid_;
        else
            return mediumLambdaSolid_;
    }

    /*!
     * \brief Returns the index of the used soil type
     *
     * The soil, can be chosen as runtime parameter:
     * 1: coarse,
     * 2: medium,
     * 3: fine
     *
     * \param globalPos The global position
     */
    const unsigned checkSoilType(const GlobalPosition &globalPos) const
    {
        return soilType_;
    }

    /*!
     * \brief This method allows the generation of a statistical field
     *        for the intrinsic permeability using GStat
     *
     * Because gstat is not open source and has to be installed manually,
     * the content of this function is deactivated (commented) by
     * default. If you have gstat installed, please uncomment
     * the lines between the curly brackets.
     *
     * \param gridView The GridView which is used by the problem
     */
    void loadIntrinsicPermeability(const GridView& gridView)
    {
//        // only load random field, if soilType is set to 0
//        if (soilType_ != 0)
//            return;
//
//        const unsigned size = gridView.size(dim);
//        permeability_.resize(size, 0.0);
//        materialLawParams_.resize(size);
//
//        bool create = true;
//        if (ParameterTree::tree().hasKey("SpatialParams.GenerateNewPermeability"))
//            create = GET_RUNTIME_PARAM_FROM_GROUP(TypeTag, bool , SpatialParams, GenerateNewPermeability);
//
//        std::string gStatControlFileName("gstatControl_2D.txt");
//        if (ParameterTree::tree().hasKey("SpatialParams.GStatControlFileName"))
//        gStatControlFileName = GET_RUNTIME_PARAM_FROM_GROUP(TypeTag, std::string, SpatialParams,
//                        GStatControlFileName);
//
//
//        std::string gStatInputFileName("gstatInput.txt");
//        if (ParameterTree::tree().hasKey("SpatialParams.GStatInputFileName"))
//                gStatInputFileName = GET_RUNTIME_PARAM_FROM_GROUP(TypeTag, std::string, SpatialParams,
//                        GStatInputFileName);
//
//        std::string permeabilityFileName("permeab.dat");
//        if (ParameterTree::tree().hasKey("SpatialParams.PermeabilityInputFileName"))
//            permeabilityFileName = GET_RUNTIME_PARAM_FROM_GROUP(TypeTag, std::string, SpatialParams,
//                        PermeabilityInputFileName);
//
//        // create random permeability object
//        GstatRandomPermeability<GridView, Scalar> randomPermeability(gridView,
//                                                                     create,
//                                                                     permeabilityFileName.c_str(),
//                                                                     gStatControlFileName.c_str(),
//                                                                     gStatInputFileName.c_str());
//
//        Scalar totalVolume = 0;
//        Scalar meanPermeability = 0;
//
////        for (unsigned n=0; n<numVertices; ++n)
////        {
////            //TODO the following can be done in a loop
////        }
//        // iterate over elements
//        ElementIterator eItEnd = gridView.template end<0> ();
//        for (ElementIterator eIt = gridView.template begin<0> (); eIt
//                != eItEnd; ++eIt)
//        {
//            Scalar perm = randomPermeability.K(*eIt)[0][0];
//            permeability_[indexSet_.index(*(eIt->template subEntity<dim> (3)))]
//                          = perm;
//
//            Scalar volume = eIt->geometry().volume();
//            totalVolume += volume;
//
//            meanPermeability += volume / perm;
//        }
//
//        for (ElementIterator eIt = gridView.template begin<0> (); eIt
//                != eItEnd; ++eIt)
//        {
//            Scalar perm = randomPermeability.K(*eIt)[0][0];
//            permeability_[indexSet_.index(*(eIt->template subEntity<dim> (1)))]
//                          = perm;
//
//            Scalar volume = eIt->geometry().volume();
//            totalVolume += volume;
//
//            meanPermeability += volume / perm;
//        }
//
//        for (ElementIterator eIt = gridView.template begin<0> (); eIt
//                != eItEnd; ++eIt)
//        {
//            Scalar perm = randomPermeability.K(*eIt)[0][0];
//            permeability_[indexSet_.index(*(eIt->template subEntity<dim> (2)))]
//                          = perm;
//
//            Scalar volume = eIt->geometry().volume();
//            totalVolume += volume;
//
//            meanPermeability += volume / perm;
//        }
//
//        for (ElementIterator eIt = gridView.template begin<0> (); eIt
//                != eItEnd; ++eIt)
//        {
//            Scalar perm = randomPermeability.K(*eIt)[0][0];
//            permeability_[indexSet_.index(*(eIt->template subEntity<dim> (0)))]
//                          = perm;
//
//            Scalar volume = eIt->geometry().volume();
//            totalVolume += volume;
//
//            meanPermeability += volume / perm;
//        }
//
//        meanPermeability /= totalVolume;
//        meanPermeability = 1.0 / meanPermeability;
//
//        //Iterate over elements
//        VertexIterator vItEnd = gridView.template end<dim> ();
//        for (VertexIterator vIt = gridView.template begin<dim> (); vIt
//                != vItEnd; ++vIt)
//        {
//            int vIdxGlobal = indexSet_.index(*vIt);
//
//            Scalar elementEntryPressure = 1./GET_RUNTIME_PARAM_FROM_GROUP(TypeTag, Scalar, SpatialParams.Medium, VgAlpha2);
//            elementEntryPressure *= sqrt(meanPermeability / permeability_[vIdxGlobal]);
//
//            materialLawParams_[vIdxGlobal].setVgAlpha(1./elementEntryPressure);
//            materialLawParams_[vIdxGlobal].setVgn(GET_RUNTIME_PARAM_FROM_GROUP(TypeTag, Scalar, SpatialParams.Medium, VgN2));
//            materialLawParams_[vIdxGlobal].setSwr(GET_RUNTIME_PARAM_FROM_GROUP(TypeTag, Scalar, SpatialParams.Medium, Swr2));
//            materialLawParams_[vIdxGlobal].setSnr(GET_RUNTIME_PARAM_FROM_GROUP(TypeTag, Scalar, SpatialParams.Medium, Snr2));
//        }
//
//        //        if (create)
//        //        {
//        Dune::VTKWriter<GridView> vtkwriter(gridView);
//        vtkwriter.addVertexData(permeability_, "absolute permeability");
//        PermeabilityType logPerm(size);
//        for (unsigned i = 0; i < size; i++)
//            logPerm[i] = log10(permeability_[i]);
//        vtkwriter.addVertexData(logPerm, "logarithm of permeability");
//        vtkwriter.write("permeability", Dune::VTK::OutputType::ascii);
    }

    /*!
     * \brief This is called from the coupled problem and creates
     *        a gnuplot output of the Pc-Sw curve
     *
     * If this function should be used, uncomment the lines between
     * the curly brackets.
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
//                    1.0 - GET_RUNTIME_PARAM_FROM_GROUP(TypeTag, Scalar, SpatialParams.Coarse, Snr1));
//        }
//        if (GET_RUNTIME_PARAM_FROM_GROUP(TypeTag, bool, SpatialParams.Medium, PlotMaterialLaw2))
//        {
//            PlotFluidMatrixLaw<MaterialLaw> mediumFluidMatrixLaw_;
//            mediumFluidMatrixLaw_.plotpC(mediumParams_,
//                    GET_RUNTIME_PARAM_FROM_GROUP(TypeTag, Scalar, SpatialParams.Medium, Swr2),
//                    1.0 - GET_RUNTIME_PARAM_FROM_GROUP(TypeTag, Scalar, SpatialParams.Medium, Snr2));
//        }
//        if (GET_RUNTIME_PARAM_FROM_GROUP(TypeTag, bool, SpatialParams.Fine, PlotMaterialLaw3))
//        {
//            PlotFluidMatrixLaw<MaterialLaw> fineFluidMatrixLaw_;
//            fineFluidMatrixLaw_.plotpC(fineParams_,
//                    GET_RUNTIME_PARAM_FROM_GROUP(TypeTag, Scalar, SpatialParams.Fine, Swr3),
//                    1.0 - GET_RUNTIME_PARAM_FROM_GROUP(TypeTag, Scalar, SpatialParams.Fine, Snr3));
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
//    MaterialLawParamsVector materialLawParams_;
//    PermeabilityType permeability_;
    const IndexSet& indexSet_;

    MaterialLawParams coarseParams_;
    MaterialLawParams mediumParams_;
    MaterialLawParams fineParams_;
};

} // end namespace Dumux

#endif // DUMUX_TWOCSTOKES_2P2C_SPATIALPARAMS_HH
