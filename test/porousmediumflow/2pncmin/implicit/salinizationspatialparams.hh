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
 * \ingroup TwoPNCMinTests
 * \brief Spatial parameters for the dissolution problem
 * where water is injected in a for flushing precipitated salt clogging a gas reservoir.
 */
#ifndef DUMUX_INJECTION_SPATIAL_PARAMETERS_HH
#define DUMUX_INJECTION_SPATIAL_PARAMETERS_HH

#include <dumux/io/gnuplotinterface.hh>
#include <dumux/io/plotmateriallaw.hh>
#include <dumux/material/spatialparams/fv.hh>
#include <dumux/material/fluidmatrixinteractions/2p/linearmaterial.hh>
#include <dumux/material/fluidmatrixinteractions/2p/regularizedvangenuchten.hh>
#include <dumux/material/fluidmatrixinteractions/2p/efftoabslaw.hh>
#include <dumux/material/fluidmatrixinteractions/porosityprecipitation.hh>
#include <dumux/material/fluidmatrixinteractions/permeabilitykozenycarman.hh>
#include <random>

namespace Dumux {

/*!
 * \ingroup TwoPNCMinTests
 * \brief Spatial parameters for the dissolution problem
 * where water is injected in a for flushing precipitated salt clogging a gas reservoir.
 */
template<class FVGridGeometry, class Scalar>
class DissolutionSpatialParams
: public FVSpatialParams<FVGridGeometry, Scalar,
                         DissolutionSpatialParams<FVGridGeometry, Scalar>>
{
    using GridView = typename FVGridGeometry::GridView;
    using FVElementGeometry = typename FVGridGeometry::LocalView;
    using SubControlVolume = typename FVElementGeometry::SubControlVolume;
    using Element = typename GridView::template Codim<0>::Entity;

    using ParentType = FVSpatialParams<FVGridGeometry, Scalar,
                                       DissolutionSpatialParams<FVGridGeometry, Scalar>>;

    using EffectiveLaw = RegularizedVanGenuchten<Scalar>;

    using GlobalPosition = typename SubControlVolume::GlobalPosition;

public:
    // type used for the permeability (i.e. tensor or scalar)
    using PermeabilityType = Scalar;
    //! export the material law type used
    using MaterialLaw = EffToAbsLaw<EffectiveLaw>;
    using MaterialLawParams = typename MaterialLaw::Params;

    DissolutionSpatialParams(std::shared_ptr<const FVGridGeometry> fvGridGeometry)
    : ParentType(fvGridGeometry)
    {
        solubilityLimit_       = getParam<Scalar>("SpatialParams.SolubilityLimit", 0.26);
        referencePorosity_     = getParam<Scalar>("SpatialParams.referencePorosity", 0.11);
        referencePermeability_ = getParam<Scalar>("SpatialParams.referencePermeability", 2.23e-14);
        sdevPorosity_ = getParam<Scalar>("SpatialParams.PorosityStandardDeviation", 1e-14);

        irreducibleLiqSat_     = getParam<Scalar>("SpatialParams.IrreducibleLiqSat", 0.2);
        irreducibleGasSat_     = getParam<Scalar>("SpatialParams.IrreducibleGasSat", 1e-3);
        vgAlpha_               = getParam<Scalar>("SpatialParams.VGAlpha", 1.5);
        vgn_                   = getParam<Scalar>("SpatialParams.VGn", 4.0);

        plotFluidMatrixInteractions_ = getParam<bool>("Output.PlotFluidMatrixInteractions");





        unsigned seed = 1; //std::chrono::system_clock::now().time_since_epoch().count();
        std::default_random_engine generator (seed);
        std::normal_distribution<double> distribution (referencePorosity_,sdevPorosity_);

        const auto numElems = this->fvGridGeometry().elementMapper().size();

        porosity_.resize(numElems, 0.0);
        permeability_.resize(numElems, 0.0);
        materialParams_.resize(numElems);

//         for (const auto& element : elements(fvGridGeometry->gridView()))

        std::cout << "numElems: " << numElems << std::endl;

        for (int eIdx = numElems-1; eIdx >= 0; --eIdx)
        {
            std::cout << "at eIdx : " << eIdx << std::endl;
//             const int eIdx = this->fvGridGeometry().elementMapper().index(element);
            //Van Genuchen parameters
            // residual saturations
            materialParams_[eIdx].setSwr(irreducibleLiqSat_);
            materialParams_[eIdx].setSnr(irreducibleGasSat_);

            //Van Genuchen parameters
            materialParams_[eIdx].setVgAlpha(vgAlpha_ );
            materialParams_[eIdx].setVgn(vgn_);


//             Scalar numberOfElements = this->fvGridGeometry().elementMapper().size();
            porosity_[eIdx] = distribution(generator);

            //avoid negative permeability values
            if(porosity_[eIdx] < 0)
            { porosity_[eIdx] = 1e-5; }

            // just valid if no precipitant in the beginning
            permeability_[eIdx] = permLaw_.evaluatePermeability(referencePermeability_, referencePorosity_, porosity_[eIdx]);

            std::cout << porosity_[eIdx] << " " << permeability_[eIdx] << std::endl;
            materialParams_[eIdx].setLeverettFactor(pow(referencePermeability_/permeability_[eIdx]*porosity_[eIdx]/referencePorosity_, 0.5));
        }
    }

//     /*!
//      * \brief This is called from the problem and creates a gnuplot output
//      *        of e.g the pc-Sw curve
//      */
//     void plotMaterialLaw()
//     {
//         PlotMaterialLaw<Scalar, MaterialLaw> plotMaterialLaw;
//         GnuplotInterface<Scalar> gnuplot(plotFluidMatrixInteractions_);
//         gnuplot.setOpenPlotWindow(plotFluidMatrixInteractions_);
//         plotMaterialLaw.addpcswcurve(gnuplot, materialParams_, irreducibleLiqSat_, 1-irreducibleGasSat_);
//         gnuplot.setOption("set xrange [0:1]");
//         gnuplot.setOption("set yrange [0:15000]");
//         gnuplot.setOption("set label \"residual\\nsaturation\" at 0.1,100000 center");
//         gnuplot.plot("pc-Sw");
//
//         gnuplot.resetAll();
//         gnuplot.setOption("set xrange [0:1]");
//         gnuplot.setOption("set yrange [0:1]");
//         plotMaterialLaw.addkrcurves(gnuplot, materialParams_, irreducibleLiqSat_, 1/*-irreducibleGasSat_*/, "fine");
//         gnuplot.plot("kr");
//     }

    /*!
     *  \brief Define the minimum porosity \f$[-]\f$ distribution
     *
     *  \param element The finite element
     *  \param scv The sub-control volume
     */
    Scalar minimalPorosity(const Element& element, const SubControlVolume &scv) const
    { return 1e-5; }


        /*!
     * \brief Returns the volume fraction of the inert component with index compIdx \f$[-]\f$
     *
     * \param element The current element
     * \param scv The sub-control volume inside the element.
     * \param elemSol The element solution
     * \param compIdx The solid component index
     * \return solid volume fraction
     */
    template<class SolidSystem, class ElementSolution>
    Scalar inertVolumeFraction(const Element& element,
                               const SubControlVolume& scv,
                               const ElementSolution& elemSol,
                               int compIdx) const
    {
        unsigned int eIdx = this->fvGridGeometry().elementMapper().index(element);
        return 1.0-porosity_[eIdx];

    }

    /*!
     *  \brief Define the reference porosity \f$[-]\f$ distribution.
     *  This is the porosity of the porous medium without any of the
     *  considered solid phases.
     *
     *  \param element The finite element
     *  \param scv The sub-control volume
     */
    Scalar referencePorosity(const Element& element, const SubControlVolume &scv) const
    {
        unsigned int eIdx = this->fvGridGeometry().elementMapper().index(element);
        return porosity_[eIdx];

    }
//     { return referencePorosity_; }

    /*! Intrinsic permeability tensor K \f$[m^2]\f$ depending
     *  on the position in the domain
     *
     *  \param element The finite volume element
     *  \param scv The sub-control volume
     *
     *  Solution dependent permeability function
     */
    template<class ElementSolution>
    PermeabilityType permeability(const Element& element,
                                  const SubControlVolume& scv,
                                  const ElementSolution& elemSol) const
    {
        auto priVars = evalSolution(element, element.geometry(), elemSol, scv.center());

        Scalar sumPrecipitates = 0.0;
        sumPrecipitates += priVars[3 /*numComp*/];

         using std::max;
         unsigned int eIdx = this->fvGridGeometry().elementMapper().index(element);

         const auto poro =  max(/*minPoro*/1e-5, porosity_[eIdx] - sumPrecipitates);

         return permLaw_.evaluatePermeability(referencePermeability_, referencePorosity_, poro);
    }

//     Scalar solidity(const SubControlVolume &scv) const
//     { return 1.0 - porosityAtPos(scv.center()); }

    Scalar solubilityLimit() const
    { return solubilityLimit_; }

    Scalar theta(const SubControlVolume &scv) const
    { return 10.0; }

//     // return the brooks-corey context depending on the position
//     const MaterialLawParams& materialLawParamsAtPos(const GlobalPosition& globalPos) const
//     { return materialParams_; }

    template<class ElementSolution>
    const MaterialLawParams& materialLawParams(const Element& element,
                                               const SubControlVolume& scv,
                                               const ElementSolution& elemSol) const
    {
        unsigned int eIdx = this->fvGridGeometry().elementMapper().index(element);
        return materialParams_[eIdx];
    }



    // define which phase is to be considered as the wetting phase
    template<class FluidSystem>
    int wettingPhaseAtPos(const GlobalPosition& globalPos) const
    { return FluidSystem::H2OIdx; }

private:

    std::vector<MaterialLawParams> materialParams_;

    PermeabilityKozenyCarman<PermeabilityType> permLaw_;

    std::vector<Scalar> porosity_;
    Scalar solubilityLimit_;
    Scalar referencePorosity_;
    Scalar referencePermeability_;
    Scalar sdevPorosity_;
    std::vector<PermeabilityType> permeability_;
    Scalar irreducibleLiqSat_;
    Scalar irreducibleGasSat_;
    Scalar vgAlpha_;
    Scalar vgn_ ;

    bool plotFluidMatrixInteractions_;
};

} // end namespace Dumux

#endif
