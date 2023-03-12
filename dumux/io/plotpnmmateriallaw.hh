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
 * \brief Interface for plotting the PNM fluid-matrix-interaction laws
 */
#ifndef DUMUX_PLOT_PNM_FLUID_MATRIX_LAW_HH
#define DUMUX_PLOT_PNM_FLUID_MATRIX_LAW_HH

#include <dune/istl/bvector.hh>
#include <dumux/io/gnuplotinterface.hh>
#include <dumux/porenetwork/2p/fluxvariablescache.hh>

namespace Dumux
{

/*!
 *\brief Interface for plotting the PNM fluid-matrix-interaction laws
 *
 */
template<class MaterialLaw, class Scalar>
class PlotLocalRules
{
    using MaterialLawParams = typename MaterialLaw::Params<Scalar>;

    template<class VV>
    class PseudoElemVolVars
    {
    public:

        using VolumeVariables = VV;

        template<class ElementSolutionVector, class Problem, class Element, class Scv>
        PseudoElemVolVars(const ElementSolutionVector& elemSol,
                          const Problem& problem,
                          const Element& element,
                          const Scv& scv)
        {
            volVars_[0].update(elemSol, problem, element, scv);
            volVars_[1].update(elemSol, problem, element, scv);
        }

         const VolumeVariables& operator [](std::size_t scvIdx) const
         { return volVars_[scvIdx]; }

         template<class SubControlVolume, typename std::enable_if_t<!std::is_integral<SubControlVolume>::value, int> = 0>
         const VolumeVariables& operator [](const SubControlVolume& scv) const
         { return volVars_[scv.indexInElement()]; }

         VolumeVariables& operator [](std::size_t scvIdx)
         { return volVars_[scvIdx]; }

         template<class SubControlVolume, typename std::enable_if_t<!std::is_integral<SubControlVolume>::value, int> = 0>
         VolumeVariables& operator [](const SubControlVolume& scv)
         { return volVars_[scv.indexInElement()]; }

    private:
        std::array<VolumeVariables, 2> volVars_;
    };

    template<class PV>
    class PseudoElemSol
    {
    public:
        using PrimaryVariables = PV;

        PseudoElemSol()
        {
            priVars_.resize(2);
        }

        //! bracket operator access
        template<typename IndexType>
        PrimaryVariables& operator [](IndexType i)
        { return priVars_[i]; }

        //! bracket operator access
        template<typename IndexType>
        const PrimaryVariables& operator [](IndexType i) const
        { return priVars_[i]; }

    private:
        Dune::BlockVector<PrimaryVariables> priVars_;
    };

public:
    //! Constructor
    PlotLocalRules()
    : numIntervals_(1000)
    { }

    /*!
     * \brief Plot the capillary pressure-saturation curve for the water-gas interphase
     *
     */
    template<class VolumeVariables, class Problem, class Element, class FVElementGeometry>
    void plotPcOverSw(const Problem& problem,
                      const Element& element,
                      const FVElementGeometry& fvGeometry,
                      const typename FVElementGeometry::SubControlVolume& scv,
                      const Scalar lowerSat = 0.0,
                      const Scalar upperSat = 1.0,
                      const std::string plotName = "")
    {
        PseudoElemSol<typename VolumeVariables::PrimaryVariables> elemSol;
        elemSol[0][VolumeVariables::Indices::pressureIdx] = 1e5;

        std::vector<Scalar> sw(numIntervals_+1);
        std::vector<Scalar> pc(numIntervals_+1);
        const Scalar satInterval = upperSat - lowerSat;

        for (int i = 0; i <= numIntervals_; i++)
        {
            sw[i] = lowerSat + satInterval * Scalar(i) / Scalar(numIntervals_);
            elemSol[0][VolumeVariables::Indices::saturationIdx] = sw[i];
            elemSol[1] = elemSol[0];
            VolumeVariables volVars;
            volVars.update(elemSol, problem, element, scv);
            pc[i] = volVars.capillaryPressure();
        }

        if (pc[0] < pc.back())
            std::reverse(pc.begin(), pc.end());

        const auto pcMinMax = std::minmax_element(pc.begin(), pc.end());

        // extend the x range by 10% on each side
        gnuplotPc_.setXRange(lowerSat - 0.1*satInterval, upperSat + 0.1*satInterval);
        gnuplotPc_.setYRange(*pcMinMax.first, *pcMinMax.second);
        gnuplotPc_.setXlabel("wetting phase saturation [-]");
        gnuplotPc_.setYlabel("capillary pressure [Pa]");
        gnuplotPc_.addDataSetToPlot(sw, pc, plotName + "_pore_" + std::to_string(scv.dofIndex()));
        gnuplotPc_.plot("pc-Sw");
    }


    // template<class VolumeVariables, class Problem, class ElementFluxVarsCache, class Element, class FVElementGeometry>
    // void plotTransmissibilities(const Problem& problem,
    //                             const Element& element,
    //                             const ElementFluxVarsCache& elementFluxVarsCache,
    //                             const FVElementGeometry& fvGeometry,
    //                             const Scalar lowerSat = 0.0,
    //                             const Scalar upperSat = 1.0,
    //                             const std::string plotName = "")
    // {
    //     PseudoElemSol<typename VolumeVariables::PrimaryVariables> elemSol;
    //     elemSol[0][VolumeVariables::Indices::pressureIdx] = 1e5;

    //     std::vector<Scalar> sw(numIntervals_+1);
    //     std::vector<Scalar> pc(numIntervals_+1);
    //     std::vector<Scalar> gwInvasion(numIntervals_+1);
    //     std::vector<Scalar> gwSnapOff(numIntervals_+1);
    //     std::vector<Scalar> gnInvasion(numIntervals_+1);
    //     std::vector<Scalar> gnSnapOff(numIntervals_+1);
    //     const Scalar satInterval = upperSat - lowerSat;

    //     const auto& scvf = fvGeometry.scvf(0);

    //     // take the smaller pore
    //     const auto& scv = [&]()
    //     {
    //         const auto& scv0 = fvGeometry.scv(0);
    //         const auto& scv1 = fvGeometry.scv(1);

    //         if (fvGeometry.gridGeometry().poreInscribedRadius(scv0.dofIndex()) < fvGeometry.gridGeometry().poreInscribedRadius(scv1.dofIndex()))
    //             return scv0;
    //         else
    //             return scv1;
    //     }();

    //     for (int i = 0; i <= numIntervals_; i++)
    //     {
    //         sw[i] = lowerSat + satInterval * Scalar(i) / Scalar(numIntervals_);
    //         elemSol[0][VolumeVariables::Indices::saturationIdx] = sw[i];
    //         elemSol[1] = elemSol[0];

    //         PseudoElemVolVars<VolumeVariables> elemVolVars(elemSol, problem, element, scv);
    //         // const auto localElemFluxVarsCache = localView(elementFluxVarsCache).bindElement(element, fvGeometry, elemVolVars);
    //         const auto& fluxVarsCache = elementFluxVarsCache[scvf];

    //         // consider invasion first
    //         pc[i] = std::max(elemVolVars[0].capillaryPressure(), elemVolVars[1].capillaryPressure());
    //         const Scalar pcEntry = problem.spatialParams().pcEntry(element, elemVolVars);
    //         bool invaded = pc[i] > pcEntry ? true : false;
    //         fluxVarsCache.update(problem, element, fvGeometry, elemVolVars, scvf, invaded);
    //         gwInvasion[i] = fluxVarsCache.transmissibility(fluxVarsCache.wPhaseIdx()) * elemVolVars[0].mobility(fluxVarsCache.wPhaseIdx());
    //         gnInvasion[i] = fluxVarsCache.transmissibility(fluxVarsCache.nPhaseIdx()) * elemVolVars[0].mobility(fluxVarsCache.nPhaseIdx());

    //         // consider snap-off
    //         const Scalar pcSnapoff = problem.spatialParams().pcSnapoff(element, elemVolVars);
    //         invaded = pc[i] > pcSnapoff ? true : false;
    //         fluxVarsCache.update(problem, element, fvGeometry, elemVolVars, scvf, invaded);
    //         gwSnapOff[i] = fluxVarsCache.transmissibility(fluxVarsCache.wPhaseIdx()) * elemVolVars[0].mobility(fluxVarsCache.wPhaseIdx());
    //         gnSnapOff[i] = fluxVarsCache.transmissibility(fluxVarsCache.nPhaseIdx()) * elemVolVars[0].mobility(fluxVarsCache.nPhaseIdx());
    //     }

    //     // make sure to always plot over S_w
    //     if (pc[0] < pc.back())
    //     {
    //         std::reverse(gwInvasion.begin(), gwInvasion.end());
    //         std::reverse(gwSnapOff.begin(), gwSnapOff.end());
    //         std::reverse(gnInvasion.begin(), gnInvasion.end());
    //         std::reverse(gnSnapOff.begin(), gnSnapOff.end());
    //     }

    //     // plot
    //     const auto eIdx = problem.gridGeometry().elementMapper().index(element);
    //     gnuplotTranssibilities_.setXRange(lowerSat, upperSat);
    //     gnuplotTranssibilities_.setXlabel("wetting phase saturation [-]");
    //     gnuplotTranssibilities_.setYlabel("throat conductivity");
    //     gnuplotTranssibilities_.setOption("set logscale y");
    //     gnuplotTranssibilities_.addDataSetToPlot(sw, gwInvasion, plotName + "_gwInvasion" + "_throat_" + std::to_string(eIdx));
    //     gnuplotTranssibilities_.addDataSetToPlot(sw, gwSnapOff, plotName + "_gwSnapOff"  + "_throat_" + std::to_string(eIdx), "with lines dashtype 2");
    //     gnuplotTranssibilities_.addDataSetToPlot(sw, gnInvasion, plotName + "_gnInvasion"  + "_throat_" + std::to_string(eIdx));
    //     gnuplotTranssibilities_.addDataSetToPlot(sw, gnSnapOff, plotName + "_gnSnapOff"  + "_throat_" + std::to_string(eIdx), "with lines dashtype 2");
    //     gnuplotTranssibilities_.plot("g-Sw");
    // }


private:
    std::size_t numIntervals_;
    GnuplotInterface<Scalar> gnuplotPc_;
    GnuplotInterface<Scalar> gnuplotTranssibilities_;

};
} // end of namespace

#endif // DUMUX_PLOT_PNM_FLUID_MATRIX_LAW_HH
