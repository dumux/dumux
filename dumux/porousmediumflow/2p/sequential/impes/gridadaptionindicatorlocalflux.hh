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
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the           *
 *   GNU General Public License for more details.                            *
 *                                                                           *
 *   You should have received a copy of the GNU General Public License       *
 *   along with this program.  If not, see <http://www.gnu.org/licenses/>.   *
 *****************************************************************************/
/*!
 * \file
 * \ingroup SequentialTwoPModel
 * \brief  Class defining a standard, saturation dependent indicator for grid adaption
 */
#ifndef DUMUX_GRIDADAPTIONINDICATOR2PLOCALFLUX_HH
#define DUMUX_GRIDADAPTIONINDICATOR2PLOCALFLUX_HH

#include <dumux/porousmediumflow/sequential/impetproperties.hh>
#include <dumux/porousmediumflow/2p/sequential/properties.hh>

namespace Dumux
{

/*!
 * \brief  Class defining a standard, saturation dependent indicator for grid adaption
 * \ingroup SequentialTwoPModel
 *
 * \tparam TypeTag The problem TypeTag
 */
template<class TypeTag>
class GridAdaptionIndicator2PLocalFlux
{
private:
    using Problem = typename GET_PROP_TYPE(TypeTag, Problem);
    using GridView = typename GET_PROP_TYPE(TypeTag, GridView);
    using Scalar = typename GET_PROP_TYPE(TypeTag, Scalar);
    using Element = typename GridView::Traits::template Codim<0>::Entity;

    using Grid = typename GET_PROP_TYPE(TypeTag, Grid);
    using IndexSet = typename Grid::LevelGridView::IndexSet;

    using SolutionTypes = typename GET_PROP(TypeTag, SolutionTypes);
    using ScalarSolutionType = typename SolutionTypes::ScalarSolution;

    using PrimaryVariables = typename GET_PROP_TYPE(TypeTag, PrimaryVariables);
    using BoundaryTypes = typename GET_PROP_TYPE(TypeTag, BoundaryTypes);
    using CellData = typename GET_PROP_TYPE(TypeTag, CellData);

    using Indices = typename GET_PROP_TYPE(TypeTag, Indices);

    enum
    {
        dim = GridView::dimension, dimWorld = GridView::dimensionworld
    };

    enum
    {
        Sw = Indices::saturationW,
        Sn = Indices::saturationNw,
        eqIdxSat = Indices::satEqIdx
    };
    enum
    {
        wPhaseIdx = Indices::wPhaseIdx,
        nPhaseIdx = Indices::nPhaseIdx,
        numPhases = GET_PROP_VALUE(TypeTag, NumPhases)
    };

    using ReferenceElements = Dune::ReferenceElements<Scalar, dim>;
    using GlobalPosition = Dune::FieldVector<Scalar, dimWorld>;
    using DimVector = Dune::FieldVector<Scalar, dim>;
    using DimMatrix = Dune::FieldMatrix<Scalar, dim, dim>;
    using SetVector = Dune::FieldVector<Scalar, 2>;

    struct SetField {
        Scalar indicator;
        Scalar volume;
        int index;
        SetField()
            :indicator(0),volume(0), index(0)
        {}
    };

    struct Comparison {
        bool operator() (const SetField& lhs, const SetField& rhs) const
        {return lhs.indicator<rhs.indicator;}
    };

    using RangeSet = std::set<SetField, Comparison>;

public:
    /*!
     * \brief Calculates the indicator used for refinement/coarsening for each grid cell.
     *
     * This standard indicator is based on the saturation gradient.
     */
    void calculateIndicator()
    {

        int size = problem_.variables().cellDataGlobal().size();

        // prepare an indicator for refinement
        indicatorVector_.resize(size);

        if (useSatInd_ || usePercentileSat_)
            indicatorVectorSat_.resize(size);
        if (useFluxInd_ || usePercentileFlux_)
            indicatorVectorFlux_.resize(size);


        RangeSet satRange;
        RangeSet fluxRange;

        SetField satVec;
        SetField fluxVec;

        Scalar totalVolume = 0;
        Scalar totalVolumeSat = 0;

        // 1) calculate Indicator -> min, maxvalues
        // Schleife über alle Leaf-Elemente
        for (const auto& element : elements(problem_.gridView()))
        {
            // Bestimme maximale und minimale Sättigung
            // Index des aktuellen Leaf-Elements
            int globalIdxI = problem_.variables().index(element);

            indicatorVector_[globalIdxI] = 0.5 * (refineBound_ + coarsenBound_);
            if (useSatInd_ || usePercentileSat_)
                indicatorVectorSat_[globalIdxI] = -1e100;
            if (useFluxInd_ || usePercentileFlux_)
                indicatorVectorFlux_[globalIdxI] = -1e100;

            const CellData& cellDataI = problem_.variables().cellData(globalIdxI);

            bool isSpecialCell = false;

            if (!checkPhysicalRange_(cellDataI))
            {
                indicatorVector_[globalIdxI] = refineBound_ + 1.0;
                isSpecialCell = true;
            }

            Scalar volume = element.geometry().volume();
            totalVolume += volume;

            if (refineAtSource_)
            {
                PrimaryVariables source(0.0);
                problem_.sourceAtPos(source, element.geometry().center());
                for (int i = 0; i < 2; i++)
                {
                    using std::abs;
                    if (abs(source[i]) > 1e-10)
                    {
                        indicatorVector_[globalIdxI] = refineBound_ + 1.0;
                        isSpecialCell = true;
                        break;
                    }
                }
            }

            if (isSpecialCell)
            {
                continue;
            }

            Scalar satI = 0.0;
            Scalar satW = 0.0;
            switch (saturationType_)
            {
            case Sw:
                satI = cellDataI.saturation(wPhaseIdx);
                satW = satI;
                break;
            case Sn:
                satI = cellDataI.saturation(nPhaseIdx);
                satW = 1.0-satI;
                break;
            }

            const typename Element::Geometry& geometry = element.geometry();
            // get corresponding reference element
            using ReferenceElements = Dune::ReferenceElements<Scalar, dim>;
            const Dune::ReferenceElement< Scalar , dim > & refElement =
                    ReferenceElements::general( geometry.type() );
            const int numberOfFaces=refElement.size(1);

            std::vector<Scalar> flux(numberOfFaces,0);

            // Calculate refinement indicator on all cells
            for (const auto& intersection : intersections(problem_.gridView(), element))
            {
                if (isSpecialCell)
                {
                    break;
                }

                if (intersection.neighbor())
                {
                const CellData& cellDataJ = problem_.variables().cellData(problem_.variables().index(intersection.outside()));
                if (!checkPhysicalRange_(cellDataJ))
                {
                    indicatorVector_[globalIdxI] = refineBound_ + 1.0;
                    isSpecialCell = true;
                    break;
                }
                }

                int idxInInside = intersection.indexInInside();

                if (useFluxInd_ || usePercentileFlux_)
                {
                    int isIndex = intersection.indexInInside();
                    flux[isIndex] += (intersection.centerUnitOuterNormal()
                                      * cellDataI.fluxData().velocityTotal(idxInInside)) * intersection.geometry().volume();

                    //Scalar velNorm = cellDataI.fluxData().velocityTotal(idxInInside).two_norm();
                    //using std::max;
                    //indicatorVectorFlux_[globalIdxI] = max(velNorm, indicatorVectorFlux_[globalIdxI]);
                }

                // exit, if it is not a neighbor
                if (intersection.boundary())
                {
                    BoundaryTypes bcTypes;
                    problem_.boundaryTypes(bcTypes, intersection);

                    for (int i = 0; i < 2; i++)
                    {
                        if (bcTypes.isNeumann(i))
                        {
                            PrimaryVariables flux(0.0);
                            problem_.neumann(flux, intersection);

                            bool fluxBound = false;
                            for (int j = 0; j < 2; j++)
                            {
                                using std::abs;
                                if (abs(flux[j]) > 1e-10)
                                {
                                    if (refineAtFluxBC_)
                                    {
                                        indicatorVector_[globalIdxI] = refineBound_ + 1.0;
                                        isSpecialCell = true;
                                        fluxBound = true;
                                        break;
                                    }
                                }
                            }
                            if (fluxBound)
                                break;
                        }
                        else if (bcTypes.isDirichlet(i))
                        {
                            if (refineAtDirichletBC_)
                            {
                                indicatorVector_[globalIdxI] = refineBound_ + 1.0;
                                isSpecialCell = true;
                                break;
                            }
                        }
                    }
                    if (useSatInd_ || usePercentileSat_)
                    {
                        Scalar satJ = satI;
                        if (bcTypes.isDirichlet(eqIdxSat))
                        {
                            PrimaryVariables sat(0.0);
                            problem_.dirichlet(sat, intersection);
                            satJ = sat[eqIdxSat];
                        }
                        using std::abs;
                        Scalar localdelta = abs(satI - satJ);
                        using std::max;
                        indicatorVectorSat_[globalIdxI] = max(indicatorVectorSat_[globalIdxI], localdelta);
                    }
                }
                else
                {
                    if (useSatInd_ || usePercentileSat_)
                    {
                        // get neighbors
                        int globalIdxJ = problem_.variables().index(intersection.outside());

                        Scalar satJ = 0.;
                        switch (saturationType_)
                        {
                        case Sw:
                            satJ = problem_.variables().cellData(globalIdxJ).saturation(wPhaseIdx);
                            break;
                        case Sn:
                            satJ = problem_.variables().cellData(globalIdxJ).saturation(nPhaseIdx);
                            break;
                        }

                        using std::abs;
                        Scalar localdelta = abs(satI - satJ);
                        using std::max;
                        indicatorVectorSat_[globalIdxI] = max(indicatorVectorSat_[globalIdxI], localdelta);
                    }
                }
            }

            if (isSpecialCell)
            {
                continue;
            }

            if (useFluxInd_ || usePercentileFlux_)
            {
                // calculate velocity on reference element as the Raviart-Thomas-0
                // interpolant of the fluxes
                Dune::FieldVector<Scalar, dim> refVelocity;
                // simplices
                if (refElement.type().isSimplex()) {
                    for (int dimIdx = 0; dimIdx < dim; dimIdx++)
                    {
                        refVelocity[dimIdx] = -flux[dim - 1 - dimIdx];
                        for (int fIdx = 0; fIdx < dim + 1; fIdx++)
                        {
                            refVelocity[dimIdx] += flux[fIdx]/(dim + 1);
                        }
                    }
                }
                // cubes
                else if (refElement.type().isCube()){
                    for (int i = 0; i < dim; i++)
                        refVelocity[i] = 0.5 * (flux[2*i + 1] - flux[2*i]);
                }
                // 3D prism and pyramids
                else {
                    DUNE_THROW(Dune::NotImplemented, "velocity output for prism/pyramid not implemented");
                }

                const DimVector& localPos = refElement.position(0, 0);

                // get the transposed Jacobian of the element mapping
                const DimMatrix jacobianT = geometry.jacobianTransposed(localPos);

                // calculate the element velocity by the Piola transformation
                DimVector elementVelocity(0);
                jacobianT.umtv(refVelocity, elementVelocity);
                elementVelocity /= geometry.integrationElement(localPos);

                Scalar velNorm = elementVelocity.two_norm();
                using std::max;
                indicatorVectorFlux_[globalIdxI] = max(velNorm, indicatorVectorFlux_[globalIdxI]);
            }


            if (useFluxInd_ || usePercentileFlux_)
            {
                fluxVec.indicator = indicatorVectorFlux_[globalIdxI];
                fluxVec.volume = volume;
                fluxVec.index = globalIdxI;
                fluxRange.insert(fluxVec);
            }

            if (useSatInd_ || usePercentileSat_)
            {
                satVec.indicator = indicatorVectorSat_[globalIdxI];
                if (satVec.indicator > 1e-8)
                {
                satVec.volume = volume;
                satVec.index = globalIdxI;
                satRange.insert(satVec);

                    totalVolumeSat += volume;
                }
            }
        }

        if (useSatInd_ || usePercentileSat_)
        {
            int range = satRange.size();
            if (range > 0)
            {
                typename RangeSet::iterator it = satRange.begin();
                Scalar minLocalDelta = (*it).indicator;

                typename RangeSet::reverse_iterator rIt = satRange.rbegin();
                Scalar maxLocalDelta = (*rIt).indicator;

                if (maxLocalDelta > 0.)
                {
                    Scalar globalDeltaDelta = maxLocalDelta - minLocalDelta;
                    for (int i = 0; i < size; i++)
                    {
                        indicatorVectorSat_[i] -= minLocalDelta;
                        indicatorVectorSat_[i] /= globalDeltaDelta;
                    }
                }

                if (usePercentileSat_)
                {
                    using std::max;
                    Scalar lowerBound = max(0.0,coarsenPercentileSat_ * totalVolumeSat);

                    Scalar accumulatedVolume = 0;
                    while (accumulatedVolume <= lowerBound && it != satRange.end())
                    {
                        indicatorVector_[(*it).index] = coarsenBound_-1;
                        accumulatedVolume += (*it).volume;
                        ++it;
                    }
                }
            }
        }

        if (useFluxInd_  || usePercentileFlux_)
        {
            int range = fluxRange.size();
            if (range > 0)
            {
                typename RangeSet::iterator it = fluxRange.begin();
                Scalar minFlux = (*it).indicator;

                typename RangeSet::reverse_iterator rIt = fluxRange.rbegin();
                Scalar maxFlux = (*rIt).indicator;

                if (maxFlux > 0.)
                {
                    Scalar globalDeltaDelta = maxFlux - minFlux;
                    for (int i = 0; i < size; i++)
                    {
                        indicatorVectorFlux_[i] -= minFlux;
                        indicatorVectorFlux_[i] /= globalDeltaDelta;
                    }
                }

                if (usePercentileFlux_)
                {
                    using std::max;
                    Scalar lowerBound = max(0.0,coarsenPercentileFlux_ * totalVolume);
                    Scalar accumulatedVolume = 0;
                    while (accumulatedVolume <= lowerBound && it != fluxRange.end())
                    {
                        indicatorVector_[(*it).index] = coarsenBound_-1;
                        accumulatedVolume += (*it).volume;
                        ++it;
                    }
                }
            }
        }

        if (useSatInd_ || usePercentileSat_)
        {
            int range = satRange.size();
            if (range > 0)
            {
                if (usePercentileSat_)
                {
                    typename RangeSet::reverse_iterator rIt = satRange.rbegin();

                    using std::max;
                    Scalar upperBound = max(0.0, refinePercentileSat_ * totalVolumeSat);

                    Scalar accumulatedVolume = 0;
                    while (accumulatedVolume <= upperBound && (*rIt).indicator > 1e-6 && rIt != satRange.rend())
                    {
                        indicatorVector_[(*rIt).index] = refineBound_+1;
                        accumulatedVolume += (*rIt).volume;
                        ++rIt;
                    }
                }


            }
        }

        if (useFluxInd_ || usePercentileFlux_)
        {
            int range = fluxRange.size();
            if (range > 0)
            {
                typename RangeSet::reverse_iterator rIt = fluxRange.rbegin();

                if (usePercentileFlux_)
                {
                    using std::min;
                    Scalar upperBound = min(totalVolume, refinePercentileFlux_ * totalVolume);
                    Scalar accumulatedVolume = 0;
                    while (accumulatedVolume <= upperBound && rIt != fluxRange.rend())
                    {
                        indicatorVector_[(*rIt).index] = refineBound_+1;
                        accumulatedVolume += (*rIt).volume;
                        ++rIt;
                    }
                }
            }
        }

        for (int idx = 0; idx < size; idx++)
        {
            if (useSatInd_ && indicatorVectorSat_[idx] > refineThresholdSat_)
            {
                indicatorVector_[idx] = refineBound_+1;
            }
            else if (useFluxInd_  && indicatorVectorFlux_[idx] > refineThresholdFlux_)
            {
                indicatorVector_[idx] = refineBound_+1;
            }
            else if (useSatInd_ && indicatorVectorSat_[idx] < coarsenThresholdSat_ && indicatorVector_[idx] < refineBound_)
            {
                indicatorVector_[idx] = coarsenBound_-1;
            }
            else if (useFluxInd_  && indicatorVectorFlux_[idx] < coarsenThresholdFlux_ && indicatorVector_[idx] < refineBound_)
            {
                indicatorVector_[idx] = coarsenBound_-1;
            }
        }
    }

    /*!
     * \brief Check if pressure data is in some physical range
     *
     * Returns true if the cell pressure is in some range
     *
     *  \param cellData A grid element
     */
    bool checkPhysicalRange_(const CellData& cellData)
    {
        for (int j = 0; j < 2; j++)
        {
            if (cellData.pressure(j) < lowerPressureBound_ || cellData.pressure(j) > upperPressureBound_)
                return false;
        }
        return true;
    }

    /*!
     * \brief Set a lower and upper pressure constraint
     *
     *  \param lowerPressureBound lower pressure value
     *  \param upperPressureBound upper pressure value
     */
    void setPressureBounds(Scalar lowerPressureBound, Scalar upperPressureBound)
    {
        lowerPressureBound_ = lowerPressureBound;
        upperPressureBound_ = upperPressureBound;
    }

    /*!
     * \brief Indicator function for marking of grid cells for refinement
     *
     * Returns true if an element should be refined.
     *
     *  \param element A grid element
     */
    bool refine(const Element& element)
    {
        int idx = problem_.elementMapper().index(element);
        return (indicatorVector_[idx] > refineBound_);
    }

    /*!
     * \brief Indicator function for marking of grid cells for coarsening
     *
     * Returns true if an element should be coarsened.
     *
     *  \param element A grid element
     */
    bool coarsen(const Element& element)
    {
        int idx = problem_.elementMapper().index(element);
        return (indicatorVector_[idx] < coarsenBound_);
    }

    /*!
     * \brief Initializes the adaption indicator class
     */
    void init()
    {
        int size = problem_.gridView().size(0);

        indicatorVector_.resize(size, -1e100);
        if (useSatInd_)
            indicatorVectorSat_.resize(size, -1e100);
        if (useFluxInd_)
            indicatorVectorFlux_.resize(size, -1e100);
    }

    /*!
     * \brief Function for changing the indicatorVector values for refinement
     *
     *  \param idx Index of cell which will be refined
     */
    void setIndicatorToRefine(int idx)
    {
        indicatorVector_[idx] = refineBound_+1;
    }

    /*!
     * \brief Function for changing the indicatorVector values for coarsening
     *
     *  \param idx Index of cell which will be coarsen
     */
    void setIndicatorToCoarse(int idx)
    {
        indicatorVector_[idx] = coarsenBound_ - 1;
    }

    /*!
     * \brief Constructs a GridAdaptionIndicator instance
     *
     *  This standard indicator is based on the saturation gradient.
     *  It checks the local gradient compared to the maximum global gradient.
     *  The indicator is compared locally to a refinement/coarsening threshold to decide whether
     *  a cell should be marked for refinement or coarsening or should not be adapted.
     *
     * \param problem The problem object
     */
    GridAdaptionIndicator2PLocalFlux (Problem& problem):
        problem_(problem), lowerPressureBound_(0.0), upperPressureBound_(std::numeric_limits<Scalar>::max())
    {
        refineThresholdSat_ = getParam<Scalar>("GridAdapt.RefineThresholdSat", 0.8);
        coarsenThresholdSat_ = getParam<Scalar>("GridAdapt.CoarsenThresholdSat", 0.2);
        refineBound_ = 2.0;
        coarsenBound_ = 1.0;
        refineThresholdFlux_ = getParam<Scalar>("GridAdapt.RefineThresholdFlux", 0.8);
        coarsenThresholdFlux_ = getParam<Scalar>("GridAdapt.CoarsenThresholdFlux", 0.2);
        refinePercentileFlux_ = getParam<Scalar>("GridAdapt.RefinePercentileFlux", 0.8);
        coarsenPercentileFlux_ = getParam<Scalar>("GridAdapt.CoarsenPercentileFlux", 0.2);
        refinePercentileSat_ = getParam<Scalar>("GridAdapt.RefinePercentileSat", 0.8);
        coarsenPercentileSat_ = getParam<Scalar>("GridAdapt.CoarsenPercentileSat", 0.2);
        refineAtDirichletBC_ = getParam<bool>("GridAdapt.RefineAtDirichletBC");
        refineAtFluxBC_ = getParam<bool>("GridAdapt.RefineAtFluxBC");
        refineAtSource_ = getParam<bool>("GridAdapt.RefineAtSource");

        useSatInd_ = (refineThresholdSat_ < 1.0 || coarsenThresholdSat_ > 0.0);
        useFluxInd_ = (refineThresholdFlux_ < 1.0 || coarsenThresholdFlux_ > 0.0);
        usePercentileSat_ = (refinePercentileSat_ > 0.0 || coarsenPercentileSat_ > 0.0);
        usePercentileFlux_ = (refinePercentileFlux_ > 0.0 || coarsenPercentileFlux_ > 0.0);

        std::cout<<"--------------------------------------------\n";
        std::cout<<"Used adaption indicators:\n";
        if (useSatInd_)
            std::cout<<"    - delta S\n";
        if (useFluxInd_)
            std::cout<<"    - total flux\n";
        std::cout<<"\n Use percentiles:\n";
        if (usePercentileSat_)
            std::cout<<"    - delta S\n";
        if (usePercentileFlux_)
            std::cout<<"    - total flux\n";
        std::cout<<"--------------------------------------------\n";
    }

private:
    Problem& problem_;
    Scalar lowerPressureBound_;
    Scalar upperPressureBound_;
    bool useSatInd_;
    bool useFluxInd_;
    bool usePercentileSat_;
    bool usePercentileFlux_;
    Scalar refineThresholdSat_;
    Scalar coarsenThresholdSat_;
    Scalar refineBound_;
    Scalar coarsenBound_;
    Scalar refineThresholdFlux_;
    Scalar coarsenThresholdFlux_;
    Scalar refinePercentileFlux_;
    Scalar coarsenPercentileFlux_;
    Scalar refinePercentileSat_;
    Scalar coarsenPercentileSat_;
    bool refineAtDirichletBC_;
    bool refineAtFluxBC_;
    bool refineAtSource_;
    std::vector<Scalar> indicatorVector_;
    std::vector<Scalar> indicatorVectorSat_;
    std::vector<Scalar> indicatorVectorFlux_;
    static const int saturationType_ = GET_PROP_VALUE(TypeTag, SaturationFormulation);
    static const bool compressibility_ = GET_PROP_VALUE(TypeTag, EnableCompressibility);

};
}

#endif
