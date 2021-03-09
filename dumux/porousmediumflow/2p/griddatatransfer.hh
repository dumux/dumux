// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*****************************************************************************
 *   See the file COPYING for full copying permissions.                      *
 *                                                                           *
 *   This program is free software: you can redistribute it and/or modify    *
 *   it under the terms of the GNU General Public License as published by    *
 *   the Free Software Foundation, either version 3 of the License, or       *
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
 * \ingroup TwoPModel
 * \brief Performs the transfer of data on a grid from before to after adaptation.
 */

#ifndef DUMUX_TWOP_GRIDDATA_TRANSFER_HH
#define DUMUX_TWOP_GRIDDATA_TRANSFER_HH

#include <memory>

#include <dune/grid/common/partitionset.hh>
#include <dune/grid/utility/persistentcontainer.hh>

#include <dumux/common/properties.hh>

#include <dumux/discretization/method.hh>
#include <dumux/discretization/extrusion.hh>
#include <dumux/discretization/elementsolution.hh>
#include <dumux/porousmediumflow/2p/formulation.hh>
#include <dumux/adaptive/griddatatransfer.hh>

namespace Dumux {

/*!
 * \ingroup TwoPModel
 * \brief Class performing the transfer of data on a grid from before to after adaptation.
 */
template<class TypeTag>
class TwoPGridDataTransfer : public GridDataTransfer
{
    using Grid = GetPropType<TypeTag, Properties::Grid>;
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using Problem = GetPropType<TypeTag, Properties::Problem>;
    using GridGeometry = GetPropType<TypeTag, Properties::GridGeometry>;
    using FVElementGeometry = typename GridGeometry::LocalView;
    using SubControlVolume = typename FVElementGeometry::SubControlVolume;
    using Extrusion = Extrusion_t<GridGeometry>;
    using GridVariables = GetPropType<TypeTag, Properties::GridVariables>;
    using VolumeVariables = GetPropType<TypeTag, Properties::VolumeVariables>;
    using PrimaryVariables = GetPropType<TypeTag, Properties::PrimaryVariables>;
    using SolutionVector = GetPropType<TypeTag, Properties::SolutionVector>;
    using Element = typename Grid::template Codim<0>::Entity;
    using ElementSolution = std::decay_t<decltype(elementSolution(std::declval<Element>(),
                                                                  std::declval<SolutionVector>(),
                                                                  std::declval<GridGeometry>()))>;
    using FluidSystem = GetPropType<TypeTag, Properties::FluidSystem>;
    using ModelTraits = GetPropType<TypeTag, Properties::ModelTraits>;
    using Indices = typename ModelTraits::Indices;

    struct AdaptedValues
    {
        AdaptedValues() : associatedMass(0.0) {}
        ElementSolution u;
        int count = 0;
        PrimaryVariables associatedMass;
        bool wasLeaf = false;
    };

    using PersistentContainer = Dune::PersistentContainer<Grid, AdaptedValues>;

    static constexpr int dim = Grid::dimension;
    static constexpr int dimWorld = Grid::dimensionworld;
    static constexpr bool isBox = GetPropType<TypeTag, Properties::GridGeometry>::discMethod == DiscretizationMethod::box;

    // saturation primary variable index
    enum { saturationIdx = Indices::saturationIdx };

    // phase indices
    enum
    {
        phase0Idx = FluidSystem::phase0Idx,
        phase1Idx = FluidSystem::phase1Idx,
    };

    // formulations
    static constexpr auto p0s1 = TwoPFormulation::p0s1;
    static constexpr auto p1s0 = TwoPFormulation::p1s0;

    // the formulation that is actually used
    static constexpr auto formulation = ModelTraits::priVarFormulation();

    // This won't work (mass conservative) for compressible fluids
    static_assert(!FluidSystem::isCompressible(phase0Idx)
                  && !FluidSystem::isCompressible(phase1Idx),
                  "This adaption helper is only mass conservative for incompressible fluids!");

    // check if the used formulation is implemented here
    static_assert(formulation == p0s1 || formulation == p1s0, "Chosen formulation not known to the TwoPGridDataTransfer");

public:
    /*!
     * \brief Constructor
     *
     * \param problem The DuMuX problem to be solved
     * \param gridGeometry The finite volume grid geometry
     * \param gridVariables The secondary variables on the grid
     * \param sol The solution (primary variables) on the grid
     */
    TwoPGridDataTransfer(std::shared_ptr<const Problem> problem,
                         std::shared_ptr<GridGeometry> gridGeometry,
                         std::shared_ptr<const GridVariables> gridVariables,
                         SolutionVector& sol)
    : GridDataTransfer()
    , problem_(problem)
    , gridGeometry_(gridGeometry)
    , gridVariables_(gridVariables)
    , sol_(sol)
    , adaptionMap_(gridGeometry->gridView().grid(), 0)
    {}

    /*!
     * \brief Stores primary variables and additional data
     *
     * To reconstruct the solution in father elements, problem properties might
     * need to be accessed. From upper level on downwards, the old solution is stored
     * into a container object, before the grid is adapted. Father elements hold averaged
     * information from the son cells for the case of the sons being coarsened.
     */
    void store() override
    {
        adaptionMap_.resize();

        const auto& grid = gridGeometry_->gridView().grid();
        for (auto level = grid.maxLevel(); level >= 0; level--)
        {
            for (const auto& element : elements(grid.levelGridView(level)))
            {
                // get map entry
                auto& adaptedValues = adaptionMap_[element];

                // put values in the map for leaf elements
                if (element.isLeaf())
                {
                    auto fvGeometry = localView(*gridGeometry_);
                    fvGeometry.bindElement(element);

                    // store current element solution
                    adaptedValues.u = ElementSolution(element, sol_, *gridGeometry_);

                    // compute mass in the scvs
                    for (const auto& scv : scvs(fvGeometry))
                    {
                        VolumeVariables volVars;
                        volVars.update(adaptedValues.u, *problem_, element, scv);

                        const auto poreVolume = Extrusion::volume(scv)*volVars.porosity();
                        adaptedValues.associatedMass[phase1Idx] += poreVolume * volVars.density(phase1Idx) * volVars.saturation(phase1Idx);
                        adaptedValues.associatedMass[phase0Idx] += poreVolume * volVars.density(phase0Idx) * volVars.saturation(phase0Idx);
                    }

                    // leaf elements always start with count = 1
                    adaptedValues.count = 1;
                    adaptedValues.wasLeaf = true;
                }
                // Average in father elements
                if (element.level() > 0)
                {
                    auto& adaptedValuesFather = adaptionMap_[element.father()];
                    // For some grids the father element is identical to the son element.
                    // In that case averaging is not necessary.
                    if(&adaptedValues != &adaptedValuesFather)
                        storeAdaptionValues(adaptedValues, adaptedValuesFather);
                }

                // The vertices of the non-leaf elements exist on the leaf as well
                // This element solution constructor uses the vertex mapper to obtain
                // the privars at the vertices, thus, this works for non-leaf elements!
                if(isBox && !element.isLeaf())
                    adaptedValues.u = ElementSolution(element, sol_, *gridGeometry_);
            }
        }
    }

    /*!
     * \brief Reconstruct missing primary variables (where elements are created/deleted)
     *
     * To reconstruct the solution in father elements, problem properties might
     * need to be accessed.
     * Starting from the lowest level, the old solution is mapped on the new grid:
     * Where coarsened, new cells get information from old father element.
     * Where refined, a new solution is reconstructed from the old father cell,
     * and then a new son is created. That is then stored into the general data
     * structure (AdaptedValues).
     */
    void reconstruct() override
    {
        // resize stuff (grid might have changed)
        adaptionMap_.resize();
        gridGeometry_->update();
        sol_.resize(gridGeometry_->numDofs());

        // vectors storing the mass associated with each vertex, when using the box method
        std::vector<Scalar> massCoeff;
        std::vector<Scalar> associatedMass;

        if(isBox)
        {
            massCoeff.resize(gridGeometry_->numDofs(), 0.0);
            associatedMass.resize(gridGeometry_->numDofs(), 0.0);
        }

        // iterate over leaf and reconstruct the solution
        for (const auto& element : elements(gridGeometry_->gridView().grid().leafGridView(), Dune::Partitions::interior))
        {
            if (!element.isNew())
            {
                const auto& adaptedValues = adaptionMap_[element];

                auto fvGeometry = localView(*gridGeometry_);
                fvGeometry.bindElement(element);

                // obtain element solution from map (divide by count!)
                auto elemSol = adaptedValues.u;
                if (!isBox)
                    elemSol[0] /= adaptedValues.count;

                const auto elementVolume = Extrusion::volume(element.geometry());
                for (const auto& scv : scvs(fvGeometry))
                {
                    VolumeVariables volVars;
                    volVars.update(elemSol, *problem_, element, scv);

                    // write solution at dof in current solution vector
                    sol_[scv.dofIndex()] = elemSol[scv.localDofIndex()];

                    const auto dofIdxGlobal = scv.dofIndex();
                    // For cc schemes, overwrite the saturation by a mass conservative one here
                    if (!isBox)
                    {
                        // only recalculate the saturations if element hasn't been leaf before adaptation
                        if (!adaptedValues.wasLeaf)
                        {
                            if (formulation == p0s1)
                            {
                                sol_[dofIdxGlobal][saturationIdx] = adaptedValues.associatedMass[phase1Idx];
                                sol_[dofIdxGlobal][saturationIdx] /= elementVolume * volVars.density(phase1Idx) * volVars.porosity();
                            }
                            else if (formulation == p1s0)
                            {
                                sol_[dofIdxGlobal][saturationIdx] = adaptedValues.associatedMass[phase0Idx];
                                sol_[dofIdxGlobal][saturationIdx] /= elementVolume * volVars.density(phase0Idx) * volVars.porosity();
                            }
                        }
                    }

                    // For the box scheme, add mass & mass coefficient to container (saturations are recalculated at the end)
                    else
                    {
                        const auto scvVolume = Extrusion::volume(scv);
                        if (formulation == p0s1)
                        {
                            massCoeff[dofIdxGlobal] += scvVolume * volVars.density(phase1Idx) * volVars.porosity();
                            associatedMass[dofIdxGlobal] += scvVolume / elementVolume * adaptedValues.associatedMass[phase1Idx];
                        }
                        else if (formulation == p1s0)
                        {
                            massCoeff[dofIdxGlobal] += scvVolume * volVars.density(phase0Idx) * volVars.porosity();
                            associatedMass[dofIdxGlobal] += scvVolume / elementVolume * adaptedValues.associatedMass[phase0Idx];
                        }
                    }
                }
            }
            else
            {
                // value is not in map, interpolate from father element
                assert(element.hasFather() && "new element does not have a father element!");

                // find the ancestor element that existed on the old grid already
                auto fatherElement = element.father();
                while(fatherElement.isNew() && fatherElement.level() > 0)
                    fatherElement = fatherElement.father();

                if(!isBox)
                {
                    const auto& adaptedValuesFather = adaptionMap_[fatherElement];

                    // obtain the mass contained in father
                    Scalar massFather = 0.0;
                    if (formulation == p0s1)
                        massFather = adaptedValuesFather.associatedMass[phase1Idx];
                    else if (formulation == p1s0)
                        massFather = adaptedValuesFather.associatedMass[phase0Idx];

                    // obtain the element solution through the father
                    auto elemSolSon = adaptedValuesFather.u;
                    elemSolSon[0] /= adaptedValuesFather.count;

                    auto fvGeometry = localView(*gridGeometry_);
                    fvGeometry.bindElement(element);

                    for (const auto& scv : scvs(fvGeometry))
                    {
                        VolumeVariables volVars;
                        volVars.update(elemSolSon, *problem_, element, scv);

                        // store constructed values of son in the current solution
                        sol_[scv.dofIndex()] = elemSolSon[0];

                        // overwrite the saturation by a mass conservative one here
                        Scalar massCoeffSon = 0.0;
                        if (formulation == p0s1)
                            massCoeffSon = Extrusion::volume(scv) * volVars.density(phase1Idx) * volVars.porosity();
                        else if (formulation == p1s0)
                            massCoeffSon = Extrusion::volume(scv) * volVars.density(phase0Idx) * volVars.porosity();
                        sol_[scv.dofIndex()][saturationIdx] =
                            ( Extrusion::volume(scv)/Extrusion::volume(fatherElement.geometry())*massFather )/massCoeffSon;
                    }
                }
                else
                {
                    auto& adaptedValuesFather = adaptionMap_[fatherElement];

                    auto fvGeometry = localView(*gridGeometry_);
                    fvGeometry.bindElement(element);

                    // interpolate solution in the father to the vertices of the new son
                    ElementSolution elemSolSon(element, sol_, *gridGeometry_);
                    const auto fatherGeometry = fatherElement.geometry();
                    for (const auto& scv : scvs(fvGeometry))
                        elemSolSon[scv.localDofIndex()] = evalSolution(fatherElement,
                                                                        fatherGeometry,
                                                                        adaptedValuesFather.u,
                                                                        scv.dofPosition());

                    // compute mass & mass coeffients for the scvs (saturations are recalculated at the end)
                    const auto fatherElementVolume = Extrusion::volume(fatherGeometry);
                    for (const auto& scv : scvs(fvGeometry))
                    {
                        VolumeVariables volVars;
                        volVars.update(elemSolSon, *problem_, element, scv);

                        const auto dofIdxGlobal = scv.dofIndex();
                        const auto scvVolume = Extrusion::volume(scv);
                        if (formulation == p0s1)
                        {
                            massCoeff[dofIdxGlobal] += scvVolume * volVars.density(phase1Idx) * volVars.porosity();
                            associatedMass[dofIdxGlobal] += scvVolume / fatherElementVolume * adaptedValuesFather.associatedMass[phase1Idx];
                        }
                        else if (formulation == p1s0)
                        {
                            massCoeff[dofIdxGlobal] += scvVolume * volVars.density(phase0Idx) * volVars.porosity();
                            associatedMass[dofIdxGlobal] += scvVolume / fatherElementVolume * adaptedValuesFather.associatedMass[phase0Idx];
                        }

                        // store constructed (pressure) values of son in the current solution (saturation comes later)
                        sol_[dofIdxGlobal] = elemSolSon[scv.localDofIndex()];
                    }
                }
            }
        }

        if(isBox)
        {
            for(std::size_t dofIdxGlobal = 0; dofIdxGlobal < gridGeometry_->numDofs(); dofIdxGlobal++)
                sol_[dofIdxGlobal][saturationIdx] = associatedMass[dofIdxGlobal] / massCoeff[dofIdxGlobal];
        }

        // reset entries in adaptation map
        adaptionMap_.resize( typename PersistentContainer::Value() );
        adaptionMap_.shrinkToFit();
        adaptionMap_.fill( typename PersistentContainer::Value() );

//! TODO: fix adaptive simulations in parallel
//#if HAVE_MPI
//        // communicate ghost data
//        using SolutionTypes = typename GetProp<TypeTag, SolutionTypes>;
//        using ElementMapper = typename SolutionTypes::ElementMapper;
//        using DataHandle = VectorExchange<ElementMapper, std::vector<CellData> >;
//        DataHandle dataHandle(problem.elementMapper(), this->cellDataGlobal());
//        problem.gridView().template communicate<DataHandle>(dataHandle,
//                                                            Dune::InteriorBorder_All_Interface,
//                                                            Dune::ForwardCommunication);
//#endif
    }

  private:

    /*!
     * \brief Stores sons entries into father element for averaging
     *
     * Sum up the adaptedValues (sons values) into father element. We store from leaf
     * upwards, so sons are stored first, then cells on the next leaf (=fathers)
     * can be averaged.
     *
     * \param adaptedValues Container for model-specific values to be adapted
     * \param adaptedValuesFather Values to be adapted of father cell
     */
    static void storeAdaptionValues(AdaptedValues& adaptedValues,
                                    AdaptedValues& adaptedValuesFather)
    {
        // Add associated mass of the child to the one of the father
        adaptedValuesFather.associatedMass += adaptedValues.associatedMass;

        if(!isBox)
        {
            // add the child's primary variables to the ones of father
            // we have to divide the child's ones in case it was composed
            // of several children as well!
            auto values = adaptedValues.u[0];
            values /= adaptedValues.count;
            adaptedValuesFather.u[0] += values;

            // keep track of the number of children that composed this father
            adaptedValuesFather.count += 1;

            // A father element is never leaf
            adaptedValuesFather.wasLeaf = false;
        }
        else
        {
            // For the box scheme, scaling of primary variables by count is obsolete
            // Thus, we always want count = 1
            adaptedValuesFather.count = 1;

            // A father element is never leaf
            adaptedValuesFather.wasLeaf = false;
        }
    }

    std::shared_ptr<const Problem> problem_;
    std::shared_ptr<GridGeometry> gridGeometry_;
    std::shared_ptr<const GridVariables> gridVariables_;
    SolutionVector& sol_;
    PersistentContainer adaptionMap_;
};

} // end namespace Dumux

#endif
