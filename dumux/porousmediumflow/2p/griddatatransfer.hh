// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup TwoPModel
 * \brief Performs the transfer of data on a grid from before to after adaptation.
 */

#ifndef DUMUX_TWOP_GRIDDATA_TRANSFER_HH
#define DUMUX_TWOP_GRIDDATA_TRANSFER_HH

#include <memory>

#include <dune/common/fvector.hh>
#include <dune/grid/common/capabilities.hh>
#include <dune/grid/common/partitionset.hh>
#include <dune/grid/utility/persistentcontainer.hh>

#include <dumux/common/properties.hh>

#include <dumux/discretization/method.hh>
#include <dumux/discretization/extrusion.hh>
#include <dumux/discretization/elementsolution.hh>
#include <dumux/geometry/volume.hh>
#include <dumux/parallel/vectorcommdatahandle.hh>
#include <dumux/porousmediumflow/2p/formulation.hh>
#include <dumux/adaptive/griddatatransfer.hh>

namespace Dumux {

/*!
 * \ingroup TwoPModel
 * \brief Class performing the transfer of data on a grid from before to after adaptation.
 */
template<class TypeTag>
class TwoPGridDataTransfer : public GridDataTransfer<GetPropType<TypeTag, Properties::Grid>>
{
    using Grid = GetPropType<TypeTag, Properties::Grid>;
    using ParentType = GridDataTransfer<Grid>;
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using Problem = GetPropType<TypeTag, Properties::Problem>;
    using GridGeometry = GetPropType<TypeTag, Properties::GridGeometry>;
    using GridView = typename GridGeometry::GridView;
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
    static constexpr bool isBox = GetPropType<TypeTag, Properties::GridGeometry>::discMethod == DiscretizationMethods::box;

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
    : ParentType()
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
    void store(const Grid& grid) override
    {
        adaptionMap_.resize();

        for (auto level = grid.maxLevel(); level >= 0; level--)
        {
            auto fvGeometry = localView(*gridGeometry_);
            for (const auto& element : elements(grid.levelGridView(level)))
            {
                // get map entry
                auto& adaptedValues = adaptionMap_[element];

                // put values in the map for leaf elements
                if (element.isLeaf())
                {
                    fvGeometry.bindElement(element);

                    // store current element solution
                    adaptedValues.u = elementSolution(element, sol_, *gridGeometry_);

                    // compute mass in the scvs
                    for (const auto& scv : scvs(fvGeometry))
                    {
                        VolumeVariables volVars;
                        volVars.update(adaptedValues.u, *problem_, element, scv);

                        const auto poreVolume = Extrusion::volume(fvGeometry, scv)*volVars.porosity();
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
                    adaptedValues.u = elementSolution(element, sol_, *gridGeometry_);
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
    void reconstruct(const Grid& grid) override
    {
        gridGeometry_->update(grid.leafGridView());
        reconstruct_();
    }

  private:

    void reconstruct_()
    {
        // resize stuff (grid might have changed)
        adaptionMap_.resize();
        sol_.resize(gridGeometry_->numDofs());

        // For the box method the mass is lumped to the vertices. For each vertex we store
        // the associated mass (component 0) and the corresponding mass coefficient
        // (component 1) in a single block, so that the partial contributions at process
        // boundaries can be summed up in one communication below.
        using VertexMass = Dune::FieldVector<Scalar, 2>;
        std::vector<VertexMass> vertexMass;

        if(isBox)
            vertexMass.assign(gridGeometry_->numDofs(), VertexMass(0.0));

        // iterate over the interior leaf elements and reconstruct the solution;
        // data on ghost entities is obtained by communication further below
        auto fvGeometry = localView(*gridGeometry_);
        for (const auto& element : elements(gridGeometry_->gridView(), Dune::Partitions::interior))
        {
            if (!element.isNew())
            {
                const auto& adaptedValues = adaptionMap_[element];
                fvGeometry.bindElement(element);

                // obtain element solution from map (divide by count!)
                auto elemSol = adaptedValues.u;
                if (!isBox)
                    elemSol[0] /= adaptedValues.count;

                const auto elementVolume = volume(element.geometry(), Extrusion{});
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
                        const auto scvVolume = Extrusion::volume(fvGeometry, scv);
                        if (formulation == p0s1)
                        {
                            vertexMass[dofIdxGlobal][0] += scvVolume / elementVolume * adaptedValues.associatedMass[phase1Idx];
                            vertexMass[dofIdxGlobal][1] += scvVolume * volVars.density(phase1Idx) * volVars.porosity();
                        }
                        else if (formulation == p1s0)
                        {
                            vertexMass[dofIdxGlobal][0] += scvVolume / elementVolume * adaptedValues.associatedMass[phase0Idx];
                            vertexMass[dofIdxGlobal][1] += scvVolume * volVars.density(phase0Idx) * volVars.porosity();
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
                while (fatherElement.isNew() && fatherElement.level() > 0)
                    fatherElement = fatherElement.father();

                if (!isBox)
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
                            massCoeffSon = Extrusion::volume(fvGeometry, scv) * volVars.density(phase1Idx) * volVars.porosity();
                        else if (formulation == p1s0)
                            massCoeffSon = Extrusion::volume(fvGeometry, scv) * volVars.density(phase0Idx) * volVars.porosity();
                        sol_[scv.dofIndex()][saturationIdx] =
                            ( Extrusion::volume(fvGeometry, scv)/volume(fatherElement.geometry(), Extrusion{})*massFather )/massCoeffSon;
                    }
                }
                else
                {
                    auto& adaptedValuesFather = adaptionMap_[fatherElement];

                    fvGeometry.bindElement(element);

                    // interpolate solution in the father to the vertices of the new son
                    ElementSolution elemSolSon(element, sol_, *gridGeometry_);
                    const auto fatherGeometry = fatherElement.geometry();
                    for (const auto& scv : scvs(fvGeometry))
                        elemSolSon[scv.localDofIndex()] = evalSolution(
                            fatherElement, fatherGeometry, adaptedValuesFather.u, scv.dofPosition()
                        );

                    // compute mass & mass coefficients for the scvs (saturations are recalculated at the end)
                    const auto fatherElementVolume = volume(fatherGeometry, Extrusion{});
                    for (const auto& scv : scvs(fvGeometry))
                    {
                        VolumeVariables volVars;
                        volVars.update(elemSolSon, *problem_, element, scv);

                        const auto dofIdxGlobal = scv.dofIndex();
                        const auto scvVolume = Extrusion::volume(fvGeometry, scv);
                        if (formulation == p0s1)
                        {
                            vertexMass[dofIdxGlobal][0] += scvVolume / fatherElementVolume * adaptedValuesFather.associatedMass[phase1Idx];
                            vertexMass[dofIdxGlobal][1] += scvVolume * volVars.density(phase1Idx) * volVars.porosity();
                        }
                        else if (formulation == p1s0)
                        {
                            vertexMass[dofIdxGlobal][0] += scvVolume / fatherElementVolume * adaptedValuesFather.associatedMass[phase0Idx];
                            vertexMass[dofIdxGlobal][1] += scvVolume * volVars.density(phase0Idx) * volVars.porosity();
                        }

                        // store constructed (pressure) values of son in the current solution (saturation comes later)
                        sol_[dofIdxGlobal] = elemSolSon[scv.localDofIndex()];
                    }
                }
            }
        }

        if(isBox)
        {
#if HAVE_MPI
            // At process-boundary vertices each process only accumulated the mass of its
            // interior elements. Sum up the partial contributions (associated mass and mass
            // coefficient together) to obtain the total mass lumped to the vertex before dividing.
            if constexpr (Dune::Capabilities::canCommunicate<Grid, dim>::v)
            {
                const auto& gridView = gridGeometry_->gridView();
                if (gridView.comm().size() > 1)
                {
                    using VertexMapper = std::decay_t<decltype(gridGeometry_->vertexMapper())>;
                    VectorCommDataHandleSum<VertexMapper, std::vector<VertexMass>, dim> massHandle(
                        gridGeometry_->vertexMapper(), vertexMass
                    );
                    gridView.communicate(
                        massHandle, Dune::InteriorBorder_InteriorBorder_Interface, Dune::ForwardCommunication
                    );
                }
            }
#endif
            // ghost vertices received no interior contribution (mass coefficient == 0); their
            // values are obtained by communication of the solution below
            for(std::size_t dofIdxGlobal = 0; dofIdxGlobal < gridGeometry_->numDofs(); dofIdxGlobal++)
                if (vertexMass[dofIdxGlobal][1] > 0.0)
                    sol_[dofIdxGlobal][saturationIdx] = vertexMass[dofIdxGlobal][0] / vertexMass[dofIdxGlobal][1];
        }

#if HAVE_MPI
        // The reconstruction above only produces valid data on interior/border entities.
        // Overwrite the values on ghost entities with those of their owner so that fluxes
        // across the process boundary are evaluated with consistent neighbor values.
        static constexpr int dofCodim = isBox ? dim : 0;
        if constexpr (Dune::Capabilities::canCommunicate<Grid, dofCodim>::v)
        {
            const auto& gridView = gridGeometry_->gridView();
            if (gridView.comm().size() > 1)
            {
                using DofMapper = std::decay_t<decltype(gridGeometry_->dofMapper())>;
                VectorCommDataHandleEqual<DofMapper, SolutionVector, dofCodim> dataHandle(
                    gridGeometry_->dofMapper(), sol_
                );
                gridView.communicate(
                    dataHandle, Dune::InteriorBorder_All_Interface, Dune::ForwardCommunication
                );
            }
        }
#endif

        // reset entries in adaptation map
        adaptionMap_.resize( typename PersistentContainer::Value() );
        adaptionMap_.shrinkToFit();
        adaptionMap_.fill( typename PersistentContainer::Value() );
    }

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
