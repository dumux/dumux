// // -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup PoreNetworkModels
 * \copydoc Dumux::PoreNetwork::VtkOutputModule
 */

#ifndef DUMUX_PNM_VTK_OUTPUT_MODULE_HH
#define DUMUX_PNM_VTK_OUTPUT_MODULE_HH

#include <dumux/io/vtkoutputmodule.hh>
#include "velocityoutput.hh"

namespace Dumux::PoreNetwork {

/*!
 * \ingroup PoreNetworkModels
 * \brief Adds vtk output fields specific to pore-network models.
 */
template<class GridVariables, class FluxVariables, class SolutionVector>
class VtkOutputModule : public Dumux::VtkOutputModule<GridVariables, SolutionVector>
{
    using ParentType = Dumux::VtkOutputModule<GridVariables, SolutionVector>;
    using Scalar = typename GridVariables::Scalar;
    using FluxVarsCache = typename GridVariables::GridFluxVariablesCache::FluxVariablesCache;

    struct ThroatFluxDataInfo { std::function<Scalar(const FluxVariables&, const FluxVarsCache&)> get; std::string name; };

public:

    //! The constructor
    VtkOutputModule(const GridVariables& gridVariables,
                    const SolutionVector& sol,
                    const std::string& name,
                    const std::string& paramGroup = "",
                    Dune::VTK::DataMode dm = Dune::VTK::conforming,
                    bool verbose = true)
    : ParentType(gridVariables, sol, name, paramGroup, dm, verbose)
    {
        if constexpr (GridVariables::VolumeVariables::numFluidPhases() >= 1)
        {
            // enable velocity output per default
            using VelocityOutput = VelocityOutput<GridVariables, FluxVariables>;
            this->addVelocityOutput(std::make_shared<VelocityOutput>(gridVariables));
        }
    }

    //! Output a scalar flux variable related to pore throats. This is basically a wrapper for the ParentType's addField method.
    //! \param f A function taking a Problem, FluxVariables and FluxVarsCache object and returning the desired scalar
    //! \param name The name of the vtk field
    void addFluxVariable(std::function<Scalar(const FluxVariables&, const FluxVarsCache&)>&& f, const std::string& name)
    {
        throatFluxDataInfo_.push_back(ThroatFluxDataInfo{f, name});
        const auto numElems = problem().gridGeometry().gridView().size(0);
        throatFluxData_.push_back(std::vector<Scalar>(numElems));
        ParentType::addField(throatFluxData_.back(), throatFluxDataInfo_.back().name, Vtk::FieldType::element);
    }

    //! Gather and process all required data and write them to a vtk file.
    void write(double time, Dune::VTK::OutputType type = Dune::VTK::ascii)
    {
        const auto gridView = problem().gridGeometry().gridView();
        const auto numElems = gridView.size(0);

        // resize all fields
        for (auto& data : throatFluxData_)
            data.resize(numElems);

        auto fvElementGeometry = localView(problem().gridGeometry());
        auto elemVolVars = localView(this->gridVariables().curGridVolVars());
        auto elemFluxVarsCache = localView(this->gridVariables().gridFluxVarsCache());
        // iterate over all elements
        for (const auto& element : elements(gridView, Dune::Partitions::interior))
        {
            const auto eIdx = problem().gridGeometry().elementMapper().index(element);

            // make sure FVElementGeometry & vol vars are bound to the element
            fvElementGeometry.bindElement(element);
            elemVolVars.bind(element, fvElementGeometry, this->sol());
            elemFluxVarsCache.bind(element, fvElementGeometry, elemVolVars);

            // treat the throat flux related data
            std::size_t dataIdx = 0;
            for (auto&& scvf : scvfs(fvElementGeometry))
            {
                if (!scvf.boundary())
                {
                    FluxVariables fluxVars;
                    fluxVars.init(problem(), element, fvElementGeometry, elemVolVars, scvf, elemFluxVarsCache);

                    dataIdx = 0;
                    for(auto& data : throatFluxData_)
                        data[eIdx] = throatFluxDataInfo_[dataIdx++].get(fluxVars, elemFluxVarsCache[scvf]);
                }
            }
        }

        // call the ParentType's write method to write out all data
        ParentType::write(time, type);

        // empty the data containers in order to save some memory
        auto clearAndShrink = [] (auto& data)
        {
            data.clear();
            data.shrink_to_fit();
        };

        for (auto& data : throatFluxData_)
            clearAndShrink(data);
    }

    //! Return a reference to the problem
    const auto& problem() const { return ParentType::problem(); }

private:
    std::vector<ThroatFluxDataInfo> throatFluxDataInfo_;
    std::list<std::vector<Scalar>> throatFluxData_;
};

} //namespace Dumux::PoreNetwork

#endif
