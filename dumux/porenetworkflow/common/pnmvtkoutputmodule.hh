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
 *
 * \copydoc Dumux::PNMVtkOutputModule
 */

#ifndef DUMUX_PNM_VTK_OUTPUT_MODULE_HH
#define DUMUX_PNM_VTK_OUTPUT_MODULE_HH

#include <utility>

#include <dumux/common/properties.hh> // TODO remove soon
#include <dumux/io/vtkoutputmodule.hh>
#include "velocityoutput.hh"

namespace Dumux
{

/*!
 * \ingroup PoreNetworkFlow
 * \brief Adds vtk output fields specific to pore-network models. TODO remove deprecated soon
 */
template<bool deprecated, class GridVariables, class FluxVariables, class SolutionVector>
class PNMVtkOutputModuleImpl : public VtkOutputModule<GridVariables, SolutionVector>
{
    using ParentType = VtkOutputModule<GridVariables, SolutionVector>;
    using Scalar = typename GridVariables::Scalar;
    using FluxVarsCache = typename GridVariables::GridFluxVariablesCache::FluxVariablesCache;

    struct ThroatFluxDataInfo { std::function<Scalar(const FluxVariables&, const FluxVarsCache&)> get; std::string name; };

public:

    //! The constructor
    PNMVtkOutputModuleImpl(const GridVariables& gridVariables,
                           const SolutionVector& sol,
                           const std::string& name,
                           const std::string& paramGroup = "",
                           Dune::VTK::DataMode dm = Dune::VTK::conforming,
                           bool verbose = true)
    : ParentType(gridVariables, sol, name, paramGroup, dm, verbose)
    {
        // enable velocity output per default
        using VelocityOutput = PNMVelocityOutput<GridVariables, FluxVariables>;
        this->addVelocityOutput(std::make_shared<VelocityOutput>(gridVariables));

        // TODO remove soon
        if constexpr(deprecated)
            deprecationMessage_();
    }

    //! Output a scalar flux variable related to pore throats. This is basically a wrapper for the ParentType's addField method.
    //! \param f A function taking a Problem, FluxVariables and FluxVarsCache object and returning the desired scalar
    //! \param name The name of the vtk field
    void addFluxVariable(std::function<Scalar(const FluxVariables&, const FluxVarsCache&)>&& f, const std::string& name)
    {
        throatFluxDataInfo_.push_back(ThroatFluxDataInfo{f, name});
        const auto numElems = problem().gridGeometry().gridView().size(0);
        throatFluxData_.push_back(std::vector<Scalar>(numElems));
        ParentType::addField(throatFluxData_.back(), throatFluxDataInfo_.back().name, ParentType::FieldType::element);
    }

    //! Gather and process all required data and write them to a vtk file.
    void write(double time, Dune::VTK::OutputType type = Dune::VTK::ascii)
    {
        const auto gridView = problem().gridGeometry().gridView();
        const auto numElems = gridView.size(0);

        // resize all fields
        for (auto& data : throatFluxData_)
            data.resize(numElems);

        // iterate over all elements
        for (const auto& element : elements(gridView, Dune::Partitions::interior))
        {
            // make sure FVElementGeometry & vol vars are bound to the element
            auto fvElementGeometry = localView(problem().gridGeometry());
            fvElementGeometry.bindElement(element);

            const auto eIdx = problem().gridGeometry().elementMapper().index(element);

            auto elemVolVars = localView(this->gridVariables().curGridVolVars());
            auto elemFluxVarsCache = localView(this->gridVariables().gridFluxVarsCache());

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

    [[deprecated("Use PNMVtkOutputModule<GridVariables, FluxVariables, SolutionVector> instead. Will be removed soon.")]]
    void deprecationMessage_() const {}

    std::vector<ThroatFluxDataInfo> throatFluxDataInfo_;
    std::list<std::vector<Scalar>> throatFluxData_;
};

// TODO remove Detail and alias soon
namespace Detail {

template<class T>
struct TypeTagHelper
{
    using GridVariables = GetPropType<T, Properties::GridVariables>;
    using FluxVariables = GetPropType<T, Properties::FluxVariables>;
    using SolutionVector = GetPropType<T, Properties::SolutionVector>;
};

template<int i>
struct ImplHelper
{
    template<class ...T>
    using type = PNMVtkOutputModuleImpl<false, T...>;
};

template<>
struct  ImplHelper<1>
{
    template<class T>
    using type = PNMVtkOutputModuleImpl<true, typename Detail::TypeTagHelper<T>::GridVariables, typename Detail::TypeTagHelper<T>::FluxVariables, typename Detail::TypeTagHelper<T>::SolutionVector>;
};

}

template<class ...T>
using PNMVtkOutputModule = typename Detail::ImplHelper<sizeof...(T)>::template type<T...>;

} //namespace Dumux


#endif // DUMUX_PNM_VTK_OUTPUT_MODULE_HH
