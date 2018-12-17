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
 * \ingroup MultiDomain
 * \brief Multidomain wrapper for multiple vtk output modules
 */
#ifndef DUMUX_MULTIDOMAIN_VTK_OUTPUT_MODULE_HH
#define DUMUX_MULTIDOMAIN_VTK_OUTPUT_MODULE_HH

#include <tuple>
#include <memory>
#include <utility>

#include <dune/common/hybridutilities.hh>
#include <dune/common/indices.hh>

#include <dumux/common/typetraits/utility.hh>
#include <dumux/io/vtkoutputmodule.hh>

namespace Dumux {

/*!
 * \ingroup MultiDomain
 * \brief A multidomain wrapper for multiple vtk output modules
 * \tparam MDTraits The multidomain traits
 * \tparam Module An output module class template that takes GridVariables and SolutionVector as arugments
 */
template<class MDTraits, template<class GV, class S> class Module = Dumux::VtkOutputModule>
class MultiDomainVtkOutputModule
{
    using MDSolutionVector = typename MDTraits::SolutionVector;
    static constexpr std::size_t numSubDomains = MDTraits::numSubDomains;

    template<std::size_t i>
    using GridVariables = typename MDTraits::template SubDomain<i>::GridVariables;

    using MDGridVars = typename MDTraits::template TupleOfSharedPtrConst<GridVariables>;

    template<std::size_t i>
    using SolutionVector = typename MDTraits::template SubDomain<i>::SolutionVector;

    template<std::size_t i>
    using VtkOutputModule = Module<GridVariables<i>, SolutionVector<i>>;

    using VtkOutputModuleTuple = typename MDTraits::template TupleOfSharedPtr<VtkOutputModule>;

public:
    //! export base types of the stored type
    template<std::size_t i>
    using Type = VtkOutputModule<i>;

    //! export pointer types the stored type
    template<std::size_t i>
    using PtrType = std::shared_ptr<Type<i>>;

    /*!
     * \brief The default constructor
     */
    MultiDomainVtkOutputModule() = default;

    /*!
     * \brief Contruct the vtk output modules
     * \param gridVars a tuple of grid variables
     * \param sol the multidomain solution vector
     * \param name the base name for the vtk output
     */
    MultiDomainVtkOutputModule(MDGridVars&& gridVars, const MDSolutionVector& sol,
                               const std::array<std::string, numSubDomains>& name)
    {
        using namespace Dune::Hybrid;
        forEach(std::make_index_sequence<numSubDomains>{}, [&](auto&& id)
        {
            constexpr auto i = std::decay_t<decltype(id)>::value;
            elementAt(vtkOutputModule_, id) = std::make_shared<Type<i>>(*std::get<i>(gridVars), sol[id], name[id]);
        });
    }

    //! initialized all vtkoutput modules with the models default output fields
    void initDefaultOutputFields()
    {
        using namespace Dune::Hybrid;
        forEach(std::make_index_sequence<numSubDomains>{}, [&](auto&& id)
        {
            constexpr auto i = std::decay_t<decltype(id)>::value;
            MDTraits::template SubDomain<i>::IOFields::initOutputModule(*elementAt(vtkOutputModule_, id));
        });
    }

    //! Write the data for this timestep to file for all output moduless
    void write(double t, Dune::VTK::OutputType type = Dune::VTK::ascii)
    {
        using namespace Dune::Hybrid;
        forEach(std::make_index_sequence<numSubDomains>{}, [&](auto&& id)
        {
            elementAt(vtkOutputModule_, id)->write(t, type);
        });
    }

    //! return the output module for domain with index i
    template<std::size_t i>
    const Type<i>& operator[] (Dune::index_constant<i> id) const
    { return *Dune::Hybrid::elementAt(vtkOutputModule_, id); }

    //! return the output module for domain with index i
    template<std::size_t i>
    Type<i>& operator[] (Dune::index_constant<i> id)
    { return *Dune::Hybrid::elementAt(vtkOutputModule_, id); }

    //! return the vtkoutput module for domain with index i
    template<std::size_t i>
    PtrType<i> get(Dune::index_constant<i> id = Dune::index_constant<i>{})
    { return Dune::Hybrid::elementAt(vtkOutputModule_, id); }

    //! set the pointer for sub domain i
    template<std::size_t i>
    void set(PtrType<i> p, Dune::index_constant<i> id = Dune::index_constant<i>{})
    { Dune::Hybrid::elementAt(vtkOutputModule_, id) = p; }

private:

    //! a tuple of pointes to all vtk output modules
    typename MDTraits::template Tuple<PtrType> vtkOutputModule_;
};

} // end namespace Dumux

#endif
