// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/****************************************************************************
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
 * \brief TODO: docme!
 */

#ifndef DUMUX_PYTHON_COMMON_FVPROBLEM_HH
#define DUMUX_PYTHON_COMMON_FVPROBLEM_HH

#include <string>
#include <memory>
#include <tuple>

#include <dune/common/fvector.hh>
#include <dune/common/exceptions.hh>
#include <dune/python/pybind11/pybind11.h>

#include <dumux/common/boundarytypes.hh>
#include <dumux/discretization/method.hh>
#include <dumux/python/common/boundarytypes.hh>
#include <dumux/python/common/fvspatialparams.hh>

namespace Dumux::Python {

/*!
 * \ingroup Common
 * \brief A C++ wrapper for a Python problem
 */
template<class GridGeometry_,  class SpatialParams_, class PrimaryVariables, bool enableInternalDirichletConstraints_>
class FVProblem
{
public:
    using GridGeometry = GridGeometry_;
    using SpatialParams = SpatialParams_;
    using Scalar = typename GridGeometry::GridView::ctype;
    using NumEqVector = Dune::FieldVector<Scalar, PrimaryVariables::dimension>;
    using Element = typename GridGeometry::GridView::template Codim<0>::Entity;
    using FVElementGeometry = typename GridGeometry::LocalView;
    using SubControlVolume = typename GridGeometry::SubControlVolume;
    using SubControlVolumeFace = typename GridGeometry::SubControlVolumeFace;
    using GlobalPosition = typename Element::Geometry::GlobalCoordinate;

    static constexpr bool isBox = GridGeometry::discMethod == DiscretizationMethods::box;
    static constexpr std::size_t numEq = static_cast<std::size_t>(PrimaryVariables::dimension);
    using BoundaryTypes = Dumux::BoundaryTypes<PrimaryVariables::dimension>;

    FVProblem(std::shared_ptr<const GridGeometry> gridGeometry,
              std::shared_ptr<SpatialParams> spatialParams,
              pybind11::object pyProblem)
    : gridGeometry_(gridGeometry)
    , pyProblem_(pyProblem)
    , name_("python_problem")
    , paramGroup_("")
    , spatialParams_(spatialParams)
    {
        if (pybind11::hasattr(pyProblem_, "name"))
            name_ = pyProblem.attr("name")().template cast<std::string>();

        if (pybind11::hasattr(pyProblem_, "paramGroup"))
            paramGroup_ = pyProblem.attr("paramGroup")().template cast<std::string>();
    }

    FVProblem(std::shared_ptr<const GridGeometry> gridGeometry,
              pybind11::object pyProblem)
    : FVProblem(gridGeometry, std::make_shared<SpatialParams>(gridGeometry), pyProblem)
    {}

    const std::string& name() const
    { return name_; }

    const std::string& paramGroup() const
    { return paramGroup_; }

    BoundaryTypes boundaryTypes(const Element &element,
                                const SubControlVolume &scv) const
    {
        if constexpr (!isBox)
            DUNE_THROW(Dune::InvalidStateException, "boundaryTypes(..., scv) called for cell-centered method.");
        else
        {
            if (pybind11::hasattr(pyProblem_, "boundaryTypes"))
                return pyProblem_.attr("boundaryTypes")(element, scv).template cast<BoundaryTypes>();
            else
                return pyProblem_.attr("boundaryTypesAtPos")(scv.dofPosition()).template cast<BoundaryTypes>();
        }
    }

    BoundaryTypes boundaryTypes(const Element &element,
                                const SubControlVolumeFace &scvf) const
    {
        if constexpr (isBox)
            DUNE_THROW(Dune::InvalidStateException, "boundaryTypes(..., scvf) called for box method.");
        else
        {
            if (pybind11::hasattr(pyProblem_, "boundaryTypes"))
                return pyProblem_.attr("boundaryTypes")(element, scvf).template cast<BoundaryTypes>();
            else
                return pyProblem_.attr("boundaryTypesAtPos")(scvf.ipGlobal()).template cast<BoundaryTypes>();
        }
    }

    PrimaryVariables dirichlet(const Element &element,
                               const SubControlVolume &scv) const
    {
        if constexpr (!isBox)
            DUNE_THROW(Dune::InvalidStateException, "dirichlet(scv) called for cell-centered method.");
        else
        {
            if (pybind11::hasattr(pyProblem_, "dirichlet"))
                return pyProblem_.attr("dirichlet")(element, scv).template cast<PrimaryVariables>();
            else
                return pyProblem_.attr("dirichletAtPos")(scv.dofPosition()).template cast<PrimaryVariables>();
        }
    }

    PrimaryVariables dirichlet(const Element &element,
                               const SubControlVolumeFace &scvf) const
    {
        if constexpr (isBox)
            DUNE_THROW(Dune::InvalidStateException, "dirichlet(scvf) called for box method.");
        else
        {
            if (pybind11::hasattr(pyProblem_, "dirichlet"))
                return pyProblem_.attr("dirichlet")(element, scvf).template cast<PrimaryVariables>();
            else
                return pyProblem_.attr("dirichletAtPos")(scvf.ipGlobal()).template cast<PrimaryVariables>();
        }
    }

    template<class ElementVolumeVariables, class ElementFluxVariablesCache>
    NumEqVector neumann(const Element& element,
                        const FVElementGeometry& fvGeometry,
                        const ElementVolumeVariables& elemVolVars,
                        const ElementFluxVariablesCache& elemFluxVarsCache,
                        const SubControlVolumeFace& scvf) const
    {
        if (pybind11::hasattr(pyProblem_, "neumann"))
            return pyProblem_.attr("neumann")(element, fvGeometry, scvf).template cast<NumEqVector>();
        else
            return pyProblem_.attr("neumannAtPos")(scvf.ipGlobal()).template cast<NumEqVector>();
    }

    template<class ElementVolumeVariables>
    NumEqVector source(const Element &element,
                       const FVElementGeometry& fvGeometry,
                       const ElementVolumeVariables& elemVolVars,
                       const SubControlVolume &scv) const
    {
        if (pybind11::hasattr(pyProblem_, "source"))
            return pyProblem_.attr("source")(element, fvGeometry, scv).template cast<NumEqVector>();
        else
            return pyProblem_.attr("sourceAtPos")(scv.dofPosition()).template cast<NumEqVector>();
    }

    NumEqVector sourceAtPos(const GlobalPosition &globalPos) const
    {
        return pyProblem_.attr("sourceAtPos")(globalPos).template cast<NumEqVector>();
    }

    template<class ElementVolumeVariables>
    NumEqVector scvPointSources(const Element& element,
                                const FVElementGeometry& fvGeometry,
                                const ElementVolumeVariables& elemVolVars,
                                const SubControlVolume& scv) const
    {
        if (pybind11::hasattr(pyProblem_, "scvPointSources"))
            return pyProblem_.attr("scvPointSources")(element, fvGeometry, scv).template cast<NumEqVector>();
        else
            return NumEqVector(0.0);
    }

    template<class Entity>
    PrimaryVariables initial(const Entity& entity) const
    {
        return pyProblem_.attr("initial")(entity).template cast<PrimaryVariables>();
    }

    static constexpr bool enableInternalDirichletConstraints()
    { return enableInternalDirichletConstraints_; }

    /*!
     * \brief Add source term derivative to the Jacobian
     * \note Only needed in case of analytic differentiation and solution dependent sources
     */
    template<class MatrixBlock, class VolumeVariables>
    void addSourceDerivatives(MatrixBlock& block,
                              const Element& element,
                              const FVElementGeometry& fvGeometry,
                              const VolumeVariables& volVars,
                              const SubControlVolume& scv) const
    {
        if (pybind11::hasattr(pyProblem_, "addSourceDerivatives"))
            pyProblem_.attr("addSourceDerivatives")(block, element, fvGeometry, scv);
    }

    [[deprecated("extrusionFactorAtPos() should now be defined in the spatial params. This interface will be removed after 3.5.")]]
    Scalar extrusionFactorAtPos(const GlobalPosition &globalPos, double defaultValue = 1.0) const
    { return 1.0; }

    template<class ElementSolution>
    [[deprecated("extrusionFactor() should now be defined in the spatial params. This interface will be removed after 3.5.")]]
    Scalar extrusionFactor(const Element& element,
                           const SubControlVolume& scv,
                           const ElementSolution& elemSol,
                           double defaultValue = 1.0) const
    { return this->extrusionFactorAtPos(scv.center()); }

    [[deprecated("temperatureAtPos() should now be defined in the spatial params. This interface will be removed after 3.5.")]]
    Scalar temperatureAtPos(const GlobalPosition &globalPos, int defaultValue = 1) const
    { return 293.0; }

    [[deprecated("temperature() should now be defined in the spatial params. This interface will be removed after 3.5.")]]
    Scalar temperature(int defaultValue = 1) const
    { return this->temperature(GlobalPosition(0.0)); }

    const GridGeometry& gridGeometry() const
    { return *gridGeometry_; }

    //! Return a reference to the underlying spatial parameters
    const SpatialParams& spatialParams() const
    { return *spatialParams_; }

    //! Return a reference to the underlying spatial parameters
    SpatialParams& spatialParams()
    { return *spatialParams_; }

private:
    std::shared_ptr<const GridGeometry> gridGeometry_;
    pybind11::object pyProblem_;
    std::string name_;
    std::string paramGroup_;
    std::shared_ptr<SpatialParams> spatialParams_;
};

// Python wrapper for the above FVProblem C++ class
template<class Problem, class... options>
void registerFVProblem(pybind11::handle scope, pybind11::class_<Problem, options...> cls)
{
    using pybind11::operator""_a;
    using namespace Dune::Python;

    using GridGeometry = typename Problem::GridGeometry;
    using SpatialParams = typename Problem::SpatialParams;
    cls.def(pybind11::init([](std::shared_ptr<const GridGeometry> gridGeometry,
                              std::shared_ptr<SpatialParams> spatialParams,
                              pybind11::object p){
        return std::make_shared<Problem>(gridGeometry, spatialParams, p);
    }));
    cls.def(pybind11::init([](std::shared_ptr<const GridGeometry> gridGeometry,
                              pybind11::object p){
        return std::make_shared<Problem>(gridGeometry, p);
    }));

    cls.def_property_readonly("name", &Problem::name);
    cls.def_property_readonly("numEq", [](Problem&){ return Problem::numEq; });
    cls.def_property_readonly("gridGeometry", &Problem::gridGeometry);

    using GridView = typename GridGeometry::GridView;
    using Element = typename GridView::template Codim<0>::Entity;
    using Vertex = typename GridView::template Codim<GridView::dimension>::Entity;

    if constexpr (Problem::isBox)
    {
        using SCV = typename Problem::SubControlVolume;
        cls.def("boundaryTypes", pybind11::overload_cast<const Element&, const SCV&>(&Problem::boundaryTypes, pybind11::const_), "element"_a, "scv"_a);
        cls.def("dirichlet", pybind11::overload_cast<const Element&, const SCV&>(&Problem::dirichlet, pybind11::const_), "element"_a, "scv"_a);
    }
    else
    {
        using SCVF = typename Problem::SubControlVolumeFace;
        cls.def("boundaryTypes", pybind11::overload_cast<const Element&, const SCVF&>(&Problem::boundaryTypes, pybind11::const_), "element"_a, "scvf"_a);
        cls.def("dirichlet", pybind11::overload_cast<const Element&, const SCVF&>(&Problem::dirichlet, pybind11::const_), "element"_a, "scvf"_a);
    }

    cls.def("neumann", &Problem::template neumann<decltype(std::ignore), decltype(std::ignore)>);
    cls.def("source", &Problem::template source<decltype(std::ignore)>);
    cls.def("sourceAtPos", &Problem::sourceAtPos);
    cls.def("initial", &Problem::template initial<Element>);
    cls.def("initial", &Problem::template initial<Vertex>);
}

} // end namespace Dumux::Python

#endif
