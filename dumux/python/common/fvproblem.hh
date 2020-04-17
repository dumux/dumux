#ifndef DUMUX_PYTHON_COMMON_FVPROBLEM_HH
#define DUMUX_PYTHON_COMMON_FVPROBLEM_HH

#include <string>
#include <memory>
#include <tuple>

#include <iostream> // for debug stuff below
#include <dune/python/pybind11/iostream.h> // for debug stuff below

#include <dune/common/fvector.hh>
#include <dune/common/exceptions.hh>
#include <dune/python/pybind11/pybind11.h>

#include <dumux/common/boundarytypes.hh>
#include <dumux/python/common/boundarytypes.hh>

namespace Dumux::Python {

/*!
 * \ingroup Common
 * \brief A C++ wrapper for a Python problem
 */
template<class GridGeometry_, class PrimaryVariables>
class FVProblem
{
public:
    using GridGeometry = GridGeometry_;
    using Scalar = typename PrimaryVariables::value_type;
    using NumEqVector = Dune::FieldVector<Scalar, PrimaryVariables::dimension>;
    using Element = typename GridGeometry::GridView::template Codim<0>::Entity;
    using FVElementGeometry = typename GridGeometry::LocalView;
    using SubControlVolume = typename GridGeometry::SubControlVolume;
    using SubControlVolumeFace = typename GridGeometry::SubControlVolumeFace;

    static constexpr bool isBox = GridGeometry::discMethod == DiscretizationMethod::box;
    static constexpr std::size_t numEq = static_cast<std::size_t>(PrimaryVariables::dimension);
    using BoundaryTypes = Dumux::BoundaryTypes<PrimaryVariables::dimension>;

    FVProblem(std::shared_ptr<const GridGeometry> gridGeometry, pybind11::object pyProblem)
    : gridGeometry_(gridGeometry), pyProblem_(pyProblem)
    {}

    std::string name() const
    {
        return pyProblem_.attr("name").template cast<std::string>();
    }

    BoundaryTypes boundaryTypes(const Element &element,
                                const SubControlVolume &scv) const
    {
        if constexpr (!isBox)
            DUNE_THROW(Dune::InvalidStateException, "boundaryTypes(..., scv) called for cell-centered method.");
        else
            return pyProblem_.attr("boundaryTypes")(element, scv).template cast<BoundaryTypes>();
    }

    BoundaryTypes boundaryTypes(const Element &element,
                                const SubControlVolumeFace &scvf) const
    {
        if constexpr (isBox)
            DUNE_THROW(Dune::InvalidStateException, "boundaryTypes(..., scvf) called for box method.");
        else
            return pyProblem_.attr("boundaryTypes")(element, scvf).template cast<BoundaryTypes>();
    }

    PrimaryVariables dirichlet(const Element &element,
                               const SubControlVolume &scv) const
    {
        if constexpr (!isBox)
            DUNE_THROW(Dune::InvalidStateException, "dirichlet(scvf) called for cell-centered method.");
        else
            return pyProblem_.attr("dirichlet")(element, scv).template cast<PrimaryVariables>();
    }

    PrimaryVariables dirichlet(const Element &element,
                               const SubControlVolumeFace &scvf) const
    {
        if constexpr (isBox)
            DUNE_THROW(Dune::InvalidStateException, "dirichlet(scvf) called for box method.");
        else
            return pyProblem_.attr("dirichlet")(element, scvf).template cast<PrimaryVariables>();
    }

    template<class ElementVolumeVariables, class ElementFluxVariablesCache>
    NumEqVector neumann(const Element& element,
                        const FVElementGeometry& fvGeometry,
                        const ElementVolumeVariables& elemVolVars,
                        const ElementFluxVariablesCache& elemFluxVarsCache,
                        const SubControlVolumeFace& scvf) const
    {
        return pyProblem_.attr("neumann")(element, fvGeometry, scvf).template cast<NumEqVector>();
    }

    template<class ElementVolumeVariables>
    NumEqVector source(const Element &element,
                       const FVElementGeometry& fvGeometry,
                       const ElementVolumeVariables& elemVolVars,
                       const SubControlVolume &scv) const
    {
        return pyProblem_.attr("source")(element, fvGeometry, scv).template cast<NumEqVector>();
    }

    template<class Entity>
    PrimaryVariables initial(const Entity& entity) const
    {
        return pyProblem_.attr("intial")(entity).template cast<PrimaryVariables>();
    }

    template<class ElementSolution>
    Scalar extrusionFactor(const Element& element,
                           const SubControlVolume& scv,
                           const ElementSolution& elemSol) const
    {
        return pyProblem_.attr("extrusionFactor")(element, scv).template cast<Scalar>();
    }

    const GridGeometry& gridGeometry() const
    { return *gridGeometry_; }

private:
    std::shared_ptr<const GridGeometry> gridGeometry_;
    pybind11::object pyProblem_;
};

// Python wrapper for the above FVProblem C++ class
template<class Problem, class... options>
void registerFVProblem(pybind11::handle scope, pybind11::class_<Problem, options...> cls)
{
    using pybind11::operator""_a;
    using namespace Dune::Python;

    using GridGeometry = typename Problem::GridGeometry;
    cls.def(pybind11::init([](std::shared_ptr<const GridGeometry> gridGeometry, pybind11::object p){
        return std::make_shared<Problem>(gridGeometry, p);
    }));

    cls.def_property_readonly("name", &Problem::name);
    cls.def_property_readonly("numEq", [](Problem&){ return Problem::numEq; });

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
    cls.def("initial", &Problem::template initial<Element>);
    cls.def("initial", &Problem::template initial<Vertex>);
    cls.def("extrusionFactor", &Problem::template extrusionFactor<decltype(std::ignore)>);
    cls.def("gridGeometry", &Problem::gridGeometry);
}


/////////////////////////////////////////////////////////
/// Some debugging/testing  stuff
///////////////////////////////////////////////////////////
template<class Problem>
void printBoundaryStuff(const Problem& problem)
{
    const auto& gg = problem.gridGeometry();
    for (const auto& element : elements(gg.gridView()))
    {
        auto fvGeometry = localView(gg);
        fvGeometry.bindElement(element);

        for (const auto& scv : scvs(fvGeometry))
        {
            const auto boundaryTypes = problem.boundaryTypes(element, scv);
            std::cout << "-- scv at " << scv.dofPosition() << ": "
                      << " isNeumann: " << std::boolalpha << boundaryTypes.hasNeumann()
                      << " isDirichlet: " << boundaryTypes.hasDirichlet() << std::endl;

            if (boundaryTypes.hasDirichlet())
                std::cout << "  -- Dirichlet values: " << problem.dirichlet(element, scv) << std::endl;
        }
    }
}

template<class P>
class PrintBoundaryStuff
{
public:
    using Problem = P;

    PrintBoundaryStuff(std::shared_ptr<const Problem> problem)
    : problem_(problem)
    {}

    void print()
    {
        printBoundaryStuff(*problem_);
    }

private:
    std::shared_ptr<const Problem> problem_;
};

template<class T, class... options>
void registerPrintBoundaryStuff(pybind11::handle scope, pybind11::class_<T, options...> cls)
{
    using Problem = typename T::Problem;
    cls.def(pybind11::init([](std::shared_ptr<const Problem> problem){
        return std::make_unique<T>(problem);
    }));
    cls.def("print", [](T& self){
        pybind11::scoped_ostream_redirect stream(std::cout, pybind11::module::import("sys").attr("stdout"));
        self.print();
    });
}

} // end namespace Dumux::Python

#endif
