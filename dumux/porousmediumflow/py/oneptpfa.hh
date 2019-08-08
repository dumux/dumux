#ifndef _DUMUX_POROUSMEDIUMFLOW_PY_ONEPTPFA_HH_
#define _DUMUX_POROUSMEDIUMFLOW_PY_ONEPTPFA_HH_

#include <dune/python/pybind11/pybind11.h>
#include <dumux/porousmediumflow/properties.hh>

#include <dumux/discretization/cctpfa.hh>
#include <dumux/porousmediumflow/1p/model.hh>
#include <dumux/porousmediumflow/1p/incompressiblelocalresidual.hh>

#include <dune/grid/yaspgrid.hh>
#include <dumux/porousmediumflow/py/onepproblem.hh>


namespace Dumux {
namespace Properties {
  namespace TTag {
    struct OnePTpfa {
      using InheritsFrom = std::tuple<OneP, CCTpfaModel>;
    };
  }

// Set the problem type
template<class TypeTag>
struct Problem<TypeTag, TTag::OnePTpfa> { using type = OnePTestProblem<TypeTag>; };

// set the spatial params
template<class TypeTag>
struct SpatialParams<TypeTag, TTag::OnePTpfa>
{
    using FVGridGeometry = GetPropType<TypeTag, Properties::FVGridGeometry>;
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using type = OnePTestSpatialParams<FVGridGeometry, Scalar>;
};

// use the incompressible local residual (provides analytic jacobian)
template<class TypeTag>
struct LocalResidual<TypeTag, TTag::OnePTpfa> { using type = OnePIncompressibleLocalResidual<TypeTag>; };

// the fluid system
template<class TypeTag>
struct FluidSystem<TypeTag, TTag::OnePTpfa>
{
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using type = FluidSystems::OnePLiquid<Scalar, Components::SimpleH2O<Scalar> >;
};

// Enable caching
template<class TypeTag>
struct EnableGridVolumeVariablesCache<TypeTag, TTag::OnePTpfa> { static constexpr bool value = false; };
template<class TypeTag>
struct EnableGridFluxVariablesCache<TypeTag, TTag::OnePTpfa> { static constexpr bool value = false; };
template<class TypeTag>
struct EnableFVGridGeometryCache<TypeTag, TTag::OnePTpfa> { static constexpr bool value = false; };


//! PythonTypeTag
namespace TTag {
  template< class GridView >
  struct PythonTypeTag {
    using InheritsFrom = std::tuple<OnePTpfa>;
    using Grid = typename GridView::Grid;
  };
}

template<class TypeTag, class GridView>
struct Grid<TypeTag, TTag::PythonTypeTag<GridView>> { using type = typename TTag::template PythonTypeTag<GridView>::Grid; };

} // end namespace Properties

template< class GridView >
using TempTypeTag = typename Properties::TTag::PythonTypeTag<GridView>;

    template< class GridView, class... options >
    void registerOnePTpfa ( pybind11::module scope,
                            pybind11::class_< TempTypeTag<GridView>, options... > cls )
    {
      cls.def( pybind11::init( []()
      {
        return TempTypeTag<GridView>();
      } ) );
    }
}

#endif
