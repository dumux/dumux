#ifndef _DUMUX_POROUSMEDIUMFLOW_PY_PROBLEM_HH_
#define _DUMUX_POROUSMEDIUMFLOW_PY_PROBLEM_HH_

#include <dune/python/pybind11/pybind11.h>
#include <dumux/porousmediumflow/py/onepproblem.hh>


namespace Dumux
{
    template< class Impl, class... options >
    void registerOnePProblem ( pybind11::class_< Impl, options... > cls )
    {
      // cls.def( "temperature", []( Impl &self ) { return self.temperature(); } );
    }

    template< class TypeTag, class... options >
    void registerOnePProblem ( pybind11::module scope,
                               pybind11::class_< OnePTestProblem< TypeTag >, options... > cls )
    {
      registerOnePProblem( cls );

      using GridView = typename GetPropType<TypeTag, Properties::Grid>::LeafGridView;
      cls.def( pybind11::init( [](const GridView& gridView)
      {
        using FVGridGeometry = GetPropType<TypeTag, Properties::FVGridGeometry>;
        auto fvGridGeometry = std::make_shared<FVGridGeometry>( gridView );
        return OnePTestProblem< TypeTag >( fvGridGeometry );
      }
      ) );
    }
}

#endif
