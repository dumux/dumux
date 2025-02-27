# Gmsh with ALUGrid

If you are using Gmsh to create more complex grids and are using "Physical groups" to mark different boundary segments, you can normally access these inside your definition of the `neumann` function (solution dependent von Neumann boundaries) in your problem file by using something like the following:

```c++
NumEqVector neumann(const Element& element,
                    const FVElementGeometry& fvGeometry,
                    const ElementVolumeVariables& elemVolVars,
                    const SubControlVolumeFace& scvf) const
{
    ...
    const auto boundaryId = gridData_->getBoundaryDomainMarker(scvf.boundaryFlag());
    ...
}
```
where `gridData_` should be a `GridData` instance and a member of your problem class. (Add the grid data as argument in your problem constructor, it can be obtained via the grid manager: `gridManager.getGridData()`.)

Additionally, you have to specify to use boundary segments and domain markers in your `params.input` file:

```ini
[Grid]
File = mygrid.msh
DomainMarkers = 1
BoundarySegments = 1
```

Unfortunately `ALUGrid` uses boundary flags differently, therefore you have to adapt your `FVGridGeometry` to use the class `BoundarySegmentIndexFlag` instead of the normal `BoundaryFlag` class.
This can be done in your property settings, where you need to add the following:

```c++
#include <dumux/common/boundaryflag.hh>

namespace Properties {

// other property settings
...

// custom grid geometry traits that enable you to use another boundary flag class
template<class GridView>
struct MyGridGeometryTraits : public BoxDefaultGridGeometryTraits<GridView>
{
    struct MyScvfTraits : public BoxDefaultScvfGeometryTraits<GridView>
    {
        // use BoundarySegmentIndexFlag
        using BoundaryFlag = BoundarySegmentIndexFlag;
    };

    using SubControlVolumeFace = BoxSubControlVolumeFace<GridView, MyScvfTraits>;
};

// custom FVGridGeometry that uses your specified boundary flag implementation
template<class TypeTag>
struct FVGridGeometry<TypeTag, TTag::MyProblem>
{
private:
    static constexpr bool enableCache = getPropValue<TypeTag, Properties::EnableFVGridGeometryCache>();
    using GridView = GetPropType<TypeTag, Properties::GridView>;
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
public:
    using type = BoxFVGridGeometry<Scalar, GridView, enableCache,
                                   MyGridGeometryTraits<GridView>>;
};
...
} // end namespace properties
```
This was how to do it when using the Box discretization, for CCTpfa it works similarly:

```c++
namespace Properties {

// other property settings
...

template<class GridView>
struct MyGridGeometryTraits : public CCTpfaDefaultGridGeometryTraits<GridView>
{
    struct MyScvfTraits : public CCTpfaDefaultScvfGeometryTraits<GridView>
    {
        // use BoundarySegmentIndexFlag
        using BoundaryFlag = BoundarySegmentIndexFlag;
    };

    using SubControlVolumeFace = CCTpfaSubControlVolumeFace<GridView, MyScvfTraits>;
};

template<class TypeTag>
struct FVGridGeometry<TypeTag, TTag::MyProblem>
{
private:
    static constexpr bool enableCache = getPropValue<TypeTag, Properties::EnableFVGridGeometryCache>();
    using GridView = GetPropType<TypeTag, Properties::GridView>;
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
public:
    using type = CCTpfaFVGridGeometry<GridView, enableCache,
                                   MyGridGeometryTraits<GridView>>;
};
...
} // end namespace properties
```

Then you can access the boundary marker/physical group in the same way as shown above