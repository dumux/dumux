# Style guide

## File system

* Try to order your new header into the existing directory structure
* Headers are named like the classes they contain (usually one class per file, exception: closely tied helper classes / functions)
* Headers are names lower case only (using underscores only if absolutely necessary for readability)
* Folder names are lower case only
* Tests should be called after model and discretization scheme using underscores and lower case only, e.g.

__Example:__
```
test_1p_tpfa
test_2p2c_box_infiltration
```

## C++

### Indent

* The default indent is 4 spaces (no tabs, not 2 spaces)

### Formatting

* Curly brackets that open or close a scope are on a separate line (except for [namespaces](#Namespaces))
* There should be a space between `if`,`else if`,`switch` and the condition
* `if`,`else if`,`switch` can omit brackets if the expression is only one line

```c++
// comment for if block, space between if and (enableGravity)
if (enableGravity)
{
    Scalar b = 1.0;
    return b*gravity_[dimWorld-1];
}

// comment for else-block
else
    return 1.0; // ok, one-liner
```

### Namespaces

* Open curly brackets on the same line
* Do not indent the code inside the namespace
* Comment closing curly brackets uniquely
__Example:__

```c++
namespace Dumux {
namespace Properties {

bool here = true;
    bool nothere = false; // not like this

} // end namespace Properties
} // end namespace Dumux
```

### Includes

* Space between `#include` and path
* C++ standard library includes first, then Dune, then others, then DuMu<sup>x</sup>
* Always use project relative paths
    - exception: the other header is in the same folder and closely related (`#include "volumevariables.hh"`)

```c++
#include <type_traits>
#include <dune/common/fvector.hh>
#include <dumux/common/exceptions.hh>
```

### Documentation

* Prefer long-style doxygen documentation of the form
* at least containing a `\brief`

```c++
/*
 * \brief Base class for all linear solvers
 */
```

* document all parameters of functions and what it returns
* optionally document which error the function might throw

### Classes

TODO

#### Member functions

TODO

#### Private data member
* Use underscore postfix

```c++
...
private:
    int steps_;
...
```

### Functions

TODO

### Property system

* Prefer class templates with regular template arguments over class templates with a `TypeTag` as template argument

TODO
