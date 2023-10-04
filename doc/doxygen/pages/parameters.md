# Runtime Parameters

Runtime simulation parameters can be parsed to the program via a parameter file, via the command line,
or can be initialized programmatically. We discuss all three approaches. Moreover, most parameters have default values as explained below.
Runtime parameters are a configuration mechanism at program runtime (avoiding recompilation).

## Initializing the parameter tree

The parameter tree is initialized by Dumux::Parameters::init.
This constructs a parameter tree [singleton](https://en.wikipedia.org/wiki/Singleton_pattern)
from which parameters can be retrieved via accessor functions.
The parameter tree can stores key-value pairs of type string.
Values can also be other parameter trees (subtrees, groups).
A simple program only initializing the parameter tree looks like this

```cpp
#include <dumux/common/parameters.hh>
int main(int argc, char** argv)
{
    Dumux::Parameters::init(argc, argv);
    return 0;
}
```
We can also specify a parameter file (default: `params.input`) to be parsed
(the expected INI file format is described below):

```cpp
#include <dumux/common/parameters.hh>
int main(int argc, char** argv)
{
    Dumux::Parameters::init(argc, argv, "params.input");
    return 0;
}
```

Default parameters can explicitly be set upon initialization
```cpp
#include <dune/common/parametertree.hh>
#include <dumux/common/parameters.hh>
int main(int argc, char** argv)
{
    Dumux::Parameters::init(argc, argv, [](Dune::ParameterTree& p){
        p["key"] = "value";
        p["group.key"] = "value2";
        ...
    });
    return 0;
}
```

The following variant omits reading command-line arguments and only
parsing the specified input parameter file `params.input`:

```cpp
#include <dumux/common/parameters.hh>
int main(int argc, char** argv)
{
    Dumux::Parameters::init("params.input");
    return 0;
}
```

## Reading Runtime Parameters

Runtime parameters can be read from the parameter tree with the functions
Dumux::getParam (converts string to requested type)

```cpp
paramname_ = getParam<TYPE>("GROUPNAME.PARAMNAME", default);
```

where the specified default value is expected to be convertible to the type `TYPE`.
Some examples are
```cpp
  bool enableGravity = getParam<bool>("Problem.EnableGravity", true);
  auto upperRight = getParam<Dune::FieldVector<double, 3>>("FreeFlow.Grid.UpperRight");
```

Not specifying a default parameter causes the function to
throws a `Dumux::ParameterException` if parameter doesn't exist in the parameter tree.

The function Dumux::getParamFromGroup traverses the parameter tree
```cpp
paramname_ = getParamFromGroup<TYPE>("GROUPNAME", "PARAMNAME", default);
```
for example
```cpp
bool enableGravity = getParamFromGroup<bool>("FreeFlow", "Problem.Gravity");
```
first looks for the key `FreeFlow.Problem.Gravity` and then looks for the key `Problem.Gravity`.
This function is useful when configuring multiple simulation components or multi-domain problem
via the single parameter tree.

Reading parameters from the tree can be a very slow operation. Therefore, we recommend to
read parameters in constructors of high-level classes
and in particular never read parameters in functions called for all elements.

## Checking existence of runtime parameters

The existence of a parameter in the parameter tree can be queried with the function Dumux::hasParam.
```cpp
if (hasParam("GROUPNAME.PARAMNAME"))
    // do something with parameter
```
There exists also a variant using hierarchical group lookup as described above
for Dumux::getParamFromGroup using the function Dumux::hasParamInGroup:

```cpp
if (hasParamInGroup("GROUPNAME","PARAMNAME"))
    // do something with parameter
```

## Parameter tree logs

There is a bookkeeping mechanism keeping track of used and unused parameters throughout the lifetime
of the program. It is useful to print a parameter report using Dumux::Parameters::print at the end
of the main file. This reports unused parameters and is great, for instance, for detecting typos in
configuration files.

```cpp
#include <dumux/common/parameters.hh>
int main(int argc, char** argv)
{
    Dumux::Parameters::init(argc, argv);
    ...
    Dumux::Parameters::print(); // print report
    return 0;
}
```

## Parameter input file

The parameter files are expected to use the
Dune INI syntax (consists of `[Group]` and `Key = Value` pairs).
An example is given below:

```ini
[Grid]
LowerLeft = 0 0
UpperRight = 60 40
Cells = 24 16

[Problem]
Name = test

[FreeFlow.Problem]
Name = test_ff
```

## Command-line arguments

All parameter that can be specified via the parameter file can also
be overwritten via the command line using the following exemplary form:
```sh
./executable params.input -Key Value -Key2 Value2 -Key3 "a b c"
```

The following could be used in combination with the parameter file above:
```sh
./executable -Grid.Refinement 2
./executable -Grid.Refinement 2 -Grid.Cells "10 10"
./executable -FreeFlow.Problem.Name test_ff_customrun
```

As first argument a parameter file can be specified:
```sh
./executable params_alt.input -Grid.Refinement 2
```

## Parameter precedence

Parameters are parsed into the parameter tree in the following precedence:

1. Command-line arguments overwrite
2. Input file arguments overwrite
3. User-default arguments overwrite
4. Global defaults
