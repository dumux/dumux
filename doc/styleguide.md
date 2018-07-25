# Style guide

## General formatting

* Use 4 spaces indent (no tabs, not 2 spaces)
* _Trailing whitespace_: source files may not contain trailing whitespace to reduce the amount of noise in diffs and during merges.
* In contrast to the remainder of the coding style guidelines, these code formatting rules are (partially) enforced automatically with a pre-commit hook. Due to the distributed nature of git, this hook can only check your commits once they arrive in the central repository, so it is important to make your local git repository check your commits as well. The dunecontrol script will automatically install such a pre-commit hook for you.


## C++

### Documentation

* Please document freely what each part of your code _should_ do. Document assumptions.
* All comments/documentation in English.
* We proclaim the Doc-Me dogma, which means whatever you do, please document it at least with

    ```c++
    //! \todo Please doc me!
    ```

* We use doxygen to generate documentation from the source code

    ```c++
    int lineOfCode = 1; // Short comment on line of code 1 that will _not_ show in doxygen
    int lineOfCode = 2; //!< Short comment on line of code 2 that will show in doxygen
    //! Short comment on the following line (line of code 3) that will show in doxygen
    int lineOfCode = 3;
    /*!
     * Longer comment on line of code 4
     * with several lines of length
     * that will show in doxygen
     */
    int lineOfCode = 4;
    ```

* Files always contain the following documentation header before the headerguard

    ```c++
    /*!
     * \file
     * \ingroup GroupName
     * \brief A short description of the file.
     */
    #ifndef DUMUX_HEADERGUARD_HH
    #define DUMUX_HEADERGUARD_HH
    ```

    where `GroupName` is a doxygen module group, click [here](https://git.iws.uni-stuttgart.de/dumux-repositories/dumux/blob/master/doc/doxygen/modules.txt) for an overview of existing groups.

* Each class should be documented using the following style

    ```c++
    /*!
     * \ingroup GroupName
     * \brief A short description of the class.
     *
     * Optional detailed description of the class
     * over several lines.
     *
     * \tparam T very short description of the template parameter T
     */
    template<class T>
    class MyClass {};
    ```

* Each free function and class member function should be documented using the following style

    ```c++
    /*!
     * \brief A short description of the function.
     *
     * Optional detailed description of the function
     * over several lines.
     *
     * \tparam paramType Very short description of paramType
     * \param paramName Very short description of paramName
     * \return returnName Very short description of returnName
     */
    template<typename paramType>
    paramType functionName(const paramType& paramName)
    {
      ...
      return returnName
    }
    ```

* Also document non-obvious function parameters at call site

    ```c++
    localResidual.eval(/*includeBoundaries=*/true);
    ```

* Possible exceptions thrown by a function

### Naming schemes

* avoid abbreviation, i.e. `saturation` instead of `s`


#### Variables / Functions

* CamelCase starting with a lower case letter, each new word starting with a capital letter
* no underscores (except private data members, see below)
* exception: If and only if a single letter that represents an
             abbreviation or index is followed by a single letter (abbreviation or index),
             CamelCase is __not__ used.
* examples
    - `pw` but `pressureW`, because "pressure" is a word
    - `srnw` but `sReg`, because "Reg" is not an abbreviation of a single letter
    - `pcgw` but `dTauDPi`, because "Tau" and "Pi" are words and longer than a letter
    - `CaCO3`, because we write chemical formulas in their chemically sensible way

* private data members end with an underscore

#### Typenames / Classes / Aliases / Namespaces

* same rules as for variables, except the first letter is capital


#### Filenames / Folders

* lower case letters (filenames, but not foldernames may contain underscores)
* Header files get the suffix `.hh`, implementation files the suffix `.cc`
* Every header file contains a unique header guard. The name should mimic the folder structure, and contain the filename,
  i.e. for a file `common/myfile.hh` it should be

    ```c++
    #ifndef DUMUX_COMMON_MYFILE_HH
    #define DUMUX_COMMON_MYFILE_HH
    ```

#### Macros
* The use of preprocessor macros is strongly discouraged. If you have to use them for whatever reason, please use capital letters only.


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

    ```c++
    namespace Dumux {
    namespace Properties {

    bool here = true;
        bool nothere = false; // not like this

    } // end namespace Properties
    } // end namespace Dumux
    ```

* Use a `Detail` namespace for hiding implementation details, e.g. for template meta programming

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

### Property system

* Prefer class templates with regular template arguments over class templates with a `TypeTag` as template argument

### Exception

* The use of exceptions for error handling is encouraged
* There is a variety of DuMux and Dune-specific exceptions you can throw
* All derive (possibly indirectly) from the class `Dune::Exception` in dune-common.



## Files and folders

* Try to order your new header into the existing directory structure
* Headers are named like the classes they contain (usually one class per file, exception: closely tied helper classes / functions)
* Headers are named lower case only (using underscores only if absolutely necessary for readability)
* Folder names are lower case only
* Tests should be called after model and discretization scheme using underscores and lower case only, e.g.

    ```
    test_1p_tpfa
    test_2p2c_box_infiltration
    ```

## CMake

* Use named arguments only
* TODO

## Python

* TODO
