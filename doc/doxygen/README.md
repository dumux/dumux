# Building the documentation

This section explain how to build the documentation pages that you are reading now.
The documentation is built with the tool [Doxygen](https://www.doxygen.nl/index.html).
Doxygen automatically creates documentation from documented source code and from
separate documentation pages which we write in Markdown
(see [Doxygen Manual](https://www.doxygen.nl/manual/markdown.html) for supported Markdown flavor).
For building the documentation, you need have

* doxygen (version >= 1.9.1)

installed on your system. To build the documentation go to the build folder of dumux.
If you have following the [installation instructions](installation.md) and use the default build
folder this will be `dumux/build-cmake`. Run

    make doc

This will build the documentation as a html page. Building the full documentation can take
several minutes. To view the documentation
open `dumux/build-cmake/doc/doxygen/html/index.html`.

## Quick guide for developing the Doxygen documentation

* Doxygen is configured in `Doxylocal`
* The main html page is configured in the `header.html` template
* Doxygen pages are in the subfolder `pages/`
* We don't add the layout type `pages` in doxygen layout (see `DoxygenDumuxLayout.xml`)
* Add pages and subpages manually into the sidebar (see `DoxygenDumuxLayout.xml`)
* TOCs in markdown files are created with the `[TOC]` command
* Modules are documented in `groups/` in several files
* The use of `@addtogroup` can be used to combine group documentation from several files (use this if the description is long)
* The base style `doxygen awesome` is described in the [Doxygen Awesome Documentation](https://jothepro.github.io/doxygen-awesome-css/)

## Building the doxygen documentation partially

For developing the page documentation it can be helpful to disable the automatic code documentation for a much faster partial build
that only contains the page documentation and group without the code documentation. To this end, comment the
line `@top_srcdir@/dumux \` in `Doxylocal`.

## Instructions for code documentation

Please follow these instructions when commenting your code for the Doxygen class documentation.
When checking doxygen, the doxygen should be built and there should be no `doxyerr.log`. In case this file exists
it contains errors and warning pointing you to the places that have to be fixed. The GitLab CI pipelines of dumux
automatically check whether the documentation can be built without errors and otherwise declines the merge request.

The Doxygen for a **file** should look like this:

@verbatim
/*!
 * \file
 * \ingroup Common
 * \brief Manages the handling of time dependent problems
 */
@endverbatim

It should always contain the `\file` first.
The `\ingroup` gives a Group that is as **precise** as possible and is part of the groups described in a mardown file in `groups/`.
The `\brief` is a **short comment**  on what happens in the file.
Alternatively a `\copybrief` could be used. **Make sure the copybrief links to the correct section and is not ambiguous!**
Additional text could be added similar to the function description.

A class could look like this:

@verbatim
/*!
 * \ingroup Common
 * \brief Manages the handling of time dependent problems.
 *
 * This class facilitates the time management of the simulation.
 * It doesn't manage [...]
 * [...]
 * [...] index starting at 0.
 *
 * \note Time and time step sizes are in units of seconds
 */
@endverbatim

The class should always have the proper `\ingroup` just like the file.
The `\brief` should always contain a **short description**, and only in very very rare occasions a class can have a `\copybrief`.
In the end this is also documentation for users that do not build doxygen and `\copybrief`s from another file are most likely not very useful.
The main part should be an in depth explanation of what is done. This should contain math-expressions where applicable.
They can look like the following:

@verbatim
```
 * This class implements a spline \f$s(x)\f$ for which, given \f$n\f$ sampling
 * points \f$x_1, \dots, x_n\f$, the following conditions hold
 *\f{ align*}{
   s(x_i)   & = y_i \quad \forall i \in \{1, \dots, n \}\\
   s'(x_1)  & = m_1 \\
   s'(x_n)  & = m_n
  \f}
 *
```
@endverbatim

For a **function** in a file the Doxygen can look like this:

@verbatim
/*!
 * \brief Set the current simulated time and the time step index.
 *
 * \param t The time \f$\mathrm{[s]}\f$ which should be jumped to
 * \param stepIdx The new time step index
 */
 void setTime(Scalar t, int stepIdx)
 { time_ = t; timeStepIdx_ = stepIdx; }
@endverbatim

The function has a **short comment** with the `\brief`.
If there are **function-arguments** that are not self-explanatory, they should be described using `\param`.
**Always explain all params or no params at all!** Otherwise Doxygen will throw an error.
Template parameters are documented with `\tparam`.
Additional Doxygen-commands that might be useful are `\note` for giving an important note/hint
on what the function does as well as `\return` which specifies the return value (if applicable).
