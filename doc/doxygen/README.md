Please follow these instructions when commenting your code for the Doxygen class documentation.
When checking doxygen, the doxygen should be built and the doxyerr.log should have as little lines/errors as possible.

The Doxygen for a **file** should look like this:
```
/*!
 * \file
 * \ingroup Common
 * \brief Manages the handling of time dependent problems
 */
```
It should always contain the `\file` first.
The `\ingroup` gives a Group that is as **precise** as possible and is part of the `modules.txt`
The `\brief` is a **short comment**  on what happens in the file. Alternatively a `\copybrief` could be used. **Make sure the copybrief links to the correct section and is not ambigous!**
Additional text could be added similar to the function description.

A class could look like this:

```
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
```
The class should always have the proper `\ingroup` just like the file.
The `\brief` should always contain a **short description**, and only in very very rare occasions a class can have a `\copybrief`. In the end this is also documentation for users that do not build doxygen and `\copybriefs` from another file are most likely not very useful.
The main part should be an in depth explanation of what is done. This should contain math-expressions where applicable. They can look like the following:
```
 * This class implements a spline \f$s(x)\f$ for which, given \f$n\f$ sampling
 * points \f$x_1, \dots, x_n\f$, the following conditions hold
 *\f{align*}{
   s(x_i)   & = y_i \quad \forall i \in \{1, \dots, n \}\\
   s'(x_1)  & = m_1 \\
   s'(x_n)  & = m_n
   \f}
*
```

For a **function** in a file the Doxygen can look like this:

```
/*!
 * \brief Set the current simulated time and the time step index.
 *
 * \param t The time \f$\mathrm{[s]}\f$ which should be jumped to
 * \param stepIdx The new time step index
 */
 void setTime(Scalar t, int stepIdx)
     { time_ = t; timeStepIdx_ = stepIdx; }
```
The function has a **short comment** with the `\brief`.
If there are **function-arguments** that are not self-explanatory, they should be described using `\param`.
**Always explain all params or no params at all!** Otherwise Doxygen will throw an error.
Template parameters are documented with `\tparam`.
Additional Doxygen-commands that might be useful are `\note` for giving an important note/hint on what the function does as well as `\return` which specifies the return value (if applicable).

April 2020 note: `\copydoc`,`\copybrief` and `\copydetails`, as well as `@copydoc`,`@copybrief` and `@copydetails` do not work with filenames in the current doxygen version.
The bug is reported and has been fixed in doxygen pull request #7693. (It will work again in doxygen 1_8_18, probably at the end of 2020.)
