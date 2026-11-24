# SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
# SPDX-License-Identifier: GPL-3.0-or-later

"""
Checks on a decorated Python class, run before any C++ is generated

A problem or spatial params class is handed to C++, which looks up its callbacks by name while
it assembles. A name that is not there, or one whose signature does not fit, is therefore found
in the middle of an element loop, several C++ frames away from the class that is missing it.
These checks move that report back to the line where the class is decorated.
"""

import inspect


class CallbackError(Exception):
    """A decorated class cannot serve as a problem or as spatial parameters"""


def _parameterCount(callback):
    """How many arguments a callback must be called with, or None if it takes any number"""

    try:
        parameters = inspect.signature(callback).parameters.values()
    except (TypeError, ValueError):  # not introspectable, e.g. a C++ or built-in callable
        return None

    if any(p.kind is inspect.Parameter.VAR_POSITIONAL for p in parameters):
        return None

    return sum(
        1
        for p in parameters
        if p.kind in (inspect.Parameter.POSITIONAL_ONLY, inspect.Parameter.POSITIONAL_OR_KEYWORD)
        and p.default is inspect.Parameter.empty
    )


def checkNumEq(cls):
    """The number of equations decides the C++ types, so it has to be there and be a count"""

    numEq = getattr(cls, "numEq", None)
    if numEq is None:
        raise CallbackError(
            f"{cls.__name__} has no numEq. A problem states how many equations it has as a class "
            f"attribute, e.g. 'numEq = 1', because the C++ types are built from it."
        )
    if not isinstance(numEq, int) or isinstance(numEq, bool) or numEq < 1:
        raise CallbackError(f"{cls.__name__}.numEq is {numEq!r}, expected a positive integer")
    return numEq


def checkCallbacks(cls, required, arities):
    """Check that one of each required group of callbacks is there and takes what C++ passes it

    ``required`` is a list of alternatives: each entry is the group of names of which the class
    has to define at least one. ``arities`` maps a callback name to the number of arguments C++
    calls it with, ``self`` excluded.
    """
    for alternatives in required:
        if not any(hasattr(cls, name) for name in alternatives):
            expected = " or ".join(
                f"{name}(self, {', '.join(_ARGUMENTS[name])})" for name in alternatives
            )
            raise CallbackError(
                f"{cls.__name__} defines none of {list(alternatives)}. Define {expected}."
            )

    for name, arity in arities.items():
        callback = getattr(cls, name, None)
        if callback is None or not callable(callback):
            continue
        count = _parameterCount(callback)
        if count is not None and count != arity + 1:
            raise CallbackError(
                f"{cls.__name__}.{name}{inspect.signature(callback)} cannot be called the way "
                f"DuMux calls it: {name}(self, {', '.join(_ARGUMENTS[name])})"
            )


#: what DuMux passes each callback, used to spell out the expected signature in a message
_ARGUMENTS = {
    "boundaryTypes": ("element", "scvOrScvf"),
    "boundaryTypesAtPos": ("globalPos",),
    "dirichlet": ("element", "scvOrScvf"),
    "dirichletAtPos": ("globalPos",),
    "neumann": ("element", "fvGeometry", "scvf"),
    "neumannAtPos": ("globalPos",),
    "source": ("element", "fvGeometry", "scv"),
    "sourceAtPos": ("globalPos",),
    "initial": ("entity",),
    "permeability": ("element", "scv", "elemSol"),
    "permeabilityAtPos": ("globalPos",),
    "porosity": ("element", "scv", "elemSol"),
    "porosityAtPos": ("globalPos",),
}

_PROBLEM_ARITIES = {name: len(args) for name, args in _ARGUMENTS.items()}


def validateProblem(cls):
    """Check a class decorated as a finite volume problem"""

    checkNumEq(cls)
    checkCallbacks(cls, [("boundaryTypes", "boundaryTypesAtPos")], _PROBLEM_ARITIES)


def validateSpatialParams(cls):
    """Check a class decorated as one-phase spatial parameters"""

    checkCallbacks(cls, [("permeability", "permeabilityAtPos")], _PROBLEM_ARITIES)
