# SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
# SPDX-License-Identifier: GPL-3.0-or-later

"""The DuMux common module, containing classes and functions needed for most simulations

DuMux is
* short for Dune for Multi-{Phase, Component, Scale, Physics, …} flow and transport in porous media
* a free and open-source simulator for flow and transport processes in porous media
* a research code written in C++
* based on Dune (Distributed and Unified Numerics Environment)
* a Dune user module in the Dune environment

https://dumux.org/
"""

from dune.generator.generator import SimpleGenerator
from dune.common.hashit import hashIt

from dumux.wrapping import cppWrapperCreator, cppWrapperClassAlias
from dumux.common.properties import Model, Property
from dumux.common.boundarytypes import BoundaryTypes
from dumux.common.fvproblem import FVProblem
from dumux.common.fvspatialparams import FVSpatialParams

from ._common import *

# pylint cannot read the names out of a compiled module without importing it
from ._common import getParam as _getParamAsString  # pylint: disable=no-name-in-module
from ._common import initParameters as _initParameters  # pylint: disable=no-name-in-module


def _toParameterString(value):
    """Spell a Python value the way DuMux's parameter parser reads it back.

    A parameter tree holds text. ``str()`` alone is not enough for the two types whose Python
    spelling differs from DuMux's: ``True`` is not ``true``, and ``[10, 10]`` is not ``10 10``.
    """
    if isinstance(value, bool):
        return "true" if value else "false"
    if isinstance(value, (list, tuple)):
        return " ".join(_toParameterString(entry) for entry in value)
    return str(value)


def initParameters(argv=None, params=None, inputFileName=""):
    """Fill the parameter tree from the command line, an input file and a dict of values.

    Every DuMux run reads its runtime parameters from this tree, so it has to be filled before
    anything that reads one is constructed. Pass ``argv=sys.argv[:1]`` to keep the command line
    out of it.
    """
    _initParameters(
        list(argv) if argv is not None else [],
        {key: _toParameterString(value) for key, value in (params or {}).items()},
        inputFileName,
    )


def getParam(name, default=None):
    """A runtime parameter, read back as the type of ``default``.

    The tree holds text, so without a default this hands back the text as written and raises if
    the parameter was never set. With one, the text is read as that type — a bool comes back a
    bool rather than the non-empty (and therefore true) string ``"false"``.
    """
    if default is None:
        return _getParamAsString(name)

    text = _getParamAsString(name, _toParameterString(default))
    if isinstance(default, (list, tuple)):
        prototype = default[0] if len(default) else ""
        return type(default)(_fromParameterString(name, entry, prototype) for entry in text.split())
    return _fromParameterString(name, text, default)


def _fromParameterString(name, text, prototype):
    """Read ``text`` back as the type of ``prototype``."""
    if isinstance(prototype, bool):
        if text.lower() not in ("true", "false", "1", "0"):
            raise ValueError(f"parameter {name} is {text!r}, which is not a truth value")
        return text.lower() in ("true", "1")
    return type(prototype)(text)
