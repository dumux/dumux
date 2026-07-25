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

# register dumux with dune-py; dumux is a namespace package and so has no __init__.py
try:
    from dune.packagemetadata import registerExternalModule
    import pathlib

    registerExternalModule(
        moduleName="dumux",
        modulePath=str(pathlib.Path(__file__).parent.parent.resolve()),
    )
except ImportError:
    pass

from dumux.common.properties import Model, Property
from dumux.common.boundarytypes import BoundaryTypes
from dumux.common.fvproblem import FVProblem
from dumux.common.fvspatialparams import FVSpatialParams
from dumux.wrapping import cppWrapperCreator, cppWrapperClassAlias

from ._common import *
