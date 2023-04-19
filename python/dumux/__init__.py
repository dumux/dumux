# SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
# SPDX-License-Identifier: GPL-3.0-or-later

"""The DuMux Python-bindings for using DuMux entirely from Python

DuMux is
* short for Dune for Multi-{Phase, Component, Scale, Physics, …} flow and transport in porous media
* a free and open-source simulator for flow and transport processes in porous media
* a research code written in C++
* based on Dune (Distributed and Unified Numerics Environment)
* a Dune user module in the Dune environment

https://dumux.org/
"""

try:
    from dune.packagemetadata import registerExternalModule
    import pathlib

    # register dumux to be recognized by dune-py (code generation module)
    # as a module of the dune univers
    registerExternalModule(
        moduleName="dumux",
        modulePath=str(pathlib.Path(__file__).parent.resolve()),
    )
except ImportError:
    pass
