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
    from dune.common import registerExternalModule

    # register dumux to be recognized by dune-py (code generation module)
    # as a module of the dune univers
    registerExternalModule("dumux")
except ImportError:
    pass
