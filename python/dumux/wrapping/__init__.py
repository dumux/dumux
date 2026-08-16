# SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
# SPDX-License-Identifier: GPL-3.0-or-later

"""Helper classes and function related Python <-> C++ interaction"""

import functools
from typing import Callable, Type

# dumux is a namespace package and so has no __init__.py in which to make itself known to
# dune-py. Every subpackage that generates C++ imports this module, so registering here covers
# each of them whichever one an application enters through.
try:
    from dune.packagemetadata import registerExternalModule
    import pathlib

    registerExternalModule(
        moduleName="dumux",
        modulePath=str(pathlib.Path(__file__).parent.parent.resolve()),
    )
except ImportError:
    pass


def cppWrapperCreator(creator: Callable) -> Callable:
    """
    Decorator for creator functions that return a C++ type with Python bindings
    resulting from C++ code generation and just-in-time compilation
    """

    def makeCreator(aliasClass: Type) -> Callable:
        @functools.wraps(creator)
        def _makeCreator(*args, **kwargs):
            return creator(*args, **kwargs)

        # make the creator assume the name of the alias class
        _makeCreator.__name__ = aliasClass.__name__
        return _makeCreator

    return makeCreator


def cppWrapperClassAlias(creator: Callable) -> Callable:
    """
    Decorator for a class alias corresponding to a creator function
    This makes the creator function appear like a class constructor

    Args:
        creator (Callable): The corresponding creator function the
                            decorated class is an alias for
    Returns:
        Callable: The creator function with the alias name
    """

    def makeCreator(aliasClass: Type) -> Callable:
        return creator(aliasClass)

    return makeCreator
