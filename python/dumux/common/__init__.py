# SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
# SPDX-License-Identifier: GPL-3.0-or-later

"""
The DuMux common module
containing classes and functions needed for most simulations
"""

from dune.generator.generator import SimpleGenerator
from dune.common.hashit import hashIt

from dumux.common.properties import Model, Property
from dumux.common.boundarytypes import BoundaryTypes
from dumux.common.fvproblem import FVProblem
from dumux.common.fvspatialparams import FVSpatialParams
from dumux.wrapping import cppWrapperCreator, cppWrapperClassAlias

from ._common import *
