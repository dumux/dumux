# SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
# SPDX-License-Identifier: GPL-3.0-or-later

add_subdirectory(bowl)
add_subdirectory(dambreak)
add_subdirectory(poiseuilleflow)
add_subdirectory(roughchannel)

set(SHALLOWWATER_TESTS
    test_shallowwater_bowl
    test_shallowwater_dambreak
    test_shallowwater_poiseuilleflow
    test_shallowwater_poiseuilleflow_unstructured
    test_shallowwater_poiseuilleflow_vertical
    test_shallowwater_roughchannel
    test_shallowwater_roughchannel_limitednikuradse
)

add_custom_target(test_shallowwater
    DEPENDS ${SHALLOWWATER_TESTS}
)
