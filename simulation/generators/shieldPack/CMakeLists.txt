set(TSHIELD_ROOT_PATH ${CMAKE_CURRENT_SOURCE_DIR})
include(cmake/TShieldConfig.cmake)

add_subdirectory (hadgen)
add_subdirectory (shield)
add_subdirectory (THadgen)
add_subdirectory (TShield)
add_subdirectory (TShieldGenerator)
