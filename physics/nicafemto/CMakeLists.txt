if("${PROJECT_NAME}" STREQUAL "MPDROOT")
	set(CMAKE_SOURCE_FAIRMCSTACK "${CMAKE_SOURCE_DIR}/mcstack")
else()
	set(CMAKE_SOURCE_FAIRMCSTACK "${CMAKE_SOURCE_DIR}/examples/common/mcstack")
endif()
set(CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/external/nicafemto/cmake/modules"  ${CMAKE_MODULE_PATH})
include(NicaFemto)
SetHeaderVariables(${CMAKE_SOURCE_DIR}/external/nicafemto)
set(NICAFEMTO_MPDCUTS 
    ${CMAKE_SOURCE_DIR}/physics/nicafemto/cuts
    ${CMAKE_SOURCE_DIR}/physics/nicafemto/cuts/paircuts/common
    ${CMAKE_SOURCE_DIR}/physics/nicafemto/cuts/paircuts/common/tpc_cylinder
    ${CMAKE_SOURCE_DIR}/physics/nicafemto/cuts/paircuts/full
    ${CMAKE_SOURCE_DIR}/physics/nicafemto/cuts/paircuts/full/tpc_modular
    ${CMAKE_SOURCE_DIR}/physics/nicafemto/cuts/paircuts/minidst
    ${CMAKE_SOURCE_DIR}/physics/nicafemto/cuts/trackcuts
    ${CMAKE_SOURCE_DIR}/physics/nicafemto/cuts/cutmonitors
)

set(NICAFEMTO_MPDFORMAT
    ${CMAKE_SOURCE_DIR}/physics/nicafemto/format/combined
    ${CMAKE_SOURCE_DIR}/physics/nicafemto/format/minidst
    ${CMAKE_SOURCE_DIR}/physics/nicafemto/format/hbtformat
    ${CMAKE_SOURCE_DIR}/physics/nicafemto/format/mpdevent
    ${CMAKE_SOURCE_DIR}/physics/nicafemto/format/mpdmc
    ${CMAKE_SOURCE_DIR}/physics/nicafemto/format/unigen
    ${CMAKE_SOURCE_DIR}/physics/nicafemto/format/
)
set (NICAFEMTO_MPDHELPERS
    ${CMAKE_SOURCE_DIR}/physics/nicafemto/nica_helpers
)
add_subdirectory(cuts)
add_subdirectory(format)
add_subdirectory(nica_helpers)
add_subdirectory(tasks)