enable_language (C CXX Fortran)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS}")
set(INCLUDE_DIRECTORIES ${ROOT_INCLUDE_DIR} ${CMAKE_CURRENT_SOURCE_DIR}/inc ${TSHIELD_ROOT_PATH}/shield/inc ${TSHIELD_ROOT_PATH}/hadgen/inc ${TSHIELD_ROOT_PATH}/THadgen/inc)
include_directories(${INCLUDE_DIRECTORIES})
link_directories(${ROOT_LIBRARY_DIR})

configure_file(inc/TShield.h inc/TShield.h COPYONLY)
configure_file(inc/LinkDef.h inc/LinkDef.h COPYONLY)
set(G__Shield_HEADERS inc/TShield.h)
set(G__Shield_LINKDEF inc/LinkDef.h)
set(G__Shield_DICTIONARY ${CMAKE_CURRENT_BINARY_DIR}/G__Shield.cxx)
ROOT_GENERATE_DICTIONARY("${G__Shield_HEADERS}" "${G__Shield_LINKDEF}" "${G__Shield_DICTIONARY}" "${INCLUDE_DIRECTORIES}")

set(TShield_SRCS
    src/TShield.cxx
    src/TShieldCallbacks.cxx
    src/TShieldOpt.cxx
    src/TShieldGeometry.cxx
    src/TShieldGeometryBool.cxx
    src/TShieldGeometryConvert.cxx
    src/TShieldGeometryOperators.cxx
    src/TShieldGeometryPrint.cxx

    ${G__Shield_DICTIONARY}
)

add_library(TShield SHARED ${TShield_SRCS})
target_link_libraries(TShield SHIELD HADGEN ${ROOT_LIBRARIES})

