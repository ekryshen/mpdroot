# Check if cmake has the required version
CMAKE_MINIMUM_REQUIRED(VERSION 3.11.0 FATAL_ERROR)

cmake_policy(SET CMP0074 NEW) # allow to have _ROOT suffix in variables containing program installation roots

include(cmake/CMakeListsDefaults.cmake) # preliminary check of input variables

enable_language(C CXX Fortran)

# This set of commands is needed at the end of ROOTMacros.cmake when building library (around line 335)
execute_process(COMMAND ${FAIRROOT_ROOT}/bin/fairroot-config --major_version OUTPUT_VARIABLE FAIRROOT_MAJOR_VERSION)
execute_process(COMMAND ${FAIRROOT_ROOT}/bin/fairroot-config --minor_version OUTPUT_VARIABLE FAIRROOT_MINOR_VERSION)
execute_process(COMMAND ${FAIRROOT_ROOT}/bin/fairroot-config --patch_version OUTPUT_VARIABLE FAIRROOT_PATCH_VERSION)
string(STRIP ${FAIRROOT_MAJOR_VERSION} FAIRROOT_MAJOR_VERSION)
string(STRIP ${FAIRROOT_MINOR_VERSION} FAIRROOT_MINOR_VERSION)
string(STRIP ${FAIRROOT_PATCH_VERSION} FAIRROOT_PATCH_VERSION)
SET(FAIRROOT_VERSION "${FAIRROOT_MAJOR_VERSION}.${FAIRROOT_MINOR_VERSION}.${FAIRROOT_PATCH_VERSION}")
SET(FAIRROOT_LIBRARY_PROPERTIES ${FAIRROOT_LIBRARY_PROPERTIES}
    VERSION "${FAIRROOT_VERSION}"
    SOVERSION "${FAIRROOT_MAJOR_VERSION}"
    SUFFIX ".so"
)

list(PREPEND CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/cmake/modules")
list(APPEND CMAKE_MODULE_PATH "$ENV{FAIRROOT_ROOT}/share/fairbase/cmake/modules")
set(FairRoot_DIR ${FAIRROOT_ROOT}) # needed by ROOTMacros.cmake

include(FairMacros) # needed by find_package(ROOT)
find_package(ROOTMpd 0.0.0 REQUIRED) # 0.0.0 - minimal requested version of ROOT - bug in FindRoot.cmake by FairRoot
# ensure that ROOT has FFTW3 support
Execute_Process(COMMAND ${ROOT_CONFIG_EXECUTABLE} --has-fftw3 OUTPUT_VARIABLE ROOT_HAS_FFTW3)
String(STRIP ${ROOT_HAS_FFTW3} ROOT_HAS_FFTW3)
if(NOT ROOT_HAS_FFTW3)
  message(FATAL_ERROR "${BoldRed}\nROOT was not built with FFTW3 support. Rebuild ROOT and try again.${ColourReset}\n")  
endif()

find_package(FairRoot REQUIRED)
find_package(FairLogger REQUIRED)
find_package(FMT REQUIRED)
find_package(Pythia8 REQUIRED)
find_package(Boost REQUIRED)
find_package(LibXml2 REQUIRED)
find_package(Geant3 REQUIRED)
find_package(Eigen3 REQUIRED)
find_package(VMC QUIET)

#set(BASE_INCLUDE_DIRECTORIES ${BASE_INCLUDE_DIRECTORIES} ${SIMPATH}/include/root ${SIMPATH}/include/vmc)
#list(APPEND BASE_INCLUDE_DIR ${ROOT_INCLUDE_DIR} ${FAIRROOT_INCLUDE_DIR} ${FairLogger_INCDIR} "${FMT_ROOT}/include")
list(APPEND BASE_INCLUDE_DIR ${ROOT_INCLUDE_DIR} ${FAIRROOT_INCLUDE_DIR} ${FairLogger_INCDIR} ${FMT_INCLUDE_DIRS}
                             ${PYTHIA8_INCLUDE_DIR} ${Boost_INCLUDE_DIRS}
                             ${LIBXML2_INCLUDE_DIRS} ${Geant3_INCLUDE_DIRS} ${Eigen3_INCLUDE_DIRS}
                             )

if(VMC_FOUND)
  list(APPEND BASE_INCLUDE_DIR ${VMC_INCLUDE_DIRS})
endif()

list(APPEND BASE_LIBRARY_DIR ${ROOT_LIB_DIR} ${FAIRROOT_LIB_DIR} ${FairLogger_LIBDIR} ${LIBXML2_LIBRARIES})

# START this is for compatibility with legacy, remove later
set(BASE_INCLUDE_DIRECTORIES ${BASE_INCLUDE_DIR})
set(BASE_LIBRARY_DIRECTORIES ${BASE_LIBRARY_DIR})
# END this is for compatibility with legacy, remove later

set(CMAKE_INSTALL_LIBDIR "lib")
set(LIBRARY_OUTPUT_PATH "${CMAKE_BINARY_DIR}/${CMAKE_INSTALL_LIBDIR}")

# Recurse into the given subdirectories. This does not actually
# cause another cmake executable to run. The same process will walk through
# the project's entire directory structure.
# LEVEL 1
add_subdirectory (core/mpdBase) # INDEPENDENT
add_subdirectory (core/mpdDst) # Base
add_subdirectory (core/mpdField) # INDEPENDENT
add_subdirectory (core/mpdPassive) # INDEPENDENT
add_subdirectory (core/mpdPid) # INDEPENDENT
add_subdirectory (detectors/bmd) # INDEPENDENT
add_subdirectory (detectors/emc) # INDEPENDENT
add_subdirectory (detectors/etof) # INDEPENDENT
add_subdirectory (detectors/ffd) # INDEPENDENT
add_subdirectory (detectors/mcord) # INDEPENDENT
add_subdirectory (detectors/sts) # INDEPENDENT
add_subdirectory (detectors/tof) # INDEPENDENT
add_subdirectory (detectors/zdc) # INDEPENDENT
add_subdirectory (simulation/generators/genFactory) # INDEPENDENT
add_subdirectory (simulation/generators/mpdGen)
add_subdirectory (simulation/generators/mpdGeneralGenerator) # INDEPENDENT
add_subdirectory (simulation/generators/unigenFormat) # INDEPENDENT

Execute_Process(COMMAND ${ROOT_CONFIG_EXECUTABLE} --has-pgsql OUTPUT_VARIABLE ROOT_HAS_PGSQL)
String(STRIP ${ROOT_HAS_PGSQL} ROOT_HAS_PGSQL)
if(ROOT_HAS_PGSQL)
  message("${BoldWhite}\nROOT was built with PGSQL support. Database part will be built.${ColourReset}\n")  
  add_subdirectory (tools/database) # INDEPENDENT
else()
  message("${BoldRed}\nROOT was not built with PGSQL support. Database part will not work.${ColourReset}\n")  
endif()

# LEVEL 2
add_subdirectory (detectors/tpc) # tof
add_subdirectory (reconstruction/tracking/kalman) # mpdfield
add_subdirectory (physics) # mpdbase mpddst
add_subdirectory (simulation/mcDst)
add_subdirectory (simulation/mcStack) # MpdGen
# LEVEL 3
add_subdirectory (tools/eventDisplay) # emc xml2 TODO - remove dependencies on root configuration
add_subdirectory (reconstruction/tracking/lheTrack) # mpdbase kalman

INSTALL(DIRECTORY gconfig/ DESTINATION gconfig)
INSTALL(DIRECTORY input/ DESTINATION input)
INSTALL(DIRECTORY geometry/ DESTINATION geometry)
INSTALL(DIRECTORY macros/ DESTINATION macros)
INSTALL(FILES config/env.sh DESTINATION config)
INSTALL(FILES config/rootmanager.dat DESTINATION config)
