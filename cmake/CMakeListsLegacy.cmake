# The name of our project is "MPDROOT".  CMakeLists files in this project can
# refer to the root source directory of the project as ${MPDROOT_SOURCE_DIR}
# or as ${CMAKE_SOURCE_DIR} and to the root binary directory of the project as
# ${MPDROOT_BINARY_DIR} or ${CMAKE_BINARY_DIR}.

# Check if cmake has the required version
CMAKE_MINIMUM_REQUIRED(VERSION 3.0 FATAL_ERROR)
enable_language(C CXX Fortran)

### CMP0025   Compiler id for Apple Clang is now AppleClang.
### CMP0042   MACOSX_RPATH is enabled by default.

foreach(p
  CMP0025 # CMake 3.0
  CMP0028 # double colon for imported and alias targets
  CMP0042 # CMake 3.0
  CMP0054 # Only interpret ``if()`` arguments as variables or keywords when unquoted.
  )
  if(POLICY ${p})
  cmake_policy(SET ${p} NEW)
  endif()
endforeach()

#In case you need Fortran
set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -std=legacy")

Option(USE_PATH_INFO "Information from PATH and LD_LIBRARY_PATH are used." ON)
Option(ALIBUILD "Flag if we are building with AliBuild." OFF)

#Check if necessary environment variables are set
#If not stop execution
IF (NOT ALIBUILD)
 IF(NOT DEFINED ENV{FAIRROOTPATH})
  MESSAGE(FATAL_ERROR "You did not define the environment variable FAIRROOTPATH which is needed to find FairRoot. Please set this variable and execute cmake again.")
 ENDIF(NOT DEFINED ENV{FAIRROOTPATH})
 IF(NOT DEFINED ENV{SIMPATH})
  MESSAGE(FATAL_ERROR "You did not define the environment variable SIMPATH which is nedded to find the external packages. Please set this variable and execute cmake again.") 
 ENDIF(NOT DEFINED ENV{SIMPATH})
 IF(NOT DEFINED ENV{ROOTSYS})
  MESSAGE(FATAL_ERROR "You did not define the environment variable ROOTSYS which is nedded to find Root. Please set this variable and execute cmake again.")
 ENDIF(NOT DEFINED ENV{ROOTSYS})

 SET(SIMPATH $ENV{SIMPATH})
 SET(FAIRROOTPATH $ENV{FAIRROOTPATH})
ELSE (NOT ALIBUILD)
  SET(FAIRROOTPATH $ENV{FAIRROOT_ROOT})
  SET(GEANT3_PATH $ENV{GEANT3_ROOT})
  SET(BOOST_ROOT $ENV{BOOST_ROOT})
  SET(GSL_DIR $ENV{GSL_ROOT})
ENDIF (NOT ALIBUILD)

# where to look first for cmake modules, before ${CMAKE_ROOT}/Modules/ 
# is checked
set(CMAKE_MODULE_PATH "${FAIRROOTPATH}/share/fairbase/cmake/modules" ${CMAKE_MODULE_PATH})
set(CMAKE_MODULE_PATH "${FAIRROOTPATH}/share/fairbase/cmake/modules_old" ${CMAKE_MODULE_PATH})
set(CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/cmake/modulesLegacy" ${CMAKE_MODULE_PATH})

Set(CheckSrcDir "${FAIRROOTPATH}/share/fairbase/cmake/checks")

find_package(FairRoot REQUIRED)
# Taken from FairSoftLib.cmake start
set(CMAKE_CONFIGURATION_TYPES "Debug" "Release" "RelWithDebInfo")
## this warnings disable compilation due to the errors in code don't use until you know you need them set(_warnings "-Wshadow -Wall -Wextra -Wpedantic")
set(CMAKE_C_FLAGS_DEBUG                "-Og -g ${_warnings}")
set(CMAKE_C_FLAGS_RELEASE              "-O2 -DNDEBUG")
set(CMAKE_C_FLAGS_RELWITHDEBINFO       "-O2 -g ${_warnings} -DNDEBUG")
set(CMAKE_C_FLAGS_RELWITHDEBUGINFO       "-O2 -g ${_warnings} -DNDEBUG")
set(CMAKE_CXX_FLAGS_DEBUG              "-Og -g ${_warnings}")
set(CMAKE_CXX_FLAGS_RELEASE            "-O2 -DNDEBUG")
set(CMAKE_CXX_FLAGS_RELWITHDEBINFO     "-O2 -g ${_warnings} -DNDEBUG")
set(CMAKE_CXX_FLAGS_RELWITHDEBUGINFO     "-O2 -g ${_warnings} -DNDEBUG")
set(CMAKE_Fortran_FLAGS_DEBUG          "-Og -g ${_warnings}")
set(CMAKE_Fortran_FLAGS_RELEASE        "-O2 -DNDEBUG")
set(CMAKE_Fortran_FLAGS_RELWITHDEBINFO "-O2 -g ${_warnings} -DNDEBUG")
set(CMAKE_Fortran_FLAGS_RELWITHDEBUGINFO "-O2 -g ${_warnings} -DNDEBUG")
## JANO - unset if previously set before flags    unset(_warnings)
# Taken from FairSoftLib.cmake end

# Load some basic macros which are needed later on
include(FairMacros)
include(WriteConfigFile)
include(CTest)
include(CheckCompiler)
include(CheckFortran)
include(InstallNicaFemto)
include(CreateMpdConfig)

# Set the build type: None, Debug, Release, RelWithDebInfo or MinSizeRel 
If(NOT CMAKE_BUILD_TYPE)
  Message("Set BuildType RelWithDebInfo")
  set(CMAKE_BUILD_TYPE RelWithDebInfo)
EndIf(NOT CMAKE_BUILD_TYPE)

#Check the compiler and set the compile and link flags
Check_Compiler()
#custom flags for build type if you so desire
#set(CMAKE_CXX_FLAGS_RELWITHDEBINFO "-O3 -g ") # -Wall -Wextra  -std=c++17
#tune for local cpu: -march=native/rocketlake/znver2 -mtune=znver2

set(LIBRARY_OUTPUT_PATH "${CMAKE_BINARY_DIR}/lib")
set(EXECUTABLE_OUTPUT_PATH "${CMAKE_BINARY_DIR}/bin")
set(INCLUDE_OUTPUT_DIRECTORY "${PROJECT_BINARY_DIR}/include")
Set(VMCWORKDIR ${CMAKE_SOURCE_DIR})

If(USE_PATH_INFO)
  Set(PATH "$PATH")
  If (APPLE)
    Set(LD_LIBRARY_PATH $ENV{DYLD_LIBRARY_PATH})
  Else (APPLE)
    Set(LD_LIBRARY_PATH $ENV{LD_LIBRARY_PATH})
  EndIf (APPLE)
Else(USE_PATH_INFO)
  STRING(REGEX MATCHALL "[^:]+" PATH $ENV{PATH})
EndIf(USE_PATH_INFO)

# Check if the user wants to build the project in the source
# directory and if the install directory is different from the build
# directory
CHECK_OUT_OF_SOURCE_BUILD()
CHECK_INSTALL_DIRECTORY()

# Check if we are on an UNIX system. If not stop with an error
# message
IF(NOT UNIX)
  MESSAGE(FATAL_ERROR "You're not on an UNIX system. The project was up to now only tested on UNIX systems, so we break here. If you want to go on please edit the CMakeLists.txt in the source directory.")
ENDIF(NOT UNIX)

# Check if the external packages are installed into a separate install
# directory
CHECK_EXTERNAL_PACKAGE_INSTALL_DIR()

# searches for needed packages
# REQUIRED means that cmake will stop if this packages are not found
# For example the framework can run without GEANT4, but ROOT is
# mandatory
find_package(ROOT 6.10.08 REQUIRED)
find_package(Pythia8)
find_package(Pythia6)
find_package(GEANT3 REQUIRED)
find_package(GEANT4)
find_package(GEANT4DATA)
find_package(GEANT4VMC)
find_package(CLHEP REQUIRED)
find_package(SSE)
find_package(PLUTO)
find_package(GENERATORS REQUIRED)
find_package(XML2 REQUIRED)
find_package(FFTW REQUIRED)
find_package(TDAQ)
find_package(Eigen3 REQUIRED)

#find_package(HEPMC)
#find_package(CUDA)
#find_package(IWYU)
#find_package(ZeroMQ)
#find_package(Protobuf)
#find_package(DDS)

IF (NOT ALIBUILD)
Set(Boost_NO_SYSTEM_PATHS TRUE)
Set(Boost_NO_BOOST_CMAKE TRUE)
If(${ROOT_LIBRARY_DIR} MATCHES /lib/root)
  set(BOOST_ROOT ${SIMPATH})
  set(GSL_DIR ${SIMPATH})
Else(${ROOT_LIBRARY_DIR} MATCHES /lib/root)
  set(BOOST_ROOT ${SIMPATH})
#  set(GSL_DIR ${SIMPATH}/lib/root)
EndIf(${ROOT_LIBRARY_DIR} MATCHES /lib/root)
ENDIF (NOT ALIBUILD)
Message("-- Looking for Boost ...")
# If an older version of boost is found both of the variables below are
# cached and in a second cmake run, a good boost version is found even 
# if the version is to old. 
# To overcome this problem both variables are cleared before checking
# for boost.
Unset(Boost_INCLUDE_DIR CACHE)
Unset(Boost_LIBRARY_DIRS CACHE)
find_package(Boost 1.41)
If (Boost_FOUND)
  Set(Boost_Avail 1)
Else (Boost_FOUND)
  Set(Boost_Avail 0)
EndIf (Boost_FOUND)

find_package(GSL REQUIRED)

# set a variable which should be used in all CMakeLists.txt
# Defines all basic include directories from fairbase
SetBasicVariables()
set(BASE_INCLUDE_DIRECTORIES ${BASE_INCLUDE_DIRECTORIES} ${SIMPATH}/include/root ${SIMPATH}/include/vmc)
include_directories(${SIMPATH}/include/root  ${SIMPATH}/include/vmc ${Eigen3_INCLUDE_DIRS})
find_package(FairLogger)
SET(PATH ${EXECUTABLE_OUTPUT_PATH} ${PATH})

# Set the library version in the main CMakeLists.txt
SET(FAIRROOT_MAJOR_VERSION 0)
SET(FAIRROOT_MINOR_VERSION 0)
SET(FAIRROOT_PATCH_VERSION 0)
SET(FAIRROOT_VERSION "${FAIRROOT_MAJOR_VERSION}.${FAIRROOT_MINOR_VERSION}.${FAIRROOT_PATCH_VERSION}")
SET(FAIRROOT_LIBRARY_PROPERTIES ${FAIRROOT_LIBRARY_PROPERTIES}
    VERSION "${FAIRROOT_VERSION}"
    SOVERSION "${FAIRROOT_MAJOR_VERSION}"
    SUFFIX ".so"
)

Generate_Version_Info()

# Recurse into the given subdirectories. This does not actually
# cause another cmake executable to run. The same process will walk through
# the project's entire directory structure.
add_subdirectory (core/mpdPassive)
add_subdirectory (core/mpdBase)
add_subdirectory (core/mpdField)
add_subdirectory (reconstruction/tracking/kalman) #MpdBase MpdField
add_subdirectory (simulation/generators/genFactory) # INDEPENDENT
add_subdirectory (simulation/generators/mpdGen)
add_subdirectory (simulation/generators/mpdGeneralGenerator) # INDEPENDENT
add_subdirectory (simulation/generators/unigenFormat) # INDEPENDENT
add_subdirectory (simulation/generators/shieldPack)
add_subdirectory (simulation/mcDst)
add_subdirectory (simulation/mcStack) # MpdGen
add_subdirectory (core/mpdPid) # MpdMCStack MpdBase
add_subdirectory (reconstruction/tracking/lheTrack) # MpdBase Kalman Sts Tof
add_subdirectory (core/mpdDst) # MpdBase LHETrack MpdPid

add_subdirectory (detectors/tof) # MpdMCStack Kalman
add_subdirectory (detectors/tpc) # Tof
add_subdirectory (detectors/etof)
add_subdirectory (detectors/zdc)
add_subdirectory (detectors/emc) # Cluster
add_subdirectory (detectors/ffd)
add_subdirectory (detectors/sts)
add_subdirectory (detectors/bmd)
add_subdirectory (detectors/mcord)
add_subdirectory (physics) #MpdBase MpdMCStack Kalman MpdPid LHETrack
add_subdirectory (tools/eventDisplay) #Emc

add_subdirectory (macro)
#add_subdirectory (detectors/bbc)

if(EXISTS "${CMAKE_SOURCE_DIR}/macro/nica_scheduler/CMakeLists.txt")
 add_subdirectory (macro/nica_scheduler)
endif()
InstallNicaFemto()

Option(BUILD_DOXYGEN "Build Doxygen" OFF)
if(BUILD_DOXYGEN)
  MESSAGE(STATUS "*** Building the Doxygen documentaion ***")
  ADD_SUBDIRECTORY(doxygen)
endif(BUILD_DOXYGEN)

WRITE_CONFIG_FILE(config.sh)

SET(VMCWORKDIR ${CMAKE_INSTALL_PREFIX})
SET(ROOT_INCLUDE_PATH ${CMAKE_INSTALL_PREFIX}/include ${ROOT_INCLUDE_PATH})
WRITE_CONFIG_FILE(install_config.sh)

# Summary ######################################################################
if(CMAKE_CXX_FLAGS)
  message(STATUS "  ")
  message(STATUS "  ${Cyan}GLOBAL CXX FLAGS${CR}  ${BGreen}${CMAKE_CXX_FLAGS}${CR}")
endif()
if(CMAKE_CONFIGURATION_TYPES)
  message(STATUS "  ")
  message(STATUS "  ${Cyan}BUILD TYPE         CXX FLAGS${CR}")
  string(TOUPPER "${CMAKE_BUILD_TYPE}" selected_type)
  foreach(type IN LISTS CMAKE_CONFIGURATION_TYPES)
   string(TOUPPER "${type}" type_upper)
   if(type_upper STREQUAL selected_type)
     pad("${type}" 18 " " type_padded)
     message(STATUS "${BGreen}* ${type_padded}${CMAKE_CXX_FLAGS_${type_upper}}${CR}")
   else()
     pad("${type}" 18 " " type_padded)
     message(STATUS "  ${BWhite}${type_padded}${CR}${CMAKE_CXX_FLAGS_${type_upper}}")
   endif()
   unset(type_padded)
   unset(type_upper)
 endforeach()
message(STATUS "  ")
message(STATUS "  (Change the build type with ${BMagenta}-DCMAKE_BUILD_TYPE=...${CR})")
endif()
################################################################################
