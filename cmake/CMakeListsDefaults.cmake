if(NOT WIN32) # some colors, taken from https://stackoverflow.com/questions/18968979/how-to-get-colorized-output-with-cmake
  string(ASCII 27 Esc)
  set(ColourReset "${Esc}[m")
  set(ColourBold  "${Esc}[1m")
  set(Red         "${Esc}[31m")
  set(Green       "${Esc}[32m")
  set(Blue        "${Esc}[34m")
  set(Yellow      "${Esc}[33m")
  set(White       "${Esc}[37m")
  set(BoldRed     "${Esc}[1;31m")
  set(BoldGreen   "${Esc}[1;32m")
  set(BoldBlue    "${Esc}[1;34m")
  set(BoldYellow  "${Esc}[1;33m")
  set(BoldWhite   "${Esc}[1;37m")
endif()

set(PTH_FOUND "")
set(PTH_MISSING "")
set(PTH_ADDEDNDUM "")

macro(checkPath defaultName envName description)
  if (NOT "${${defaultName}}" STREQUAL "") # I have already set variable with path
    list(APPEND PTH_FOUND "${description}: ${${defaultName}}")
  elseif (NOT "$ENV{${envName}}" STREQUAL "") # I have global variable of expected name
    set(${defaultName} $ENV{${envName}})
    list(APPEND PTH_FOUND "${description}: ${${defaultName}}")
  else () # I don't have global variable of given name
    list(APPEND PTH_MISSING "${description}: -D${defaultName} (variable ${envName})")
    string(APPEND PTH_ADDENDUM "\nexport ${defaultName}=<path to ${envName}>")
  endif()
endmacro()

if ("$ENV{MPDROOT}" STREQUAL "") # user did not set MPDROOT variable
  message(FATAL_ERROR "${BoldRed}\n$MPDROOT not set.${ColourReset}\n${BoldYellow}Build without setting environment variable holding installation (and execution) directory is not allowed.${ColourReset}\n")
endif()

# force CMAKE_INSTALL_PREFIX to MPDROOT
IF(CMAKE_INSTALL_PREFIX_INITIALIZED_TO_DEFAULT)
  SET(CMAKE_INSTALL_PREFIX $ENV{MPDROOT} CACHE PATH "" FORCE)
ENDIF(CMAKE_INSTALL_PREFIX_INITIALIZED_TO_DEFAULT)

checkPath(CMAKE_INSTALL_PREFIX MPDROOT "MPDROOT")
checkPath(FAIRROOT_ROOT FAIRROOT_ROOT "FairRoot root")
checkPath(FAIRLOGGER_ROOT FAIRLOGGER_ROOT "FairLogger root")
checkPath(GEANT3_ROOT GEANT3_ROOT "Geant3 root")
checkPath(BOOST_ROOT BOOST_ROOT "Boost root")
checkPath(EIGEN3_ROOT EIGEN3_ROOT "Eigen3 root")

if ((NOT ("${PC_LIBXML_INCLUDEDIR}" STREQUAL "")) AND (NOT ("${PC_LIBXML_LIBDIR}" STREQUAL ""))) # if user has provided LIBXML path, don't look for LIBXML2_ROOT
  list(APPEND PTH_FOUND "libXML2 include dir: ${PC_LIBXML_INCLUDEDIR}")
  list(APPEND PTH_FOUND "libXML2 library dir: ${PC_LIBXML_LIBDIR}")
else()
  checkPath(LIBXML2_ROOT LIBXML2_ROOT "libXML2 root")
endif()
checkPath(FMT_ROOT FMT_ROOT "FMT root")

if (PTH_MISSING) # Some dependencies are missing. I will write to the screen, what I have found and what I am missing
  if (PTH_FOUND)
    message("${Green}\nPaths (packages) found:${ColourReset}")
    foreach (_PTH_FOUND ${PTH_FOUND})
      message("${BoldGreen}${_PTH_FOUND}${ColourReset}")
    endforeach()
  endif()

  if (PTH_MISSING)
    message("${Red}\nPaths (packages) missing:${ColourReset}")
    foreach (_PTH_MISSING ${PTH_MISSING})
      message("${BoldRed}${_PTH_MISSING}${ColourReset}")
    endforeach()
    message("${Yellow}\nCreate the missing variables:${ColourReset}")
    message("${BoldYellow}${PTH_ADDENDUM}\n${ColourReset}")
    message("${Yellow}and re-run your cmake command${ColourReset}\n")

    # go through input commandline parameters and build proper commandline string. Not necessary, just polite.
    get_cmake_property(CACHE_VARS CACHE_VARIABLES)
    foreach(CACHE_VAR ${CACHE_VARS})
      get_property(CACHE_VAR_HELPSTRING CACHE ${CACHE_VAR} PROPERTY HELPSTRING)
      if(CACHE_VAR_HELPSTRING STREQUAL "No help, variable specified on the command line.")
        get_property(CACHE_VAR_TYPE CACHE ${CACHE_VAR} PROPERTY TYPE)
        if(CACHE_VAR_TYPE STREQUAL "UNINITIALIZED")
          set(CACHE_VAR_TYPE)
        else()
          set(CACHE_VAR_TYPE :${CACHE_VAR_TYPE})
        endif()
        set(CMAKE_ARGS "${CMAKE_ARGS} -D${CACHE_VAR}${CACHE_VAR_TYPE}=\"${${CACHE_VAR}}\"")
      endif()
    endforeach()
    message("${BoldYellow}cmake${CMAKE_ARGS} ${CMAKE_SOURCE_DIR}${ColourReset}\n")
    message(FATAL_ERROR "") # terminate function
  endif()
else() # I have found all dependencies, no text to the screen
  if ("${ENV${Geant3_DIR}}" STREQUAL "") # I have not set Geant3_DIR
    set(ENV{Geant3_DIR} ${GEANT3_ROOT})
  endif()
  set(ENV{PKG_CONFIG_PATH} "$ENV{PKG_CONFIG_PATH}:${LIBXML2_ROOT}/lib/pkgconfig")
endif()

if (CMAKE_INSTALL_PREFIX STREQUAL CMAKE_SOURCE_DIR) # user wants to install MPDROOT into sources directory
  message(FATAL_ERROR "${BoldRed}\n$MPDROOT may not point to the source directory.${ColourReset}\n${BoldYellow}In-source install is not allowed.${ColourReset}\n")
endif()
