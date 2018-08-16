# Find FairLogger installation
# Once done this will define
#  FAIRLOGGER_FOUND - system has found FairLogger installation
#  FairLogger_INCDIR - FairLogger include directory
#  FairLogger_LIBDIR - FairLogger library directory

find_path(FairLogger_INCDIR NAMES fairlogger/Logger.h PATHS
  ${FAIRLOGGER_ROOT}/include
  $ENV{FAIRLOGGER_ROOT}/include
  ${FAIRROOTPATH}/include
  NO_DEFAULT_PATH
)

find_path(FairLogger_LIBDIR NAMES libFairLogger.so PATHS
  ${FAIRLOGGER_ROOT}/lib
  $ENV{FAIRLOGGER_ROOT}/lib
  ${FAIRROOTPATH}/lib
  NO_DEFAULT_PATH
)

if (FairLogger_INCDIR AND FairLogger_LIBDIR)
  set(FAIRLOGGER_FOUND TRUE)
  SET(LD_LIBRARY_PATH ${LD_LIBRARY_PATH} ${FairLogger_LIBDIR})
  if (NOT FAIRLOGGER_FIND_QUIETLY)
    message(STATUS "Looking for FairLogger... found at ${FairLogger_LIBDIR}")
  endif (NOT FAIRLOGGER_FIND_QUIETLY)
else (FairLogger_INCDIR AND FairLogger_LIBDIR)
  set(FAIRLOGGER_FOUND FALSE)
  if (FAIRLOGGER_FIND_REQUIRED)
    message(FATAL_ERROR "Looking for FairLogger... not found")
  else (FAIRLOGGER_FIND_REQUIRED)
    message(STATUS "Looking for FairLogger... not found")
  endif (FAIRLOGGER_FIND_REQUIRED)
endif (FairLogger_INCDIR AND FairLogger_LIBDIR)
