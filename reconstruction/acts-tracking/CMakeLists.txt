# Create a library called "libKalman" which includes the source files given in
# the array .
# The extension is already found.  Any number of sources could be listed here.

function (GET_ACTS_DIR)
  string (FIND $ENV{PATH} "Acts" ActsFound)
  if (ActsFound STREQUAL "-1")
    set (ACTS_PATH "" PARENT_SCOPE)
    return()
  endif()
  if (UNIX)
    SET (_path_separator "/")
    string(REGEX REPLACE "^(.*[:])?([^:]*Acts[^:]+).*" "\\2" ACTS_PATH $ENV{PATH})
  else()
    SET (_path_separator "\\")
    string(REGEX REPLACE "^(.*[;])?([^;]*Acts[^;]+).*" "\\2" ACTS_PATH $ENV{PATH})
  endif()
  SET (_path_separator ${_path_separator} PARENT_SCOPE)
  set (ACTS_PATH ${ACTS_PATH}${_path_separator}..${_path_separator}  PARENT_SCOPE)
endfunction()

GET_ACTS_DIR()

set(INCLUDE_DIRECTORIES
  ${BASE_INCLUDE_DIRECTORIES}
  ${CMAKE_CURRENT_LIST_DIR}
  ${CMAKE_SOURCE_DIR}/core/mpdBase
  ${CMAKE_SOURCE_DIR}/core/mpdField
  ${CMAKE_SOURCE_DIR}/detectors/tpc
  ${CMAKE_SOURCE_DIR}/detectors/tpc/geometry
  ${CMAKE_SOURCE_DIR}/reconstruction/acts-tracking
  ${CMAKE_SOURCE_DIR}/reconstruction/tracking/kalman
  ${CMAKE_SOURCE_DIR}/reconstruction/tracking/lheTrack
  ${CMAKE_SOURCE_DIR}/simulation/mcStack
  ${ACTS_PATH}/include
)

Set(SYSTEM_INCLUDE_DIRECTORIES
  ${ROOT_INCLUDE_PATH}
)

include_directories(${INCLUDE_DIRECTORIES})

set(LINK_DIRECTORIES
  ${ACTS_PATH}/lib64
  ${Boost_INCLUDE_DIRS}${_path_separator}..${_path_separator}lib
)
set (_path_separator)

link_directories(${LINK_DIRECTORIES})

# List of source files
set(SRCS
  ActsExamples/TGeoDetector/TGeoDetector.cxx
  ActsExamples/TGeoDetector/TGeoITkModuleSplitter.cxx
  MpdTpcConfig.cxx
  MpdTpcDetector.cxx
  MpdTpcDigitization.cxx
  MpdTpcMagneticField.cxx
  MpdTpcRunner.cxx
  MpdTpcSpacePointMaking.cxx
  MpdTpcTrackEstimation.cxx
  MpdTpcTrackFinding.cxx
  MpdTpcTrackSeeding.cxx
  MpdTpcTracker.cxx
  PlotUtils.cxx
)

Set(HEADERS
  ActsExamples/Detector/IBaseDetector.hpp
  ActsExamples/EventData/GeometryContainers.hpp
  ActsExamples/EventData/Index.hpp
  ActsExamples/EventData/IndexSourceLink.hpp
  ActsExamples/EventData/Measurement.hpp
  ActsExamples/EventData/ProtoTrack.hpp
  ActsExamples/EventData/SimSpacePoint.hpp
  ActsExamples/EventData/Track.hpp
  ActsExamples/EventData/Trajectories.hpp
  ActsExamples/TGeoDetector/JsonTGeoDetectorConfig.hpp
  ActsExamples/TGeoDetector/TGeoDetector.hpp
  ActsExamples/TGeoDetector/TGeoITkModuleSplitter.hpp
  ActsExamples/Utilities/GroupBy.hpp
  ActsExamples/Utilities/Options.hpp
  ActsExamples/Utilities/OptionsFwd.hpp
  ActsExamples/Utilities/Range.hpp
  MpdTpcAlgorithm.h
  MpdTpcConfig.h
  MpdTpcContext.h
  MpdTpcDetector.h
  MpdTpcDigitization.h
  MpdTpcEventData.h
  MpdTpcEventStorage.h
  MpdTpcInputHit.h
  MpdTpcMagneticField.h
  MpdTpcRunner.h
  MpdTpcSpacePointMaking.h
  MpdTpcTrackEstimation.h
  MpdTpcTrackFinding.h
  MpdTpcTrackSeeding.h
  MpdTpcTracker.h
  PlotUtils.h
)

Set(LINKDEF MpdTpcTrackerLinkDef.h)
Set(LIBRARY_NAME MpdTpcActsTracker)
Set(DEPENDENCIES Kalman MpdField ActsCore ActsFatras ActsPluginTGeo boost_program_options)

GENERATE_LIBRARY()
