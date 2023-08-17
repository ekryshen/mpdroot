# Create a library called "libKalman" which includes the source files given in
# the array .
# The extension is already found.  Any number of sources could be listed here.

function (GET_ACTS_DIR)
  string (FIND $ENV{PATH} "Acts" ActsFound)
  if (ActsFound STREQUAL "-1")
    set (ACTS_PATH "" PARENT_SCOPE)
    return()
  endif()
  string(REGEX REPLACE "^(.*[:])?([^:]*Acts[^:]+).*" "\\2" ACTS_PATH $ENV{PATH})
  set (ACTS_PATH ${ACTS_PATH}/..  PARENT_SCOPE)
endfunction()

GET_ACTS_DIR()

set(INCLUDE_DIRECTORIES
  ${ACTS_PATH}/include
  ${BASE_INCLUDE_DIRECTORIES}
  ${CMAKE_CURRENT_LIST_DIR}
  ${CMAKE_SOURCE_DIR}/core/mpdBase
# ${CMAKE_SOURCE_DIR}/core/mpdField
  ${CMAKE_SOURCE_DIR}/detectors/tpc
  ${CMAKE_SOURCE_DIR}/detectors/tpc/clusterHitFinder/mlem
  ${CMAKE_SOURCE_DIR}/detectors/tpc/fairTpc
  ${CMAKE_SOURCE_DIR}/reconstruction/acts-tracking
  ${CMAKE_SOURCE_DIR}/reconstruction/tracking/kalman
# ${CMAKE_SOURCE_DIR}/reconstruction/tracking/lheTrack
  ${CMAKE_SOURCE_DIR}/simulation/mcStack
)

Set(SYSTEM_INCLUDE_DIRECTORIES
  ${ROOT_INCLUDE_PATH}
)

include_directories(${INCLUDE_DIRECTORIES})

set(LINK_DIRECTORIES
  ${ACTS_PATH}/lib64
)

link_directories(${LINK_DIRECTORIES})

# List of source files
set(SRCS
  ActsExamples/Io/Performance/CKFPerformanceWriter.cxx
  ActsExamples/TGeoDetector/TGeoDetector.cxx
  ActsExamples/TGeoDetector/TGeoITkModuleSplitter.cxx
  DumpUtils.cxx
  MpdTpcConfig.cxx
  MpdTpcDetector.cxx
  MpdTpcDigitization.cxx
  MpdTpcMagneticField.cxx
  MpdTpcParticleSelector.cxx
  MpdTpcRunner.cxx
  MpdTpcSpacePointMaking.cxx
  MpdTpcTrack.cxx
  MpdTpcTrackEstimation.cxx
  MpdTpcTrackFinding.cxx
  MpdTpcTrackSeeding.cxx
  MpdTpcTracker.cxx
  PlotUtils.cxx
  PerfUtils.cxx
)

Set(HEADERS
  ActsExamples/Detector/IBaseDetector.hpp
  ActsExamples/EventData/GeometryContainers.hpp
  ActsExamples/EventData/Index.hpp
  ActsExamples/EventData/IndexSourceLink.hpp
  ActsExamples/EventData/Measurement.hpp
  ActsExamples/EventData/ProtoTrack.hpp
  ActsExamples/EventData/SimSpacePoint.hpp
  ActsExamples/EventData/SimParticle.hpp
  ActsExamples/EventData/Track.hpp
  ActsExamples/EventData/Trajectories.hpp
  ActsExamples/Framework/AlgorithmContext.hpp
  ActsExamples/Framework/BareAlgorithm.hpp
  ActsExamples/Framework/IAlgorithm.hpp
  ActsExamples/Framework/IWriter.hpp
  ActsExamples/Framework/ProcessCode.hpp
  ActsExamples/Framework/WhiteBoard.hpp
  ActsExamples/Framework/WriterT.hpp
  ActsExamples/Io/Performance/CKFPerformanceWriter.hpp
  ActsExamples/TGeoDetector/JsonTGeoDetectorConfig.hpp
  ActsExamples/TGeoDetector/TGeoDetector.hpp
  ActsExamples/TGeoDetector/TGeoITkModuleSplitter.hpp
  ActsExamples/TruthTracking/TruthSeedSelector.hpp
  ActsExamples/Utilities/GroupBy.hpp
  ActsExamples/Utilities/Helpers.hpp
  ActsExamples/Utilities/Options.hpp
  ActsExamples/Utilities/OptionsFwd.hpp
  ActsExamples/Utilities/Paths.hpp
  ActsExamples/Utilities/Range.hpp
  ActsExamples/Validation/DuplicationPlotTool.hpp
  ActsExamples/Validation/EffPlotTool.hpp
  ActsExamples/Validation/FakeRatePlotTool.hpp
  ActsExamples/Validation/ResPlotTool.hpp
  ActsExamples/Validation/TrackClassification.hpp
  ActsExamples/Validation/TrackSummaryPlotTool.hpp
  DumpUtils.h
  MpdTpcAlgorithm.h
  MpdTpcConfig.h
  MpdTpcDetector.h
  MpdTpcDigitization.h
  MpdTpcEventData.h
  MpdTpcInputHit.h
  MpdTpcMagneticField.h
  MpdTpcParticleSelector.h
  MpdTpcRunner.h
  MpdTpcSpacePointMaking.h
  MpdTpcTrack.h
  MpdTpcTrackEstimation.h
  MpdTpcTrackFinding.h
  MpdTpcTrackSeeding.h
  MpdTpcTracker.h
  PlotUtils.h
  PerfUtils.h
)

Set(LINKDEF MpdTpcTrackerLinkDef.h)
Set(LIBRARY_NAME MpdTpcActsTracker)
#                Kalman MpdField                                                       boost_program_options
Set(DEPENDENCIES                 ActsCore ActsFatras ActsPluginTGeo ActsExamplesCommon)

GENERATE_LIBRARY()
