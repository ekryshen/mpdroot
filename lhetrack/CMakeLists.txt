# Create a library called "libLHETrack" 

set(INCLUDE_DIRECTORIES
 ${BASE_INCLUDE_DIRECTORIES}
  ${CMAKE_SOURCE_DIR}/mcstack
  ${CMAKE_SOURCE_DIR}/mpdbase
  ${CMAKE_SOURCE_DIR}/mpdfield 
  ${CMAKE_SOURCE_DIR}/detectors/tpc
  ${CMAKE_SOURCE_DIR}/kalman
  ${CMAKE_SOURCE_DIR}/detectors/sts
  ${CMAKE_SOURCE_DIR}/lhetrack
)

Set(SYSTEM_INCLUDE_DIRECTORIES
 ${ROOT_INCLUDE_DIR}
)

include_directories(${INCLUDE_DIRECTORIES})

set(LINK_DIRECTORIES
 ${ROOT_LIBRARY_DIR}
 ${FAIRROOT_LIBRARY_DIR}
)

link_directories(${LINK_DIRECTORIES})

# List of source files
set(SRCS
#TpcLheTrackFinder.cxx
#TpcLheSegments.cxx
#TpcLheHit.cxx
#TpcLheTrackCuts.cxx      
#TpcLheCMPoint.cxx
#TpcLheTrack.cxx          
#TpcLheCMTrack.cxx        
#TpcLheTrackFitter.cxx
#TpcLheHitsMaker.cxx      
#TpcLhePoint.cxx 
#TpcLheKalmanFilter.cxx
#TpcLheKalmanTrack.cxx
MpdTpcKalmanFilter.cxx
MpdItsKalmanTrack.cxx
MpdEctKalmanTrack.cxx
#MpdEctTrackFinderTof.cxx
#MpdEctTrackFinderTofTpc.cxx
#MpdEctTrackFinderTpc.cxx
#MpdEctTrackFollow2Tpc.cxx
#MpdSftTrackFinderTpc.cxx
#MpdTpcTrackFollow2Sft.cxx
MpdKfPrimaryVertexFinder.cxx
MpdTrackFinderIts5spd.cxx
MpdVertexZfinder.cxx
MpdTpcDedxTask.cxx
#MpdKfV0Fitter.cxx
MpdMotherFitterPart.cxx
#MpdMotherFitterTrac.cxx
MpdParticle.cxx
MpdVector.cxx
MpdVectorFinder.cxx
MpdItsToTpcMatching.cxx
)

# List of header files
set(HEADERS
#TpcLheTrackFinder.h
#TpcLheSegments.h
#TpcLheHit.h
#TpcLheTrackCuts.h      
#TpcLheCMPoint.h
#TpcLheTrack.h          
#TpcLheCMTrack.h        
#TpcLheTrackFitter.h
#TpcLheHitsMaker.h      
#TpcLhePoint.h 
#TpcLheKalmanFilter.h
#TpcLheKalmanTrack.h
MpdTpcKalmanFilter.h
MpdItsKalmanTrack.h
MpdEctKalmanTrack.h
#MpdEctTrackFinderTof.h
#MpdEctTrackFinderTofTpc.h
#MpdEctTrackFinderTpc.h
#MpdEctTrackFollow2Tpc.h
#MpdSftTrackFinderTpc.h
#MpdTpcTrackFollow2Sft.h
MpdKfPrimaryVertexFinder.h
MpdTrackFinderIts5spd.h
MpdVertexZfinder.h
MpdTpcDedxTask.h
#MpdKfV0Fitter.h
MpdMotherFitterPart.h
#MpdMotherFitterTrac.h
#lhe.h
MpdParticle.h
MpdVector.h
MpdVectorFinder.h
MpdItsToTpcMatching.h
)

Set(LINKDEF LHETrackLinkDef.h)
Set(LIBRARY_NAME LHETrack)
Set(DEPENDENCIES Kalman MpdBase)

GENERATE_LIBRARY()
