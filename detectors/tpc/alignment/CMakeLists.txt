# Create a library called "libtpcAlignment"

set(INCLUDE_DIRECTORIES
  ${BASE_INCLUDE_DIRECTORIES}
  ${CMAKE_CURRENT_LIST_DIR}
  ${CMAKE_CURRENT_LIST_DIR}/muons
  ${CMAKE_CURRENT_LIST_DIR}/laserRays
  ${CMAKE_SOURCE_DIR}/detectors/tpc/
  ${CMAKE_SOURCE_DIR}/detectors/tpc/geometry
)

set(SYSTEM_INCLUDE_DIRECTORIES
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
  #muons/TpcMissAlignment.cxx
  muons/TpcAlignment.cxx
  muons/Mille.cxx
  muons/DerLMpdTpc.cxx
  muons/DerHMpdTpc.cxx
  muons/RecoMiniMC.cxx
  #muons/TpcMiniMC.cxx
  laserRays/CorrectionMatrix.cxx
  laserRays/Debug.cxx
  laserRays/EigenCalculator.cxx
  laserRays/EigenCoefficients.cxx
  laserRays/FileHelper.cxx
  laserRays/Linalg.cxx
  laserRays/Point.cxx
  laserRays/Rays.cxx
  laserRays/Runner.cxx
  laserRays/SystemOfLinearEquations.cxx
  laserRays/TpcAlignmentTaskLaserRays.cxx
  laserRays/TpcCorrectionTaskLaserRays.cxx
  laserRays/Track.cxx
)

Set(HEADERS
  #muons/TpcMissAlignment.h
  muons/TpcAlignment.h
  muons/Mille.h
  muons/DerLMpdTpc.h
  muons/DerHMpdTpc.h
  muons/RecoMiniMC.h
  #muons/TpcMiniMC.h
  laserRays/Constants.h
  laserRays/CorrectionMatrix.h
  laserRays/Debug.h
  laserRays/EigenCalculator.h
  laserRays/EigenCoefficients.h
  laserRays/Enums.h
  laserRays/FileHelper.h
  laserRays/Linalg.h
  laserRays/Point.h
  laserRays/Rays.h
  laserRays/Runner.h
  laserRays/SystemOfLinearEquations.h
  laserRays/TpcAlignmentTaskLaserRays.h
  laserRays/TpcCorrectionTaskLaserRays.h
  laserRays/Track.h
  laserRays/Types.h
)

# see alignment/muons/macro/ReadMe.md
if(ALIGNMENT_XI2)
  set (SRCS ${SRCS} muons/TpcMiniMC.cxx)
  set (HEADERS ${HEADERS} muons/TpcMiniMC.h)
else()
  set (SRCS ${SRCS} muons/TpcMissAlignment.cxx)
  set (HEADERS ${HEADERS} muons/TpcMissAlignment.h)
endif(ALIGNMENT_XI2)

Set(LINKDEF tpcAlignmentLinkDef.h)
Set(LIBRARY_NAME tpcAlignment)
Set(DEPENDENCIES Core Base tpc MpdField MpdBase Minuit2)

GENERATE_LIBRARY()
