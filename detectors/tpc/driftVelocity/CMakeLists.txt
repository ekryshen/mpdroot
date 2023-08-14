# Create a library called "libtpcPid"

set(INCLUDE_DIRECTORIES
 ${BASE_INCLUDE_DIRECTORIES}
 ${CMAKE_CURRENT_LIST_DIR}
 ${CMAKE_SOURCE_DIR}/detectors/tpc
 ${CMAKE_SOURCE_DIR}/detectors/tpc/clusterHitFinder/mlem
 ${CMAKE_SOURCE_DIR}/detectors/tpc/digitizer/AZlt
 ${CMAKE_SOURCE_DIR}/detectors/tpc/fairTpc
 ${CMAKE_SOURCE_DIR}/reconstruction/tracking/kalman
 ${CMAKE_SOURCE_DIR}/simulation/mcStack
 ${CMAKE_SOURCE_DIR}/tools/database/uniDb
 ${CMAKE_SOURCE_DIR}/tools/database/uniCommon
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
 curve_fit.cxx
 TpcLaserGridVelocityDbUtils.cxx
 TpcLaserGridVelocityMakeMap.cxx
 TpcLaserGridVelocityAdjByMap.cxx
 TpcLaserGridPreprocess.cxx
)

Set(HEADERS
 curve_fit.h
 TpcLaserGridVelocityDbUtils.h
 TpcLaserGridConsts.h
 TpcLaserGridVelocityMakeMap.h
 TpcLaserGridVelocityAdjByMap.h
 TpcLaserGridPreprocess.h
)

Set(LINKDEF tpcDriftVelLinkDef.h)
Set(LIBRARY_NAME tpcDriftVelocity)
Set(DEPENDENCIES MathMore tpcClusterHitFinder UniDb)

GENERATE_LIBRARY()
