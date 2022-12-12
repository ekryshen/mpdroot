# Create a library called "libEventDisplay" 

Set(INCLUDE_DIRECTORIES
  ${BASE_INCLUDE_DIRECTORIES}
  ${CMAKE_CURRENT_LIST_DIR}
  ${CMAKE_CURRENT_LIST_DIR}/mpd
  ${CMAKE_SOURCE_DIR}/config
)

Set(SYSTEM_INCLUDE_DIRECTORIES
  ${ROOT_INCLUDE_DIR} 
  ${Boost_INCLUDE_DIRS}
  ${GEANT3_INCLUDE_DIR}
  ${XML2_INCLUDE_DIR}
)

set(LINK_DIRECTORIES
  ${ROOT_LIBRARY_DIR}
  ${XML2_LIBRARY_DIR}
  ${FAIRROOT_LIBRARY_DIR}
)

set(SRCS
  MpdEventManager.cxx
  MpdMCTracks.cxx
  MpdMCStack.cxx
  MpdEventManagerEditor.cxx
  MpdMCTracksEditor.cxx
  MpdBoxSet.cxx
  MpdBoxSetDraw.cxx
  MpdBoxSetEditor.cxx
  MpdPointSetDraw.cxx
  MpdMCPointDraw.cxx
  MpdHitDraw.cxx
  MpdXMLNode.cxx
  MpdHitPointSetDraw.cxx
  MpdWebScreenshots.cxx
  mpd/MpdGlobalTrackDraw.cxx
  mpd/MpdEmcTowerDraw.cxx
  mpd/MpdZdcTowerDraw.cxx
)

Set(INCLUDE_DIRECTORIES
  ${INCLUDE_DIRECTORIES}
  ${CMAKE_SOURCE_DIR}/core/mpdBase
  ${CMAKE_SOURCE_DIR}/core/mpdDst
  ${CMAKE_SOURCE_DIR}/detectors/emc
  ${CMAKE_SOURCE_DIR}/detectors/tpc
  ${CMAKE_SOURCE_DIR}/detectors/tpc/clusterHitFinder
  ${CMAKE_SOURCE_DIR}/detectors/zdc
  ${CMAKE_SOURCE_DIR}/reconstruction/tracking/kalman
  ${CMAKE_SOURCE_DIR}/reconstruction/tracking/lheTrack
  ${CMAKE_SOURCE_DIR}/simulation/mcStack
)

# TDAQ not implemented yet
# if (TDAQ_FOUND)
#   Set(SYSTEM_INCLUDE_DIRECTORIES
#     ${SYSTEM_INCLUDE_DIRECTORIES} 
#     ${TDAQ_INCLUDE_PATH}
#   )
#   Set(LINK_DIRECTORIES
#     ${LINK_DIRECTORIES} 
#     ${TDAQ_LIBRARY_PATH}
#   )
# endif (TDAQ_FOUND)

Include_Directories(${INCLUDE_DIRECTORIES})
Include_Directories(SYSTEM ${SYSTEM_INCLUDE_DIRECTORIES})

link_directories(${LINK_DIRECTORIES})

Set(HEADERS)

# TDAQ not implemented yet
# if (TDAQ_FOUND)
#   Set(LINKDEF mpd/EventMDLinkDef.h)
# else (TDAQ_FOUND)
  Set(LINKDEF mpd/EventMLinkDef.h)
# endif (TDAQ_FOUND)

if (DEFINED LIBXML2_LIBRARIES) # if we have found XML2 library using find_package, then use this found libxml2
  SET(MY_XML2 ${LIBXML2_LIBRARIES})
else()
  SET(MY_XML2 xml2)
endif()

Set(DEPENDENCIES Emc Geane Base Eve Geom EG Physics Ged Gui Graf3d Tree Core ${MY_XML2})

# TDAQ not implemented yet
# if (TDAQ_FOUND)
#   Set(DEPENDENCIES emon emon-dal cmdline ipc owl omniORB4 omnithread ers ${DEPENDENCIES})
# endif (TDAQ_FOUND)

Set(LIBRARY_NAME EventDisplay)

GENERATE_LIBRARY()

INSTALL(FILES config/coloring.xsd DESTINATION etc)
INSTALL(FILES config/eventDisplay.xml DESTINATION etc)
