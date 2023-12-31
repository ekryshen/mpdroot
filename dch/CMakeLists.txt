# Create a library called "libBase" which includes the source files given in
# the array .
# The extension is already found.  Any number of sources could be listed here.

set(INCLUDE_DIRECTORIES
${ROOT_INCLUDE_DIR} 
${CMAKE_SOURCE_DIR}/geobase 
${CMAKE_SOURCE_DIR}/parbase 
${CMAKE_SOURCE_DIR}/base
${CMAKE_SOURCE_DIR}/cbmbase 
${CMAKE_SOURCE_DIR}/mcstack 
${CMAKE_SOURCE_DIR}/dch
${CMAKE_SOURCE_DIR}/fairtools
)

include_directories(${INCLUDE_DIRECTORIES})

set(LINK_DIRECTORIES
${ROOT_LIBRARY_DIR}
)
 
link_directories(${LINK_DIRECTORIES})

set(Dch_SRCS
MpdDchGeo.cxx 
MpdDchGeoPar.cxx 
MpdDch.cxx 
MpdDchContFact.cxx 
MpdDchPoint.cxx
MpdDchHit.cxx 
MpdDchHitProducer.cxx 
)

# fill list of header files from list of source files
# by exchanging the file extension 
CHANGE_FILE_EXTENSION(*.cxx *.h Dch_HEADERS "${Dch_SRCS}")

set(Dch_LINKDEF MpdDchLinkDef.h)
set(Dch_DICTIONARY ${CMAKE_CURRENT_BINARY_DIR}/MpdDchDict.cxx) 

ROOT_GENERATE_DICTIONARY("${Dch_HEADERS}" "${Dch_LINKDEF}" "${Dch_DICTIONARY}" "${INCLUDE_DIRECTORIES}")


set(Dch_SRCS ${Dch_SRCS} ${Dch_DICTIONARY})


add_library(Dch SHARED ${Dch_SRCS})
target_link_libraries(Dch ${ROOT_LIBRARIES})
set_target_properties(Dch PROPERTIES VERSION 0.0.1 SOVERSION 0)

################ install ###################
install(TARGETS Dch DESTINATION ${CMAKE_BINARY_DIR}/lib)
