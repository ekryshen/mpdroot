set(INCLUDE_DIRECTORIES
${BASE_INCLUDE_DIRECTORIES}
${CMAKE_CURRENT_LIST_DIR}
${CMAKE_SOURCE_DIR}/tools/database/uniCommon
${CMAKE_SOURCE_DIR}/tools/database/uniDb
${CMAKE_SOURCE_DIR}/decoder/Devices
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
MpdTDCDigit.cxx
MpdSyncDigit.cxx
MpdTTVXSDigit.cxx
MpdRawDataDecoder.cxx
MpdDecoderImpl.cxx
MpdEventHeader.cxx
MpdDevice.cxx
MpdFileSource.cxx
)

Set(HEADERS
MpdTDCDigit.h
MpdSyncDigit.h
MpdTTVXSDigit.h
MpdRawDataDecoder.h
MpdDecoderImpl.h
MpdEventHeader.h
TLVBlock.h
MpdDevice.h
MpdFileSource.h
)

Set(LINKDEF DecoderLinkDef.h)
Set(LIBRARY_NAME Decoder)
Set(DEPENDENCIES Core UniCommon UniDb)

GENERATE_LIBRARY()
