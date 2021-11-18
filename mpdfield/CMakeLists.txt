# Create a library called "libMpdField" which includes the source files given in
# the array .
# The extension is already found.  Any number of sources could be listed here.

set(INCLUDE_DIRECTORIES
 ${BASE_INCLUDE_DIRECTORIES}
 ${CMAKE_SOURCE_DIR}/mpdfield 
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
MpdConstField.cxx       
MpdFieldContFact.cxx 
MpdFieldMap.cxx      
MpdFieldMapData.cxx  
MpdFieldPar.cxx      
MpdFieldCreator.cxx  
MpdRegion.cxx        
MpdMapPar.cxx        
MpdMultiFieldPar.cxx 
MpdConstPar.cxx    
MpdMultiField.cxx  
MpdFieldMapSym3.cxx  
)

Set(HEADERS)
Set(LINKDEF FieldLinkDef.h)
Set(LIBRARY_NAME MpdField)
Set(DEPENDENCIES)

GENERATE_LIBRARY()
