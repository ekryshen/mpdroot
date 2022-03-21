set(PACKAGE THadgen) # this line comes after last add_subdirectory
# Create a library called "libTHadgen"
set(INCLUDE_DIRECTORIES
 ${CMAKE_CURRENT_SOURCE_DIR}
 ${TSHIELD_ROOT_PATH}/hadgen
)

Set(SYSTEM_INCLUDE_DIRECTORIES
 ${ROOT_INCLUDE_DIR}
)

include_directories(${INCLUDE_DIRECTORIES})
Include_Directories(SYSTEM ${SYSTEM_INCLUDE_DIRECTORIES})

link_directories(
 ${ROOT_LIBRARY_DIR}
 ${FAIRROOT_LIBRARY_DIR}
)

# List of source files
set(SRCS
    THadgen.cxx
    THadgenRoutines.cxx
)

set(HEADERS THadgen.h)

Set(LINKDEF ${PACKAGE}LinkDef.h)
Set(LIBRARY_NAME ${PACKAGE})
Set(DEPENDENCIES HADGEN Core Base)

GENERATE_LIBRARY()
unset(PACKAGE) # not really necessary, it will be overwritten by following CMakeLists.txt
