# Create a library called "libHADGEN"
add_subdirectory (parse)

include_directories(${CMAKE_CURRENT_BINARY_DIR} ${CMAKE_CURRENT_BINARY_DIR}/parse)

Set(LINK_DIRECTORIES
 ${ROOT_LIBRARY_DIR}
 ${FAIRROOT_LIBRARY_DIR}
)
 
Link_Directories(${LINK_DIRECTORIES})

# List of source files
set(SRCS
    hadgen_common_blocks.c
    hadgen_options.c
    hadgen_prev.c
    hadgen_tabnuc.c
    hadgen_time.c
    hadgen_setdky.c

    HistAGT.f
    Qgsmn1.f
    Qgsmn2.f
    casevp.f
    coales.f
    deex.f
    files.f
    genagt.f
    gqaar.f
    gqamn.f
    gqstn.f
    hadgen_so.f
    preco.f
    ranlux.f
)

# List of header files
set(HEADERS
    hadgen.h
    hadgen_common_blocks.h
)

Set(LIBRARY_NAME HADGEN)
Set(DEPENDENCIES Core Base)

GENERATE_LIBRARY()

add_dependencies(HADGEN parse_run)
