enable_language (C Fortran)

include_directories(inc ${TSHIELD_ROOT_PATH}/hadgen/inc)
set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -w")

Set(LINK_DIRECTORIES
 ${ROOT_LIBRARY_DIR}
 ${FAIRROOT_LIBRARY_DIR}
)
 
Link_Directories(${LINK_DIRECTORIES})

set(SRCS
    src/shield_decode_tree.c 
    src/shield_geo.c 
    src/shield_init.c 
    src/shield_iso.c 
    src/shield_media.c 
    src/shield_setget.c 
     
    src/Fragments.f
    src/Inspar.f
    src/Loent28.f
    src/Medium.f
    src/Microd.f
    src/Shield_so.f
    src/Tree.f
    src/decay1.f
    src/elscat.f
    src/gemca.f
    src/files.f
    src/gentree.f
    src/lndata28.f
    src/range.f
    src/Output.f
)

Set(LIBRARY_NAME SHIELD)
Set(DEPENDENCIES HADGEN)

GENERATE_LIBRARY()
