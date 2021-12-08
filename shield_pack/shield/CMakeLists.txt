include_directories(${CMAKE_CURRENT_SOURCE_DIR} ${TSHIELD_ROOT_PATH}/hadgen)

Set(LINK_DIRECTORIES
 ${ROOT_LIBRARY_DIR}
 ${FAIRROOT_LIBRARY_DIR}
)
 
Link_Directories(${LINK_DIRECTORIES})

set(SRCS
    shield_decode_tree.c 
    shield_geo.c 
    shield_init.c 
    shield_iso.c 
    shield_media.c 
    shield_setget.c 
     
    Fragments.f
    Inspar.f
    Loent28.f
    Medium.f
    Microd.f
    Shield_so.f
    Tree.f
    decay1.f
    elscat.f
    gemca.f
    files.f
    gentree.f
    lndata28.f
    range.f
    Output.f
)

Set(LIBRARY_NAME SHIELD)
Set(DEPENDENCIES HADGEN)

GENERATE_LIBRARY()
