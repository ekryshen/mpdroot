enable_language (C Fortran)

link_directories(${ROOT_LIBRARY_DIR})

# List of source files
set(SRCS
    src/fparse.f
    src/main.c
)

add_executable(parse ${SRCS})

#Copy files for run generated executable.
configure_file(TABNUC.DAT TABNUC.DAT COPYONLY)
configure_file(ATAB.DAT   ATAB.DAT   COPYONLY)

#Run executable after build (always outdated target)
add_custom_target(parse_run DEPENDS parse TABNUC.DAT ATAB.DAT VERBATIM)
add_custom_command(
    TARGET parse_run
    POST_BUILD
    COMMAND parse TABNUC.DAT
    DEPENDS parse TABNUC.DAT ATAB.DAT
    COMMENT "" #Parsing TABNUC.DAT and ATAB.DAT for HADGEN
    VERBATIM
)
