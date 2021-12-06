MACRO(InstallNicaFemto)
set(NICAFEMTO_SUBDIR_BUILD ON)
set(NICAFEMTO_HASH c8ee254fa293e5676f224e8d12f3c2c209093f8e) 
if(NOT EXISTS "${CMAKE_SOURCE_DIR}/external/nicafemto/CMakeLists.txt")
    execute_process(COMMAND git clone -b master https://git.jinr.ru/nica/nicafemto.git ${CMAKE_SOURCE_DIR}/external/nicafemto)
    execute_process(COMMAND git reset --hard  ${NICAFEMTO_HASH} WORKING_DIRECTORY  ${CMAKE_SOURCE_DIR}/external/nicafemto)
endif()
execute_process(COMMAND git fetch WORKING_DIRECTORY  ${CMAKE_SOURCE_DIR}/external/nicafemto)
    execute_process(COMMAND git reset --hard  ${NICAFEMTO_HASH} WORKING_DIRECTORY  ${CMAKE_SOURCE_DIR}/external/nicafemto)
        set(NICAFEMTO_CORE "${CMAKE_SOURCE_DIR}/external/nicafemto")
        if(NOT DEFINED NICAFEMTO_SUBDIR_BUILD)
            MESSAGE(FATAL_ERROR "SUBIDR NOT FOUND")
        endif()
        add_subdirectory (${CMAKE_SOURCE_DIR}/external/nicafemto)
        add_subdirectory (${CMAKE_SOURCE_DIR}/physics/pwg3/femtoscopy/nicafemto)
ENDMACRO(InstallNicaFemto)
