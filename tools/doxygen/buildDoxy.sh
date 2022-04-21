#!/bin/bash
echo "INPUT = \ " > inputDirectories.txt
find ../../ -type d -exec bash -O dotglob -c '  # search for directories containing source files
    for dirpath do
        ok=false
        seen_files=false
        set -- "$dirpath"/*
        for name do
            [ -d "$name" ] && continue  # skip dirs
            seen_files=true
            case "${name##*/}" in
                *.cpp|*.cxx|*.h) ok=true; break
            esac
        done

        if [[ $seen_files == false ]] # empty directory
        then
          continue
        fi

        if [[ $ok == false ]] # no source file in directory
        then
          continue
        fi

        if [[ "$dirpath" = *"/legacy"* ]] # skip legacy folders
        then
          continue
        fi

        if [[ "$dirpath" = "../../build"* ]] # skip build folder
        then
          continue
        fi

        if [[ "$dirpath" = "../../external"* ]] # skip external (nicafemto) folder
        then
          continue
        fi

        if [[ "$dirpath" = "../../macro"* ]] # skip directory with macros
        then
          continue
        fi

        if [[ "$dirpath" = "../../nicafemto"* ]] # we have no nicafemto support for now
        then
          continue
        fi

        echo "        $dirpath/ \ " >> inputDirectories.txt  # directory containing source files
    done' bash {} +
echo "" >> inputDirectories.txt
#doxygen MPDrootDoxy
#rm inputDirectories.txt
