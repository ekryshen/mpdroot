#bin/bash
. ~/mpdroot/build/config.sh
valgrind \
--tool=memcheck \
--leak-check=full \
--log-file=val_runMC.log \
--suppressions=/opt/fairsoft/tools/root/etc/valgrind-root.supp \
root.exe -l -b -q $VMCWORKDIR/macro/mpd/runMC.C

