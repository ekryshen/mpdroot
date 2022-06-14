add_test(runMC ${ROOTSYS}/bin/root -b -q -l "${PROJECT_SOURCE_DIR}/macro/mpd/runMC.C(\"auau.09gev.mbias.98k.ftn14\",\"$VMCWORKDIR/macro/mpd/evetest.root\",0,2)")
SET_TESTS_PROPERTIES (runMC PROPERTIES TIMEOUT "200")
SET_TESTS_PROPERTIES (runMC PROPERTIES PASS_REGULAR_EXPRESSION "Macro finished successfully.")

add_test(reco ${ROOTSYS}/bin/root -b -q -l ${PROJECT_SOURCE_DIR}/macro/mpd/reco.C)
SET_TESTS_PROPERTIES (reco PROPERTIES TIMEOUT "200")
SET_TESTS_PROPERTIES (reco PROPERTIES PASS_REGULAR_EXPRESSION "Macro finished successfully.")
