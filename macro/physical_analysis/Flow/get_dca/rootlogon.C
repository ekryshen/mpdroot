/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id: rootlogon.C,v 1.4 2007/05/04 11:18:51 ivana Exp $ */

/// By Laurent Aphecetche

{
  //cout << "Loading MPD libraries ..." << endl;
  //gROOT->LoadMacro("${ALICE_ROOT}/MUON/loadlibs.C");
  //gInterpreter->ProcessLine("loadlibs()");
    
  cout << "Setting include path ..." << endl;
  TString includePath = "-I${ROOT_INCLUDE_DIR} ";
  includePath        += "-I${Boost_INCLUDE_DIRS} ";
  includePath        += "-I${VMCWORKDIR}/core/mpdBase ";
  includePath        += "-I${VMCWORKDIR}/core/mpdPid ";
  includePath        += "-I${VMCWORKDIR}/detectors/etof ";
  includePath        += "-I${VMCWORKDIR}/detectors/tof ";
  includePath        += "-I${VMCWORKDIR}/detectors/tpc ";
  includePath        += "-I${VMCWORKDIR}/kalman ";
  includePath        += "-I${VMCWORKDIR}/lhetrack ";
  includePath        += "-I${VMCWORKDIR}/simulation/mcStack ";
  /*
  includePath        += "-I${ALICE_ROOT}/RAW ";
  includePath        += "-I${ALICE_ROOT}/FASTSIM ";
  includePath        += "-I${ALICE_ROOT}/EVGEN ";
  includePath        += "-I${ALICE_ROOT}/SHUTTLE/TestShuttle ";
  includePath        += "-I${ALICE_ROOT}/ITS ";
  includePath        += "-I${ALICE_ROOT}/MUON ";
  includePath        += "-I${ALICE_ROOT}/MUON/mapping";
  */
  gSystem->SetIncludePath(includePath.Data());
}
