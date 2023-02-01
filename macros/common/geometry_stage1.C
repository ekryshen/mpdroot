void geometry_stage1(FairRunSim *fRun) {
  // Set Material file Name
  fRun->SetMaterials("media.geo");

  // Create and add detectors
  //-------------------------
  FairModule *Cave = new FairCave("CAVE");
  Cave->SetGeometryFileName("cave.geo");
  fRun->AddModule(Cave);

  FairModule *Pipe = new FairPipe("PIPE");
  Pipe->SetGeometryFileName("pipe_v3.root");
  fRun->AddModule(Pipe);

  FairModule *Magnet = new FairMagnet("MAGNET");
  Magnet->SetGeometryFileName("magnet_v5.root");
  fRun->AddModule(Magnet);

  FairModule *Cradle = new FairMagnet("CRADLE");
  Cradle->SetGeometryFileName("cradle_v1.root");
  fRun->AddModule(Cradle);

  FairDetector *Ffd = new MpdFfd("FFD", kTRUE);
  Ffd->SetGeometryFileName("FFD_v8.root");
  fRun->AddModule(Ffd);

  FairDetector *Tpc = new TpcDetector("TPC", kTRUE);
  Tpc->SetGeometryFileName("tpc_v9.root");
  fRun->AddModule(Tpc);

  FairDetector *Tof = new MpdTof("TOF", kTRUE);
  Tof->SetGeometryFileName("tof_v8_3.root");
  fRun->AddModule(Tof);

  FairDetector *Emc = new MpdEmcKI("EMC", kTRUE);
  Emc->SetGeometryFileName("emc_v3.root");
  fRun->AddModule(Emc);

  FairDetector *Zdc = new MpdZdc("ZDC", kTRUE);
  Zdc->SetGeometryFileName(
      "zdc_oldnames_7sect_v1_no_overlaps_w_pipe_magnet.root");
  fRun->AddModule(Zdc);

  /* unsupported as of 12.2022
  FairDetector *mcord = new MpdMcord("MCORD",kTRUE);
  mcord->SetGeometryFileName("mcord_v3.root");
  fRun->AddModule(mcord);
  */
}
