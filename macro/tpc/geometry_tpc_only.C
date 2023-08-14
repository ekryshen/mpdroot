void geometry_tpc_only(FairRunSim *fRun) {
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

  FairDetector *Tpc = new TpcDetector("TPC", kTRUE);
  Tpc->SetGeometryFileName("tpc_v10.root");
  fRun->AddModule(Tpc);
}
