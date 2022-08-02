void DecayConfig()
{
  // External decayer configuration file

  MpdDecayerPyt8* decayer = MpdDecayerPyt8::Instance();
  gMC->SetExternalDecayer(decayer);
  //gMC->SetUserDecay(3122); // user decay of Lambda-hyperon
  decayer->AddMotherPdg(+3122); // user decay of polarized Lambda-hyperon
  decayer->AddMotherPdg(-3122); // user decay of polarized antiLambda-hyperon
}
