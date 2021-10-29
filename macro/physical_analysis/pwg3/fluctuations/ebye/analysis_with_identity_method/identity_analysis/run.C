void run()
{
    gROOT->LoadMacro("src/TIdentityFunctions.cc+g");
    gROOT->LoadMacro("src/TIdentity2D.cc");
    gROOT->ProcessLine(".x analyse_identity_dEdx.C");
}
