#include <TFile.h>
#include <TH1D.h>
#include <TH1F.h>
#include <TH2D.h>
#include <TProfile.h>
#include <TList.h>
#include <TString.h>
#include <TMath.h>
#include "TSystem.h"
#include "TROOT.h"
#include "TFileMerger.h"

#include <iostream>

using namespace std;

void mergeScript(TString inDir, TString outDir)
{
    // load libraries
    gSystem->Load("libCore.so");
    gSystem->Load("libGeom.so");
    gSystem->Load("libVMC.so");
    gSystem->Load("libPhysics.so");

    TString string = ".! ls "+ inDir +" > files.txt";
    gROOT->ProcessLine(string);
    TFileMerger merger;
    
    int n;
    n=0;
    
    ifstream base("files.txt");
    char *str = new char [1000];
    while (!base.eof())
    {
        base.getline(str, 1000, '\n');
        cout << str << endl;
        if ( TString(str).Length() < 3 ) break;
        merger.AddFile(str);
        n++;
    }
    
    cout << "Start merging..." << endl;
    merger.OutputFile(outDir);
    merger.Merge();
    cout << "Merging finished." << endl;
}
