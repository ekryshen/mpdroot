// =========================================================================================================================================================================
// This code creates tree that contains sim and rec info about events. Then this tree will be used in the Unfolding procedure to build detector Response Matrix.

// Author: Daria Prokhorova (SPbSU), 10.2020
// =========================================================================================================================================================================

#include <Rtypes.h>
#include <TString.h>
#include <TChain.h>
#include <TFile.h>
#include <TH1.h>
#include <TH2.h>
#include <TVector3.h>
#include <TSystem.h>

#include "TTree.h"
#include "TFile.h"

#include "TRandom.h"
#include "TString.h"

#include "TLorentzVector.h"
#include "TDatabasePDG.h"
#include "TParticlePDG.h"
#include "TMath.h"
#include "TString.h"
#include "TCanvas.h"

#include <iomanip>
#include <iostream>
using std::cout;
using std::endl;
using namespace std;

#include "TRandom3.h"
#include "TH1D.h"
#include "TH2D.h"

#include <string>
#include <TStyle.h>
#include <TSystem.h>
#include <TCanvas.h>
#include <TFile.h>
#include <TProfile.h>
#include <TList.h>
#include <TString.h>
#include <TMath.h>
#include <TF1.h>
#include <TLegend.h>
#include <TROOT.h>
#include <TChain.h>
#include <vector>
#include <iostream>
#include <fstream>
#include <TNtuple.h>
#include <TBranch.h>
#include <TKey.h>
#include <TMarker.h>
#include <TGraph.h>
#include <TCanvas.h>
#include <TGraphErrors.h>
#include <TMultiGraph.h>
#include <TLegend.h>
#include <TAxis.h>
#include <TGaxis.h>
#include "TLatex.h"
#include "TImage.h"
#include <TLine.h>
#include <vector>

using namespace std;

R__ADD_INCLUDE_PATH($VMCWORKDIR)
#include "/*add_your_path*/mpdroot/macro/mpd/mpdloadlibs.C"

#include "/*add_your_path*/mpdroot/mpddst/MpdMiniEvent/MpdMiniDstReader.h"
#include "/*add_your_path*/mpdroot/mpddst/MpdMiniEvent/MpdMiniDst.h"
#include "/*add_your_path*/fairSoft/mpdroot/mpddst/MpdMiniEvent/MpdMiniEvent.h"
#include "/*add_your_path*/mpdroot/mpddst/MpdMiniEvent/MpdMiniTrack.h"
#include "/*add_your_path*/mpdroot/mpddst/MpdMiniEvent/MpdMiniMcEvent.h"
#include "/*add_your_path*/mpdroot/mpddst/MpdMiniEvent/MpdMiniMcTrack.h"

R__LOAD_LIBRARY(/*add_your_path*/mpdroot/mpddst/MpdMiniEvent/libMpdMiniDst)

int NumbOne_winds = 10; // number of pseudorapidity intervals for N, pT and PTN unfolding
int Numb_winds = 10; // number of pseudorapidity intervals for NFNB unfolding
const float mPion = 0.137;
Double_t wind_shift = 0.1; // step of pseudorapidity interval
Double_t wind_width = 0.1; // width of pseudorapidity interval
Double_t eta_max = 1.; // maximum of eta

void fill_miss(Double_t eta_sim, Double_t pt_sim, Double_t charge_sim,
               vector<Double_t> &tempAll0, vector<Double_t> &tempPlus0, vector<Double_t> &tempMinus0, vector<Double_t> &tempAll0_rec, vector<Double_t> &tempPlus0_rec, vector<Double_t> &tempMinus0_rec,
               vector<Double_t> &tempAll1, vector<Double_t> &tempPlus1, vector<Double_t> &tempMinus1, vector<Double_t> &tempAll1_rec, vector<Double_t> &tempPlus1_rec, vector<Double_t> &tempMinus1_rec,
               vector<Double_t> &tempAll2, vector<Double_t> &tempPlus2, vector<Double_t> &tempMinus2, vector<Double_t> &tempAll2_rec, vector<Double_t> &tempPlus2_rec, vector<Double_t> &tempMinus2_rec,
               vector<Double_t> &tempAll3, vector<Double_t> &tempPlus3, vector<Double_t> &tempMinus3, vector<Double_t> &tempAll3_rec, vector<Double_t> &tempPlus3_rec, vector<Double_t> &tempMinus3_rec,
               vector<Double_t> &tempAll4, vector<Double_t> &tempPlus4, vector<Double_t> &tempMinus4, vector<Double_t> &tempAll4_rec, vector<Double_t> &tempPlus4_rec, vector<Double_t> &tempMinus4_rec,
               vector<Double_t> &tempAll5, vector<Double_t> &tempPlus5, vector<Double_t> &tempMinus5, vector<Double_t> &tempAll5_rec, vector<Double_t> &tempPlus5_rec, vector<Double_t> &tempMinus5_rec,
               vector<Double_t> &tempAll6, vector<Double_t> &tempPlus6, vector<Double_t> &tempMinus6, vector<Double_t> &tempAll6_rec, vector<Double_t> &tempPlus6_rec, vector<Double_t> &tempMinus6_rec,
               vector<Double_t> &tempAll7, vector<Double_t> &tempPlus7, vector<Double_t> &tempMinus7, vector<Double_t> &tempAll7_rec, vector<Double_t> &tempPlus7_rec, vector<Double_t> &tempMinus7_rec,
               vector<Double_t> &tempAll8, vector<Double_t> &tempPlus8, vector<Double_t> &tempMinus8, vector<Double_t> &tempAll8_rec, vector<Double_t> &tempPlus8_rec, vector<Double_t> &tempMinus8_rec,
               vector<Double_t> &tempAll9, vector<Double_t> &tempPlus9, vector<Double_t> &tempMinus9, vector<Double_t> &tempAll9_rec, vector<Double_t> &tempPlus9_rec, vector<Double_t> &tempMinus9_rec)
{
    for(int i = 0; i<NumbOne_winds;i++)
    {
        if((eta_sim<(eta_max-i*wind_shift))&&(eta_sim>(-eta_max+i*wind_shift)))
        {
            if(i==0)
            {
                tempAll0.push_back(pt_sim);
                tempAll0_rec.push_back(-1);
                if(charge_sim>0)
                {
                    tempPlus0.push_back(pt_sim);
                    tempPlus0_rec.push_back(-1);
                }
                if(charge_sim<0)
                {
                    tempMinus0.push_back(pt_sim);
                    tempMinus0_rec.push_back(-1);
                }
            }
            if(i==1)
            {
                tempAll1.push_back(pt_sim);
                tempAll1_rec.push_back(-1);
                if(charge_sim>0)
                {
                    tempPlus1.push_back(pt_sim);
                    tempPlus1_rec.push_back(-1);
                }
                if(charge_sim<0)
                {
                    tempMinus1.push_back(pt_sim);
                    tempMinus1_rec.push_back(-1);
                }
            }
            if(i==2)
            {
                tempAll2.push_back(pt_sim);
                tempAll2_rec.push_back(-1);
                if(charge_sim>0)
                {
                    tempPlus2.push_back(pt_sim);
                    tempPlus2_rec.push_back(-1);
                }
                if(charge_sim<0)
                {
                    tempMinus2.push_back(pt_sim);
                    tempMinus2_rec.push_back(-1);
                }
            }
            if(i==3)
            {
                tempAll3.push_back(pt_sim);
                tempAll3_rec.push_back(-1);
                if(charge_sim>0)
                {
                    tempPlus3.push_back(pt_sim);
                    tempPlus3_rec.push_back(-1);
                }
                if(charge_sim<0)
                {
                    tempMinus3.push_back(pt_sim);
                    tempMinus3_rec.push_back(-1);
                }
            }
            if(i==4)
            {
                tempAll4.push_back(pt_sim);
                tempAll4_rec.push_back(-1);
                if(charge_sim>0)
                {
                    tempPlus4.push_back(pt_sim);
                    tempPlus4_rec.push_back(-1);
                }
                if(charge_sim<0)
                {
                    tempMinus4.push_back(pt_sim);
                    tempMinus4_rec.push_back(-1);
                }
            }
            if(i==5)
            {
                tempAll5.push_back(pt_sim);
                tempAll5_rec.push_back(-1);
                if(charge_sim>0)
                {
                    tempPlus5.push_back(pt_sim);
                    tempPlus5_rec.push_back(-1);
                }
                if(charge_sim<0)
                {
                    tempMinus5.push_back(pt_sim);
                    tempMinus5_rec.push_back(-1);
                }
            }
            if(i==6)
            {
                tempAll6.push_back(pt_sim);
                tempAll6_rec.push_back(-1);
                if(charge_sim>0)
                {
                    tempPlus6.push_back(pt_sim);
                    tempPlus6_rec.push_back(-1);
                }
                if(charge_sim<0)
                {
                    tempMinus6.push_back(pt_sim);
                    tempMinus6_rec.push_back(-1);
                }
            }
            if(i==7)
            {
                tempAll7.push_back(pt_sim);
                tempAll7_rec.push_back(-1);
                if(charge_sim>0)
                {
                    tempPlus7.push_back(pt_sim);
                    tempPlus7_rec.push_back(-1);
                }
                if(charge_sim<0)
                {
                    tempMinus7.push_back(pt_sim);
                    tempMinus7_rec.push_back(-1);
                }
            }
            if(i==8)
            {
                tempAll8.push_back(pt_sim);
                tempAll8_rec.push_back(-1);
                if(charge_sim>0)
                {
                    tempPlus8.push_back(pt_sim);
                    tempPlus8_rec.push_back(-1);
                }
                if(charge_sim<0)
                {
                    tempMinus8.push_back(pt_sim);
                    tempMinus8_rec.push_back(-1);
                }
            }
            if(i==9)
            {
                tempAll9.push_back(pt_sim);
                tempAll9_rec.push_back(-1);
                if(charge_sim>0)
                {
                    tempPlus9.push_back(pt_sim);
                    tempPlus9_rec.push_back(-1);
                }
                if(charge_sim<0)
                {
                    tempMinus9.push_back(pt_sim);
                    tempMinus9_rec.push_back(-1);
                }
            }
        }
    }
}

void fill_fake(Double_t eta_rec, Double_t pt_rec, Double_t charge_rec,
               vector<Double_t> &tempAll0, vector<Double_t> &tempPlus0, vector<Double_t> &tempMinus0, vector<Double_t> &tempAll0_rec, vector<Double_t> &tempPlus0_rec, vector<Double_t> &tempMinus0_rec,
               vector<Double_t> &tempAll1, vector<Double_t> &tempPlus1, vector<Double_t> &tempMinus1, vector<Double_t> &tempAll1_rec, vector<Double_t> &tempPlus1_rec, vector<Double_t> &tempMinus1_rec,
               vector<Double_t> &tempAll2, vector<Double_t> &tempPlus2, vector<Double_t> &tempMinus2, vector<Double_t> &tempAll2_rec, vector<Double_t> &tempPlus2_rec, vector<Double_t> &tempMinus2_rec,
               vector<Double_t> &tempAll3, vector<Double_t> &tempPlus3, vector<Double_t> &tempMinus3, vector<Double_t> &tempAll3_rec, vector<Double_t> &tempPlus3_rec, vector<Double_t> &tempMinus3_rec,
               vector<Double_t> &tempAll4, vector<Double_t> &tempPlus4, vector<Double_t> &tempMinus4, vector<Double_t> &tempAll4_rec, vector<Double_t> &tempPlus4_rec, vector<Double_t> &tempMinus4_rec,
               vector<Double_t> &tempAll5, vector<Double_t> &tempPlus5, vector<Double_t> &tempMinus5, vector<Double_t> &tempAll5_rec, vector<Double_t> &tempPlus5_rec, vector<Double_t> &tempMinus5_rec,
               vector<Double_t> &tempAll6, vector<Double_t> &tempPlus6, vector<Double_t> &tempMinus6, vector<Double_t> &tempAll6_rec, vector<Double_t> &tempPlus6_rec, vector<Double_t> &tempMinus6_rec,
               vector<Double_t> &tempAll7, vector<Double_t> &tempPlus7, vector<Double_t> &tempMinus7, vector<Double_t> &tempAll7_rec, vector<Double_t> &tempPlus7_rec, vector<Double_t> &tempMinus7_rec,
               vector<Double_t> &tempAll8, vector<Double_t> &tempPlus8, vector<Double_t> &tempMinus8, vector<Double_t> &tempAll8_rec, vector<Double_t> &tempPlus8_rec, vector<Double_t> &tempMinus8_rec,
               vector<Double_t> &tempAll9, vector<Double_t> &tempPlus9, vector<Double_t> &tempMinus9, vector<Double_t> &tempAll9_rec, vector<Double_t> &tempPlus9_rec, vector<Double_t> &tempMinus9_rec)
{
    for(int i = 0; i<NumbOne_winds;i++)
    {
        if((eta_rec<(eta_max-i*wind_shift))&&(eta_rec>(-eta_max+i*wind_shift)))
        {
            if(i==0)
            {
                tempAll0.push_back(-1);
                tempAll0_rec.push_back(pt_rec);
                if(charge_rec>0)
                {
                    tempPlus0.push_back(-1);
                    tempPlus0_rec.push_back(pt_rec);
                }
                if(charge_rec<0)
                {
                    tempMinus0.push_back(-1);
                    tempMinus0_rec.push_back(pt_rec);
                }
            }
            if(i==1)
            {
                tempAll1.push_back(-1);
                tempAll1_rec.push_back(pt_rec);
                if(charge_rec>0)
                {
                    tempPlus1.push_back(-1);
                    tempPlus1_rec.push_back(pt_rec);
                }
                if(charge_rec<0)
                {
                    tempMinus1.push_back(-1);
                    tempMinus1_rec.push_back(pt_rec);
                }
            }
            if(i==2)
            {
                tempAll2.push_back(-1);
                tempAll2_rec.push_back(pt_rec);
                if(charge_rec>0)
                {
                    tempPlus2.push_back(-1);
                    tempPlus2_rec.push_back(pt_rec);
                }
                if(charge_rec<0)
                {
                    tempMinus2.push_back(-1);
                    tempMinus2_rec.push_back(pt_rec);
                }
            }
            if(i==3)
            {
                tempAll3.push_back(-1);
                tempAll3_rec.push_back(pt_rec);
                if(charge_rec>0)
                {
                    tempPlus3.push_back(-1);
                    tempPlus3_rec.push_back(pt_rec);
                }
                if(charge_rec<0)
                {
                    tempMinus3.push_back(-1);
                    tempMinus3_rec.push_back(pt_rec);
                }
            }
            if(i==4)
            {
                tempAll4.push_back(-1);
                tempAll4_rec.push_back(pt_rec);
                if(charge_rec>0)
                {
                    tempPlus4.push_back(-1);
                    tempPlus4_rec.push_back(pt_rec);
                }
                if(charge_rec<0)
                {
                    tempMinus4.push_back(-1);
                    tempMinus4_rec.push_back(pt_rec);
                }
            }
            if(i==5)
            {
                tempAll5.push_back(-1);
                tempAll5_rec.push_back(pt_rec);
                if(charge_rec>0)
                {
                    tempPlus5.push_back(-1);
                    tempPlus5_rec.push_back(pt_rec);
                }
                if(charge_rec<0)
                {
                    tempMinus5.push_back(-1);
                    tempMinus5_rec.push_back(pt_rec);
                }
            }
            if(i==6)
            {
                tempAll6.push_back(-1);
                tempAll6_rec.push_back(pt_rec);
                if(charge_rec>0)
                {
                    tempPlus6.push_back(-1);
                    tempPlus6_rec.push_back(pt_rec);
                }
                if(charge_rec<0)
                {
                    tempMinus6.push_back(-1);
                    tempMinus6_rec.push_back(pt_rec);
                }
            }
            if(i==7)
            {
                tempAll7.push_back(-1);
                tempAll7_rec.push_back(pt_rec);
                if(charge_rec>0)
                {
                    tempPlus7.push_back(-1);
                    tempPlus7_rec.push_back(pt_rec);
                }
                if(charge_rec<0)
                {
                    tempMinus7.push_back(-1);
                    tempMinus7_rec.push_back(pt_rec);
                }
            }
            if(i==8)
            {
                tempAll8.push_back(-1);
                tempAll8_rec.push_back(pt_rec);
                if(charge_rec>0)
                {
                    tempPlus8.push_back(-1);
                    tempPlus8_rec.push_back(pt_rec);
                }
                if(charge_rec<0)
                {
                    tempMinus8.push_back(-1);
                    tempMinus8_rec.push_back(pt_rec);
                }
            }
            if(i==9)
            {
                tempAll9.push_back(-1);
                tempAll9_rec.push_back(pt_rec);
                if(charge_rec>0)
                {
                    tempPlus9.push_back(-1);
                    tempPlus9_rec.push_back(pt_rec);
                }
                if(charge_rec<0)
                {
                    tempMinus9.push_back(-1);
                    tempMinus9_rec.push_back(pt_rec);
                }
            }
        }
    }
}

void fill_RM(Double_t eta_sim, Double_t pt_sim, Double_t charge_sim, Double_t eta_rec, Double_t pt_rec, Double_t charge_rec,
             vector<Double_t> &tempAll0, vector<Double_t> &tempPlus0, vector<Double_t> &tempMinus0, vector<Double_t> &tempAll0_rec, vector<Double_t> &tempPlus0_rec, vector<Double_t> &tempMinus0_rec,
             vector<Double_t> &tempAll1, vector<Double_t> &tempPlus1, vector<Double_t> &tempMinus1, vector<Double_t> &tempAll1_rec, vector<Double_t> &tempPlus1_rec, vector<Double_t> &tempMinus1_rec,
             vector<Double_t> &tempAll2, vector<Double_t> &tempPlus2, vector<Double_t> &tempMinus2, vector<Double_t> &tempAll2_rec, vector<Double_t> &tempPlus2_rec, vector<Double_t> &tempMinus2_rec,
             vector<Double_t> &tempAll3, vector<Double_t> &tempPlus3, vector<Double_t> &tempMinus3, vector<Double_t> &tempAll3_rec, vector<Double_t> &tempPlus3_rec, vector<Double_t> &tempMinus3_rec,
             vector<Double_t> &tempAll4, vector<Double_t> &tempPlus4, vector<Double_t> &tempMinus4, vector<Double_t> &tempAll4_rec, vector<Double_t> &tempPlus4_rec, vector<Double_t> &tempMinus4_rec,
             vector<Double_t> &tempAll5, vector<Double_t> &tempPlus5, vector<Double_t> &tempMinus5, vector<Double_t> &tempAll5_rec, vector<Double_t> &tempPlus5_rec, vector<Double_t> &tempMinus5_rec,
             vector<Double_t> &tempAll6, vector<Double_t> &tempPlus6, vector<Double_t> &tempMinus6, vector<Double_t> &tempAll6_rec, vector<Double_t> &tempPlus6_rec, vector<Double_t> &tempMinus6_rec,
             vector<Double_t> &tempAll7, vector<Double_t> &tempPlus7, vector<Double_t> &tempMinus7, vector<Double_t> &tempAll7_rec, vector<Double_t> &tempPlus7_rec, vector<Double_t> &tempMinus7_rec,
             vector<Double_t> &tempAll8, vector<Double_t> &tempPlus8, vector<Double_t> &tempMinus8, vector<Double_t> &tempAll8_rec, vector<Double_t> &tempPlus8_rec, vector<Double_t> &tempMinus8_rec,
             vector<Double_t> &tempAll9, vector<Double_t> &tempPlus9, vector<Double_t> &tempMinus9, vector<Double_t> &tempAll9_rec, vector<Double_t> &tempPlus9_rec, vector<Double_t> &tempMinus9_rec)
{
    for(int i = 0; i<NumbOne_winds;i++)
    {
        if((eta_sim<(eta_max-i*wind_shift))&&(eta_sim>(-eta_max+i*wind_shift)))
        {
            if(i==0)
            {
                tempAll0.push_back(pt_sim);
                tempAll0_rec.push_back(pt_rec);
                if(charge_sim>0)
                {
                    tempPlus0.push_back(pt_sim);
                    tempPlus0_rec.push_back(pt_rec);
                }
                if(charge_sim<0)
                {
                    tempMinus0.push_back(pt_sim);
                    tempMinus0_rec.push_back(pt_rec);
                }
            }
            if(i==1)
            {
                tempAll1.push_back(pt_sim);
                tempAll1_rec.push_back(pt_rec);
                if(charge_sim>0)
                {
                    tempPlus1.push_back(pt_sim);
                    tempPlus1_rec.push_back(pt_rec);
                }
                if(charge_sim<0)
                {
                    tempMinus1.push_back(pt_sim);
                    tempMinus1_rec.push_back(pt_rec);
                }
            }
            if(i==2)
            {
                tempAll2.push_back(pt_sim);
                tempAll2_rec.push_back(pt_rec);
                if(charge_sim>0)
                {
                    tempPlus2.push_back(pt_sim);
                    tempPlus2_rec.push_back(pt_rec);
                }
                if(charge_sim<0)
                {
                    tempMinus2.push_back(pt_sim);
                    tempMinus2_rec.push_back(pt_rec);
                }
            }
            if(i==3)
            {
                tempAll3.push_back(pt_sim);
                tempAll3_rec.push_back(pt_rec);
                if(charge_sim>0)
                {
                    tempPlus3.push_back(pt_sim);
                    tempPlus3_rec.push_back(pt_rec);
                }
                if(charge_sim<0)
                {
                    tempMinus3.push_back(pt_sim);
                    tempMinus3_rec.push_back(pt_rec);
                }
            }
            if(i==4)
            {
                tempAll4.push_back(pt_sim);
                tempAll4_rec.push_back(pt_rec);
                if(charge_sim>0)
                {
                    tempPlus4.push_back(pt_sim);
                    tempPlus4_rec.push_back(pt_rec);
                }
                if(charge_sim<0)
                {
                    tempMinus4.push_back(pt_sim);
                    tempMinus4_rec.push_back(pt_rec);
                }
            }
            if(i==5)
            {
                tempAll5.push_back(pt_sim);
                tempAll5_rec.push_back(pt_rec);
                if(charge_sim>0)
                {
                    tempPlus5.push_back(pt_sim);
                    tempPlus5_rec.push_back(pt_rec);
                }
                if(charge_sim<0)
                {
                    tempMinus5.push_back(pt_sim);
                    tempMinus5_rec.push_back(pt_rec);
                }
            }
            if(i==6)
            {
                tempAll6.push_back(pt_sim);
                tempAll6_rec.push_back(pt_rec);
                if(charge_sim>0)
                {
                    tempPlus6.push_back(pt_sim);
                    tempPlus6_rec.push_back(pt_rec);
                }
                if(charge_sim<0)
                {
                    tempMinus6.push_back(pt_sim);
                    tempMinus6_rec.push_back(pt_rec);
                }
            }
            if(i==7)
            {
                tempAll7.push_back(pt_sim);
                tempAll7_rec.push_back(pt_rec);
                if(charge_sim>0)
                {
                    tempPlus7.push_back(pt_sim);
                    tempPlus7_rec.push_back(pt_rec);
                }
                if(charge_sim<0)
                {
                    tempMinus7.push_back(pt_sim);
                    tempMinus7_rec.push_back(pt_rec);
                }
            }
            if(i==8)
            {
                tempAll8.push_back(pt_sim);
                tempAll8_rec.push_back(pt_rec);
                if(charge_sim>0)
                {
                    tempPlus8.push_back(pt_sim);
                    tempPlus8_rec.push_back(pt_rec);
                }
                if(charge_sim<0)
                {
                    tempMinus8.push_back(pt_sim);
                    tempMinus8_rec.push_back(pt_rec);
                }
            }
            if(i==9)
            {
                tempAll9.push_back(pt_sim);
                tempAll9_rec.push_back(pt_rec);
                if(charge_sim>0)
                {
                    tempPlus9.push_back(pt_sim);
                    tempPlus9_rec.push_back(pt_rec);
                }
                if(charge_sim<0)
                {
                    tempMinus9.push_back(pt_sim);
                    tempMinus9_rec.push_back(pt_rec);
                }
            }
        }
    }
}

int main(int argc, char* argv[])
{
    
    MpdMiniDstReader* miniDstReader = new MpdMiniDstReader(argv[1]);
    miniDstReader->Init();
    
    miniDstReader->SetStatus("*",0);
    miniDstReader->SetStatus("Event*",1);
    miniDstReader->SetStatus("Track*",1);
    miniDstReader->SetStatus("McEvent*",1);
    miniDstReader->SetStatus("McTrack*",1);
    
    std::cout << "Status has been set" << std::endl;
    
    if (!miniDstReader->chain()) {
        std::cout << "No chain has been found." << std::endl;
    }
    
    Long64_t eventsInTree = miniDstReader->tree()->GetEntries();
    std::cout << "eventsInTree: " << eventsInTree << std::endl;
    Long64_t events2read = miniDstReader->chain()->GetEntries();
    std::cout << "Number of events to read: " << events2read
    << std::endl;
    
    TFile *oFile = new TFile("out.root", "RECREATE");
    
    TTree *N_Pt_Tree;
    N_Pt_Tree = new TTree("N_Pt_Tree", "An example of a ROOT tree");
    
    TBranch *BranchMultF_rec;                        TBranch *BranchMultF_sim;
    TBranch *BranchMultB_rec;                        TBranch *BranchMultB_sim;
    TBranch *BranchMultPlusF_rec;                    TBranch *BranchMultPlusF_sim;
    TBranch *BranchMultPlusB_rec;                    TBranch *BranchMultPlusB_sim;
    TBranch *BranchMultMinusF_rec;                   TBranch *BranchMultMinusF_sim;
    TBranch *BranchMultMinusB_rec;                   TBranch *BranchMultMinusB_sim;
    TBranch *BranchMultOne_rec;                      TBranch *BranchMultOne_sim;
    TBranch *BranchMultOnePlus_rec;                  TBranch *BranchMultOnePlus_sim;
    TBranch *BranchMultOneMinus_rec;                 TBranch *BranchMultOneMinus_sim;
    
    TBranch *BranchPTOne_rec;                        TBranch *BranchPTOne_sim;
    TBranch *BranchPT2One_rec;                       TBranch *BranchPT2One_sim;
    TBranch *BranchPTOnePlus_rec;                    TBranch *BranchPTOnePlus_sim;
    TBranch *BranchPT2OnePlus_rec;                   TBranch *BranchPT2OnePlus_sim;
    TBranch *BranchPTOneMinus_rec;                   TBranch *BranchPTOneMinus_sim;
    TBranch *BranchPT2OneMinus_rec;                  TBranch *BranchPT2OneMinus_sim;
    
    TBranch *BranchptOne_rec;                       TBranch *BranchptOne_sim;
    TBranch *BranchptOnePlus_rec;                   TBranch *BranchptOnePlus_sim;
    TBranch *BranchptOneMinus_rec;                  TBranch *BranchptOneMinus_sim;
    
    TBranch *BranchNumbOne_winds;
    TBranch *BranchNumb_winds;
    
    int multOne_rec[NumbOne_winds],multOnePlus_rec[NumbOne_winds], multOneMinus_rec[NumbOne_winds], multF_rec[Numb_winds],multB_rec[Numb_winds], multPlusF_rec[Numb_winds],multPlusB_rec[Numb_winds], multMinusF_rec[Numb_winds], multMinusB_rec[Numb_winds];
    
    int multOne_sim[NumbOne_winds],multOnePlus_sim[NumbOne_winds], multOneMinus_sim[NumbOne_winds], multF_sim[Numb_winds],multB_sim[Numb_winds], multPlusF_sim[Numb_winds],multPlusB_sim[Numb_winds], multMinusF_sim[Numb_winds], multMinusB_sim[Numb_winds];
    
    Double_t PTOne_rec[NumbOne_winds], PT2One_rec[NumbOne_winds],  PTOnePlus_rec[NumbOne_winds],PT2OnePlus_rec[NumbOne_winds],PTOneMinus_rec[NumbOne_winds], PT2OneMinus_rec[NumbOne_winds], PTOne_sim[NumbOne_winds], PT2One_sim[NumbOne_winds], PTOnePlus_sim[NumbOne_winds], PT2OnePlus_sim[NumbOne_winds], PTOneMinus_sim[NumbOne_winds], PT2OneMinus_sim[NumbOne_winds];
    
    BranchNumbOne_winds = N_Pt_Tree->Branch("NumbOne_winds", &NumbOne_winds, "NumbOne_winds/I");
    BranchNumb_winds = N_Pt_Tree->Branch("Numb_winds", &Numb_winds, "Numb_winds/I");
    
    BranchMultF_rec = N_Pt_Tree->Branch("multF_rec", multF_rec, "multF_rec[Numb_winds]/I");
    BranchMultB_rec = N_Pt_Tree->Branch("multB_rec", multB_rec, "multB_rec[Numb_winds]/I");
    BranchMultPlusF_rec = N_Pt_Tree->Branch("multPlusF_rec", multPlusF_rec, "multPlusF_rec[Numb_winds]/I");
    BranchMultMinusF_rec = N_Pt_Tree->Branch("multMinusF_rec", multMinusF_rec, "multMinusF_rec[Numb_winds]/I");
    BranchMultPlusB_rec = N_Pt_Tree->Branch("multPlusB_rec", multPlusB_rec, "multPlusB_rec[Numb_winds]/I");
    BranchMultMinusB_rec = N_Pt_Tree->Branch("multMinusB_rec", multMinusB_rec, "multMinusB_rec[Numb_winds]/I");
    BranchMultOne_rec = N_Pt_Tree->Branch("multOne_rec", multOne_rec, "multOne_rec[NumbOne_winds]/I");
    BranchMultOnePlus_rec = N_Pt_Tree->Branch("multOnePlus_rec", multOnePlus_rec, "multOnePlus_rec[NumbOne_winds]/I");
    BranchMultOneMinus_rec = N_Pt_Tree->Branch("multOneMinus_rec", multOneMinus_rec, "multOneMinus_rec[NumbOne_winds]/I");
    
    BranchMultF_sim = N_Pt_Tree->Branch("multF_sim", multF_sim, "multF_sim[Numb_winds]/I");
    BranchMultB_sim = N_Pt_Tree->Branch("multB_sim", multB_sim, "multB_sim[Numb_winds]/I");
    BranchMultPlusF_sim = N_Pt_Tree->Branch("multPlusF_sim", multPlusF_sim, "multPlusF_sim[Numb_winds]/I");
    BranchMultMinusF_sim = N_Pt_Tree->Branch("multMinusF_sim", multMinusF_sim, "multMinusF_sim[Numb_winds]/I");
    BranchMultPlusB_sim = N_Pt_Tree->Branch("multPlusB_sim", multPlusB_sim, "multPlusB_sim[Numb_winds]/I");
    BranchMultMinusB_sim = N_Pt_Tree->Branch("multMinusB_sim", multMinusB_sim, "multMinusB_sim[Numb_winds]/I");
    BranchMultOne_sim = N_Pt_Tree->Branch("multOne_sim", multOne_sim, "multOne_sim[NumbOne_winds]/I");
    BranchMultOnePlus_sim = N_Pt_Tree->Branch("multOnePlus_sim", multOnePlus_sim, "multOnePlus_sim[NumbOne_winds]/I");
    BranchMultOneMinus_sim = N_Pt_Tree->Branch("multOneMinus_sim", multOneMinus_sim, "multOneMinus_sim[NumbOne_winds]/I");
    
    BranchPTOne_rec = N_Pt_Tree->Branch("PTOne_rec", PTOne_rec, "PTOne_rec[NumbOne_winds]/D");
    BranchPT2One_rec = N_Pt_Tree->Branch("PT2One_rec", PT2One_rec, "PT2One_rec[NumbOne_winds]/D");
    BranchPTOnePlus_rec = N_Pt_Tree->Branch("PTOnePlus_rec", PTOnePlus_rec, "PTOnePlus_rec[NumbOne_winds]/D");
    BranchPT2OnePlus_rec = N_Pt_Tree->Branch("PT2OnePlus_rec", PT2OnePlus_rec, "PT2OnePlus_rec[NumbOne_winds]/D");
    BranchPTOneMinus_rec = N_Pt_Tree->Branch("PTOneMinus_rec", PTOneMinus_rec, "PTOneMinus_rec[NumbOne_winds]/D");
    BranchPT2OneMinus_rec = N_Pt_Tree->Branch("PT2OneMinus_rec", PT2OneMinus_rec, "PT2OneMinus_rec[NumbOne_winds]/D");
    
    BranchPTOne_sim = N_Pt_Tree->Branch("PTOne_sim", PTOne_sim, "PTOne_sim[NumbOne_winds]/D");
    BranchPT2One_sim = N_Pt_Tree->Branch("PT2One_sim", PT2One_sim, "PT2One_sim[NumbOne_winds]/D");
    BranchPTOnePlus_sim = N_Pt_Tree->Branch("PTOnePlus_sim", PTOnePlus_sim, "PTOnePlus_sim[NumbOne_winds]/D");
    BranchPT2OnePlus_sim = N_Pt_Tree->Branch("PT2OnePlus_sim", PT2OnePlus_sim, "PT2OnePlus_sim[NumbOne_winds]/D");
    BranchPTOneMinus_sim = N_Pt_Tree->Branch("PTOneMinus_sim", PTOneMinus_sim, "PTOneMinus_sim[NumbOne_winds]/D");
    BranchPT2OneMinus_sim = N_Pt_Tree->Branch("PT2OneMinus_sim", PT2OneMinus_sim, "PT2OneMinus_sim[NumbOne_winds]/D");
    
    vector<vector<Double_t> > particles_pt_in_event_sim;
    vector<vector<Double_t> > particles_ptPlus_in_event_sim;
    vector<vector<Double_t> > particles_ptMinus_in_event_sim;
    
    vector<vector<Double_t> > particles_pt_in_event_rec;
    vector<vector<Double_t> > particles_ptPlus_in_event_rec;
    vector<vector<Double_t> > particles_ptMinus_in_event_rec;
    
    BranchptOne_sim = N_Pt_Tree->Branch("ptOne_sim", &particles_pt_in_event_sim);
    BranchptOnePlus_sim = N_Pt_Tree->Branch("ptOnePlus_sim", &particles_ptPlus_in_event_sim);
    BranchptOneMinus_sim = N_Pt_Tree->Branch("ptOneMinus_sim", &particles_ptMinus_in_event_sim);
    
    BranchptOne_rec = N_Pt_Tree->Branch("ptOne_rec", &particles_pt_in_event_rec);
    BranchptOnePlus_rec = N_Pt_Tree->Branch("ptOnePlus_rec", &particles_ptPlus_in_event_rec);
    BranchptOneMinus_rec = N_Pt_Tree->Branch("ptOneMinus_rec", &particles_ptMinus_in_event_rec);
    
    Double_t p_sim, px_sim, py_sim, pz_sim, pt_sim, phi_sim, eta_sim, E_sim, y_sim, teta_sim, charge_sim;
    Double_t p_rec, px_rec, py_rec, pz_rec, pt_rec, phi_rec, eta_rec, E_rec, y_rec, teta_rec, charge_rec;
    
    // loop over events
    
    for (Long64_t iEvent = 0; iEvent < events2read; iEvent++)
    {
        vector<double> goodVertexTracksPt;
        vector<double> goodVertexTracksEta;
        vector<double> goodVertexTracksPhi;
        vector<double> goodVertexTracksPx;
        vector<double> goodVertexTracksPy;
        vector<double> goodVertexTracksPz;
        vector<double> goodVertexTracksCharge;
        vector<double> goodTracksIndex;
        vector<MpdMiniTrack*> goodRecTracks;
        
        vector<Double_t>tempAll0; vector<Double_t>tempPlus0; vector<Double_t>tempMinus0; vector<Double_t>tempAll0_rec; vector<Double_t>tempPlus0_rec; vector<Double_t>tempMinus0_rec;
        vector<Double_t>tempAll1; vector<Double_t>tempPlus1; vector<Double_t>tempMinus1; vector<Double_t>tempAll1_rec; vector<Double_t>tempPlus1_rec; vector<Double_t>tempMinus1_rec;
        vector<Double_t>tempAll2; vector<Double_t>tempPlus2; vector<Double_t>tempMinus2; vector<Double_t>tempAll2_rec; vector<Double_t>tempPlus2_rec; vector<Double_t>tempMinus2_rec;
        vector<Double_t>tempAll3; vector<Double_t>tempPlus3; vector<Double_t>tempMinus3; vector<Double_t>tempAll3_rec; vector<Double_t>tempPlus3_rec; vector<Double_t>tempMinus3_rec;
        vector<Double_t>tempAll4; vector<Double_t>tempPlus4; vector<Double_t>tempMinus4; vector<Double_t>tempAll4_rec; vector<Double_t>tempPlus4_rec; vector<Double_t>tempMinus4_rec;
        vector<Double_t>tempAll5; vector<Double_t>tempPlus5; vector<Double_t>tempMinus5; vector<Double_t>tempAll5_rec; vector<Double_t>tempPlus5_rec; vector<Double_t>tempMinus5_rec;
        vector<Double_t>tempAll6; vector<Double_t>tempPlus6; vector<Double_t>tempMinus6; vector<Double_t>tempAll6_rec; vector<Double_t>tempPlus6_rec; vector<Double_t>tempMinus6_rec;
        vector<Double_t>tempAll7; vector<Double_t>tempPlus7; vector<Double_t>tempMinus7; vector<Double_t>tempAll7_rec; vector<Double_t>tempPlus7_rec; vector<Double_t>tempMinus7_rec;
        vector<Double_t>tempAll8; vector<Double_t>tempPlus8; vector<Double_t>tempMinus8; vector<Double_t>tempAll8_rec; vector<Double_t>tempPlus8_rec; vector<Double_t>tempMinus8_rec;
        vector<Double_t>tempAll9; vector<Double_t>tempPlus9; vector<Double_t>tempMinus9; vector<Double_t>tempAll9_rec; vector<Double_t>tempPlus9_rec; vector<Double_t>tempMinus9_rec;
        
        Bool_t readEvent = miniDstReader->readMiniEvent(iEvent);
        if (!readEvent) {
            std::cout << "Something went wrong! Nothing to analyze..." << std::endl;
            break;
        }
        
        MpdMiniDst *dst = miniDstReader->miniDst();
        MpdMiniMcEvent *mcevent = dst->mcEvent();
        
        // event selection
        if(mcevent->impactParameter()>4.59567)
        {
            continue;
        }
        
        MpdMiniEvent *event = dst->event();
        if (!event)
            continue;
        
        TVector3 pVtx = event->primaryVertex();
        double pVtxX = pVtx.X();
        double pVtxY = pVtx.Y();
        double pVtxZ = pVtx.Z();
        
        for(int i = 0; i<NumbOne_winds;i++)
        {
            multOne_rec[i]=0; multOnePlus_rec[i]=0; multOneMinus_rec[i]=0;
            PTOne_rec[i]=0; PT2One_rec[i]=0; PTOnePlus_rec[i]=0; PT2OnePlus_rec[i]=0; PTOneMinus_rec[i]=0; PT2OneMinus_rec[i]=0;
        }
        for(int i = 0; i<Numb_winds;i++)
        {
            multF_rec[i]=0; multB_rec[i]=0; multPlusF_rec[i]=0; multPlusB_rec[i]=0; multMinusF_rec[i]=0; multMinusB_rec[i]=0;
        }
        
        for(int i = 0; i<NumbOne_winds;i++)
        {
            multOne_sim[i]=0;multOnePlus_sim[i]=0; multOneMinus_sim[i]=0;
            PTOne_sim[i]=0; PT2One_sim[i]=0; PTOnePlus_sim[i]=0; PT2OnePlus_sim[i]=0; PTOneMinus_sim[i]=0; PT2OneMinus_sim[i]=0;
        }
        for(int i = 0; i<Numb_winds;i++)
        {
            multF_sim[i]=0;multB_sim[i]=0;multPlusF_sim[i]=0;multPlusB_sim[i]=0;multMinusF_sim[i]=0; multMinusB_sim[i]=0;
        }
        
        // event cuts
        
        if(((TMath::Abs(event->primaryVertex().Z())<20))&&(TMath::Abs(event->primaryVertex().Z())>0.01)&&(dst->numberOfTracks()>2))
        {
            // loop over rec tracks

            for (Int_t iTrk = 0; iTrk < dst->numberOfTracks(); iTrk++) {
                
                MpdMiniTrack *miniTrack = dst->track(iTrk);
                
                if (!miniTrack) continue;
                
                bool primary=miniTrack->isPrimary();
                if (!primary)
                {
                    continue;
                }
                
                if (miniTrack->pMom().Perp() > 2 || miniTrack->pMom().Perp() < 0.15)
                {
                    continue;
                }
                
                if (miniTrack->gDCAx(pVtxX)*miniTrack->gDCAx(pVtxX) + miniTrack->gDCAy(pVtxY)*miniTrack->gDCAy(pVtxY)+miniTrack->gDCAz(pVtxZ)*miniTrack->gDCAz(pVtxZ)>5)
                {
                    continue;
                }
                
                if(miniTrack->nHits()<30)
                {
                    continue;
                }
                
                if((miniTrack->chi2()/miniTrack->nHits())>=5)
                {
                    continue;
                }
                
                p_rec = 0; px_rec=0; py_rec=0; pz_rec=0; pt_rec=0; phi_rec=0; eta_rec=0; E_rec=0; y_rec=0; teta_rec=0; charge_rec=0;
                
                p_rec = miniTrack->pMom().Mag();
                px_rec =  miniTrack->pMom().Px();
                py_rec =  miniTrack->pMom().Py();
                pz_rec =  miniTrack->pMom().Pz();
                pt_rec = miniTrack->pMom().Perp();
                phi_rec = TMath::ATan2(py_rec,px_rec);
                eta_rec = 0.5 * log((sqrt(px_rec*px_rec + py_rec*py_rec + pz_rec*pz_rec) + pz_rec) / (sqrt(px_rec*px_rec + py_rec*py_rec + pz_rec*pz_rec) - pz_rec));
                E_rec = sqrt(pow(mPion, 2) + pow(px_rec, 2) + pow(py_rec, 2) + pow(pz_rec, 2));
                y_rec = 0.5* log((E_rec + pz_rec) / (E_rec - pz_rec));
                teta_rec = 2 * atan(exp(-eta_rec));
                charge_rec = miniTrack->charge();
                
                if((eta_rec>1)||(eta_rec<-1))
                {
                    continue;
                }
                
                if(miniTrack->hasMcTrack())
                {
                    goodVertexTracksPt.push_back(pt_rec);
                    goodVertexTracksEta.push_back(eta_rec);
                    goodVertexTracksPhi.push_back(phi_rec);
                    goodVertexTracksPx.push_back(px_rec);
                    goodVertexTracksPy.push_back(py_rec);
                    goodVertexTracksPz.push_back(pz_rec);
                    goodVertexTracksCharge.push_back(charge_rec);
                    goodTracksIndex.push_back(iTrk);
                    goodRecTracks.push_back(miniTrack);
                }
                else
                {
                    // we filled fake rec tracks
                    fill_fake(eta_rec,pt_rec,charge_rec,
                              tempAll0, tempPlus0, tempMinus0, tempAll0_rec, tempPlus0_rec, tempMinus0_rec,
                              tempAll1, tempPlus1, tempMinus1, tempAll1_rec, tempPlus1_rec, tempMinus1_rec,
                              tempAll2, tempPlus2, tempMinus2, tempAll2_rec, tempPlus2_rec, tempMinus2_rec,
                              tempAll3, tempPlus3, tempMinus3, tempAll3_rec, tempPlus3_rec, tempMinus3_rec,
                              tempAll4, tempPlus4, tempMinus4, tempAll4_rec, tempPlus4_rec, tempMinus4_rec,
                              tempAll5, tempPlus5, tempMinus5, tempAll5_rec, tempPlus5_rec, tempMinus5_rec,
                              tempAll6, tempPlus6, tempMinus6, tempAll6_rec, tempPlus6_rec, tempMinus6_rec,
                              tempAll7, tempPlus7, tempMinus7, tempAll7_rec, tempPlus7_rec, tempMinus7_rec,
                              tempAll8, tempPlus8, tempMinus8, tempAll8_rec, tempPlus8_rec, tempMinus8_rec,
                              tempAll9, tempPlus9, tempMinus9, tempAll9_rec, tempPlus9_rec, tempMinus9_rec);
                }
                
                //            for(int i = 0; i<Numb_winds;i++)
                //            {
                //
                //                if((eta_rec<(eta_max-i*wind_shift))&&(eta_rec>eta_max-wind_width-i*wind_shift))
                //                {
                //                    multF_rec[i]++;
                //                    if(charge_rec>0){multPlusF_rec[i]++;}
                //                    if(charge_rec<0){multMinusF_rec[i]++;}
                //                }
                //
                //                if((eta_rec>(-eta_max+i*wind_shift))&&(eta_rec<-eta_max+wind_width+i*wind_shift))
                //                {
                //                    multB_rec[i]++;
                //                    if(charge_rec>0){multPlusB_rec[i]++;}
                //                    if(charge_rec<0){multMinusB_rec[i]++;}
                //                }
                //            }
                //
                //            for(int i = 0; i<NumbOne_winds;i++)
                //            {
                //                if((eta_rec<(eta_max-i*wind_shift))&&(eta_rec>(-eta_max+i*wind_shift)))
                //                {
                //                    multOne_rec[i]++;
                //                    PTOne_rec[i]+=pt_rec;
                //                    PT2One_rec[i]+=pow(pt_rec,2);
                //
                //                    if(charge_rec>0)
                //                    {
                //                        multOnePlus_rec[i]++;
                //                        PTOnePlus_rec[i]+=pt_rec;
                //                        PT2OnePlus_rec[i]+=pow(pt_rec,2);
                //                    }
                //                    if(charge_rec<0)
                //                    {
                //                        multOneMinus_rec[i]++;
                //                        PTOneMinus_rec[i]+=pt_rec;
                //                        PT2OneMinus_rec[i]+=pow(pt_rec,2);
                //                    }
                //                }
                //            }
            } // end of rec track loop
            
            bool found_at_least_one_match = false;
            // loop over mc tracks
            for (Int_t iTrk = 0; iTrk < dst->numberOfMcTracks(); iTrk++) {
                
                MpdMiniMcTrack* miniMcTrack = dst->mcTrack(iTrk);
                
                if (!miniMcTrack)
                    continue;
                
                p_sim=0; px_sim=0; py_sim=0; pz_sim=0; pt_sim=0; phi_sim=0; eta_sim=0; E_sim=0; y_sim=0; teta_sim=0; charge_sim=0;
                
                p_sim = miniMcTrack->p().Mag();
                px_sim =  miniMcTrack->p().X();
                py_sim =  miniMcTrack->p().Y();
                pz_sim =  miniMcTrack->p().Z();
                pt_sim = sqrt(px_sim*px_sim+py_sim*py_sim);
                phi_sim = TMath::ATan2(py_sim,px_sim);
                eta_sim = 0.5 * log((sqrt(px_sim*px_sim + py_sim*py_sim + pz_sim*pz_sim) + pz_sim) / (sqrt(px_sim*px_sim + py_sim*py_sim + pz_sim*pz_sim) - pz_sim));
                E_sim = sqrt(pow(mPion, 2) + pow(px_sim, 2) + pow(py_sim, 2) + pow(pz_sim, 2));
                y_sim = 0.5* log((E_sim + pz_sim) / (E_sim - pz_sim));
                teta_sim = 2 * atan(exp(-eta_sim));
                
                // sim track cuts
                
                bool generator=miniMcTrack->isFromGenerator();
                if(!generator)
                {
                    continue;
                }
                
                Int_t pdgId_MC = miniMcTrack->pdgId();
                if (TDatabasePDG::Instance()->GetParticle(pdgId_MC)==NULL) {
                    continue;
                }
                
                charge_sim = (miniMcTrack->charge())/3.0;
                if (charge_sim==0) {
                    continue;
                }
                
                if (pt_sim >2 || pt_sim <0.15)
                {
                    continue;
                }
                
                if((eta_sim>1)||(eta_sim<-1))
                {
                    continue;
                }
                
                for(int i = 0; i<Numb_winds;i++)
                {
                    if((eta_sim<(eta_max-i*wind_shift))&&(eta_sim>eta_max-wind_width-i*wind_shift))
                    {
                        multF_sim[i]++;
                        if(charge_sim>0){multPlusF_sim[i]++;}
                        if(charge_sim<0){multMinusF_sim[i]++;}
                    }
                    
                    if((eta_sim>(-eta_max+i*wind_shift))&&(eta_sim<-eta_max+wind_width+i*wind_shift))
                    {
                        multB_sim[i]++;
                        if(charge_sim>0){multPlusB_sim[i]++;}
                        if(charge_sim<0){multMinusB_sim[i]++;}
                    }
                }
                for(int i = 0; i<NumbOne_winds;i++)
                {
                    if((eta_sim<(eta_max-i*wind_shift))&&(eta_sim>(-eta_max+i*wind_shift)))
                    {
                        multOne_sim[i]++;
                        PTOne_sim[i]+=pt_sim;
                        PT2One_sim[i]+=pow(pt_sim,2);
                        
                        if(charge_sim>0)
                        {
                            multOnePlus_sim[i]++;
                            PTOnePlus_sim[i]+=pt_sim;
                            PT2OnePlus_sim[i]+=pow(pt_sim,2);
                        }
                        if(charge_sim<0)
                        {
                            multOneMinus_sim[i]++;
                            PTOneMinus_sim[i]+=pt_sim;
                            PT2OneMinus_sim[i]+=pow(pt_sim,2);
                        }
                    }
                }
                
                vector <UShort_t> miniTrackIdx = miniMcTrack->recoTrackIds();
                
                if(miniTrackIdx.size() != 0)
                {
                    vector <MpdMiniTrack*> miniTracks;
                    for (auto it : miniTrackIdx)
                        miniTracks.push_back(dst->track(it));
                    
                    double best_diff_px(0), best_diff_py(0), best_diff_pz(0);
                    int best_index=0;
                    double best_p2(0);
                    
                    if (goodTracksIndex.size()!=0)
                    {
                        //there are good rec tracks in this event
                        
                        for (auto it : miniTrackIdx)
                        {
                            if(dst->track(it)->mcTrackIndex()==miniMcTrack->id())
                            {
                                // ids are the same
                                if(std::find(goodTracksIndex.begin(), goodTracksIndex.end(), it) != goodTracksIndex.end())
                                {
                                    //matched rec track is in the list of good ones
                                    double tmp_p_rec_found(0),tmp_px_rec_found(0), tmp_py_rec_found(0), tmp_pz_rec_found(0); int tmp_charge_rec_found(0);
                                    double diff_px(0), diff_py(0), diff_pz(0);
                                    double tmp_p2(0);
                                    
                                    int index = std::distance(goodTracksIndex.begin(), find(goodTracksIndex.begin(), goodTracksIndex.end(), it));
                                    
                                    tmp_px_rec_found =goodVertexTracksPx.at(index);
                                    tmp_py_rec_found =goodVertexTracksPy.at(index);
                                    tmp_pz_rec_found =goodVertexTracksPz.at(index);
                                    tmp_charge_rec_found=goodVertexTracksCharge.at(index);
                                    
                                    if(tmp_charge_rec_found == charge_sim)
                                    {
                                        diff_px = abs(tmp_px_rec_found-px_sim);
                                        diff_py = abs(tmp_py_rec_found-py_sim);
                                        diff_pz = abs(tmp_pz_rec_found-pz_sim);
                                        
                                        tmp_p2=diff_px*diff_px+diff_py*diff_py+diff_pz*diff_pz;
                                        
                                        if(best_p2==0)
                                        {
                                            best_p2=tmp_p2;
                                            best_index = index;
                                        }
                                        else
                                        {
                                            if(tmp_p2<best_p2)
                                            {
                                                best_p2=tmp_p2;
                                                best_index = index;
                                            }
                                        }
                                    }
                                }
                                else
                                {
                                    // it is matched to sim track but it is not in the list of good rec tracks -- it is missed
                                    fill_miss(eta_sim,pt_sim,charge_sim,
                                              tempAll0, tempPlus0, tempMinus0, tempAll0_rec, tempPlus0_rec, tempMinus0_rec,
                                              tempAll1, tempPlus1, tempMinus1, tempAll1_rec, tempPlus1_rec, tempMinus1_rec,
                                              tempAll2, tempPlus2, tempMinus2, tempAll2_rec, tempPlus2_rec, tempMinus2_rec,
                                              tempAll3, tempPlus3, tempMinus3, tempAll3_rec, tempPlus3_rec, tempMinus3_rec,
                                              tempAll4, tempPlus4, tempMinus4, tempAll4_rec, tempPlus4_rec, tempMinus4_rec,
                                              tempAll5, tempPlus5, tempMinus5, tempAll5_rec, tempPlus5_rec, tempMinus5_rec,
                                              tempAll6, tempPlus6, tempMinus6, tempAll6_rec, tempPlus6_rec, tempMinus6_rec,
                                              tempAll7, tempPlus7, tempMinus7, tempAll7_rec, tempPlus7_rec, tempMinus7_rec,
                                              tempAll8, tempPlus8, tempMinus8, tempAll8_rec, tempPlus8_rec, tempMinus8_rec,
                                              tempAll9, tempPlus9, tempMinus9, tempAll9_rec, tempPlus9_rec, tempMinus9_rec);
                                }
                            }//current mc track has the same id as mc track matched to the checked rec track
                            else
                            {   // ids are different
                                fill_miss(eta_sim,pt_sim,charge_sim,
                                          tempAll0, tempPlus0, tempMinus0, tempAll0_rec, tempPlus0_rec, tempMinus0_rec,
                                          tempAll1, tempPlus1, tempMinus1, tempAll1_rec, tempPlus1_rec, tempMinus1_rec,
                                          tempAll2, tempPlus2, tempMinus2, tempAll2_rec, tempPlus2_rec, tempMinus2_rec,
                                          tempAll3, tempPlus3, tempMinus3, tempAll3_rec, tempPlus3_rec, tempMinus3_rec,
                                          tempAll4, tempPlus4, tempMinus4, tempAll4_rec, tempPlus4_rec, tempMinus4_rec,
                                          tempAll5, tempPlus5, tempMinus5, tempAll5_rec, tempPlus5_rec, tempMinus5_rec,
                                          tempAll6, tempPlus6, tempMinus6, tempAll6_rec, tempPlus6_rec, tempMinus6_rec,
                                          tempAll7, tempPlus7, tempMinus7, tempAll7_rec, tempPlus7_rec, tempMinus7_rec,
                                          tempAll8, tempPlus8, tempMinus8, tempAll8_rec, tempPlus8_rec, tempMinus8_rec,
                                          tempAll9, tempPlus9, tempMinus9, tempAll9_rec, tempPlus9_rec, tempMinus9_rec);
                            } // end of miss
                        } // loop over all found rec to this sim
                        
                        if(best_p2!=0)
                        {
                            found_at_least_one_match=true;
                            
                            double pt_rec_found(0), eta_rec_found(0), phi_rec_found(0), charge_rec_found(0);
                            
                            pt_rec_found = goodVertexTracksPt.at(best_index);
                            eta_rec_found = goodVertexTracksEta.at(best_index);
                            phi_rec_found = goodVertexTracksPhi.at(best_index);
                            charge_rec_found=goodVertexTracksCharge.at(best_index);
                            
                            for(int i = 0; i<Numb_winds;i++)
                            {
                                if((eta_rec_found<(eta_max-i*wind_shift))&&(eta_rec_found>eta_max-wind_width-i*wind_shift)&&(eta_sim<(eta_max-i*wind_shift))&&(eta_sim>eta_max-wind_width-i*wind_shift))
                                {
                                    multF_rec[i]++;
                                    if(charge_rec_found>0){multPlusF_rec[i]++;}
                                    if(charge_rec_found<0){multMinusF_rec[i]++;}
                                }
                                
                                if((eta_rec_found>(-eta_max+i*wind_shift))&&(eta_rec_found<-eta_max+wind_width+i*wind_shift)&&(eta_sim>(-eta_max+i*wind_shift))&&(eta_sim<-eta_max+wind_width+i*wind_shift))
                                {
                                    multB_rec[i]++;
                                    if(charge_rec_found>0){multPlusB_rec[i]++;}
                                    if(charge_rec_found<0){multMinusB_rec[i]++;}
                                }
                            }
                            
                            for(int i = 0; i<NumbOne_winds;i++)
                            {
                                if((eta_rec_found<(eta_max-i*wind_shift))&&(eta_rec_found>(-eta_max+i*wind_shift))&&(eta_sim<(eta_max-i*wind_shift))&&(eta_sim>(-eta_max+i*wind_shift)))
                                {
                                    multOne_rec[i]++;
                                    PTOne_rec[i]+=pt_rec_found;
                                    PT2One_rec[i]+=pow(pt_rec_found,2);
                                    
                                    if(charge_rec_found>0)
                                    {
                                        multOnePlus_rec[i]++;
                                        PTOnePlus_rec[i]+=pt_rec_found;
                                        PT2OnePlus_rec[i]+=pow(pt_rec_found,2);
                                    }
                                    if(charge_rec_found<0)
                                    {
                                        multOneMinus_rec[i]++;
                                        PTOneMinus_rec[i]+=pt_rec_found;
                                        PT2OneMinus_rec[i]+=pow(pt_rec_found,2);
                                    }
                                }
                            }
                            
                            // here fill RM
                            fill_RM(eta_sim,pt_sim,charge_sim,eta_rec_found,pt_rec_found,charge_sim,
                                    tempAll0, tempPlus0, tempMinus0, tempAll0_rec, tempPlus0_rec, tempMinus0_rec,
                                    tempAll1, tempPlus1, tempMinus1, tempAll1_rec, tempPlus1_rec, tempMinus1_rec,
                                    tempAll2, tempPlus2, tempMinus2, tempAll2_rec, tempPlus2_rec, tempMinus2_rec,
                                    tempAll3, tempPlus3, tempMinus3, tempAll3_rec, tempPlus3_rec, tempMinus3_rec,
                                    tempAll4, tempPlus4, tempMinus4, tempAll4_rec, tempPlus4_rec, tempMinus4_rec,
                                    tempAll5, tempPlus5, tempMinus5, tempAll5_rec, tempPlus5_rec, tempMinus5_rec,
                                    tempAll6, tempPlus6, tempMinus6, tempAll6_rec, tempPlus6_rec, tempMinus6_rec,
                                    tempAll7, tempPlus7, tempMinus7, tempAll7_rec, tempPlus7_rec, tempMinus7_rec,
                                    tempAll8, tempPlus8, tempMinus8, tempAll8_rec, tempPlus8_rec, tempMinus8_rec,
                                    tempAll9, tempPlus9, tempMinus9, tempAll9_rec, tempPlus9_rec, tempMinus9_rec);
                            
                            goodVertexTracksPt.erase(goodVertexTracksPt.begin() + best_index);
                            goodVertexTracksEta.erase(goodVertexTracksEta.begin() + best_index);
                            goodVertexTracksPhi.erase(goodVertexTracksPhi.begin() + best_index);
                            goodVertexTracksPx.erase(goodVertexTracksPx.begin() + best_index);
                            goodVertexTracksPy.erase(goodVertexTracksPy.begin() + best_index);
                            goodVertexTracksPz.erase(goodVertexTracksPz.begin() + best_index);
                            goodVertexTracksCharge.erase(goodVertexTracksCharge.begin() + best_index);
                            goodTracksIndex.erase(goodTracksIndex.begin() + best_index);
                            goodRecTracks.erase(goodRecTracks.begin() + best_index);
                            
                        } // found match
                    } // good rec tracks list is not empty
                    else
                    {
                        // this event doesn't contain good tracks - they are missed by cuts
                        fill_miss(eta_sim,pt_sim,charge_sim,
                                  tempAll0, tempPlus0, tempMinus0, tempAll0_rec, tempPlus0_rec, tempMinus0_rec,
                                  tempAll1, tempPlus1, tempMinus1, tempAll1_rec, tempPlus1_rec, tempMinus1_rec,
                                  tempAll2, tempPlus2, tempMinus2, tempAll2_rec, tempPlus2_rec, tempMinus2_rec,
                                  tempAll3, tempPlus3, tempMinus3, tempAll3_rec, tempPlus3_rec, tempMinus3_rec,
                                  tempAll4, tempPlus4, tempMinus4, tempAll4_rec, tempPlus4_rec, tempMinus4_rec,
                                  tempAll5, tempPlus5, tempMinus5, tempAll5_rec, tempPlus5_rec, tempMinus5_rec,
                                  tempAll6, tempPlus6, tempMinus6, tempAll6_rec, tempPlus6_rec, tempMinus6_rec,
                                  tempAll7, tempPlus7, tempMinus7, tempAll7_rec, tempPlus7_rec, tempMinus7_rec,
                                  tempAll8, tempPlus8, tempMinus8, tempAll8_rec, tempPlus8_rec, tempMinus8_rec,
                                  tempAll9, tempPlus9, tempMinus9, tempAll9_rec, tempPlus9_rec, tempMinus9_rec);
                    } // end of missed of good tracks in event
                }// this sim has at least one rec track
                else
                {
                    fill_miss(eta_sim,pt_sim,charge_sim,
                              tempAll0, tempPlus0, tempMinus0, tempAll0_rec, tempPlus0_rec, tempMinus0_rec,
                              tempAll1, tempPlus1, tempMinus1, tempAll1_rec, tempPlus1_rec, tempMinus1_rec,
                              tempAll2, tempPlus2, tempMinus2, tempAll2_rec, tempPlus2_rec, tempMinus2_rec,
                              tempAll3, tempPlus3, tempMinus3, tempAll3_rec, tempPlus3_rec, tempMinus3_rec,
                              tempAll4, tempPlus4, tempMinus4, tempAll4_rec, tempPlus4_rec, tempMinus4_rec,
                              tempAll5, tempPlus5, tempMinus5, tempAll5_rec, tempPlus5_rec, tempMinus5_rec,
                              tempAll6, tempPlus6, tempMinus6, tempAll6_rec, tempPlus6_rec, tempMinus6_rec,
                              tempAll7, tempPlus7, tempMinus7, tempAll7_rec, tempPlus7_rec, tempMinus7_rec,
                              tempAll8, tempPlus8, tempMinus8, tempAll8_rec, tempPlus8_rec, tempMinus8_rec,
                              tempAll9, tempPlus9, tempMinus9, tempAll9_rec, tempPlus9_rec, tempMinus9_rec);
                } // end of else for sim track with no rec
            } // end over sim track loop
            
            if(goodVertexTracksPt.size()!=0)
            {
                // fill addition fakes
                for(std::vector<double>::iterator it = goodVertexTracksPt.begin(); it != goodVertexTracksPt.end(); ++it)
                {
                    int pt_fake_position = std::distance(goodVertexTracksPt.begin(), find(goodVertexTracksPt.begin(), goodVertexTracksPt.end(), *it));
                    double pt_rec_found(0);
                    pt_rec_found = goodVertexTracksPt.at(pt_fake_position);
                    
                    fill_fake(eta_rec,pt_rec_found,charge_rec,
                              tempAll0, tempPlus0, tempMinus0, tempAll0_rec, tempPlus0_rec, tempMinus0_rec,
                              tempAll1, tempPlus1, tempMinus1, tempAll1_rec, tempPlus1_rec, tempMinus1_rec,
                              tempAll2, tempPlus2, tempMinus2, tempAll2_rec, tempPlus2_rec, tempMinus2_rec,
                              tempAll3, tempPlus3, tempMinus3, tempAll3_rec, tempPlus3_rec, tempMinus3_rec,
                              tempAll4, tempPlus4, tempMinus4, tempAll4_rec, tempPlus4_rec, tempMinus4_rec,
                              tempAll5, tempPlus5, tempMinus5, tempAll5_rec, tempPlus5_rec, tempMinus5_rec,
                              tempAll6, tempPlus6, tempMinus6, tempAll6_rec, tempPlus6_rec, tempMinus6_rec,
                              tempAll7, tempPlus7, tempMinus7, tempAll7_rec, tempPlus7_rec, tempMinus7_rec,
                              tempAll8, tempPlus8, tempMinus8, tempAll8_rec, tempPlus8_rec, tempMinus8_rec,
                              tempAll9, tempPlus9, tempMinus9, tempAll9_rec, tempPlus9_rec, tempMinus9_rec);
                }
            } // end of rec that were not matched to sim
            
            // matching is done
            
            particles_pt_in_event_sim.push_back(tempAll0);
            particles_pt_in_event_sim.push_back(tempAll1); particles_pt_in_event_sim.push_back(tempAll2); particles_pt_in_event_sim.push_back(tempAll3);
            particles_pt_in_event_sim.push_back(tempAll4); particles_pt_in_event_sim.push_back(tempAll5); particles_pt_in_event_sim.push_back(tempAll6);
            particles_pt_in_event_sim.push_back(tempAll7); particles_pt_in_event_sim.push_back(tempAll8); particles_pt_in_event_sim.push_back(tempAll9);
            
            particles_ptPlus_in_event_sim.push_back(tempPlus0);
            particles_ptPlus_in_event_sim.push_back(tempPlus1); particles_ptPlus_in_event_sim.push_back(tempPlus2); particles_ptPlus_in_event_sim.push_back(tempPlus3);
            particles_ptPlus_in_event_sim.push_back(tempPlus4); particles_ptPlus_in_event_sim.push_back(tempPlus5); particles_ptPlus_in_event_sim.push_back(tempPlus6);
            particles_ptPlus_in_event_sim.push_back(tempPlus7); particles_ptPlus_in_event_sim.push_back(tempPlus8); particles_ptPlus_in_event_sim.push_back(tempPlus9);
            
            particles_ptMinus_in_event_sim.push_back(tempMinus0);
            particles_ptMinus_in_event_sim.push_back(tempMinus1); particles_ptMinus_in_event_sim.push_back(tempMinus2); particles_ptMinus_in_event_sim.push_back(tempMinus3);
            particles_ptMinus_in_event_sim.push_back(tempMinus4); particles_ptMinus_in_event_sim.push_back(tempMinus5); particles_ptMinus_in_event_sim.push_back(tempMinus6);
            particles_ptMinus_in_event_sim.push_back(tempMinus7); particles_ptMinus_in_event_sim.push_back(tempMinus8); particles_ptMinus_in_event_sim.push_back(tempMinus9);
            
            particles_pt_in_event_rec.push_back(tempAll0_rec);
            particles_pt_in_event_rec.push_back(tempAll1_rec); particles_pt_in_event_rec.push_back(tempAll2_rec); particles_pt_in_event_rec.push_back(tempAll3_rec);
            particles_pt_in_event_rec.push_back(tempAll4_rec); particles_pt_in_event_rec.push_back(tempAll5_rec); particles_pt_in_event_rec.push_back(tempAll6_rec);
            particles_pt_in_event_rec.push_back(tempAll7_rec); particles_pt_in_event_rec.push_back(tempAll8_rec); particles_pt_in_event_rec.push_back(tempAll9_rec);
            
            particles_ptPlus_in_event_rec.push_back(tempPlus0_rec);
            particles_ptPlus_in_event_rec.push_back(tempPlus1_rec); particles_ptPlus_in_event_rec.push_back(tempPlus2_rec); particles_ptPlus_in_event_rec.push_back(tempPlus3_rec);
            particles_ptPlus_in_event_rec.push_back(tempPlus4_rec); particles_ptPlus_in_event_rec.push_back(tempPlus5_rec); particles_ptPlus_in_event_rec.push_back(tempPlus6_rec);
            particles_ptPlus_in_event_rec.push_back(tempPlus7_rec); particles_ptPlus_in_event_rec.push_back(tempPlus8_rec); particles_ptPlus_in_event_rec.push_back(tempPlus9_rec);
            
            particles_ptMinus_in_event_rec.push_back(tempMinus0_rec);
            particles_ptMinus_in_event_rec.push_back(tempMinus1_rec); particles_ptMinus_in_event_rec.push_back(tempMinus2_rec); particles_ptMinus_in_event_rec.push_back(tempMinus3_rec);
            particles_ptMinus_in_event_rec.push_back(tempMinus4_rec); particles_ptMinus_in_event_rec.push_back(tempMinus5_rec); particles_ptMinus_in_event_rec.push_back(tempMinus6_rec);
            particles_ptMinus_in_event_rec.push_back(tempMinus7_rec); particles_ptMinus_in_event_rec.push_back(tempMinus8_rec); particles_ptMinus_in_event_rec.push_back(tempMinus9_rec);
            
            tempAll0.clear(); tempPlus0.clear(); tempMinus0.clear(); tempAll0_rec.clear(); tempPlus0_rec.clear(); tempMinus0_rec.clear();
            tempAll1.clear(); tempPlus1.clear(); tempMinus1.clear(); tempAll1_rec.clear(); tempPlus1_rec.clear(); tempMinus1_rec.clear();
            tempAll2.clear(); tempPlus2.clear(); tempMinus2.clear(); tempAll2_rec.clear(); tempPlus2_rec.clear(); tempMinus2_rec.clear();
            tempAll3.clear(); tempPlus3.clear(); tempMinus3.clear(); tempAll3_rec.clear(); tempPlus3_rec.clear(); tempMinus3_rec.clear();
            tempAll4.clear(); tempPlus4.clear(); tempMinus4.clear(); tempAll4_rec.clear(); tempPlus4_rec.clear(); tempMinus4_rec.clear();
            tempAll5.clear(); tempPlus5.clear(); tempMinus5.clear(); tempAll5_rec.clear(); tempPlus5_rec.clear(); tempMinus5_rec.clear();
            tempAll6.clear(); tempPlus6.clear(); tempMinus6.clear(); tempAll6_rec.clear(); tempPlus6_rec.clear(); tempMinus6_rec.clear();
            tempAll7.clear(); tempPlus7.clear(); tempMinus7.clear(); tempAll7_rec.clear(); tempPlus7_rec.clear(); tempMinus7_rec.clear();
            tempAll8.clear(); tempPlus8.clear(); tempMinus8.clear(); tempAll8_rec.clear(); tempPlus8_rec.clear(); tempMinus8_rec.clear();
            tempAll9.clear(); tempPlus9.clear(); tempMinus9.clear(); tempAll9_rec.clear(); tempPlus9_rec.clear(); tempMinus9_rec.clear();
            
            N_Pt_Tree->Fill();
            
            particles_pt_in_event_sim.clear();
            particles_ptPlus_in_event_sim.clear();
            particles_ptMinus_in_event_sim.clear();
            
            particles_pt_in_event_rec.clear();
            particles_ptPlus_in_event_rec.clear();
            particles_ptMinus_in_event_rec.clear();
            
        } // event passed rec event cuts
        else
        {
            // missed events
            for(int i = 0; i<Numb_winds;i++)
            {
                if((eta_rec<(eta_max-i*wind_shift))&&(eta_rec>eta_max-wind_width-i*wind_shift))
                {
                    multF_rec[i]=-1;
                    multPlusF_rec[i]=-1;
                    multMinusF_rec[i]=-1;
                }
                
                if((eta_rec>(-eta_max+i*wind_shift))&&(eta_rec<-eta_max+wind_width+i*wind_shift))
                {
                    multB_rec[i]=-1;
                    multPlusB_rec[i]=-1;
                    multMinusB_rec[i]=-1;
                }
            }
            
            for(int i = 0; i<NumbOne_winds;i++)
            {
                if((eta_rec<(eta_max-i*wind_shift))&&(eta_rec>(-eta_max+i*wind_shift)))
                {
                    multOne_rec[i]=-1;
                    PTOne_rec[i]=-1;
                    PT2One_rec[i]=-1;
                    
                    multOnePlus_rec[i]=-1;
                    PTOnePlus_rec[i]=-1;
                    PT2OnePlus_rec[i]=-1;
                    
                    multOneMinus_rec[i]=-1;
                    PTOneMinus_rec[i]=-1;
                    PT2OneMinus_rec[i]=-1;
                }
            }
        }
    } //for(Long64_t iEvent=0; iEvent<events2read; iEvent++)
    
    cout<<"file writing"<<endl;
    
    oFile->cd();
    N_Pt_Tree->Write();
    oFile->Close();
    
    cout<<"job is done"<<endl;
    
    miniDstReader->Finish();
    
    return 0;
}
