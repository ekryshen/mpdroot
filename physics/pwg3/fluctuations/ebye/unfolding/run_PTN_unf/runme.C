#include "TROOT.h"
#include "TSystem.h"

void runme(int n_subsamples, int n_unf_iterations, int n_bins){
    
    gSystem->Load("/*put_your_path*/RooUnfold/libRooUnfold");
    gROOT->ProcessLine(".L final_unfolding_PTN.cc+");
    std::cout << "code compiled" << std::endl;
    
    char cmd[50];
    sprintf(cmd,"final_unfolding_PTN(%i,%i,%i)", n_subsamples, n_unf_iterations, n_bins);
    gROOT->ProcessLine(cmd);
  
}
