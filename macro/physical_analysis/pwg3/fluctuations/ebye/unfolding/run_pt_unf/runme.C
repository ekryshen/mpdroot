#include "TROOT.h"
#include "TSystem.h"

void runme(int n_subsamples, int n_unf_iterations){
    
    gSystem->Load("/*put_your_path*/RooUnfold/libRooUnfold");
    gROOT->ProcessLine(".L final_unfolding_pt.cc+");
    std::cout << "code compiled" << std::endl;
    
    char cmd[50];
    sprintf(cmd,"final_unfolding_pt(%i,%i)", n_subsamples, n_unf_iterations);
    gROOT->ProcessLine(cmd);
  
}
