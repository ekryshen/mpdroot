//------------------------------------------------------------------------------------------------------------------------
/// \class MpdFileSource
///
/// \brief
/// \author Sergei Lobastov (LHE, JINR, Dubna)
//------------------------------------------------------------------------------------------------------------------------
#include<iostream>

#include <TROOT.h>

#include "FairRootManager.h"
#include "MpdEventHeader.h"
#include "MpdFileSource.h"
using namespace std;

ClassImp(MpdFileSource);
//------------------------------------------------------------------------------------------------------------------------
MpdFileSource::MpdFileSource(const char *flnm)
   : FairFileSource(TString(flnm))
{

}
//------------------------------------------------------------------------------------------------------------------------
Bool_t		MpdFileSource::Init()
{
	if(IsInitialized)
    	{
        	LOG(debug)<<"MpdFileSource already initialized.";
        	return kTRUE;
    	}

    	if(!fInChain)
    	{
        	fInChain = new TChain(FairRootManager::GetTreeName(), "/cbmroot");
        	LOG(debug)<<"[MpdFileSource::Init()] chain created.";

        	FairRootManager::Instance()->SetInChain(fInChain);
    	}
    
	fInChain->Add(fRootFile->GetName());

    // Get The list of branches from the input file and add it to the
    // actual list of existing branches.
    // Add this list of branches also to the map of input trees, which
    // stores the information which branches belong to which input tree.
    // There is at least one primary input tree, but there can be many
    // additional friend trees.
    // This information is needed to add new files to the correct friend
    // tree. With this information it is also possible to check if the
    // input files which are added to the input chain all have the same
    // branch structure. Without this check it is possible to add trees
    // with a different branch structure but the same tree name. ROOT
    // probably only checks if the name of the tree is the same.

	TString chainName = fInputTitle;
	fInputLevel.push_back(chainName);
	fCheckInputBranches[chainName] = new list<TString>;

	TObjArray* fBranchList = fInChain->GetListOfBranches();
	// if no any branches then exit (some errors occured)
	if(fBranchList == nullptr)        return kFALSE;

	size_t Nentries = fBranchList->GetEntries();
	LOG(debug)<<"Entries in the chain "<<Nentries;

	auto ppObj = new TObject*[Nentries];
	for (int i = 0; i < Nentries; i++)
	{
        	auto pBranch = (TBranch*) fBranchList->At(i);
        	TString ObjName = pBranch->GetName();

        	fCheckInputBranches[chainName]->push_back(ObjName.Data());
        	FairRootManager::Instance()->AddBranchToList(ObjName.Data());

        	ppObj[i] = nullptr;
        	//ActivateObject(&(ppObj[i]), ObjName);
        	fInChain->SetBranchAddress(ObjName, &ppObj[i]);
        	FairRootManager::Instance()->RegisterInputObject(ObjName, ppObj[i]);

        	LOG(debug)<<"Branch name "<<ObjName.Data()<<", activated="<<fInChain->GetBranchStatus(ObjName.Data());
    	}
    
	// Add all additional input files to the input chain and do a
	// consitency check
	auto temp = gFile; // Store global gFile pointer for safety reasons.
    	for (auto iter = fInputChainList.begin(), iterEnd = fInputChainList.end(); iter != iterEnd; iter++)
    	{     
        	// Temporarily open the input file to extract information which
        	// is needed to bring the friend trees in the correct order
        	auto inputFile = new TFile(*iter);
        	if (inputFile->IsZombie())  LOG(fatal)<<"Error opening the file "<<(*iter).Data()<<" which should be added to the input chain or as friend chain";
        
        	// Check if the branchlist is the same as for the first input file.
        	Bool_t isOk = CompareBranchList(inputFile, chainName);
        	if (!isOk)
        	{
            		LOG(fatal)<<"Branch structure of the input file "<<fRootFile->GetName()<<" and the file to be added "<<(*iter).Data()<<" are different.";
            		return kFALSE;
        	}
        
        	// Add the file to the input chain
        	fInChain->Add(*iter);
        
        	// Close the temporarly file.
        	inputFile->Close();   
    	}
  	gFile = temp; // Restore global gFile pointer for safety reasons.

    	fNoOfEntries = fInChain->GetEntries(); 
    	LOG(debug)<<"Entries in this Source "<<fNoOfEntries;

fRootFile->ls();

	// Get the folder structure from file which describes the input tree.
	// There are two different names possible, so check both.
	fCbmroot = dynamic_cast<TFolder*>(fRootFile->Get(FairRootManager::GetFolderName()));

	if (!fCbmroot) 
	{
         fCbmroot = dynamic_cast<TFolder*>(fRootFile->Get("cbmroot"));
          if (!fCbmroot) {
               fCbmroot = dynamic_cast<TFolder*>(fRootFile->Get("cbmout"));
             if (!fCbmroot) {
                 fCbmroot = gROOT->GetRootFolder()->AddFolder(FairRootManager::GetFolderName(), "Main Folder");
              } else {
                 fCbmroot->SetName(FairRootManager::GetFolderName());
             }
          }
	}

	fListFolder->Add(fCbmroot);

	for (Int_t i = 0; i < fListFolder->GetEntriesFast(); i++) 
	{
		auto fold = static_cast<TFolder*>(fListFolder->At(i));
		fEvtHeader = static_cast<FairEventHeader*>(fold->FindObjectAny("MpdEventHeader."));

         	if(fEvtHeader)  ActivateObject(reinterpret_cast<TObject**>(&fEvtHeader), "MpdEventHeader.");
     	}

    	AddFriendsToChain();

return kTRUE;
}
//------------------------------------------------------------------------------------------------------------------------
void 		MpdFileSource::FillEventHeader(FairEventHeader* header)
{
	auto mpdHeader = dynamic_cast<MpdEventHeader*>(header);
	if(mpdHeader)
	{
		(*mpdHeader) = (*(MpdEventHeader*)fEvtHeader);
	}

	// control shot :-)
	FairFileSource::FillEventHeader(header);
}
//------------------------------------------------------------------------------------------------------------------------
