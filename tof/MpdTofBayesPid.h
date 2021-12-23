//------------------------------------------------------------------------------------------------------------------------
#ifndef MpdTofBayesPid_HH
#define MpdTofBayesPid_HH 

//------------------------------------------------------------------------------------------------------------------------
/// \class MpdTofBayesPid
/// 
/// \brief 
/// \author Sergei Lobastov (LHE, JINR, Dubna)
// usage:
// insert MpdTofBayesPid task into FairRunAna task chain AFTER MpdTofMatching task.
// stage 1: use  auto tofPid = new MpdTofBayesPid("PDFflnm.root") to create PDFs 
// stage 2: use  auto tofPid = new MpdTofBayesPid("PDFflnm.root", iterationNmb, "Priors_prefix", "nextIter_Priors_prefix") to create next iteration Priors; 
// 		 fRun->AddTask(tofPid) -- add task to chain 
// For iterationNmb=1 is created default const Priors automatically.
// Usually 3-4 iterations are required for Priors convergence.
// Pid probabilities is copied to MpdTofMatchingData::fBayesPid[4] TRANSPARENT(don't saved to disk) public member data.
//------------------------------------------------------------------------------------------------------------------------

#include "TSystem.h"
#include <TGraph.h>
#include <TFile.h>
#include <TString.h>
#include <TList.h>
#include <TH1D.h>

#include <TObject.h>
#include "FairTask.h"

class TH2D;
class TH3D;
class TClonesArray;
//class FairMCEventHeader;

//------------------------------------------------------------------------------------------------------------------------
class MpdBayesPriors : public TObject 
{
	static const size_t		fNpdf = 4;
	static const char*		fSpeciesNames[fNpdf]; 
	static const TH1D*		fPriorOrigin;

	size_t				fIterNmb; // iteration number [1,2,3,..]		
	TString				fFlnm;	
	TList				fList;		
	TH1D*				fPriors[fNpdf] = {}; // clones of fPriorOrigin

	TString				_getFullFlnm(size_t iterNmb, const char* flnm)const;
	TString				_getTHname(size_t iterNmb, size_t pid) const;
	void				_createPriors();
	void				_setDefaultPriors();
	bool				_loadPriors();
	void				_add(TH1D *h1);

public: 
	MpdBayesPriors(size_t iterNmb, const char* flnm, bool doLoad); 

	const TH1D* 		GetPrior(size_t pid)const{ assert(pid < fNpdf); return fPriors[pid];}
	double 			GetPrior(size_t pid, double beta) const;
	void			FillPrior(size_t pid, double beta, double weight = 1.);
	void			Yields2Priors(TH1D** data);
 	void			Write();
	void			DumpPriors(const char* comment, std::ostream& os) const;


	static const TH1D*	GetPriorOrigin(){return fPriorOrigin;}	
	static Int_t		GetPriorXbins(){return fPriorOrigin->GetXaxis()->GetNbins();}
	static size_t		GetNpdf(){return fNpdf;}
	static const char*	GetSpeciesName(size_t index);

ClassDef(MpdBayesPriors, 1)
};
//------------------------------------------------------------------------------------------------------------------------
class PidProbabilitiesMatrix 
{ 
	static const size_t	fNpdf = 4;

	TH1D*			fProb[fNpdf][fNpdf] = {};// εPID[i][j] - the probability to identify a species i as a species j (clone of MpdBayesPriors::fPriorOrigin)
	TH1D*			fA_true[fNpdf] = {}; 	 // A_true - true abundances of each species (clone of MpdBayesPriors::fPriorOrigin, used for test of fProb)
	TH1D*			fA_meas[fNpdf] = {}; 	 // A_meas - measured abundances of each species (clone of MpdBayesPriors::fPriorOrigin, used for test of fProb)
	double			fPidThres = -1.; 	 // if negative -  max. probability method(1) used, else - used as probability threshold for fixed threshold method(2) 
	TString			fName, fFlnm;
	TList			fList;	

	void 		_printLine(TH1D **data, const char* comment=nullptr, std::ostream& os=std::cout);
	void		_add(TH1D *h1);

public:
	PidProbabilitiesMatrix(double thresh, const char* prefix = nullptr);
	~PidProbabilitiesMatrix();

	double 		Calc_A_meas(size_t index)const; // MUST be called after FinalizeMatrix
	void		SetThresh(double value); 	// negative or [0., 1.]
	void		FinalizeMatrices();
	void		Update(double* prob_H_S, double Pt, Int_t pdgcode);
	void 		Print(const char* comment=nullptr,  std::ostream& os=std::cout);
	void		SetName(const char* name) { fName = name; };
	const char*	GetFlnm() const { return  fFlnm.Data(); };
	const TList&	GetHistos()const{return fList;};
	TH1D*		CreatePriorClone(TString name);
};
//------------------------------------------------------------------------------------------------------------------------
class MpdTofBayesPid : public FairTask  
{
	static const size_t		fNpdf = 4;

	TString				fPriorPrefix = "TofBayesPidPriors";
	const MpdBayesPriors		*fPriors = nullptr;		// used current priors
	MpdBayesPriors			*fNextPrior = nullptr;		// calc. next iteration priors

	static const size_t		fNmatrix = 10;
	PidProbabilitiesMatrix		*fEffMatrix[fNmatrix] = {};

	TH1D				*fYield[fNpdf] = {};		// clone of MpdBayesPriors::fPriorOrigin) 
	double				fProb_H_S[fNpdf] = {};		// 
	double				fProb_S_H[fNpdf] = {}; 		//  calc. from PDFs
	TH2D				*fH_S[fNpdf] = {};		// calc. H_S distribution
	double				fPdfIntegrals[fNpdf] = {};	// calc. from PDFs
	TH3D				*hPdf[fNpdf] = {}; 		// PDFs: 3D-function of mass², P, dE/dx (all species histos MUST BE same dim)

	Bool_t				fUseMCData; 			// true - MC data need to use
	Bool_t 				fDoTest = false; 		// true - make QA tests
	Bool_t 				fDoPDFs = false; 		// false - process iteration of data, true - create PDFs only

	TString				fQAflnm, fPDFflnm;
      	TList				fQAlist, fPDFlist;

        TClonesArray 			*aTofMatchings = nullptr;
        TClonesArray 			*aTPCkfTracks = nullptr;
        TClonesArray 			*aMcTracks = nullptr;
        TClonesArray 			*aTofHits = nullptr;

	bool		CreateLoadPDFs();
	void		FillPDFs();
	bool		_getProbs(Int_t rank, double mass2, double P, double dedx, double *data);
	bool		GetProb_S_H(double mass2, double P, double dedx, Int_t& rank, double *data);
	void		DumpProbs(double *data, const char* comment = nullptr, Int_t pdg = 0, std::ostream& os = std::cout) const;

	void		CalcProb_H_S(double P);

public:
	// process iterNmb iteration (load PDFs and priors from files; at 1th iteration used default priors =1.))
  	MpdTofBayesPid(const char *PDFflnm, size_t iterNmb, const char* priorPrefix = nullptr, const char* nextPriorPrefix = nullptr,
			 const char *QAflnm = nullptr, const char *taskName = "TOF_pid",  Int_t verbose = 1, Bool_t useMCdata = true);

	// create PDFs ONLY!!! (0th iteration, MUST be MC input data)
  	MpdTofBayesPid(const char *PDFflnm); 

	MpdTofBayesPid() = delete;
  	virtual ~MpdTofBayesPid();
  
	virtual InitStatus	Init();
	virtual void		Exec(Option_t * option);
	virtual void		Finish();

	template<typename T>void 	Add(TList&, T *hist);
	static	Int_t			GetIndex(Int_t pdgcode); // return pid index (same as fPdfNames[] array) for this pdfcode, negative - if don't found
	void				SetProbMatrixThresh(size_t i, double value); // set threshould for ith PidProbabilitiesMatrix

ClassDef(MpdTofBayesPid, 1)
};
//------------------------------------------------------------------------------------------------------------------------
#endif

