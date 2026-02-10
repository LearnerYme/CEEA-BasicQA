/*	Analysis Maker - Yige Huang
 *	Generation: 1
 *	Series: Vallina
 *	Inspired by Xionghong's AnalysisMaker
**/

#ifndef __ANAMAKER__
#define __ANAMAKER__

#include <vector>
#include <string>

#include "TString.h"
#include "TFile.h"

class CeeMCTrack;
class CeeTpcPoint;
class CeeMWDCPoint;
class CeeiTOFPoint;
class CeeeTOFPoint;
class CeeZdcPoint;
class CeeTpcHit;
class CeeiTOFHit;
class CeeEtofHit;
class CeeZdcHit;
class CeeTpcTrack;
class CeeTpcVertex;
class CeeMWDC_Hit;
class CeeMWDC_Track;

class TH1F;
class TH2F;
class TChain;
class TClonesArray;

class AnaMaker {
	public:
	
	AnaMaker(std::vector<std::string> fileList, std::string outName);
	~AnaMaker();

	bool Init();
	void Make(Int_t iEv);
	void Finish();

	Int_t GetAnaEntries() { return mNev; }

	private:
	
	TFile* fOutFile;
	TChain* fChain;
	std::string fOutName;
	std::vector<std::string> files;

	// histograms
	// Event-wise
	TH1F* h1TPCE_Vz;
	TH2F* h2TPCE_Vx_Vy;
	// MC tracks
	TH1F* h1MC_p;
	TH1F* h1MC_pt;
	TH1F* h1MC_px;
	TH1F* h1MC_py;
	TH2F* h2MC_px_py;
	TH1F* h1MC_pz;
	TH1F* h1MC_eta;
	TH1F* h1MC_phi;
	TH1F* h1MC_theta;
	TH2F* h2MC_eta_pt;
	TH2F* h2MC_theta_phi;
	TH2F* h2MC_y_pt_proton;
	// TPC tracks
	TH1F* h1TPC_p;
	TH1F* h1TPC_pt;
	TH1F* h1TPC_px;
	TH1F* h1TPC_py;
	TH2F* h2TPC_px_py;
	TH1F* h1TPC_pz;
	TH1F* h1TPC_eta;
	TH1F* h1TPC_phi;
	TH1F* h1TPC_theta;
	TH2F* h2TPC_dEdx_p;
	TH2F* h2TPC_eta_pt;
	TH2F* h2TPC_theta_phi;
	TH2F* h2TPC_y_pt_proton;
	TH1F* h1TPC_sDcax;
	TH1F* h1TPC_sDcay;
	TH1F* h1TPC_sDcaz;
	TH1F* h1TPC_Dca;
	TH1F* h1TPC_nHitsFit;
	TH1F* h1TPC_nHitsDedx;
	TH1F* h1TPC_nHitsRatio;
	TH1F* h1TPC_nHitsTrackLength;
	// MWDC tracks
	TH1F* h1MWDC_p;
	TH1F* h1MWDC_pt;
	TH1F* h1MWDC_px;
	TH1F* h1MWDC_py;
	TH2F* h2MWDC_px_py;
	TH1F* h1MWDC_pz;
	TH1F* h1MWDC_eta;
	TH1F* h1MWDC_phi;
	TH1F* h1MWDC_theta;
	TH2F* h2MWDC_dEdx_p;
	TH2F* h2MWDC_eta_pt;
	TH2F* h2MWDC_theta_phi;
	TH2F* h2MWDC_y_pt_proton;
	TH1F* h1MWDC_nHit;
	TH1F* h1MWDC_nTrueHit;
	// iTOF hits
	TH1F* h1ITOF_beta;
	TH2F* h2ITOF_p_invbeta; // this p is reconstructed by mass from MC
	TH2F* h2ITOF_p_mass2; // similar with invbeta
	TH2F* h2ITOF_theta_phi; // coordinate
	TH2F* h2ITOF_perp_z; // coordinate
	// eTOF hits
	TH1F* h1ETOF_beta;
	TH2F* h2ETOF_p_invbeta; 
	TH2F* h2ETOF_p_mass2; 
	TH2F* h2ETOF_theta_phi;
	TH2F* h2ETOF_perp_z;

	// CEE objects
	CeeMCTrack* mMCTrack;
	CeeiTOFPoint* miTOFPoint;
	CeeeTOFPoint* meTOFPoint;
	CeeTpcHit* mTpcHit;
	CeeiTOFHit* miTOFHit;
	CeeEtofHit* mEtofHit;
	CeeTpcTrack* mTpcTrack;
	CeeTpcVertex* mTpcVertex;
	CeeMWDC_Track* mMWDC_Track;
	// CeeTpcPoint* mTpcPoint;
	// CeeZdcPoint* mZdcPoint;
	// CeeZdcHit* mZdcHit;
	// CeeMWDC_Hit* mMWDC_Hit;
	// CeeMWDCPoint* mMWDCPoint;

	// clone arrays
	TClonesArray* fMCTrack;
	TClonesArray* fiTOFPoint;
	TClonesArray* feTOFPoint;
	TClonesArray* fTpcHit;
	TClonesArray* fiTOFHit;
	TClonesArray* fEtofHit;
	TClonesArray* fTpcTrack;
	TClonesArray* fTpcVertex;
	TClonesArray* fMWDC_Track;
	// TClonesArray* fZdcHit;
	// TClonesArray* fZdcPoint;
	// TClonesArray* fMWDC_Hit;
	// TClonesArray* fTpcPoint;
	// TClonesArray* fMWDCPoint;

	Int_t mNev;

};

#endif
