#include <iostream>
#include <vector>

#include "TSystem.h"
#include "TSystemFile.h"
#include "TSystemDirectory.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TTree.h"
#include "TList.h"
#include "TVector3.h"
#include "TMath.h"
#include "TClonesArray.h"

#include "FairLink.h"

#include "CeeMCTrack.h"
#include "CeeTpcPoint.h"
#include "CeeMWDCPoint.h"
#include "CeeiTOFPoint.h"
#include "CeeeTOFPoint.h"
#include "CeeZdcPoint.h"
#include "CeeTpcHit.h"
#include "CeeiTOFHit.h"
#include "CeeEtofHit.h"
#include "CeeZdcHit.h"
#include "CeeTpcTrack.h"
#include "CeeTpcVertex.h"
#include "CeeMWDC_Hit.h"
#include "CeeMWDC_Track.h"

#include "AnaMaker.h"

using std::cout;
using std::endl;

AnaMaker::AnaMaker(std::vector<std::string> fileList, std::string outName) : 
// fTpcPoint(nullptr),
// fMWDCPoint(nullptr),
// fZdcPoint(nullptr),
// fZdcHit(nullptr),
// fMWDC_Hit(nullptr),
fChain(nullptr),
fMCTrack(nullptr),
fiTOFPoint(nullptr),
feTOFPoint(nullptr),
fTpcHit(nullptr),
fiTOFHit(nullptr),
fEtofHit(nullptr),
fTpcTrack(nullptr),
fTpcVertex(nullptr),
fMWDC_Track(nullptr),
mNev(0),
files(std::move(fileList)), 
fOutName(outName) {
	fMCTrack = new TClonesArray("CeeMCTrack");
	fiTOFPoint = new TClonesArray("CeeiTOFPoint");
	feTOFPoint = new TClonesArray("CeeeTOFPoint");
	fTpcHit = new TClonesArray("CeeTpcHit");
	fiTOFHit = new TClonesArray("CeeiTOFHit");
	fEtofHit = new TClonesArray("CeeEtofHit");
	fTpcTrack = new TClonesArray("CeeTpcTrack");
	fTpcVertex = new TClonesArray("CeeTpcVertex");
	fMWDC_Track = new TClonesArray("CeeMWDC_Track");
	// fTpcPoint = new TClonesArray("CeeTpcPoint");
	// fMWDCPoint = new TClonesArray("CeeMWDCPoint");
	// fZdcPoint = new TClonesArray("CeeZdcPoint");
	// fZdcHit = new TClonesArray("CeeZdcHit");
	// fMWDC_Hit = new TClonesArray("CeeMWDC_Hit");
}

AnaMaker::~AnaMaker() {
	if (fMCTrack)    { fMCTrack->Clear("C"); delete fMCTrack; }
	if (fiTOFPoint)  { fiTOFPoint->Clear("C"); delete fiTOFPoint; }
	if (feTOFPoint)  { feTOFPoint->Clear("C"); delete feTOFPoint; }
	if (fTpcHit)     { fTpcHit->Clear("C"); delete fTpcHit; }
	if (fiTOFHit)    { fiTOFHit->Clear("C"); delete fiTOFHit; }
	if (fEtofHit)    { fEtofHit->Clear("C"); delete fEtofHit; }
	if (fTpcTrack)   { fTpcTrack->Clear("C"); delete fTpcTrack; }
	if (fTpcVertex)  { fTpcVertex->Clear("C"); delete fTpcVertex; }
	if (fMWDC_Track) { fMWDC_Track->Clear("C"); delete fMWDC_Track; }
	//if (fZdcHit)     { fZdcHit->Clear("C"); delete fZdcHit; }
	//if (fZdcPoint)   { fZdcPoint->Clear("C"); delete fZdcPoint; }
	//if (fMWDC_Hit)   { fMWDC_Hit->Clear("C"); delete fMWDC_Hit; }
	//if (fTpcPoint)   { fTpcPoint->Clear("C"); delete fTpcPoint; }
	//if (fMWDCPoint)  { fMWDCPoint->Clear("C"); delete fMWDCPoint; }
}

void AnaMaker::Make(Int_t iEv) {

	if (iEv > mNev) {
		cout << "[LOG] - AnaMaker: Entry index out of range " << iEv << " / " << mNev << endl;
	}

	if (iEv % 500 == 0) {
		std::cout << "[LOG] - AnaMaker: |==> " 
			<< std::setw(8) << std::setfill(' ') << iEv << " / "
			<< std::setw(8) << std::setfill(' ') << mNev << " >==|" << endl;
	}

	fChain->GetEntry(iEv);

	// Event (vertex)
	for (int i = 0; i < fTpcVertex->GetEntries(); i++) {
		mTpcVertex = dynamic_cast<CeeTpcVertex*>(fTpcVertex->At(i));
		if (mTpcVertex) { break; } // get the first valid vertex -> I don't really understand why the structure is like this
	}
	if (!mTpcVertex) { return; }
	TVector3 vertex(mTpcVertex->GetX()[0], mTpcVertex->GetY()[0], mTpcVertex->GetZ()[0]);
	h1TPCE_Vz->Fill(vertex.Z());
	h2TPCE_Vx_Vy->Fill(vertex.X(), vertex.Y());

	// MC tracks
	for (int i = 0; i < fMCTrack->GetEntries(); i++) {
		mMCTrack = dynamic_cast<CeeMCTrack*>(fMCTrack->At(i));
		if (!mMCTrack) { continue; }
		TVector3 mom3(mMCTrack->GetPx(), mMCTrack->GetPy(), mMCTrack->GetPz());
		double E = mom3.Mag2() + 0.938272 * 0.938272; // as proton
		double rapidity = 0.5 * log((E + mom3.Z()) / (E - mom3.Z()));

		h1MC_p->Fill(mom3.Mag());
		h1MC_pt->Fill(mom3.Perp());
		h1MC_px->Fill(mom3.X());
		h1MC_py->Fill(mom3.Y());
		h1MC_pz->Fill(mom3.Z());
		h2MC_px_py->Fill(mom3.X(), mom3.Y());
		h1MC_eta->Fill(mom3.Eta());
		h1MC_phi->Fill(mom3.Phi());
		h1MC_theta->Fill(mom3.Theta());
		h2MC_eta_pt->Fill(mom3.Eta(), mom3.Perp());
		h2MC_theta_phi->Fill(mom3.Theta(), mom3.Phi());
		h2MC_y_pt_proton->Fill(rapidity, mom3.Perp());
	}

    // TPC tracks
	for (int i = 0; i < fTpcTrack->GetEntries(); i++) {
		mTpcTrack = dynamic_cast<CeeTpcTrack*>(fTpcTrack->At(i));
		if (!mTpcTrack) { continue; }
		TVector3 mom3(mTpcTrack->GetPx(), mTpcTrack->GetPy(), mTpcTrack->GetPz());
		double E = mom3.Mag2() + 0.938272 * 0.938272; // as proton
		double rapidity = 0.5 * log((E + mom3.Z()) / (E - mom3.Z()));

		h1TPC_p->Fill(mom3.Mag());
		h1TPC_pt->Fill(mom3.Perp());
		h1TPC_px->Fill(mom3.X());
		h1TPC_py->Fill(mom3.Y());
		h1TPC_pz->Fill(mom3.Z());
		h2TPC_px_py->Fill(mom3.X(), mom3.Y());
		h1TPC_eta->Fill(mom3.Eta());
		h1TPC_phi->Fill(mom3.Phi());
		h1TPC_theta->Fill(mom3.Theta());
		h2TPC_dEdx_p->Fill(mom3.Mag(), mTpcTrack->GetdEdx() / 1000.0 * 26.7);
		h2TPC_eta_pt->Fill(mom3.Eta(), mom3.Perp());
		h2TPC_theta_phi->Fill(mom3.Theta(), mom3.Phi());
		h2TPC_y_pt_proton->Fill(rapidity, mom3.Perp());
		h1TPC_sDcax->Fill(mTpcTrack->GetDcaX());
		h1TPC_sDcay->Fill(mTpcTrack->GetDcaY());
		h1TPC_sDcaz->Fill(mTpcTrack->GetDcaZ());
		h1TPC_Dca->Fill(sqrt(
			mTpcTrack->GetDcaX() * mTpcTrack->GetDcaX() + 
			mTpcTrack->GetDcaY() * mTpcTrack->GetDcaY() + 
			mTpcTrack->GetDcaZ() * mTpcTrack->GetDcaZ()
		));
		h1TPC_nHitsFit->Fill(mTpcTrack->GetnHitsFit());
		h1TPC_nHitsDedx->Fill(mTpcTrack->GetnHitsDedx());
		h1TPC_nHitsRatio->Fill(
			mTpcTrack->GetnHitsMax() == 0 ? 
			0 : mTpcTrack->GetnHitsFit() * 1.0 / mTpcTrack->GetnHitsMax()
		);
		h1TPC_nHitsTrackLength->Fill(mTpcTrack->GetFitTracklength());
	}

	// MWDC tracks
	for (int i = 0; i < fMWDC_Track->GetEntries(); i++) {
		mMWDC_Track = dynamic_cast<CeeMWDC_Track*>(fMWDC_Track->At(i));
		if (!mMWDC_Track) { continue; }
		TVector3 mom3(mMWDC_Track->p);
		double E = mom3.Mag2() + 0.938272 * 0.938272; // as proton
		double rapidity = 0.5 * log((E + mom3.Z()) / (E - mom3.Z()));

		h1MWDC_p->Fill(mom3.Mag());
		h1MWDC_pt->Fill(mom3.Perp());
		h1MWDC_px->Fill(mom3.X());
		h1MWDC_py->Fill(mom3.Y());
		h1MWDC_pz->Fill(mom3.Z());
		h2MWDC_px_py->Fill(mom3.X(), mom3.Y());
		h1MWDC_eta->Fill(mom3.Eta());
		h1MWDC_phi->Fill(mom3.Phi());
		h1MWDC_theta->Fill(mom3.Theta());
		h2MWDC_dEdx_p->Fill(mom3.Mag(), mMWDC_Track->dEdx / 1000.0 * 26.7);
		h2MWDC_eta_pt->Fill(mom3.Eta(), mom3.Perp());
		h2MWDC_theta_phi->Fill(mom3.Theta(), mom3.Phi());
		h2MWDC_y_pt_proton->Fill(rapidity, mom3.Perp());
		h1MWDC_nHit->Fill(mMWDC_Track->num_hit);
		h1MWDC_nTrueHit->Fill(mMWDC_Track->num_mchit_true);
	}

	// iTOF hits
	for (int i = 0; i < fiTOFHit->GetEntries(); i++) {
		miTOFHit = dynamic_cast<CeeiTOFHit*>(fiTOFHit->At(i));
		if (!miTOFHit) { continue; }
		Int_t link2MCTrack = miTOFHit->GetLink(0).GetIndex();
		Int_t link2iTOFPoint = miTOFHit->GetLink(1).GetIndex();
		if (link2iTOFPoint <= 0) { continue; }
		mMCTrack = dynamic_cast<CeeMCTrack*>(fMCTrack->At(link2MCTrack));
		miTOFPoint = dynamic_cast<CeeiTOFPoint*>(fiTOFPoint->At(link2iTOFPoint));
		if (!mMCTrack) { continue; }
		if (!miTOFPoint) { continue; }

		TVector3 pos_hit(miTOFHit->GetX(), miTOFHit->GetY(), miTOFHit->GetZ() + 35);
		double dl = pos_hit.Mag() * 0.01; // cm -> m
		double dt = miTOFHit->GetTime() * 29.98;
		double beta = dt == 0 ? -1 : dl / dt;
		double gamma = 1.0 / sqrt(1 - beta * beta);
		double mass = mMCTrack->GetMass();
		double p_iTOF = mass * beta * gamma;
		
		h1ITOF_beta->Fill(beta);
		h2ITOF_p_invbeta->Fill(p_iTOF, 1.0 / beta);
		h2ITOF_p_mass2->Fill(p_iTOF, mass*mass);
		h2ITOF_theta_phi->Fill(pos_hit.Theta(), pos_hit.Phi());
		h2ITOF_perp_z->Fill(pos_hit.Perp(), pos_hit.Z());
	}

	// eTOF hits
	for (int i = 0; i < fEtofHit->GetEntries(); i++) {
		mEtofHit= dynamic_cast<CeeEtofHit*>(fEtofHit->At(i));
		if (!mEtofHit) { continue; }
		Int_t link2MCTrack = mEtofHit->GetLink(0).GetIndex();
		Int_t link2eTOFPoint = mEtofHit->GetLink(1).GetIndex();
		if (link2eTOFPoint <= 0) { continue; }
		mMCTrack = dynamic_cast<CeeMCTrack*>(fMCTrack->At(link2MCTrack));
		meTOFPoint = dynamic_cast<CeeeTOFPoint*>(feTOFPoint->At(link2eTOFPoint));
		if (!mMCTrack) { continue; }
		if (!meTOFPoint) { continue; }

		TVector3 pos_hit(mEtofHit->GetX(), mEtofHit->GetY(), mEtofHit->GetZ() + 35);
		double dl = pos_hit.Mag() * 0.01; // cm -> m
		double dt = miTOFHit->GetTime() * 29.98;
		double beta = dt == 0 ? -1 : dl / dt;
		double gamma = 1.0 / sqrt(1 - beta * beta);
		double mass = mMCTrack->GetMass();
		double p_eTOF = mass * beta * gamma;
		
		h1ETOF_beta->Fill(beta);
		h2ETOF_p_invbeta->Fill(p_eTOF, 1.0 / beta);
		h2ETOF_p_mass2->Fill(p_eTOF, mass*mass);
		h2ETOF_theta_phi->Fill(pos_hit.Theta(), pos_hit.Phi());
		h2ETOF_perp_z->Fill(pos_hit.Perp(), pos_hit.Z());
	}
}

bool AnaMaker::Init() {
	fChain = new TChain("ceesim");

	for (const auto& fName : files) {
		cout << "[LOG] - AnaMaker: Add " << fName << " into the chain" << endl;
		fChain->Add(fName.c_str());
	}
	if (mNev == 0) {
		mNev = fChain->GetEntries();
		cout << "[LOG] - AnaMaker: In total " << mNev << " events in chain" << endl;
	}
	
	cout << "[LOG] - AnaMaker: Initializing tree branch address" << endl;
	// branches in the tree
	fChain->SetBranchAddress("MCTrack", &fMCTrack);
	fChain->SetBranchAddress("CeeiTOFPoint", &fiTOFPoint);
	fChain->SetBranchAddress("CeeeTOFPoint", &feTOFPoint);
	fChain->SetBranchAddress("CeeTpcHit", &fTpcHit);
	fChain->SetBranchAddress("CeeiTOFHit", &fiTOFHit);
	fChain->SetBranchAddress("CeeEtofHit", &fEtofHit);
	fChain->SetBranchAddress("CeeTpcTrack", &fTpcTrack);
	fChain->SetBranchAddress("CeeTpcVertex", &fTpcVertex);
	fChain->SetBranchAddress("MWDC_Track", &fMWDC_Track);
	//fChain->SetBranchAddress("CeeZdcPoint", &fZdcPoint);
	//fChain->SetBranchAddress("CeeZdcHit", &fZdcHit);
	// fChain->SetBranchAddress("MWDC_Hit", &fMWDC_Hit);
	//fChain->SetBranchAddress("CeeTpcPoint", &fTpcPoint);
	//fChain->SetBranchAddress("CeeMWDCPoint", &fMWDCPoint);

	// histograms
	cout << "[LOG] - AnaMaker: Initializing histograms" << endl;
	h1TPCE_Vz = new TH1F(
		"h1TPCE_Vz",
		"h1TPCE_Vz;V_{z} (cm);Counts",
		100, -40, -30
	);
	h2TPCE_Vx_Vy = new TH2F(
		"h2TPCE_Vz",
		"h2TPCE_Vz;V_{x} (cm);V_{y} (cm)",
		90, -1.5, 1.5,
		90, -1.5, 1.5
	);
	h1MC_p = new TH1F(
		"h1MC_p",
		"h1MC_p;p (GeV/c);Counts",
		100, 0, 10
	);
	h1MC_pt = new TH1F(
		"h1MC_pt",
		"h1MC_pt;p_{T} (GeV/c);Counts",
		100, 0, 5
	);
	h1MC_px = new TH1F(
		"h1MC_px",
		"h1MC_px;p_{x} (GeV/c);Counts",
		100, -2, 2
	);
	h1MC_py = new TH1F(
		"h1MC_py",
		"h1MC_py;p_{y} (GeV/c);Counts",
		100, -2, 2
	);
	h2MC_px_py = new TH2F(
		"h2MC_px_py",
		"h2MC_px_py;p_{x} (GeV/c);p_{y} (GeV/c);Counts",
		100, -2, 2,
		100, -2, 2
	);
	h1MC_pz = new TH1F(
		"h1MC_pz",
		"h1MC_pz;p_{z} (GeV/c);Counts",
		200, -5, 5
	);
	h1MC_eta = new TH1F(
		"h1MC_eta",
		"h1MC_eta;#eta;Counts",
		600, -6, 6
	);
	h1MC_phi = new TH1F(
		"h1MC_phi",
		"h1MC_phi;#phi;Counts",
		400, -4, 4
	);
	h1MC_theta = new TH1F(
		"h1MC_theta",
		"h1MC_theta;#theta;Counts",
		200, 0, 4
	);
	h2MC_eta_pt = new TH2F(
		"h2MC_eta_pt",
		"h2MC_eta_pt;#eta;p_{T} (GeV/c);Counts",
		600, -6, 6,
		500, 0, 5
	);
	h2MC_theta_phi = new TH2F(
		"h2MC_theta_phi",
		"h2MC_theta_phi;#theta;#phi;Counts",
		200, 0, 4,
		400, -4, 4
	);
	h2MC_y_pt_proton = new TH2F(
		"h2MC_y_pt_proton", 
		"h2MC_y_pt_proton;y (proton);p_{T} (GeV/c);Counts",
		400, -1, 1,
		500, 0, 2.5
	);
	h1TPC_p = new TH1F(
		"h1TPC_p",
		"h1TPC_p;p (GeV/c);Counts",
		100, 0, 10
	);
	h1TPC_pt = new TH1F(
		"h1TPC_pt",
		"h1TPC_pt;p_{T} (GeV/c);Counts",
		100, 0, 5
	);
	h1TPC_px = new TH1F(
		"h1TPC_px",
		"h1TPC_px;p_{x} (GeV/c);Counts",
		100, -2, 2
	);
	h1TPC_py = new TH1F(
		"h1TPC_py",
		"h1TPC_py;p_{y} (GeV/c);Counts",
		100, -2, 2
	);
	h2TPC_px_py = new TH2F(
		"h2TPC_px_py",
		"h2TPC_px_py;p_{x} (GeV/c);p_{y} (GeV/c);Counts",
		100, -2, 2,
		100, -2, 2
	);
	h1TPC_pz = new TH1F(
		"h1TPC_pz",
		"h1TPC_pz;p_{z} (GeV/c);Counts",
		200, -5, 5
	);
	h1TPC_eta = new TH1F(
		"h1TPC_eta",
		"h1TPC_eta;#eta;Counts",
		600, -6, 6
	);
	h1TPC_phi = new TH1F(
		"h1TPC_phi",
		"h1TPC_phi;#phi;Counts",
		400, -4, 4
	);
	h1TPC_theta = new TH1F(
		"h1TPC_theta",
		"h1TPC_theta;#theta;Counts",
		200, 0, 4
	);
	h2TPC_dEdx_p = new TH2F(
		"h2TPC_dEdx_p",
		"h2TPC_dEdx_p;p (GeV/c); dE/dx (keV/cm);Counts",
		500, 0, 5,
		300, 0, 30
	);
	h2TPC_eta_pt = new TH2F(
		"h2TPC_eta_pt",
		"h2TPC_eta_pt;#eta;p_{T} (GeV/c);Counts",
		600, -6, 6,
		500, 0, 5
	);
	h2TPC_theta_phi = new TH2F(
		"h2TPC_theta_phi",
		"h2TPC_theta_phi;#theta;#phi;Counts",
		200, 0, 4,
		400, -4, 4
	);
	h2TPC_y_pt_proton = new TH2F(
		"h2TPC_y_pt_proton", 
		"h2TPC_y_pt_proton;y (proton);p_{T} (GeV/c);Counts",
		400, -1, 1,
		500, 0, 2.5
	);
	h1TPC_sDcax = new TH1F(
		"h1TPC_sDcax",
		"h1TPC_sDcax;sDCA_{x} (cm);Counts",
		200, -5, 5
	);
	h1TPC_sDcay = new TH1F(
		"h1TPC_sDcay",
		"h1TPC_sDcay;sDCA_{y} (cm);Counts",
		200, -5, 5
	);
	h1TPC_sDcaz = new TH1F(
		"h1TPC_sDcaz",
		"h1TPC_sDcaz;sDCA_{z} (cm);Counts",
		200, -5, 5
	);
	h1TPC_Dca = new TH1F(
		"h1TPC_Dca",
		"h1TPC_Dca;DCA (cm);Counts",
		100, 0, 5
	);
	h1TPC_nHitsFit = new TH1F(
		"h1TPC_nHitsFit",
		"h1TPC_nHitsFit;nHitsFit;Counts",
		60, 0, 60
	);
	h1TPC_nHitsDedx = new TH1F(
		"h1TPC_nHitsDedx",
		"h1TPC_nHitsDedx;nHitsDedx;Counts",
		60, 0, 60
	); 
	h1TPC_nHitsRatio = new TH1F(
		"h1TPC_nHitsRatio",
		"h1TPC_nHitsRatio;nHitsRatio;Counts",
		100, 0, 1
	);
	h1TPC_nHitsTrackLength = new TH1F(
		"h1TPC_nHitsTrackLength",
		"h1TPC_nHitsTrackLength;nHitsTrackLength (cm);Counts",
		500, 0, 50
	);
	h1MWDC_p = new TH1F(
		"h1MWDC_p",
		"h1MWDC_p;p (GeV/c);Counts",
		100, 0, 10
	);
	h1MWDC_pt = new TH1F(
		"h1MWDC_pt",
		"h1MWDC_pt;p_{T} (GeV/c);Counts",
		100, 0, 5
	);
	h1MWDC_px = new TH1F(
		"h1MWDC_px",
		"h1MWDC_px;p_{x} (GeV/c);Counts",
		100, -2, 2
	);
	h1MWDC_py = new TH1F(
		"h1MWDC_py",
		"h1MWDC_py;p_{y} (GeV/c);Counts",
		100, -2, 2
	);
	h2MWDC_px_py = new TH2F(
		"h2MWDC_px_py",
		"h2MWDC_px_py;p_{x} (GeV/c);p_{y} (GeV/c);Counts",
		100, -2, 2,
		100, -2, 2
	);
	h1MWDC_pz = new TH1F(
		"h1MWDC_pz",
		"h1MWDC_pz;p_{z} (GeV/c);Counts",
		200, -5, 5
	);
	h1MWDC_eta = new TH1F(
		"h1MWDC_eta",
		"h1MWDC_eta;#eta;Counts",
		600, -6, 6
	);
	h1MWDC_phi = new TH1F(
		"h1MWDC_phi",
		"h1MWDC_phi;#phi;Counts",
		400, -4, 4
	);
	h1MWDC_theta = new TH1F(
		"h1MWDC_theta",
		"h1MWDC_theta;#theta;Counts",
		200, 0, 4
	);
	h2MWDC_dEdx_p = new TH2F(
		"h2MWDC_dEdx_p",
		"h2MWDC_dEdx_p;p (GeV/c); dE/dx (keV/cm);Counts",
		500, 0, 5,
		300, 0, 30
	);
	h2MWDC_eta_pt = new TH2F(
		"h2MWDC_eta_pt",
		"h2MWDC_eta_pt;#eta;p_{T} (GeV/c);Counts",
		600, -6, 6,
		500, 0, 5
	);
	h2MWDC_theta_phi = new TH2F(
		"h2MWDC_theta_phi",
		"h2MWDC_theta_phi;#theta;#phi;Counts",
		200, 0, 4,
		400, -4, 4
	);
	h2MWDC_y_pt_proton = new TH2F(
		"h2MWDC_y_pt_proton", 
		"h2MWDC_y_pt_proton;y (proton);p_{T} (GeV/c);Counts",
		400, -1, 1,
		500, 0, 2.5
	);
	h1MWDC_nHit = new TH1F(
		"h1MWDC_nHit",
		"h1MWDC_nHit;nHit;Counts",
		100, 0, 100
	);
	h1MWDC_nTrueHit = new TH1F(
		"h1MWDC_nTrueHit",
		"h1MWDC_nTrueHit;nHit (true);Counts",
		100, 0, 100
	);
	h1ITOF_beta = new TH1F(
		"h1ITOF_beta",
		"h1ITOF_beta;#beta;Counts",
		120, -0.1, 1.1
	);
	h2ITOF_p_invbeta = new TH2F(
		"h2ITOF_p_invbeta",
		"h2ITOF_p_invbeta;p (GeV/c);1/#beta",
		500, 0, 5,
		300, 0, 3
	);
	h2ITOF_p_mass2 = new TH2F(
		"h2ITOF_p_mass2",
		"h2ITOF_p_mass2;p (GeV/c);m^{2} (GeV^{2}/c^{4});Counts",
		500, 0, 5,
		200, 0, 2
	);
	h2ITOF_theta_phi = new TH2F(
		"h2ITOF_theta_phi",
		"h2ITOF_theta_phi;#theta;#phi;Counts",
		200, 0, 4,
		400, -4, 4
	);
	h2ITOF_perp_z = new TH2F(
		"h2ITOF_perp_z",
		"h2ITOF_perp_z;X_{#perp} (cm);X_{z} (cm)",
		400, -4, 4,
		400, -4, 4
	);
	h1ETOF_beta = new TH1F(
		"h1ETOF_beta",
		"h1ETOF_beta;#beta;Counts",
		120, -0.1, 1.1
	);
	h2ETOF_p_invbeta = new TH2F(
		"h2ETOF_p_invbeta",
		"h2ETOF_p_invbeta;p (GeV/c);1/#beta",
		500, 0, 5,
		300, 0, 3
	);
	h2ETOF_p_mass2 = new TH2F(
		"h2ETOF_p_mass2",
		"h2ETOF_p_mass2;p (GeV/c);m^{2} (GeV^{2}/c^{4});Counts",
		500, 0, 5,
		200, 0, 2
	);
	h2ETOF_theta_phi = new TH2F(
		"h2ETOF_theta_phi",
		"h2ETOF_theta_phi;#theta;#phi;Counts",
		200, 0, 4,
		400, -4, 4
	);
	h2ETOF_perp_z = new TH2F(
		"h2ETOF_perp_z",
		"h2ETOF_perp_z;X_{#perp} (cm);X_{z} (cm)",
		400, -4, 4,
		400, -4, 4
	);

	return true;
}

void AnaMaker::Finish() {
	fOutFile = new TFile(fOutName.c_str(), "recreate");
	if (!fOutFile->IsOpen()) {
		cout << "[LOG] - AnaMaker: Failed to create the output file" << endl;
		return;
	}

	h1TPCE_Vz->Write();
	h2TPCE_Vx_Vy->Write();

	h1MC_p->Write();
	h1MC_pt->Write();
	h1MC_px->Write();
	h1MC_py->Write();
	h2MC_px_py->Write();
	h1MC_pz->Write();
	h1MC_eta->Write();
	h1MC_phi->Write();
	h1MC_theta->Write();
	h2MC_eta_pt->Write();
	h2MC_theta_phi->Write();
	h2MC_y_pt_proton->Write();

	h1TPC_p->Write();
	h1TPC_pt->Write();
	h1TPC_px->Write();
	h1TPC_py->Write();
	h2TPC_px_py->Write();
	h1TPC_pz->Write();
	h1TPC_eta->Write();
	h1TPC_phi->Write();
	h1TPC_theta->Write();
	h2TPC_dEdx_p->Write();
	h2TPC_eta_pt->Write();
	h2TPC_theta_phi->Write();
	h2TPC_y_pt_proton->Write();
	h1TPC_sDcax->Write();
	h1TPC_sDcay->Write();
	h1TPC_sDcaz->Write();
	h1TPC_Dca->Write();
	h1TPC_nHitsFit->Write();
	h1TPC_nHitsDedx->Write();
	h1TPC_nHitsRatio->Write();
	h1TPC_nHitsTrackLength->Write();

	h1MWDC_p->Write();
	h1MWDC_pt->Write();
	h1MWDC_px->Write();
	h1MWDC_py->Write();
	h2MWDC_px_py->Write();
	h1MWDC_pz->Write();
	h1MWDC_eta->Write();
	h1MWDC_phi->Write();
	h1MWDC_theta->Write();
	h2MWDC_dEdx_p->Write();
	h2MWDC_eta_pt->Write();
	h2MWDC_theta_phi->Write();
	h2MWDC_y_pt_proton->Write();
	h1MWDC_nHit->Write();
	h1MWDC_nTrueHit->Write();

	h1ITOF_beta->Write();
	h2ITOF_p_invbeta->Write();
	h2ITOF_p_mass2->Write();
	h2ITOF_theta_phi->Write();
	h2ITOF_perp_z->Write();

	h1ETOF_beta->Write();
	h2ETOF_p_invbeta->Write(); 
	h2ETOF_p_mass2->Write(); 
	h2ETOF_theta_phi->Write();
	h2ETOF_perp_z->Write();

	fOutFile->Close();

	cout << "[LOG] - AnaMaker: This is the end of this program" << endl;
}
