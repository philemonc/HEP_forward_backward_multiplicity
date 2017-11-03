/*
	NOTES
	=====
	- Relative directories should be updated to work with your computer; search for "//reldir"

*/


#include "TROOT.h"
#include "TTree.h"
#include "TFile.h"
#include "TChain.h"
#include "TChainElement.h"
#include "TDirectory.h"
#include "TH1.h"
#include "TH2.h"
#inclide "TProfile.h"
#include "TLorentzVector.h"
#include "TCanvas.h"
#include "TMath.h"
#include "Math/Vector3D.h"
#include "Math/Vector4D.h"
#include "TString.h"

using namespace ROOT::Math;

#include <vector>
#include <iostream>
#include <fstream>
#include <algorithm>
#include <cmath>

using namespace std;

#ifdef __MAKECINT__
    #pragma link C++ class vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > >++;
#endif


vector<TString>* getListOfFiles(TString strfiles){

  vector<TString>* vfiles = new vector<TString>;

  if(strfiles.Contains(".root")){
    TChain chain("tree/tree","");
    chain.Add(strfiles);
    TObjArray* fileElements=chain.GetListOfFiles();
    TIter next(fileElements);
    TChainElement *chEl=0;
    while (( chEl=(TChainElement*)next() )) {
      vfiles->push_back(TString(chEl->GetTitle()));
    }
  }
  else if(strfiles.Contains(".txt")){
    ifstream txtfile;
    txtfile.open(strfiles);
    if(!txtfile) {
      cout<<"Unable to read the txt file where the rootfiles are." << endl ;
      cout << strfiles << " doesn't exist." << endl << "Aborting ...";
      exit(0);
    }
    string filename;
    while(txtfile>>filename && filename!="EOF")
      vfiles->push_back(TString(filename));
    txtfile.close();
  }
  else {
    cout << "Unknown type of input to get files. Must contain either .root or .txt extension." << endl << "Aborting ..." << endl;
    exit(0);
  }
  cout << "[getListOfFiles] Will run on " << vfiles->size() << " files" << endl;
  return vfiles;
}

void forward_backward()
{
	// Recalculate error bars when scaling/adding/dividing histograms
	TH1::SetDefaultSumw2(1);   
	
    //============================== implement cuts ==============================
    const double eta_cut = 2.;	// TODO
    const double vtx_number_cut = 1.0;
    const float vtxxysize = 0.2;
    const float vtxzsize = 10;
    const double pt_cut = 0.5;
    const int lumi_cut = 90;
    const float ptErr_pt_cut = 0.05;
    const float dz_dzErr_cut = 3;
    const float d0_d0Err_cut = 3;
    const float dof_cut = 4;
	
	const double eta_fb_upper = 0.8;
	const double eta_fb_lower = 0.2;
	

    //============================== Variables ==============================

    double fdata_evt = 0;
    double fdata_trkphi = 0;
    double fdata_trketa = 0;
    double fdata_trkd0 = 0;
    double fdata_trkpt = 0;
    double fdata_trkdz = 0;
    double fdata_trkdpt = 0;
    double fdata_dz_sigmadz, fdata_d0_sigmad0, fdata_sigmapt_pt, fdata_dz_sigmadzcalc, fdata_d0_sigmad0calc, fdata_sigmad0, fdata_sigmad0run1, fdata_d0, fdata_d0_sigmad0run1, fdata_trkdx, fdata_trkdy;
    double fdata_wx, fdata_wy,fdata_wz, fdata_vtxxysize, fdata_vtxzsize, fdata_sigmad0calc, fdata_dz, fdata_sigmadz, fdata_dz_sigmadzrun1, fdata_d0calc, fdata_dzcalc, fdata_d0calc_sigmad0, fdata_dzcalc_sigmadz; 
    double fdata_dzBS_dzErr, fdata_dzvtxBS_dzErr, fdata_dxyBS_d0Err, fdata_dxyvtxBS_d0Err;
	double fdata_numberofvtxxBS, fdata_vtxxBSvalue, fdata_vtxxBSlower, fdata_vtxxBSupper;
    double fdata_numberofvtxyBS, fdata_vtxyBSvalue, fdata_vtxyBSlower, fdata_vtxyBSupper;
    double fdata_numberofvtxzBS, fdata_vtxzBSvalue, fdata_vtxzBSlower, fdata_vtxzBSupper;
    double fdata_vtxzminusvtxz, fdata_multiplicity, fdata_forward_multiplicity, fdata_backward_multiplicity;
    double fdata_multiplicity_norm = 0;
	//double fdata_fb_multiplicity_norm = 0;
    double fdata_sqvtxx = 0;
    double fdata_sqvtxxnumber = 0;
    double fdata_sqvtxy = 0;
    double fdata_sqvtxynumber = 0;
    double fdata_sqvtxz = 0;
    double fdata_sqvtxznumber = 0;
    double fdata_numberoftrkdx = 0;
    double fdata_numberoftrkdy = 0;
    double fdata_numselectedvtxz = 0;
    double fdata_numvtxzminusvtxz = 0;
    double fdata_trkvalidhits = 0;
    double fdata_trkchi2n = 0;
    int ndata_numberofvtxx, ndata_numberofvtxy, ndata_numberofvtxz, ndata_numberofvtxxBS, ndata_numberofvtxyBS, ndata_numberofvtxzBS;
    int ndata_totaltrk;

	
	//============================== Histos for pT, eta, phi ==============================
    TH1D *data_pt_histo = new TH1D ("data pT", "Normalized data p_{T}", 200, 0, 10);
    TH1D *data_eta_histo = new TH1D("data eta", "Normalised data #eta", 50, -2.5, 2.5);
    TH1D *data_phi_histo = new TH1D ("data phi", "Normalized_data #phi", 60, -3, 3);

	
	//============================== Histos for dz and sigma_dz ==============================
    TH1D *data_dzleaf = new TH1D ("data_dz", "data dz", 400, -10, 10);
    TH1D *data_dzvtxBS_dzErr = new TH1D ("data_dzvtxBS_dzErr", "data d_{z} vtxBS/d_{z} Err", 160, -20, 20);	
    TH1D *data_dxyvtxBS_d0Err = new TH1D ("data_dxyvtxBS_d0Err", "data d_{xy} vtxBS/d_{0} Err", 160, -20, 20);	


	//============================== Histos for d0 and sigma_d0 ==============================
    TH1D *data_d0leaf = new TH1D ("data_d0_calcrun1", "data d_{0} calc run 1", 100, -1, 1); //plot of d0 using leaf value to determine how much data is cut
   
   
	//============================== Histos for pT and sigma_pT ==============================
    TH1D *data_sigmapt_pt = new TH1D ("data_sigmapt_pt", "data #sigma_{p_{T}}/p_{T}", 20, 0, 0.2);

	
	//============================== Histos for chi2n and ValidHits ==============================
    TH1D *data_validhits = new TH1D ("Tracks_vs_validhits", "Tracks vs validhits", 50, 0, 50);
    TH1D *data_chi2n = new TH1D ("Tracks_vs_chi2n", "Tracks vs #chi^{2/ndof}", 50, 0, 5);

	
	//============================== Histos for Multiplicity ==============================
    TH1D *data_normmultiplicity = new TH1D("Normalized_Multiplicity", "Normalized Multiplicity", 200, 0, 200);
	TH1D *data_multiplicity = new TH1D("Multiplicity", "Multiplicity", 200, 0, 200);
	TH1D *data_multiplicity_ZB = new TH1D("Multiplicity_ZB", "Multiplicity_ZB", 200, 0, 200);
	TH1D *data_multiplicity_HM85 = new TH1D("Multiplicity_HM85", "Multiplicity_HM85", 200, 0, 200);
	TH1D *data_vtxzminusvtxz = new TH1D ("vtxzminusvtxz", "vtxzminusvtxz", 800, -4, 4);
    TH1D *data_vtxzposn = new TH1D ("vtxzposn", "vtxzpos^{n}", 400, -20, 20);
	
	
	//============================== Histos for Forward Backward Multiplicity ==============================
	TH2D *data_fb_multiplicity = new TH2D ("FB Multiplicity", "FB Multiplicity", 50, 0, 50, 50, 0, 50);
	TH1D *data_fb_eta = new TH1D("FB eta", "FB eta ", 50, -2.5, 2.5); // for checking if eta cuts are done correctly
	// TODO TProfile *data_nBnF = new TProfile("<N_B> vs N_F", "<N_B> vs N_F", 

	
	//============================== retrieve ROOT file ==============================

	vector<TString> *datafiles = new vector<TString>();
    cout << "Getting list of files..." << endl;

    datafiles = getListOfFiles("/mnt/c/Users/Zongjin-Dell/Desktop/Multiplicity/list_of_files.txt");		//reldir
    cout << "File list stored" << endl;

        TFile *datafile;
        TTree *datatree;

    for(vector<TString>::iterator itlistdatafiles = datafiles->begin() ; itlistdatafiles != datafiles->end(); ++itlistdatafiles)
    {
        cout << "Opening new file " << *itlistdatafiles << endl;

        datafile = new TFile(*itlistdatafiles, "READ");

        cout << "Opened " << *itlistdatafiles << endl;
        datatree = (TTree*)datafile->Get("tree/tree");

	
		//============================== define variables to read TTree ==============================

		vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > > *data_tracks = 0;
		datatree->SetBranchAddress("trP4", &data_tracks);
		
		int ndata_run = 0;
		datatree->SetBranchAddress("Run", &ndata_run);
		
		int ndata_lumi = 0;
		datatree->SetBranchAddress("Lumi", &ndata_lumi);

		vector<double> *fvecdata_vtxx = 0;
		datatree->SetBranchAddress("vtxx", &fvecdata_vtxx);

		vector<double> *fvecdata_vtxy = 0;
		datatree->SetBranchAddress("vtxy", &fvecdata_vtxy);

		vector<double> *fvecdata_vtxz = 0;
		datatree->SetBranchAddress("vtxz", &fvecdata_vtxz);
		
		Double_t fvecdata_BSx = 0;
		datatree->SetBranchAddress("BSx", &fvecdata_BSx);

		Double_t fvecdata_BSy = 0;
		datatree->SetBranchAddress("BSy", &fvecdata_BSy);

		Double_t fvecdata_BSz = 0;
		datatree->SetBranchAddress("BSz", &fvecdata_BSz);

		vector<double> *fvecdata_vtxxBS = 0;
		datatree->SetBranchAddress("vtxxBS", &fvecdata_vtxxBS);

		vector<double> *fvecdata_vtxyBS = 0;
		datatree->SetBranchAddress("vtxyBS", &fvecdata_vtxyBS);

		vector<double> *fvecdata_vtxzBS = 0;
		datatree->SetBranchAddress("vtxzBS", &fvecdata_vtxzBS);

		vector<int> *nvecdata_highpurity = 0;
		datatree->SetBranchAddress("highPurity", &nvecdata_highpurity);

		vector<int> *nvecdata_vtxndof = 0;
		datatree->SetBranchAddress("vtxndof", &nvecdata_vtxndof);
		
		vector<double> *fvecdata_dzvtxBS = 0;
		datatree->SetBranchAddress("dzvtxBS", &fvecdata_dzvtxBS);	

		vector<double> *fvecdata_dxyvtxBS = 0;
		datatree->SetBranchAddress("dxyvtxBS", &fvecdata_dxyvtxBS);	

		vector<double> *fvecdata_dzerr = 0;
		datatree->SetBranchAddress("dzerr", &fvecdata_dzerr);

		vector<double> *fvecdata_d0err = 0;
		datatree->SetBranchAddress("d0err", &fvecdata_d0err);

		vector<double> *fvecdata_pterr = 0;
		datatree->SetBranchAddress("ptErr", &fvecdata_pterr);

		vector<double> *fvecdata_vtxxerr = 0;
		datatree->SetBranchAddress("vtxxErr", &fvecdata_vtxxerr);

		vector<double> *fvecdata_vtxyerr = 0;
		datatree->SetBranchAddress("vtxyErr", &fvecdata_vtxyerr);

		vector<double> *fvecdata_vtxzerr = 0;
		datatree->SetBranchAddress("vtxzErr", &fvecdata_vtxzerr);
		
		vector<double> *fvecdata_vtxxErrBS = 0;
		datatree->SetBranchAddress("vtxxErrBS", &fvecdata_vtxxErrBS);

		vector<double> *fvecdata_vtxyErrBS = 0;
		datatree->SetBranchAddress("vtxyErrBS", &fvecdata_vtxyErrBS);

		vector<double> *fvecdata_vtxzErrBS = 0;
		datatree->SetBranchAddress("vtxzErrBS", &fvecdata_vtxzErrBS);

		vector<int> *nvecdata_validhits = 0;
		datatree->SetBranchAddress("nValidHits", &nvecdata_validhits);

		vector<double> *fvecdata_trackschi2n = 0;
		datatree->SetBranchAddress("chi2n", &fvecdata_trackschi2n);
		
		
		vector<int> *triggerZB = 0;
		datatree->SetBranchAddress("trigger", &triggerZB);
		
		vector<int> *triggerMB = 0;
		datatree->SetBranchAddress("triggerMB", &triggerMB);
		
		vector<int> *triggerHM60 = 0;
		datatree->SetBranchAddress("triggerHM60", &triggerHM60);
		
		vector<int> *triggerHM85 = 0;
		datatree->SetBranchAddress("triggerHM85", &triggerHM85);	

		vector<int> *triggerHM110 = 0;
		datatree->SetBranchAddress("triggerHM110", &triggerHM110);
		
		vector<int> *triggerHM135 = 0;
		datatree->SetBranchAddress("triggerHM135", &triggerHM135);
		
		vector<int> *triggerHM160 = 0;
		datatree->SetBranchAddress("triggerHM160", &triggerHM160);	
		
		
		
		Int_t ndata_totalEvt = (Int_t)datatree->GetEntries();
		cout << "There is a total of " << ndata_totalEvt << " events." << endl;

		
		for (Int_t i = 0; i < ndata_totalEvt; ++i)
		{
			datatree->GetEntry(i);
			if (
			
			//good lumisections for runs
			((ndata_lumi >= 1 && ndata_lumi <=103) && ndata_run == 254987) ||
			(((ndata_lumi >= 1 && ndata_lumi <=120) || (ndata_lumi >= 153 && ndata_lumi <=214) ) && ndata_run == 254989) ||
			((ndata_lumi >= 1 && ndata_lumi <=28) && ndata_run == 254993) ||
			((ndata_lumi >= 4 && ndata_lumi <=414) && ndata_run == 255019) ||
			(((ndata_lumi >= 1 && ndata_lumi <=34) || (ndata_lumi >= 36 && ndata_lumi <=209) || (ndata_lumi >= 306 && ndata_lumi <=325) || (ndata_lumi >= 327 && ndata_lumi <=343) ) && ndata_run == 255029) ||
			(((ndata_lumi >= 31 && ndata_lumi <=101) || (ndata_lumi >= 125 && ndata_lumi <=231) || (ndata_lumi >= 233 && ndata_lumi <=493) || (ndata_lumi >= 586 && ndata_lumi <=1054) || (ndata_lumi >= 1096 && ndata_lumi <=1199)) && ndata_run == 255031)  //255031 lumicuts
			
			)
			
			{
				ndata_totaltrk = data_tracks->size();
				++fdata_evt;

				int vtxdof = 0;

				ndata_numberofvtxx = fvecdata_vtxx->size();
				ndata_numberofvtxy = fvecdata_vtxy->size();
				ndata_numberofvtxz = fvecdata_vtxz->size();

				//============================== Start of Vertex Loop ==============================

				for (int vtxnumber = 0; vtxnumber != ndata_numberofvtxx; ++vtxnumber)
				{
					if((*nvecdata_vtxndof)[vtxdof] > dof_cut)
					{
						if (ndata_numberofvtxz == vtx_number_cut)
						{
							data_vtxzposn->Fill((*fvecdata_vtxz)[vtxnumber]);
							++fdata_numselectedvtxz;
							fdata_multiplicity = 0;
							fdata_backward_multiplicity = 0;
							fdata_forward_multiplicity = 0;
							
							//============================== Start of Trk Loop ==============================

							for (int t = 0; t != ndata_totaltrk; ++t)
							{
								XYZTVector data_vec = (*data_tracks)[t];
								fdata_dzvtxBS_dzErr = ((*fvecdata_dzvtxBS)[t])/((*fvecdata_dzerr)[t]);
								// fdata_dxyBS_d0Err = ((*fvecdata_dxyBS)[t])/((*fvecdata_d0err)[t]);
								fdata_dxyvtxBS_d0Err = ((*fvecdata_dxyvtxBS)[t])/((*fvecdata_d0err)[t]);									
								
								//============================== d0 ==============================

								

								//============================== pT ==============================

								fdata_sigmapt_pt = (((*fvecdata_pterr)[t])/(data_vec.Pt())); //both leaf values

								//============================== Trk cuts ==============================

								if ((*nvecdata_highpurity)[t] == 1)
								{
									if (fabs(data_vec.Eta()) <= eta_cut && fabs(fdata_sigmapt_pt) < ptErr_pt_cut && fabs(fdata_dzvtxBS_dzErr) < dz_dzErr_cut && fabs(fdata_dxyvtxBS_d0Err) < d0_d0Err_cut)
									{
										data_pt_histo->Fill(data_vec.Pt());
										++fdata_trkpt;
									}

									//if (data_vec.Pt() >= pt_cut && fabs(fdata_sigmapt_pt) < ptErr_pt_cut && fabs(fdata_dz_sigmadzrun1 < dz_dzErr_cut))
									//if ((data_vec.Pt() >= pt_cut) && (fabs(fdata_sigmapt_pt) < ptErr_pt_cut))
									if (fabs(fdata_sigmapt_pt) < ptErr_pt_cut && fabs(fdata_dzvtxBS_dzErr) < dz_dzErr_cut && fabs(fdata_dxyvtxBS_d0Err) < d0_d0Err_cut && data_vec.Pt() >= pt_cut)
									{
										data_eta_histo->Fill(data_vec.Eta());
										++fdata_trketa;
									}

									if(fabs(data_vec.Eta()) <= eta_cut && fabs(fdata_dzvtxBS_dzErr) < dz_dzErr_cut && fabs(fdata_dxyvtxBS_d0Err) < d0_d0Err_cut && data_vec.Pt() >= pt_cut)
									{
										data_sigmapt_pt->Fill(fdata_sigmapt_pt);
										++fdata_trkdpt;
									}
									
									if(fabs(data_vec.Eta()) <= eta_cut && fabs(fdata_sigmapt_pt) < ptErr_pt_cut && fabs(fdata_dxyvtxBS_d0Err) < d0_d0Err_cut && data_vec.Pt() >= pt_cut)
									{
										data_dzvtxBS_dzErr->Fill(fdata_dzvtxBS_dzErr);
										++fdata_trkdz;
									}

									if(fabs(data_vec.Eta()) <= eta_cut && fabs(fdata_sigmapt_pt) < ptErr_pt_cut && fabs(fdata_dzvtxBS_dzErr) < dz_dzErr_cut && data_vec.Pt() >= pt_cut)
									{
										data_dxyvtxBS_d0Err->Fill(fdata_dxyvtxBS_d0Err);
										++fdata_trkd0;
									}
									
									if (fabs(data_vec.Eta()) <= eta_cut && fabs(fdata_sigmapt_pt) < ptErr_pt_cut && fabs(fdata_dzvtxBS_dzErr) < dz_dzErr_cut && fabs(fdata_dxyvtxBS_d0Err) < d0_d0Err_cut && data_vec.Pt() >= pt_cut)
									{
										data_phi_histo->Fill(data_vec.Phi());
										++fdata_trkphi;

										data_validhits->Fill((*nvecdata_validhits)[t]);
										++fdata_trkvalidhits;

										data_chi2n->Fill((*fvecdata_trackschi2n)[t]);
										++fdata_trkchi2n;
										
										++fdata_multiplicity;
									}
									
									if (fabs(data_vec.Eta()) >= eta_fb_lower && fabs(data_vec.Eta()) <= eta_fb_upper)
									{
										if (data_vec.Eta() > 0)
										{
											++fdata_forward_multiplicity;
										}
										else
										{
											++fdata_backward_multiplicity;
										}
										
										data_fb_eta->Fill(data_vec.Eta());
										//++fdata_fb_multiplicity_norm;
									}
								}
							}

							
							// ============================== Filling of Multiplicity histograms ==============================
							
							data_multiplicity->Fill(fdata_multiplicity);
							data_normmultiplicity->Fill(fdata_multiplicity);
							++fdata_multiplicity_norm;
							
							// ============================== Filling of Forward Backward Multiplicity histograms ==============================
									
							data_fb_multiplicity->Fill(fdata_forward_multiplicity, fdata_backward_multiplicity);
							
							
							
							//============================== Filling of triggered efficiencies ==============================
															
							if ((*triggerZB)[0] == 1)
							{
								data_multiplicity_ZB->Fill(fdata_multiplicity);
							}	
							
							if ((*triggerHM85)[0] == 1)
							{
								data_multiplicity_HM85->Fill(fdata_multiplicity);
							}	
							//============================== End of Trk Loop ==============================
						}
					}
				}
				//============================== End of Vertex Loop ==============================
			}
		}
		//============================== End of Evt Loop ==============================
	}
	//============================== End of File Loop ==============================
    
	cout << "Before plotting." << endl;

	
	//============================== Create root file, plot histograms ==============================
	
	TFile data_plot("/mnt/c/Users/Zongjin-Dell/Desktop/Multiplicity/fb.root", "recreate");		//reldir


    data_eta_histo->Scale(1./fdata_trketa);
    data_eta_histo->SetMaximum(0.04);
    data_eta_histo->GetXaxis()->SetTitle("#eta");
    data_eta_histo->GetYaxis()->SetTitleOffset(1.3);
    data_eta_histo->GetYaxis()->SetTitle("Fraction of Tracks");
    // data_eta_histo->Draw();
    data_eta_histo->Write();

    //canvas->cd(3);
    data_phi_histo->Scale(1/fdata_trkphi);
    data_phi_histo->SetMinimum(0.014);
    data_phi_histo->SetMaximum(0.024);
    data_phi_histo->GetXaxis()->SetTitle("#phi");
    data_phi_histo->GetYaxis()->SetTitleOffset(1.3);
    data_phi_histo->GetYaxis()->SetTitle("Fraction of Tracks");
    data_phi_histo->Write();

    //canvas->cd(1);
    // gPad->SetLogy();
    data_pt_histo->Scale(1/fdata_trkpt);
    data_pt_histo->GetXaxis()->SetTitle("p_{T}");
    data_pt_histo->GetYaxis()->SetTitleOffset(1.3);
    data_pt_histo->GetYaxis()->SetTitle("Fraction of Tracks");
    data_pt_histo->Write();

    
	data_dzvtxBS_dzErr->Scale(1/fdata_trkdz);
	data_dzvtxBS_dzErr->Write();
	// data_dxyBS_d0Err->Scale(1/fdata_trkdz);
	// data_dxyBS_d0Err->Write();
	data_dxyvtxBS_d0Err->Scale(1/fdata_trkdz);
	data_dxyvtxBS_d0Err->Write();
	
    data_sigmapt_pt->Scale(1/fdata_trkdpt);
    data_sigmapt_pt->SetMinimum(1E-8);
    data_sigmapt_pt->GetXaxis()->SetTitle("#sigma{p_{T}}/p_{T}");
    data_sigmapt_pt->GetYaxis()->SetTitleOffset(1.2);
    data_sigmapt_pt->GetYaxis()->SetTitle("Fraction of Tracks");
    data_sigmapt_pt->Write();

    data_validhits->Scale(1/fdata_trkvalidhits);
    data_validhits->GetXaxis()->SetTitle("# Valid Hits");
    data_validhits->GetYaxis()->SetTitleOffset(1.3);
    data_validhits->GetYaxis()->SetTitle("Fraction of Tracks");
    data_validhits->Write();

    data_chi2n->Scale(1/fdata_trkchi2n);
    data_chi2n->GetXaxis()->SetTitle("#chi^2/dof");
    data_chi2n->GetYaxis()->SetTitleOffset(1.3);
    data_chi2n->GetYaxis()->SetTitle("Fraction of Tracks");
    data_chi2n->Write();

    data_multiplicity->GetXaxis()->SetTitle("Multiplicity");
    data_multiplicity->GetYaxis()->SetTitleOffset(1.3);
    data_multiplicity->GetYaxis()->SetTitle("Number of Tracks");
    data_multiplicity->Write();
	
	
    data_multiplicity_ZB->GetXaxis()->SetTitle("Multiplicity_ZB");
    data_multiplicity_ZB->GetYaxis()->SetTitleOffset(1.3);
    data_multiplicity_ZB->GetYaxis()->SetTitle("Number of Tracks");
    data_multiplicity_ZB->Write();

    data_multiplicity_HM85->GetXaxis()->SetTitle("Multiplicity_HM85");
    data_multiplicity_HM85->GetYaxis()->SetTitleOffset(1.3);
    data_multiplicity_HM85->GetYaxis()->SetTitle("Number of Tracks");
    data_multiplicity_HM85->Write();	
	
	TH1D *data_multiplicity_HM85_efficiency = (TH1D*)data_multiplicity_HM85->Clone("data_multiplicity_HM85_efficiency");
	data_multiplicity_HM85_efficiency->Divide(data_multiplicity);
	data_multiplicity_HM85_efficiency->Write();
	
	TH1D *data_multiplicity_ZB_efficiency = (TH1D*)data_multiplicity_ZB->Clone("data_multiplicity_ZB_efficiency");
	data_multiplicity_ZB_efficiency->Divide(data_multiplicity);
	data_multiplicity_ZB_efficiency->Write();


    data_normmultiplicity->Scale(1/fdata_multiplicity_norm);
    data_normmultiplicity->GetXaxis()->SetTitle("Multiplicity");
    data_normmultiplicity->GetYaxis()->SetTitleOffset(1.3);
    data_normmultiplicity->GetYaxis()->SetTitle("Fraction of Tracks");
    data_normmultiplicity->Write();
	
	
	// Forward-backward Multiplicity Plots
	data_fb_multiplicity->GetXaxis()->SetTitle("Forward Multiplicity");
    data_fb_multiplicity->GetYaxis()->SetTitleOffset(1.3);
    data_fb_multiplicity->GetYaxis()->SetTitle("Backward Multiplicity");
    data_fb_multiplicity->Write();
	
	// For checking whether eta cuts are done properly
    data_fb_eta->GetXaxis()->SetTitle("#eta");
    data_fb_eta->GetYaxis()->SetTitleOffset(1.3);
    data_fb_eta->GetYaxis()->SetTitle("Number of Tracks");
    data_fb_eta->Write();
	

    data_plot.Write();

}

#ifndef __CINT__
int main () { forward_backward(); return 0; }  // Main program when run stand-alone
#endif
