
#include <iostream>
#include <TString.h>
#include <TCanvas.h>
#include "TBranch.h"
#include "TTree.h"
#include "TSystem.h"
#include "TH1F.h"
#include "TH1D.h"
#include "TF1.h"
#include "TMath.h"
#include "TLegend.h"
#include "TLatex.h"
#include <TH2F.h>
#include <TCut.h>
#include <TStyle.h>

#include <math.h>
#include <string>
#include <TChain.h>
#include <TFile.h>

#include "TH1.h"
#include "TROOT.h"
#include "TPad.h"
#include "TRandom.h"
#include "THStack.h"
#include "TH2.h"
#include "TLatex.h"
#include "TPaveText.h"
#include "TFrame.h"
#include "TAxis.h"
#include <iostream>
#include <fstream>
#include <iomanip>
#include <TCut.h>

#include <RooFit.h>
#include <RooRealVar.h>
#include <RooDataHist.h>
#include <RooDataSet.h>
#include <RooHistPdf.h>
#include <RooPlot.h>
#include <RooAddPdf.h>
#include <RooGlobalFunc.h>
#include <RooGaussian.h>
#include <RooPolynomial.h>
#include <RooGenericPdf.h>
#include <RooFormulaVar.h>
#include <RooChebychev.h>
#include <RooExponential.h>
#include <RooFitResult.h>
#include <RooWorkspace.h>
#include <RooMCStudy.h>
#include "RooHist.h"
#include "RooConstVar.h"
#include "RooMsgService.h"
#include "RooVoigtian.h"

using namespace RooFit;
using namespace std;

const static int eta_bins_number=10;
const double eta_min=-5;
const double eta_max=5;

const static int pt_bins_number=20;
const double pt_min=0;
const double pt_max=10;

const static int phi_bins_number=24;
const double phi_min=-M_PI;
const double phi_max=M_PI;



RooRealVar MuMu_mass("MuMu_mass","MuMu_mass",70,110);
RooRealVar SMuMu_pt("SMuMu_pt","SMuMu_pt",0);
RooRealVar SMuMu_phi("SMuMu_phi","SMuMu_phi",0);
RooRealVar SMuMu_eta("SMuMu_eta","SMuMu_eta",0);

static TString GT="";

TLatex *tlxg=new TLatex();

class FitOut{
	public:
	double mean;
	double mean_err;
	double sigma;
	double sigma_err;

	FitOut(double a, double b, double c, double d): mean(a), mean_err(b), sigma(c),sigma_err(d){}

};

/*--------------------------------------------------------------------*/
void makeNicePlotStyle(RooPlot* plot)
/*--------------------------------------------------------------------*/
{ 
  plot->GetXaxis()->CenterTitle(true);
  plot->GetYaxis()->CenterTitle(true);
  plot->GetXaxis()->SetTitleFont(42); 
  plot->GetYaxis()->SetTitleFont(42);  
  plot->GetXaxis()->SetTitleSize(0.05);
  plot->GetYaxis()->SetTitleSize(0.05);
  plot->GetXaxis()->SetTitleOffset(0.9);
  plot->GetYaxis()->SetTitleOffset(1.3);
  plot->GetXaxis()->SetLabelFont(42);
  plot->GetYaxis()->SetLabelFont(42);
  plot->GetYaxis()->SetLabelSize(.05);
  plot->GetXaxis()->SetLabelSize(.05);
}



// FitOut ZMassUnBinFit(){
//}


FitOut ZMassBinFit(RooDataHist datahist,TString s_cut="",TString s_name="nocut")
// FitOut ZMassUnBinFit(RooDataSet datahist,TString s_cut="",TString s_name="nocut")
{

	// cout<<"sumEntries() = "<<datahist.sumEntries()<<endl;
	// RooPlot *frame = SMuMu_mass.frame();
	RooPlot *frame = MuMu_mass.frame();
	datahist.plotOn(frame);
	// RooPlot* massframe=new RooPlot("massframe","Dmass",Dmass,DsDataFitRangeLow,DsDataFitRangeHigh,nbin_DmassDraw);
	//RDH.plotOn(frame);
	//frame->Draw();

	RooMsgService::instance().setGlobalKillBelow(RooFit::FATAL);
	TCanvas *c1 = new TCanvas();
	c1->Clear();

	c1->SetLeftMargin(0.15);
	c1->SetRightMargin(0.10);

	RooRealVar mean("#mu","mean",90.0, 60.0, 120.0);
	RooRealVar width("width","width",5.0, 0.0, 120.0);
	RooRealVar sigma("#sigma","sigma",5.0, 0.0, 120.0);
	//RooBreitWigner gauss("gauss","gauss",x,mean,sigma);
	RooVoigtian voigt("voigt","voigt",MuMu_mass,mean,width,sigma);

	RooRealVar lambda("#lambda", "slope", -0.01, -100., 1.);
	RooExponential expo("expo", "expo", MuMu_mass, lambda);

	RooRealVar b("N_{b}", "Number of background events",0, datahist.sumEntries()/10);
	RooRealVar s("N_{s}", "Number of signal events", 0, datahist.sumEntries());

	RooAddPdf fullModel("fullModel", "Signal + Background Model", RooArgList(voigt, expo), RooArgList(s, b));

	// auto r = fullModel.fitTo(datahist, RooFit::PrintLevel(-1), RooFit::Save(),RooFit::Range(70.,110.));
	auto r = fullModel.fitTo(datahist, RooFit::PrintLevel(-1), RooFit::Save(), NumCPU(20));
	// r = fullModel.fitTo(datahist, RooFit::PrintLevel(-1), RooFit::Save());
	// r = fullModel.fitTo(datahist, RooFit::PrintLevel(-1), RooFit::Save());

	fullModel.plotOn(frame,RooFit::LineColor(kRed));
	fullModel.plotOn(frame,RooFit::Components(expo),RooFit::LineStyle(kDashed)) ; //Other option
	fullModel.paramOn(frame,RooFit::Layout(0.65,0.90,0.90));
	frame->getAttText()->SetTextSize(0.03);

	makeNicePlotStyle(frame);

	// Redraw data on top and print / store everything
	datahist.plotOn(frame);
	frame->GetYaxis()->SetTitle("n. of events");
	frame->GetXaxis()->SetTitle("M_{#mu#mu}");
	frame->SetTitle(GT.Data());
	// TString histName = hist->GetName();
	// frame->SetName("frame"+histName);
	// frame->SetTitle(hist->GetTitle());
	frame->Draw();


	TLatex *tlx=new TLatex();
	// tlx->DrawLatexNDC(0.2,0.8,GT.Data());
	tlx->DrawLatexNDC(0.2,0.8,s_cut.Data());

	gSystem->Exec(Form("mkdir -p %s/fitplot", GT.Data()));
	c1->Print(Form("%s/fitplot/%s.pdf",GT.Data(),s_name.Data()));
	// double mean_1=1;
	// double error_1=2;

	FitOut fitRes(mean.getValV(), mean.getError(), sigma.getValV(), sigma.getError());
	return fitRes;

}
//===========================================

TH1D* FuncL0_GetTH1D_MassSpectrum(TH3D* th3d_mass_pt_phi[eta_bins_number], int idx_eta_min, int idx_eta_max, int idx_pt_min, int idx_pt_max, int idx_phi_min, int idx_phi_max)
{
	cout<<"Processing FuncL0_GetTH1D_MassSpectrum"<<endl;
	if(idx_eta_min>idx_eta_max){idx_eta_min=1; idx_eta_max=eta_bins_number;}
	
	TH3D* th3d_mass_pt_phi_merged=(TH3D*) th3d_mass_pt_phi[idx_eta_min-1]->Clone("th3d_mass_pt_phi_merged");
	for(int idx_eta=idx_eta_min-1; idx_eta<idx_eta_max-1; idx_eta++)
	{
		th3d_mass_pt_phi_merged->Add(th3d_mass_pt_phi[idx_eta],1);
		cout<<"idx_eta="<<idx_eta<<"   idx_eta_min-1="<<idx_eta_min-1<<"   idx_eta_max-1="<<idx_eta_max-1<<endl;
	}cout<<"   idx_eta_min-1="<<idx_eta_min-1<<"   idx_eta_max-1="<<idx_eta_max-1<<endl;
	TH1D* th1d_mass;
	th1d_mass=th3d_mass_pt_phi_merged->ProjectionX("th1d_mass",idx_pt_min, idx_pt_max, idx_phi_min, idx_phi_max,"d");
	return th1d_mass;
}

TString FuncGeneral_GetCutString(TString cut_variable_name,int bins_number, int idx_bin_min, int idx_bin_max, double val_min, double val_max)//
{
	double val_1=val_min+(val_max-val_min)*(idx_bin_min-1)/bins_number;
	double val_2=val_min+(val_max-val_min)*(idx_bin_max-0)/bins_number;

	TString tstring_CutName=Form("%.1f < %s < %.1f", val_1, cut_variable_name.Data(), val_2);
	if(cut_variable_name=="#phi"){tstring_CutName=Form("%.1f #phi< %s < %.1f #phi", val_1/TMath::Pi(), cut_variable_name.Data(), val_2/TMath::Pi());}
	cout<<"FuncGeneral_GetCutString"<<"  "<<tstring_CutName<<endl;
	return tstring_CutName;
}


TString FuncL1_GetCutString(int idx_eta_min, int idx_eta_max, int idx_pt_min, int idx_pt_max, int idx_phi_min, int idx_phi_max)
{
	TString tstring_Cuteta=FuncGeneral_GetCutString("#eta",eta_bins_number, idx_eta_min, idx_eta_max, eta_min, eta_max);
	TString tstring_Cutpt =FuncGeneral_GetCutString("p_{T}",  pt_bins_number , idx_pt_min , idx_pt_max , pt_min , pt_max );
	TString tstring_Cutphi=FuncGeneral_GetCutString("#phi",phi_bins_number, idx_phi_min, idx_phi_max, phi_min, phi_max);

	//TString tstring_CutName=Form("%s \n %s \n %s ", tstring_Cutphi.Data(), tstring_Cutpt.Data(), tstring_Cuteta.Data());//TString tstring_CutName=Form("#splitine{%s}{#splitine{%s}{%s}}", tstring_Cutphi.Data(), tstring_Cutpt.Data(), tstring_Cuteta.Data());
	TString tstring_CutName=Form("#splitline{%s}{#splitline{%s}{%s}}", tstring_Cuteta.Data(), tstring_Cutpt.Data(), tstring_Cutphi.Data());
	return tstring_CutName;
}


void Macro_Fitting_DiMuonValidation(TString inputfile_name="./merged_v4.root", TString GlobalTagName="103X_dataRun2_Prompt_v3")//TString tstring_inputfilename
{
	


	TH3D* th3d_mass_pt_phi[eta_bins_number];

	GT=GlobalTagName;
	gSystem->Exec(Form("mkdir -p %s",GT.Data()));
	gSystem->Exec(Form("mkdir -p %s/fitResultPlot", GT.Data()));
	
	TFile *inputfile=TFile::Open(inputfile_name.Data());
	TDirectoryFile * tdirectory=(TDirectoryFile *)inputfile->Get("myanalysis");	
	for(int idx_eta=0; idx_eta<eta_bins_number; idx_eta++)
	{
		cout<<"idx_eta="<<idx_eta<<endl;
		th3d_mass_pt_phi[idx_eta]=(TH3D*) tdirectory->Get(Form("th3d_mass_pt_phi_eta%d",idx_eta));
		cout<<Form("th3d_mass_pt_phi_eta%d",idx_eta)<<"->Entries()="<<th3d_mass_pt_phi[idx_eta]->GetEntries()<<endl;
	}


	TFile* outpufile=TFile::Open(Form("%s/output.root",GT.Data()),"recreate");
	//==================================================
	const int eta_cut_number=10;
	const int phi_cut_number=24;
	const int pt_cut_number=10;

	//==============Produce mass v.s. eta====================================
	cout<<"Start Produce mass v.s. eta"<<endl;
	outpufile->cd();
	TH1D *th1d_Mass_VS_eta=new TH1D("th1d_Mass_VS_eta","th1d_Mass_VS_eta",eta_cut_number, eta_min, eta_max);
	for(int idx_eta=1; idx_eta<eta_cut_number+1; idx_eta++)
	{
		int	idx_bin_eta_min=(idx_eta-1)*eta_bins_number/eta_cut_number +1; cout<<"idx_bin_eta_min="<<idx_bin_eta_min<<endl;
		int idx_bin_eta_max=(idx_eta-0)*eta_bins_number/eta_cut_number +0; cout<<"idx_bin_eta_max="<<idx_bin_eta_max<<endl;
		TH1D* th1d_i=FuncL0_GetTH1D_MassSpectrum(th3d_mass_pt_phi,idx_bin_eta_min,idx_bin_eta_max, 1,pt_bins_number,1,phi_bins_number);
		TString s_cut=FuncL1_GetCutString(idx_bin_eta_min,idx_bin_eta_max, 1,pt_bins_number,1,phi_bins_number);
		TString s_name=Form("eta%d",idx_eta);	

		RooPlot *massframe = MuMu_mass.frame();
		RooDataHist dh_temp("dh_temp","dh_temp",MuMu_mass,Import(*th1d_i));
		FitOut fitR=ZMassBinFit(dh_temp,s_cut,s_name);
		th1d_Mass_VS_eta->SetBinContent(idx_eta,fitR.mean);
		th1d_Mass_VS_eta->SetBinError(idx_eta,fitR.mean_err);
	}
	
	TCanvas* c_mass_vs_eta=new TCanvas();
	c_mass_vs_eta->cd();
	gStyle->SetOptStat(0);
	th1d_Mass_VS_eta->SetMarkerStyle(kFullCircle);
	th1d_Mass_VS_eta->SetMarkerColor(kRed);
	th1d_Mass_VS_eta->SetLineColor(kRed);


	th1d_Mass_VS_eta->SetMaximum(91.5);
	th1d_Mass_VS_eta->SetMinimum(90);
	th1d_Mass_VS_eta->GetXaxis()->SetTitle("#eta");	
	th1d_Mass_VS_eta->GetYaxis()->SetTitle("Mass mean (GeV)");	
	th1d_Mass_VS_eta->Draw();
	tlxg->DrawLatexNDC(0.2,0.8,Form("%s",GT.Data()));
	c_mass_vs_eta->Print(Form("%s/fitResultPlot/mass_VS_eta.pdf",GT.Data()));
	th1d_Mass_VS_eta->Write();
	cout<<"End Produce mass v.s. eta"<<endl;
	//======================================================================
	
	//==============Produce mass v.s. pt====================================
	cout<<"Start Produce mass v.s. pt"<<endl;
	outpufile->cd();
	TH1D *th1d_Mass_VS_pt=new TH1D("th1d_Mass_VS_pt","th1d_Mass_VS_pt",pt_cut_number, pt_min, pt_max);
	for(int idx_pt=1; idx_pt<pt_cut_number+1; idx_pt++)
	{
		int	idx_bin_pt_min=(idx_pt-1)*pt_bins_number/pt_cut_number +1;
		int idx_bin_pt_max=(idx_pt-0)*pt_bins_number/pt_cut_number +0;
		TH1D* th1d_i=FuncL0_GetTH1D_MassSpectrum(th3d_mass_pt_phi,1,eta_bins_number,idx_bin_pt_min,idx_bin_pt_max,1,phi_bins_number);
		TString s_cut=FuncL1_GetCutString(1,eta_bins_number,idx_bin_pt_min,idx_bin_pt_max,1,phi_bins_number);
		TString s_name=Form("pt%d",idx_pt);	

		RooPlot *massframe = MuMu_mass.frame();
		RooDataHist dh_temp("dh_temp","dh_temp",MuMu_mass,Import(*th1d_i));
		FitOut fitR=ZMassBinFit(dh_temp,s_cut,s_name);
		th1d_Mass_VS_pt->SetBinContent(idx_pt,fitR.mean);
		th1d_Mass_VS_pt->SetBinError(idx_pt,fitR.mean_err);
	}
	
	TCanvas* c_mass_vs_pt=new TCanvas();
	c_mass_vs_pt->cd();
	gStyle->SetOptStat(0);
	th1d_Mass_VS_pt->SetMarkerStyle(kFullCircle);
	th1d_Mass_VS_pt->SetMarkerColor(kRed);
	th1d_Mass_VS_pt->SetLineColor(kRed);


	th1d_Mass_VS_pt->SetMaximum(91.5);
	th1d_Mass_VS_pt->SetMinimum(90);
	th1d_Mass_VS_pt->GetXaxis()->SetTitle("p_{T}");	
	th1d_Mass_VS_pt->GetYaxis()->SetTitle("Mass mean (GeV)");	
	th1d_Mass_VS_pt->Draw();
	tlxg->DrawLatexNDC(0.2,0.8,Form("%s",GT.Data()));
	c_mass_vs_pt->Print(Form("%s/fitResultPlot/mass_VS_pt.pdf",GT.Data()));
	th1d_Mass_VS_pt->Write();
	cout<<"End Produce mass v.s. pt"<<endl;
	//=======================================================================
	//==============Produce mass v.s. phi====================================
	cout<<"Start Produce mass v.s. phi"<<endl;
	outpufile->cd();
	TH1D *th1d_Mass_VS_phi=new TH1D("th1d_Mass_VS_phi","th1d_Mass_VS_phi",phi_cut_number, phi_min, phi_max);
	for(int idx_phi=1; idx_phi<phi_cut_number+1; idx_phi++)
	{
		int	idx_bin_phi_min=(idx_phi-1)*phi_bins_number/phi_cut_number +1;
		int idx_bin_phi_max=(idx_phi-0)*phi_bins_number/phi_cut_number +0;
		TH1D* th1d_i=FuncL0_GetTH1D_MassSpectrum(th3d_mass_pt_phi,1,eta_bins_number,1,pt_bins_number,idx_bin_phi_min,idx_bin_phi_max);
		TString s_cut=FuncL1_GetCutString(1,eta_bins_number,1,pt_bins_number,idx_bin_phi_min,idx_bin_phi_max);
		TString s_name=Form("phi%d",idx_phi);	

		RooPlot *massframe = MuMu_mass.frame();
		RooDataHist dh_temp("dh_temp","dh_temp",MuMu_mass,Import(*th1d_i));
		FitOut fitR=ZMassBinFit(dh_temp,s_cut,s_name);
		th1d_Mass_VS_phi->SetBinContent(idx_phi,fitR.mean);
		th1d_Mass_VS_phi->SetBinError(idx_phi,fitR.mean_err);
	}
	
	TCanvas* c_mass_vs_phi=new TCanvas();
	c_mass_vs_phi->cd();
	gStyle->SetOptStat(0);
	th1d_Mass_VS_phi->SetMarkerStyle(kFullCircle);
	th1d_Mass_VS_phi->SetMarkerColor(kRed);
	th1d_Mass_VS_phi->SetLineColor(kRed);


	th1d_Mass_VS_phi->SetMaximum(91.5);
	th1d_Mass_VS_phi->SetMinimum(90);
	th1d_Mass_VS_phi->GetXaxis()->SetTitle("#phi");	
	th1d_Mass_VS_phi->GetYaxis()->SetTitle("Mass mean (GeV)");	
	th1d_Mass_VS_phi->Draw();
	tlxg->DrawLatexNDC(0.2,0.8,Form("%s",GT.Data()));
	c_mass_vs_phi->Print(Form("%s/fitResultPlot/mass_VS_phi.pdf",GT.Data()));
	th1d_Mass_VS_phi->Write();
	cout<<"End Produce mass v.s. phi"<<endl;
	//=======================================================================
	/**/
	


}

int main()
{
	Macro_Fitting_DiMuonValidation("./merged_v4.root", "103X_dataRun2_Prompt_v3_2");//"./merged_v4.root"
	return 0;
}
