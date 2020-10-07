
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



using namespace RooFit;
using namespace std;

static const int colors_array[]={2, 3, 4, 7, 30, 6, 9,46,36,8,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24};


RooRealVar MuMu_mass("MuMu_mass","MuMu_mass",70,110);

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
void Draw_th1d(TH1D* th1d_input, TString variable_name)
{
	TCanvas* c=new TCanvas();
	c->cd();
	gStyle->SetOptStat(0);
	th1d_input->SetMarkerStyle(kFullCircle);
	th1d_input->SetMarkerColor(kRed);
	th1d_input->SetLineColor(kRed);

	th1d_input->SetMaximum(91.5);
	th1d_input->SetMinimum(90);
	th1d_input->GetXaxis()->SetTitle(variable_name.Data());
	th1d_input->GetXaxis()->SetTitleOffset(1.2);	
	th1d_input->GetYaxis()->SetTitle("Mass mean (GeV)");	
	th1d_input->Draw();
	tlxg->DrawLatexNDC(0.2,0.8,Form("%s",GT.Data()));
	c->Print(Form("%s/fitResultPlot/mass_VS_%s.pdf",GT.Data(), variable_name.Data()));
}
//110X_dataRun2_v13  92X_dataRun2_Prompt_v11
const static int variables_number=8;
const TString tstring_variables_name[variables_number]={"CosThetaCS","DeltaEta","EtaMinus","EtaPlus","PhiCS","PhiMinus","PhiPlus","Pt"};
void Fitting_GetMassmeanVSvariables(TString inputfile_name, TString GlobalTagName)//TString tstring_inputfilename
{
	cout<<"Stage 0"<<endl;
	//================================
	
	TH2D* th2d_mass_variables[variables_number];
	
	//===============================
	cout<<"Stage 0.5"<<endl;
	TFile *inputfile=TFile::Open(inputfile_name.Data());
	TDirectoryFile * tdirectory=(TDirectoryFile *)inputfile->Get("myanalysis");	
	for(int i=0; i<variables_number; i++)
	{
		TString th2d_name=Form("th2d_mass_%s",tstring_variables_name[i].Data());
		th2d_mass_variables[i]=(TH2D*) tdirectory->Get(th2d_name);
		cout<<th2d_mass_variables[i]->GetEntries()<<endl;
	}

	//===============================
	cout<<"Stage 1"<<endl;
	GT=GlobalTagName;
	gSystem->Exec(Form("mkdir -p %s",GT.Data()));
	gSystem->Exec(Form("mkdir -p %s/fitResultPlot", GT.Data()));
	TFile* outpufile=TFile::Open(Form("%s/output.root",GT.Data()),"recreate");
	TH1D* th1d_variables_meanmass[variables_number];
	TH1D* th1d_variables_entries[variables_number];
	const int variables_rebin[variables_number]={1, 1, 1  , 1  , 1  , 1   , 1   , 5};
	for(int i=0; i<variables_number; i++)
	{	
		TString th1d_name=Form("th1d_meanmass_%s",tstring_variables_name[i].Data());

		th2d_mass_variables[i]->RebinY(variables_rebin[i]);
		th1d_variables_meanmass[i]=th2d_mass_variables[i]->ProjectionY(th1d_name,1,1,"d");
		for(int j=0; j<th1d_variables_meanmass[i]->GetNbinsX(); j++)
		{
			cout<<"th1d_variables_meanmass[i]->GetNbinsX()="<<th1d_variables_meanmass[i]->GetNbinsX()<<endl;
			cout<<"th2d_mass_variables[i]->GetNbinsY()="<<th2d_mass_variables[i]->GetNbinsY()<<endl;
			th1d_variables_meanmass[i]->SetBinContent(j,0);
			th1d_variables_meanmass[i]->SetBinError(j,0);

			TString th1d_mass_temp_name=Form("th1d_mass_%s_%d",tstring_variables_name[i].Data(),j);
			TH1D* th1d_i=th2d_mass_variables[i]->ProjectionX(th1d_mass_temp_name,j,j,"d");
			th1d_i->Write(th1d_mass_temp_name);
			TString s_cut=Form("");
			TString s_name=Form("%s_%d",tstring_variables_name[i].Data(),j);
			RooPlot *massframe = MuMu_mass.frame();
			RooDataHist dh_temp("dh_temp","dh_temp",MuMu_mass,Import(*th1d_i));
			FitOut fitR=ZMassBinFit(dh_temp,s_cut,s_name);
			th1d_variables_meanmass[i]->SetBinContent(j,fitR.mean);
			th1d_variables_meanmass[i]->SetBinError(j,fitR.mean_err);
		}
		
		Draw_th1d(th1d_variables_meanmass[i], tstring_variables_name[i]);
		th1d_variables_meanmass[i]->Write(th1d_name);

		TString th1d_name_entries=Form("th1d_entries_%s",tstring_variables_name[i].Data());
		th1d_variables_entries[i]=th2d_mass_variables[i]->ProjectionY(th1d_name_entries,0,-1,"d");
		th1d_variables_entries[i]->GetXaxis()->SetTitle(tstring_variables_name[i].Data());
		th1d_variables_entries[i]->GetYaxis()->SetTitle("Entry");
		th1d_variables_entries[i]->Write(th1d_name_entries);
	}

	outpufile->Write();
	outpufile->Close();
	delete outpufile;

}

const static int max_file_number=10;
void Draw_TH1D_forMultiRootFiles(int file_number, TString file_names[max_file_number], TString label_names[max_file_number], TString th1d_name, TString output_name)
{
	TH1D* th1d_input[max_file_number];
	TFile* file_input[max_file_number];
	for(int idx_file=0; idx_file<file_number; idx_file++)
	{
		file_input[idx_file]=TFile::Open(file_names[idx_file]);
		th1d_input[idx_file]=(TH1D*) file_input[idx_file]->Get(th1d_name);
	}

	TCanvas* c=new TCanvas();
	TLegend* lg=new TLegend(0.2,0.7,0.5,0.95);
	c->cd();
	gStyle->SetOptStat(0);
	th1d_input[0]->SetTitle("");
	for(int idx_file=0; idx_file<file_number; idx_file++)
	{
		th1d_input[idx_file]->SetMarkerColor(colors_array[idx_file]);
		th1d_input[idx_file]->SetLineColor(colors_array[idx_file]);
		th1d_input[idx_file]->Draw("same");
		lg->AddEntry(th1d_input[idx_file],label_names[idx_file]);
	}
	lg->Draw("same");
	c->SaveAs(output_name);

}

void Draw_fitResults_forMultiGT(int GT_number, TString GT_name[max_file_number])
{
	TH1D* th1d_variables_meanmass_fromfiles[variables_number];
	for(int idx_GT=0; idx_GT<GT_number; idx_GT++)
	{

	}
}

void Macro_Fitting_DiMuonValidation()
{
	Fitting_GetMassmeanVSvariables("./v6_TrackRefitter_110X_dataRun2_v13.root", "110X_dataRun2_v13");//"./merged_v4.root"
	Fitting_GetMassmeanVSvariables("./v6_TrackRefitter_92X_dataRun2_Prompt_v11.root", "92X_dataRun2_Prompt_v11");

	TString GT_names[max_file_number];
	int files_number=2;
	GT_names[0]="110X_dataRun2_v13";
	GT_names[1]="92X_dataRun2_Prompt_v11";

	TString file_names[max_file_number];
	file_names[0]=Form("./%s/output.root",GT_names[0].Data());
	file_names[1]=Form("./%s/output.root",GT_names[1].Data());

	for(int idx_variable=0; idx_variable<variables_number; idx_variable++)
	{
		TString th1d_name=Form("th1d_meanmass_%s",tstring_variables_name[idx_variable].Data());
		Draw_TH1D_forMultiRootFiles(files_number, file_names, GT_names, th1d_name, Form("meanmass_%s_GTs.pdf",tstring_variables_name[idx_variable].Data()));
		TString th1d_name_entries=Form("th1d_entries_%s",tstring_variables_name[idx_variable].Data());
		Draw_TH1D_forMultiRootFiles(files_number, file_names, GT_names, th1d_name_entries, Form("entries_%s_GTs.pdf",tstring_variables_name[idx_variable].Data()));
		
	}
	
	return 0;

}
int main()
{
	Macro_Fitting_DiMuonValidation();
	
	return 0;
}
