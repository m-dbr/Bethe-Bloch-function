#include <TH1F.h>
#include <TH2F.h>
#include <TFile.h>
#include <TSystem.h>
#include <TCanvas.h>
#include <TString.h>
#include <TStyle.h>
#include <TF1.h>
#include <TMultiGraph.h>
#include <TGraphErrors.h>
#include <TMath.h>

#include <string>
#include <time.h>
using namespace std;

TCanvas *cnv[100];
Float_t gPmin;
Float_t gPmax;
TF1* gBethe_Bloch_p;
TF1* gBethe_Bloch_pi;
TF1* gBethe_Bloch_pi_from_proton;
TF1* gBethe_Bloch_p_from_pion;
TF1* gBethe_Bloch_k;
TF1* GetGaussFitPol(TH1D* h, Float_t p, const Double_t* params);
TF1* GetGaussFitPol_pi(TH1D* h, Float_t p, const Double_t* params);

void TPstyle();
void DrawSlice(TH1D* h,TCanvas* canv);
void DyVsY2(Int_t nSlicesInit, TH2F* h2, string particle, float &par1, float &par2, float &par3, float &par4);

Double_t BetheBloch(Double_t *x, Double_t *par);

TGraphErrors* mean_pi = 0;
TGraphErrors* mean_p  = 0;

Double_t fcn1(Int_t &npar, Double_t *gin, Double_t *par, Int_t iflag);
Double_t fcn2(Int_t &npar, Double_t *gin, Double_t *par, Int_t iflag);
void fcn(Int_t &npar, Double_t *gin, Double_t &f, Double_t *par, Int_t iflag);
Double_t BetheBloch1(Double_t x, Double_t *par);
Double_t BetheBloch2(Double_t x, Double_t *par);


//void FitSlices(Int_t slic, Int_t runid, Int_t apver, Int_t nSlices, Int_t saveCalibParams);
/*
//________________________________________________________________________________________________________
int main(int argc, char** argv)
{

  FitSlices(1, 1234, 4, 15, 0);

}
*/

//___________MAIN FUNCTION_______________________________________________________________________//
void FitSlices(Int_t slic=1, Int_t runid=1234, Int_t apver=7, Int_t nSlices=15, Int_t saveCalibParams=0)
{

  gStyle->SetOptStat(0);
  gStyle->SetOptFit(1);

  gPmin=0.7;
  gPmax=1.6;

  Double_t X0=1.44; 
  Double_t X1=4; 
  Double_t miu=3.77; 
  Double_t s=1.53;
  
  ///// pion dEdx
  Float_t pmin = 0.1;
  Float_t pmax = 50.0;

  // PROTON
  Double_t mass_prot = 0.9382720799426998;
  Double_t par_p[] = {X0, X1, miu, s, mass_prot};

  gBethe_Bloch_p = new TF1("Bethe Bloch", BetheBloch, -0.5, 3, 5); 
  gBethe_Bloch_p->SetParameters(par_p);
  gBethe_Bloch_p->SetParNames("X0", "X1", "miu", "s", "mass");
  gBethe_Bloch_p->FixParameter(4, mass_prot);

//   PION
  Double_t mass_pi = 0.13957;
  Double_t par_pi[] = {X0, X1, miu, s, mass_pi};

  gBethe_Bloch_pi = new TF1("Bethe Bloch pi", BetheBloch, -0.5, 3, 5);
  gBethe_Bloch_pi->SetParameters(par_pi);
  gBethe_Bloch_pi->SetParNames("X0", "X1", "miu", "s", "mass");
  gBethe_Bloch_pi->FixParameter(4, mass_pi);

  // KAON
  Double_t mass_kaon = 0.493677;
//
//  Double_t par_kaon[] = {X0, X1, miu, s, mass_kaon};
//
//  gBethe_Bloch_k = new TF1("Bethe Bloch pi", BetheBloch, -0.5, 3, 5);
//  gBethe_Bloch_k->SetParameters(par_kaon);
//  gBethe_Bloch_k->SetParNames("X0","X1","miu","s", "mass_kaon");
//  gBethe_Bloch_k->FixParameter(4, mass_kaon);

// Take PROTONS from root file
  TString dir = "LambAnaModule";
  
  TString fineName=Form("d0ana_run%d_p5_ap%d.root", runid, apver);
  
  TFile* InputFile = TFile::Open(Form("../root_files/%s",fineName.Data()));
  if(!InputFile)return;
  
  TH2F* hr = (TH2F*)InputFile->Get(Form("%s/hProtons_log",dir.Data()));
  TH2F* hr2 = (TH2F*)InputFile->Get(Form("%s/hPbar_log",dir.Data()));

  hr->Add(hr2);

// Take PIONS from root file
  TH2F* hr_pi = (TH2F*)InputFile->Get(Form("%s/hPiPlus_log",dir.Data()));
  TH2F* hr_pi2 = (TH2F*)InputFile->Get(Form("%s/hPiMinus_log",dir.Data()));

  TString dir2 = "K0AnaModule";

  TH2F* hr_pi3 = (TH2F*)InputFile->Get(Form("%s/hPiPlus_log",dir2.Data()));
  TH2F* hr_pi4 = (TH2F*)InputFile->Get(Form("%s/hPiMinus_log",dir2.Data()));

  hr_pi->Add(hr_pi2);
  hr_pi->Add(hr_pi3);
  hr_pi->Add(hr_pi4);

  //TH2F* hrplot = (TH2F*)hr->Clone("hrplot");
  TH2F* hrplot = (TH2F*)hr_pi->Clone("hrplot");

  float par1 = 0;
  float par2 = 0;
  float par3 = 0;
  float par4 = 0;

// Slicing function for pions/protons
  string particle = "pion";
  gPmin=0.5;
  gPmax=1.2;
  DyVsY2(nSlices,hr_pi, particle, par1, par2, par3, par4);

  cout << "!!!!!!!!!!!!!!" << endl;
  cout << par1 << endl;
  cout << par2 << endl;
  cout << par3 << endl;
  cout << par4 << endl;
  cout << "!!!!!!!!!!!!!!" << endl;

  if ( slic == 0 ) {
      string particle2 = "proton";
      gPmin=1.0;
      gPmax=1.4;
      DyVsY2(nSlices,hr, particle2, par1, par2, par3, par4);
      cout << "!!!!!!!!!!!!!!" << endl;
      cout << par1 << endl;
      cout << par2 << endl;
      cout << par3 << endl;
      cout << par4 << endl;
      cout << "!!!!!!!!!!!!!!" << endl;
//      -0.717312
//  1.9083
//  4.43072
//  2.21258

  }


//    Double_t X0=1.44;
//    Double_t X1=4;
//    Double_t miu=3.77;
//    Double_t s=1.53;

  TCanvas* c1 = new TCanvas ("dEdx","dEdx",100,100,650,560);

  c1->Divide(1,1);
  c1->cd(1);

  TList* list_pi = hr_pi->GetListOfFunctions();


  TMultiGraph* mg1 = new TMultiGraph();
  
  mean_pi     = (TGraphErrors*)list_pi->FindObject(Form("Mean"));
//  TF1* Bethe_Bloch_p    = mean_p->GetFunction("Bethe_Bloch_p");
  TF1* Bethe_Bloch_pi    = (TF1*)list_pi->FindObject(Form("Bethe_Bloch_pi"));


    cout << "PAR1= " << Bethe_Bloch_pi->GetParameter(0) << endl;
    cout << "PAR2= " << Bethe_Bloch_pi->GetParameter(1) << endl;
    cout << "PAR3= " << Bethe_Bloch_pi->GetParameter(2) << endl;
    cout << "PAR4= " << Bethe_Bloch_pi->GetParameter(3) << endl;
    cout << "PAR5= " << Bethe_Bloch_pi->GetParameter(4) << endl;



  cout << "mean_pi " << mean_pi << endl;
  cout << "Bethe_Bloch_pi " << Bethe_Bloch_pi << endl;


  mean_p = 0;
  if (slic == 0 ) {
    TList* list_p = hr->GetListOfFunctions();
    mean_p = (TGraphErrors*)list_p->FindObject(Form("Mean"));
    hrplot->GetListOfFunctions()->Add(mean_p);
  }
  
  mg1->Add(mean_pi);
  if(mean_p)mg1->Add(mean_p);

  hrplot->GetXaxis()->SetTitle("log(p) [GeV/c]");
  hrplot->GetYaxis()->SetTitle("dEdx [arbitrary]");

  hrplot->GetListOfFunctions()->Add(mean_pi);
//  hrplot->GetListOfFunctions()->Add(mean_pi);

  hrplot->Draw();
  hr->Draw("same");

//  Bethe_Bloch_p->SetRange(0, 2.5);
//  Bethe_Bloch_p->SetLineColor(6);
//  Bethe_Bloch_p->Draw("same");



    // Bethe Bloch for pi with parameters from proton fit
    Double_t par_pi_from_pion[] = {par1, par2, par3, par4, mass_prot};

    cout << "!!!!!!!!!!!!!!" << endl;
    cout << par1 << endl;
    cout << par2 << endl;
    cout << par3 << endl;
    cout << par4 << endl;
    cout << "!!!!!!!!!!!!!!" << endl;

    gBethe_Bloch_p_from_pion = new TF1("Bethe Bloch pis", BetheBloch, -0.5, 3, 5);
    gBethe_Bloch_p_from_pion->SetParameters(par_pi_from_pion);
    gBethe_Bloch_p_from_pion->SetParNames("X0", "X1", "miu", "s", "mass");
//    gBethe_Bloch_pi_from_proton->FixParameter(4, mass_pi);

    gBethe_Bloch_p_from_pion->SetLineColor(1);
    gBethe_Bloch_p_from_pion->Draw("same");

    //  gBethe_Bloch_p->SetLineColor(2);
    //  gBethe_Bloch_p->SetLineStyle(2);
    //
    //  gBethe_Bloch_p->Draw("same");
    
    //
    //  TLatex text;
    //  text.SetTextFont(42);
    //  text.SetTextSize(0.5);
    //  text.SetTextColor(gBethe_Bloch_pi->GetLineColor());
    //  text.DrawLatex(0, 0.5, Form("pions from Gabor"));


    Bool_t MinuitFit = kTRUE;
    if(!MinuitFit)  return; 


//---------------------------------------------------------------------
//     Blocks for fitting procedure
//---------------------------------------------------------------------
        const Int_t npars = 4;
        TMinuit *gMinuit = new TMinuit(npars);  //initialize TMinuit with a maximum of 3 params
        gMinuit->SetFCN(fcn);

	cout<<"tu1"<<endl;

        Double_t arglist[10];
        Int_t ierflg = 0;

        arglist[0] = 1;
        gMinuit->mnexcm("SET ERR", arglist ,1,ierflg);

// Set starting values and step sizes for parameters
        //Double_t vstart[] = {1.047, 4.502 ,3.729, 1.607}; // from fit to pions
        Double_t vstart[] = {-4.641, 4.19 ,3.246, 2.145}; // from fit to protons
        Double_t step[]   = {0.1, 0.1 , 0.01, 0.01};
        gMinuit->mnparm(0, "X0",  vstart[0], step[0], 0,0,ierflg);
        gMinuit->mnparm(1, "X1",  vstart[1], step[1], 0,0,ierflg);
        gMinuit->mnparm(2, "miu", vstart[2], step[2], 0,0,ierflg);
        gMinuit->mnparm(3, "s",   vstart[3], step[3], 0,0,ierflg);


// Now ready for minimization step
        arglist[0] = 650;
        arglist[1] = 1.;
        gMinuit->mnexcm("MIGRAD", arglist ,2,ierflg);

// Print results
        Double_t amin,edm,errdef;
        Int_t nvpar,nparx,icstat;
        gMinuit->mnstat(amin,edm,errdef,nvpar,nparx,icstat);
        Int_t iuext;
        TString chnam;   // The name of the parameter
        Double_t val;    // The current (external) value of the parameter
        Double_t errl;   // The current estimate of the parameter uncertainty
        Double_t xlolim; // The lower bound (or zero if no limits)
        Double_t xuplim; // The upper bound (or zero if no limits)
        Int_t iuint;     // The internal parameter number
        Double_t currentPar[4];
        for (Int_t i=0; i< npars;i++) {
            gMinuit->mnpout(i, chnam, currentPar[i], errl, xlolim, xuplim, iuint);
	    cout<<"par="<<i<<" parval="<<currentPar[i]<<endl;
        }

	Double_t mpars_p[5] = {currentPar[0], currentPar[1], currentPar[2], currentPar[3], mass_prot};
	
	TF1* fBethe_Bloch_p = new TF1("fBethe Bloch p", BetheBloch, -0.5, 3, 5);
	fBethe_Bloch_p->SetParameters(mpars_p);
	fBethe_Bloch_p->SetParNames("X0", "X1", "miu", "s", "mass");
	
	fBethe_Bloch_p->SetLineColor(4);
	fBethe_Bloch_p->SetLineWidth(2);
	fBethe_Bloch_p->Draw("same");

	Double_t mpars_pi[5] = {currentPar[0], currentPar[1], currentPar[2], currentPar[3], mass_pi};
	
	TF1* fBethe_Bloch_pi = new TF1("fBethe Bloch pi", BetheBloch, -0.5, 3, 5);
	fBethe_Bloch_pi->SetParameters(mpars_pi);
	fBethe_Bloch_pi->SetParNames("X0", "X1", "miu", "s", "mass");
	
	fBethe_Bloch_pi->SetLineColor(6);
	fBethe_Bloch_pi->SetLineWidth(2);
	fBethe_Bloch_pi->Draw("same");

	Double_t mpars_k[5] = {currentPar[0], currentPar[1], currentPar[2], currentPar[3], mass_kaon};
	
	TF1* fBethe_Bloch_k = new TF1("fBethe Bloch k", BetheBloch, -0.5, 3, 5);
	fBethe_Bloch_k->SetParameters(mpars_k);
	fBethe_Bloch_k->SetParNames("X0", "X1", "miu", "s", "mass");
	
	fBethe_Bloch_k->SetLineColor(8);
	fBethe_Bloch_k->SetLineWidth(2);
	fBethe_Bloch_k->Draw("same");


// 1. Bethe Bloch 1: for pions
/*
        TF1 *fun_1=new TF1("Bethe Bloch 1",BetheBloch1, -0.5, 3, 4);

        fun_1->SetParameters(currentPar);
        fun_1->SetLineColor(kBlue);
        fun_1->SetLineStyle(1);
        fun_1->SetLineWidth(2);
//     mg1->Add(new TGraph(fun_1));
        fun_1->Draw("same");

        TF1 *fun_2=new TF1("Bethe Bloch 2",BetheBloch2, -0.5, 3, 4);
        fun_2->SetParameters(currentPar);
        fun_2->SetLineColor(kRed);
        fun_2->SetLineStyle(1);
        fun_2->SetLineWidth(4);
//     mg1->Add(new TGraph(fun_2));
        fun_2->Draw("same");
        TLegend *legend = new TLegend(0.52,0.95,0.95,0.85);
        legend->AddEntry(gr1,Form("%f+%fx",currentPar[0], currentPar[1]));
        legend->AddEntry(gr2,Form("%f+%fx + %fx^2",currentPar[0], currentPar[1], currentPar[2]));
        legend->Draw();
//-
        c1_1->Modified();
        c1->cd();
        c1->Modified();
//
        c1->Update();
*/

  
}

//______________________________________________________________________________
Double_t fcn1(Int_t &npar, Double_t *gin, Double_t *par, Int_t iflag)
{
    Int_t nbins = mean_pi->GetN();
    Int_t i;

//calculate chisquare
    Double_t chisq = 0;
    for (i=0;i<  nbins; i++) {
      //Double_t x = mean_pi->GetPointX(i);
      //Double_t y = mean_pi->GetPointY(i);
      Double_t x;
      Double_t y;
      mean_pi->GetPoint(i,x,y);

      cout<<"i="<<i<<" x="<<x<<" y="<<y<<"  EY="<<mean_pi->GetErrorY(i)<<"  func="<<BetheBloch1(x,par)<<endl;
      Double_t delta  = (y-BetheBloch1(x,par))/mean_pi->GetErrorY(i);
      chisq += delta*delta;
    }
    cout<<"fcn1="<<chisq<<endl;
    Int_t ii;
    cin>>ii;
    return chisq;
}

//______________________________________________________________________________
Double_t fcn2(Int_t &npar, Double_t *gin, Double_t *par, Int_t iflag)
{
    Int_t nbins = mean_p->GetN();
    Int_t i;

//calculate chisquare
    Double_t chisq = 0;
    for (i=0;i<  nbins; i++) {
      //Double_t x = mean_p->GetPointX(i);
      //Double_t y = mean_p->GetPointY(i);
	Double_t x;
	Double_t y;
	mean_p->GetPoint(i,x,y);
	
	cout<<"i="<<i<<" x="<<x<<" y="<<y<<"  EY="<<mean_p->GetErrorY(i)<<"  func="<<BetheBloch2(x,par)<<endl;
        Double_t delta  = (y-BetheBloch2(x,par))/mean_p->GetErrorY(i);
	chisq += delta*delta;
    }
    cout<<"fcn2="<<chisq<<endl;
    Int_t ii;
    cin>>ii;
    return chisq;
}
//______________________________________________________________________________
void fcn(Int_t &npar, Double_t *gin, Double_t &f, Double_t *par, Int_t iflag)
{
  Double_t f1 =  fcn1(npar,gin,par, iflag);
  Double_t f2 =  fcn2(npar,gin,par, iflag);
  f = f1 + f2;
  cout<<"fcn1="<<f1<<"  fcn2="<<f2<<"   fcn="<<f<<endl;
  return;
}
//-------------------------------------------------------------------


//_______________________________________________________________________
void DyVsY2(Int_t nSlicesInit, TH2F* h2, string particle, float &par1, float &par2, float &par3, float &par4)
{
  // this is general method to perform pid on m2 versus momentum
  // plot. The results are saved in file at path location.
  // path should point to the particular (Tof/Rich) selector
  // pid directory
      TH1D *slice[100];

      Float_t pmin = gPmin;
      Float_t pmax = gPmax;


      Float_t pstepInit = (pmax - pmin) / nSlicesInit;
      cout << "initial hbw/2 is " << pstepInit / 2. << endl;

      Float_t pstep = pstepInit;

      //Float_t pL = pmin;
      //Float_t pU;
      Int_t nSlices = 0;

//      Float_t pstep_array[] = {3.3, 3.3, 3.3, 3.3, 3.3, 3.3, 3.3, 3.3, 4, 5, 6};
//      nSlicesInit = sizeof(pstep_array) / 4;

      Float_t pL;
      Float_t pU;

      for (Int_t i = 0; i < nSlicesInit; i++) {
          //if(i>0)pL = pU;
          // if(pL>7)pstep = 2*pstepInit;
          //pU = pL + pstep_array[i];

          cout << "DyVsY2: number of slices is: " << nSlicesInit << endl;
          //pU = pmin + pstep * (i+1);
          pL = pmin + pstep * i;
          pU = pmin + pstep * (i + 1);

          Int_t binL = h2->GetXaxis()->FindBin(pL);
          Int_t binU = h2->GetXaxis()->FindBin(pU);
          Float_t p = 0.5 * (pL + pU);
          // if(p>pmax)break;
          nSlices++;
          slice[i] = (TH1D *) h2->ProjectionY(Form("slice%d", i), binL, binU);
          slice[i]->SetTitle(Form("p range: %f - %f GeV/c", pL, pU));
          slice[i]->Sumw2();
          //slice[i]->Rebin(binsToMerge);
          //slice[i]->Scale(1./(Float_t)binsToMerge);
          cout << "i=" << i << "  pL=" << pL << "   pU=" << pU << " max=" << slice[i]->GetMaximum() << endl;
          if (slice[i]->GetMaximum() < 3) {
              slice[i]->Rebin(2);
              //slice[i]->Scale(1./(Float_t)2);
          }
      }

      pmax = pU;

      cout << "Starting initial fit of slices" << endl;

      Float_t momPart[100];
      Float_t momPartErr[100];
      Float_t meanPart[100];
      Float_t meanPartErr[100];
      Float_t sigmaPart[100];
      Float_t sigmaPartProt[100];
      Float_t sigmaPartErr[100];
      Float_t betaPart[100];
      Float_t betaPartProt[100];
      Float_t betaPartErr[100];

      pL = pmin;
      pstep = pstepInit;

      for (Int_t i = 0; i < nSlicesInit; i++) {
          Float_t pL = pmin + pstep * i;
          Float_t pU = pmin + pstep * (i + 1);
          //if(i>0)pL = pU;
          // if(pL>7)pstep = 2*pstepInit;
          //pU = pL + pstep_array[i];
          Int_t binL = h2->GetXaxis()->FindBin(pL);
          Int_t binU = h2->GetXaxis()->FindBin(pU);
          Float_t p = 0.5 * (pL + pU);
          // if(p>pmax)break;

          TH1D *h = slice[i];

          // najwieksza lp zliczen:
          Int_t maxBin = h->GetMaximumBin();
          // max w binie
          Float_t cons = h->GetMaximum();
          // wartosc dla srodka binu
          Float_t mean = h->GetBinCenter(maxBin);
          Float_t sigm = 1.3;

          // ewentualne tło
          Float_t consB = h->GetBinContent(maxBin / 2);
          Float_t meanB = h->GetBinCenter(maxBin);
          Float_t sigmB = 100 * sigm;

          // do fitu parametry
          Double_t params[6];

          params[0] = cons;
          params[1] = mean;
          params[2] = sigm;
          params[3] = consB;
          params[4] = meanB;
          params[5] = sigmB;

          cout << "i=" << i << "  pL=" << pL << "   pU=" << pU << "  p=" << 0.5 * (pL + pU) << endl;
          TF1 *func;

          if (particle == "proton") {
              func = GetGaussFitPol(slice[i], 0.5 * (pL + pU), params);
          } else if (particle == "pion") {
              func = GetGaussFitPol_pi(slice[i], 0.5 * (pL + pU), params);
          }

          Bool_t drawSlices = kTRUE;
          if (drawSlices) {
              cnv[i] = new TCanvas(Form("slice%d", i), Form("slice%d", i), 5 * i, 5 * i, 400, 560);
              cnv[i]->SetTopMargin(0.02);
              cnv[i]->SetRightMargin(0.05);
              cnv[i]->SetLeftMargin(0.16);
              cnv[i]->SetBottomMargin(0.15);
              cnv[i]->Divide(1, 1);

              DrawSlice(slice[i], cnv[i]);
          }


          momPart[i] = p;
          momPartErr[i] = p * 0.005; // this should depend on field settings
          meanPart[i] = func->GetParameter(1); // to co powinien dawac Bethe Bloch
          //meanPart[i] = slice[i]->GetMean(); // to co powinien dawac Bethe Bloch
          meanPartErr[i] = TMath::Abs(func->GetParError(1));
          sigmaPart[i] = TMath::Abs(func->GetParameter(2));
          sigmaPartErr[i] = TMath::Abs(func->GetParError(2));
          betaPart[i] = p / TMath::Sqrt(p * p + 0.880354496);
          betaPartErr[i] = betaPart[i] * 0.000001;

      }


      h2->GetListOfFunctions()->Clear();

      TGraphErrors *Mean = new TGraphErrors(nSlices, momPart, meanPart, momPartErr, meanPartErr);
      TGraphErrors *Sigma = new TGraphErrors(nSlices, momPart, sigmaPart, momPartErr, sigmaPartErr);
      TGraphErrors *SigmaVsMean = new TGraphErrors(nSlices, meanPart, sigmaPart, meanPartErr, sigmaPartErr);
      TGraphErrors *SigmaVsBeta = new TGraphErrors(nSlices, betaPart, sigmaPart, betaPartErr, sigmaPartErr);
      TGraphErrors *SigmaVsBetaProt = new TGraphErrors(nSlices, betaPartProt, sigmaPartProt, betaPartErr, sigmaPartErr);

      Mean->SetName(Form("Mean"));
      Sigma->SetName(Form("Sigma"));
      SigmaVsMean->SetName("SigmaVsMean");
      SigmaVsBeta->SetName("SigmaVsBeta");
      SigmaVsBetaProt->SetName("SigmaVsBetaProt");

      Mean->SetMarkerStyle(20);
      Mean->SetMarkerSize(0.05);
      Mean->SetMarkerColor(2);
      Mean->SetLineColor(2);
      Mean->SetLineWidth(3);

      Sigma->SetMarkerStyle(20);
      Sigma->SetMarkerSize(0.05);
      Sigma->SetMarkerColor(4);
      Sigma->SetLineColor(4);
      Sigma->SetLineWidth(3);

      SigmaVsMean->SetMarkerStyle(20);
      SigmaVsMean->SetMarkerSize(0.05);
      SigmaVsMean->SetMarkerColor(4);
      SigmaVsMean->SetLineColor(4);
      SigmaVsMean->SetLineWidth(3);

      SigmaVsBeta->SetMarkerStyle(20);
      SigmaVsBeta->SetMarkerSize(0.05);
      SigmaVsBeta->SetMarkerColor(6);
      SigmaVsBeta->SetLineColor(4);
      SigmaVsBeta->SetLineWidth(3);

      SigmaVsBetaProt->SetMarkerStyle(20);
      SigmaVsBetaProt->SetMarkerSize(0.05);
      SigmaVsBetaProt->SetMarkerColor(2);
      SigmaVsBetaProt->SetLineColor(2);
      SigmaVsBetaProt->SetLineWidth(3);

      // PROTON
    TF1 *meanFit_pion = 0;

    if (particle == "pion") {
          meanFit_pion = new TF1("Bethe_Bloch_pi", BetheBloch, -0.5, 3, 5);
          meanFit_pion->SetParNames("X0", "X1", "miu", "s");
          meanFit_pion->SetLineColor(3);
          meanFit_pion->SetLineWidth(4);
          meanFit_pion->SetParameters(gBethe_Bloch_pi->GetParameter(0), gBethe_Bloch_pi->GetParameter(1),
                                      gBethe_Bloch_pi->GetParameter(2), gBethe_Bloch_pi->GetParameter(3),
                                      gBethe_Bloch_pi->GetParameter(4));
          meanFit_pion->FixParameter(4,  gBethe_Bloch_pi->GetParameter(4));
          Mean->Fit(meanFit_pion, "w", "", pmin - 0.1, pmax + 0.1);

          cout << " gBethe_Bloch_pi->GetParameter(4) = " << gBethe_Bloch_pi->GetParameter(4) << endl;
          cout << " meanFit_pion = " <<meanFit_pion->GetParameter(4) << endl;

          Mean->GetListOfFunctions()->Add(meanFit_pion);

          par1 = meanFit_pion->GetParameter(0);
          par2 = meanFit_pion->GetParameter(1);
          par3 = meanFit_pion->GetParameter(2);
          par4 = meanFit_pion->GetParameter(3);

      }

    h2->GetListOfFunctions()->Add(Mean);
    h2->GetListOfFunctions()->Add(meanFit_pion);
    

      cout << "functions added" << endl;
}


//___________________________________________________________________________
void DrawSlice(TH1D* h,TCanvas* canv)
{
  
  TString name = h->GetName();
  name.Append("_Fit");
  TF1* func = (TF1*)h->FindObject(name.Data());
  func->SetLineWidth(2);
  func->SetLineColor(2);
  cout<<"drawing function: "<<func->GetName()<<" canvas: "<<canv->GetName()<<" h: "<<h->GetName()<<endl;
  canv->cd(1);
  //h->SetAxisRange(-10,10);
  h->Draw();
  func->Draw("same");
  
}  

//______________________________________________________________
TF1* GetGaussFitPol(TH1D* h, Float_t p, const Double_t* params)
{
  // h - histograms to fit
  // params - pointer to list of parameters 
  // n - number of Gaussians to consider
  // The method create function fit it to histogram and return.

  // Rozkład Gaussa
  TString g1 = "[0]*exp(-0.5*((x-[1])/[2])^2)";

  // Polynomial pierwszy wyraz liniowy, drugi kwadratowy ale względem 1.13
  // parabola, ktrora ma wierzchołek w 1.13
  // cala funkcja to jest g1 + g2
  TString g2 = "[3] + [4]*(x-1.13) + [5]*(x-1.13)*(x-1.13)";
  TString name = h->GetName();
  name.Append("_Fit");
  cout<<"name: "<<name.Data()<<endl;
  TF1* func = new TF1(name.Data(),Form("%s + %s",g1.Data(),g2.Data()),0.7,1.8);  

  Float_t sig = 0.054;
  Float_t dd = 0.01;

  func->SetParameters(params);

  func->SetParLimits(0,1,500);

    func->SetParLimits(1,gBethe_Bloch_p->Eval(p)-3*dd, gBethe_Bloch_p->Eval(p)+3*dd);
    func->SetParLimits(2,sig-0.015,sig+0.015);
    cout<<"gBethe_Bloch_p->Eval(p) : "<<gBethe_Bloch_p->Eval(p)<<" p="<<p<<endl;


  func->FixParameter(3,0);
  func->FixParameter(4,0);
  func->FixParameter(5,0);

  //func->SetParNames("const","mean","sigma","constB","meanB","sigmaB");
  func->SetParNames("const","mean","sigma","p0","p1","p2");
  

  //limits on mean, sigma and constant  

  Float_t rmin = 0.7;
  Float_t rmax =  1.8;
 
  h->Rebin(4);
  //h->Scale(0.25);    
  
  h->Fit(func,"w0","",rmin,rmax);

  //func->SetParLimits(1,func->GetParameter(1)-1*dd, func->GetParameter(1)+1*dd);
  //func->SetParLimits(2,func->GetParameter(2)-0.005,func->GetParameter(2)+0.035);

  h->Fit(func,"w0","",rmin,rmax);
  
  h->GetListOfFunctions()->Add(func);
  
  return func;
}

//______________________________________________________________
TF1* GetGaussFitPol_pi(TH1D* h, Float_t p, const Double_t* params)
{
// h - histograms to fit
// params - pointer to list of parameters
// n - number of Gaussians to consider
// The method create function fit it to histogram and return.

// Rozkład Gaussa
TString g1 = "[0]*exp(-0.5*((x-[1])/[2])^2)";

// Polynomial pierwszy wyraz liniowy, drugi kwadratowy ale względem 1.13
// parabola, ktrora ma wierzchołek w 1.13
// cala funkcja to jest g1 + g2
TString g2 = "[3] + [4]*(x-1.13) + [5]*(x-1.13)*(x-1.13)";
TString name = h->GetName();
name.Append("_Fit");
cout<<"name: "<<name.Data()<<endl;
TF1* func = new TF1(name.Data(),Form("%s + %s",g1.Data(),g2.Data()),0.7,1.8);

Float_t sig = 0.054;
Float_t dd = 0.01;

func->SetParameters(params);

func->SetParLimits(0,1,500);

func->SetParLimits(1,gBethe_Bloch_pi->Eval(p)-3*dd, gBethe_Bloch_pi->Eval(p)+3*dd);
func->SetParLimits(2,sig-0.015,sig+0.015);
cout<<"gBethe_Bloch_p->Eval(p) : "<<gBethe_Bloch_pi->Eval(p)<<" p="<<p<<endl;


func->FixParameter(3,0);
func->FixParameter(4,0);
func->FixParameter(5,0);

//func->SetParNames("const","mean","sigma","constB","meanB","sigmaB");
func->SetParNames("const","mean","sigma","p0","p1","p2");


//limits on mean, sigma and constant

Float_t rmin = 0.7;
Float_t rmax =  1.8;

h->Rebin(4);
//h->Scale(0.25);

h->Fit(func,"w0","",rmin,rmax);

//func->SetParLimits(1,func->GetParameter(1)-1*dd, func->GetParameter(1)+1*dd);
//func->SetParLimits(2,func->GetParameter(2)-0.005,func->GetParameter(2)+0.035);

h->Fit(func,"w0","",rmin,rmax);

h->GetListOfFunctions()->Add(func);

return func;
}


//__________________________________________________________________________
Double_t BetheBloch(Double_t *x, Double_t *par) 
{

  // x + 4 paramtery Bethego-Blocha:
  Double_t log10p = x[0];
  Double_t p = TMath::Power(10,log10p);

  Double_t X0 = par[0]; // interval limits
  Double_t X1 = par[1];
  Double_t miu = par[2]; // betha*gamma value where Bethe-Bloch function gets minimum
  Double_t s = par[3]; // tam, gdzie wysyca sie Bethe Bloch dla duzych pedow

  
  // do funkcji wartości:
  Double_t mass = par[4];
  Double_t betha = p/TMath::Sqrt(p*p + mass*mass); // betha = p/sqrt(p^2 + m^2)
  Double_t gamma = 1/TMath::Sqrt(1 - betha*betha); // gamma = 1/sqrt(1 - beta^2)
  Double_t alpha = -2;
  Double_t X = TMath::Log10(betha * gamma);   // betha*gamma=p/mc
  Double_t B = miu/TMath::Sqrt(1 + miu*miu);                                  
  Double_t E0 = 1./(miu*miu);
  Double_t b = miu*miu + TMath::Log(TMath::Abs(1 - B*B));                    
  Double_t M = (X1 - X0) / ((s/E0 - b + 1)/(2 * TMath::Log(10)) - X0);          
  Double_t XA = X0 - (X0 - X1)/M;
  Double_t a = 2. * TMath::Log(10.) * TMath::Power(X1-X0,1.-M) / M;     


  Double_t delta = 0;
  if (X < X0) {
    delta = 0;
  }
  else if (X >= X0 || X < X1){
    delta = 2 * TMath::Log(10.) * (X - XA) + a * TMath::Power(X1 - X, M);
  }
  else if ( X1 <= X ){
    delta = 2 * TMath::Log(10) * (X - XA);
  }

  Double_t dEdx = E0 * TMath::Power(betha, alpha) * (b + 2*TMath::Log(gamma) - betha*betha - delta);

  return dEdx;
}

Double_t BetheBloch1(Double_t x, Double_t *par)
{
    // Bethe Bloch assuming pion mass

    // x + 4 paramtery Bethego-Blocha:
    Double_t log10p = x;
    Double_t p = TMath::Power(10,log10p);

    Double_t X0 = par[0]; // interval limits
    Double_t X1 = par[1];
    Double_t miu = par[2]; // betha*gamma value where Bethe-Bloch function gets minimum
    Double_t s = par[3]; // tam, gdzie wysyca sie Bethe Bloch dla duzych pedow


    // do funkcji wartości:
    Double_t mass = 0.13957;
    Double_t betha = p/TMath::Sqrt(p*p + mass*mass); // betha = p/sqrt(p^2 + m^2)
    Double_t gamma = 1/TMath::Sqrt(1 - betha*betha); // gamma = 1/sqrt(1 - beta^2)
    Double_t alpha = -2;
    Double_t X = TMath::Log10(betha * gamma);   // betha*gamma=p/mc
    Double_t B = miu/TMath::Sqrt(1 + miu*miu);
    Double_t E0 = 1./(miu*miu);
    Double_t b = miu*miu + TMath::Log(TMath::Abs(1 - B*B));
    Double_t M = (X1 - X0) / ((s/E0 - b + 1)/(2 * TMath::Log(10)) - X0);
    Double_t XA = X0 - (X0 - X1)/M;
    Double_t a = 2. * TMath::Log(10.) * TMath::Power(X1-X0,1.-M) / M;


    Double_t delta = 0;
    if (X < X0) {
        delta = 0;
    }
    else if (X >= X0 || X < X1){
        delta = 2 * TMath::Log(10.) * (X - XA) + a * TMath::Power(X1 - X, M);
    }
    else if ( X1 <= X ){
        delta = 2 * TMath::Log(10) * (X - XA);
    }

    Double_t dEdx = E0 * TMath::Power(betha, alpha) * (b + 2*TMath::Log(gamma) - betha*betha - delta);

    return dEdx;
}

Double_t BetheBloch2(Double_t x, Double_t *par)
{
    // Bethe Bloch assuming proton mass

    // x + 4 paramtery Bethego-Blocha:
    Double_t log10p = x;
    Double_t p = TMath::Power(10,log10p);

    Double_t X0 = par[0]; // interval limits
    Double_t X1 = par[1];
    Double_t miu = par[2]; // betha*gamma value where Bethe-Bloch function gets minimum
    Double_t s = par[3]; // tam, gdzie wysyca sie Bethe Bloch dla duzych pedow


    // do funkcji wartości:
    Double_t mass = 0.9382720799426998;
    Double_t betha = p/TMath::Sqrt(p*p + mass*mass); // betha = p/sqrt(p^2 + m^2)
    Double_t gamma = 1/TMath::Sqrt(1 - betha*betha); // gamma = 1/sqrt(1 - beta^2)
    Double_t alpha = -2;
    Double_t X = TMath::Log10(betha * gamma);   // betha*gamma=p/mc
    Double_t B = miu/TMath::Sqrt(1 + miu*miu);
    Double_t E0 = 1./(miu*miu);
    Double_t b = miu*miu + TMath::Log(TMath::Abs(1 - B*B));
    Double_t M = (X1 - X0) / ((s/E0 - b + 1)/(2 * TMath::Log(10)) - X0);
    Double_t XA = X0 - (X0 - X1)/M;
    Double_t a = 2. * TMath::Log(10.) * TMath::Power(X1-X0,1.-M) / M;


    Double_t delta = 0;
    if (X < X0) {
        delta = 0;
    }
    else if (X >= X0 || X < X1){
        delta = 2 * TMath::Log(10.) * (X - XA) + a * TMath::Power(X1 - X, M);
    }
    else if ( X1 <= X ){
        delta = 2 * TMath::Log(10) * (X - XA);
    }

    Double_t dEdx = E0 * TMath::Power(betha, alpha) * (b + 2*TMath::Log(gamma) - betha*betha - delta);

    return dEdx;
}

//___________________________________________________________________________
void TPstyle()
{
  gStyle->SetCanvasBorderMode(0);
  gStyle->SetCanvasBorderSize(2);
  gStyle->SetCanvasColor(10);
  
  gStyle->SetFrameFillColor(10);
  gStyle->SetFrameBorderSize(1);
  gStyle->SetFrameBorderMode(-1);
  gStyle->SetFrameLineWidth(2);
  gStyle->SetFrameLineColor(1);

  gStyle->SetHistFillColor(0);
  gStyle->SetHistLineWidth(2);
  gStyle->SetHistLineColor(1);

  gStyle->SetPadColor(0);
  gStyle->SetPadBorderSize(0);
  gStyle->SetPadBorderMode(-1);

  gStyle->SetTextFont(133);
  gStyle->SetTitleFont(133, "xyz");
  gStyle->SetLabelFont(133, "xyz");
  gStyle->SetLabelSize(20, "xyz");
  gStyle->SetTitleSize(20, "xyz");
  gStyle->SetTitleOffset(1, "y");  
  gStyle->SetTitleOffset(1, "x");
  gStyle->SetEndErrorSize(2);

  gStyle->SetTitleColor(1,"X");
  gStyle->SetTitleColor(1,"Y");

  gStyle->SetOptDate(0);
  gStyle->SetOptTitle(0);
  gStyle->SetOptStat(0);
  gStyle->SetOptFit(0);
}


