TCanvas *cnv[100];
Float_t gPmin;
Float_t gPmax;
TF1* gBethe_Bloch_p;
TF1* gBethe_Bloch_pi;
TF1* gBethe_Bloch_pi_from_proton;
TF1* gBethe_Bloch_k;
TF1* GetGaussFitPol(TH1D* h, Float_t p, const Double_t* params);
TF1* GetGaussFitPol_pi(TH1D* h, Float_t p, const Double_t* params);

void TPstyle();
void DrawSlice(TH1D* h,TCanvas* canv);
void DyVsY2(Int_t nSlicesInit, TH2F* h2, string particle, float &par1, float &par2, float &par3, float &par4);

Double_t BetheBloch(Double_t *x, Double_t *par);

//___________MAIN FUNCTION_______________________________________________________________________//
void FitSlices(Int_t runid=1234, Int_t apver=23, Int_t nSlices=15, Int_t saveCalibParams=0)
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
//  Double_t mass_kaon = 0.493677;
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

  TH2F* hrplot = (TH2F*)hr->Clone("hrplot");

  float par1 = 0;
  float par2 = 0;
  float par3 = 0;
  float par4 = 0;

// Slicing function for pions/protons
  string particle = "proton";
  DyVsY2(nSlices,hr, particle, par1, par2, par3, par4);

  string particle2 = "pion";
  DyVsY2(nSlices,hr_pi, particle2, par1, par2, par3, par4);

//    Double_t X0=1.44;
//    Double_t X1=4;
//    Double_t miu=3.77;
//    Double_t s=1.53;

par1 = TMath::Abs(par1); // DLACZEGO UJEMNA WARTOSC DLA X0 Z PROTONOW WYCHODZI???
cout << "!!!!!!!!!!!!!!!!!!" << endl;
cout << "X0 = " << par1 << endl;
cout << "X1 = " << par2 << endl;
cout << "miu = " << par3 << endl;
cout << "s = " << s << endl;
cout << "!!!!!!!!!!!!!!!!!!" << endl;

// TODO: dlaczego X1 wychodzi na minusie??
// TODO: piony z protonow nie robia sie
// TODO: protony z pionow cos brzydko wyszly
// TODO: poprawic slicowanie pionow bo teraz sa za wysoko
// TODO: dodac podpisy i kolory ladne



  TCanvas* c1 = new TCanvas ("dEdx","dEdx",100,100,650,560);

  c1->Divide(1,1);
  c1->cd(1);

  TList* list_p = hr->GetListOfFunctions();
  TList* list_pi = hr_pi->GetListOfFunctions();

  TMultiGraph* mg1 = new TMultiGraph();
  
  TGraphErrors* mean_p     = (TGraphErrors*)list_p->FindObject(Form("Mean"));
  TGraphErrors* mean_pi     = (TGraphErrors*)list_pi->FindObject(Form("Mean"));

  mg1->Add(mean_p);
  mg1->Add(mean_pi);


  hrplot->GetXaxis()->SetTitle("log(p) [GeV/c]");
  hrplot->GetYaxis()->SetTitle("dEdx [arbitrary]");

  hrplot->GetListOfFunctions()->Add(mean_p);
  hrplot->GetListOfFunctions()->Add(mean_pi);

  hrplot->Draw();
  hr_pi->Draw("same");

  // Bethe Bloch for pi with parameters from proton fit
    Double_t par_pi_from_proton[] = {par1, par2, par3, par4, mass_pi};

    gBethe_Bloch_pi_from_proton = new TF1("Bethe Bloch pis", BetheBloch, -0.5, 3, 5);
    gBethe_Bloch_pi_from_proton->SetParameters(par_pi_from_proton);
    gBethe_Bloch_pi_from_proton->SetParNames("X0", "X1", "miu", "s", "mass");
    gBethe_Bloch_pi_from_proton->FixParameter(4, mass_pi);

    gBethe_Bloch_pi_from_proton->SetLineColor(1);
    gBethe_Bloch_pi_from_proton->Draw("same");


  //  gBethe_Bloch_p->SetLineColor(2);
//  gBethe_Bloch_p->SetLineStyle(2);
//
//  gBethe_Bloch_p->Draw("same");
//


//
//  TLatex text;
//  text.SetTextFont(42);
//  text.SetTextSize(0.5);
//  text.SetTextColor(gBethe_Bloch_pi->GetLineColor());
//  text.DrawLatex(0, 0.5, Form("pions from Gabor"));
  
  return;  

  TCanvas* c2 = new TCanvas ("sigmaVsmean","",200,200,650,560);
  
  c2->Divide(1,1);
  c2->cd(1); 
  mg1->Draw("a");
  
}

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

      Float_t pstep_array[] = {3.3, 3.3, 3.3, 3.3, 3.3, 3.3, 3.3, 3.3, 4, 5, 6};
      nSlicesInit = sizeof(pstep_array) / 4;

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
      if (particle == "proton") {
          TF1 *meanFit_prot = new TF1("Bethe Bloch p", BetheBloch, -0.5, 3, 5);
          meanFit_prot->SetParNames("X0", "X1", "miu", "s");
          meanFit_prot->SetLineColor(3);
          meanFit_prot->SetLineWidth(4);
          meanFit_prot->SetParameters(gBethe_Bloch_p->GetParameter(0), gBethe_Bloch_p->GetParameter(1),
                                      gBethe_Bloch_p->GetParameter(2), gBethe_Bloch_p->GetParameter(3),
                                      gBethe_Bloch_p->GetParameter(4));
          meanFit_prot->FixParameter(4, gBethe_Bloch_p->GetParameter(4));
          Mean->Fit(meanFit_prot, "w", "", pmin - 0.1, pmax + 0.1);

          Mean->GetListOfFunctions()->Add(meanFit_prot);

          par1 = meanFit_prot->GetParameter(0);
          par2 = meanFit_prot->GetParameter(1);
          par3 = meanFit_prot->GetParameter(2);
          par4 = meanFit_prot->GetParameter(3);

      }

      h2->GetListOfFunctions()->Add(Mean);



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

  //  gStyle->SetStatColor(10);
  //  gStyle->SetTitleColor(10);

  gStyle->SetTextFont(133);
  gStyle->SetTitleFont(133, "xyz");
  gStyle->SetLabelFont(133, "xyz");
  gStyle->SetLabelSize(20, "xyz");
  gStyle->SetTitleSize(20, "xyz");
  gStyle->SetTitleOffset(1, "y");  
  gStyle->SetTitleOffset(1, "x");
  gStyle->SetEndErrorSize(2);

  //  gStyle->SetTitleFont(42);
  //  gStyle->SetLabelFont(42,"X");
  //  gStyle->SetLabelFont(42,"Y");
  //  gStyle->SetStatFont(42);

  //  gStyle->SetTitleOffset(1.1,"X");
  //  gStyle->SetTitleOffset(0.9,"Y");
  gStyle->SetTitleColor(1,"X");
  gStyle->SetTitleColor(1,"Y");

  gStyle->SetOptDate(0);
  gStyle->SetOptTitle(0);
  gStyle->SetOptStat(0);
  gStyle->SetOptFit(0);
}


// dostalismy pion z Gabora pion z fitow do protonow chce dostac,
// robie funkcje dla pionow wykorzystujac parametry nie z Gabora tylk oz fitow i dodatkowa linia (TLatex)
// to samo dla pionow tylko nie z parametrami z gabora tylko
// fitowac piony i rysowac kaony

// slicowanie dla pionow i rysowanie pionow