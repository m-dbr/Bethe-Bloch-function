#include "TMinuit.h"
#include "TROOT.h"
#include "TCanvas.h"
#include "TStyle.h"
#include "TGraphBentErrors.h"
#include "TAxis.h"
#include "TMultiGraph.h"
#include "TLegend.h"
// Valeri Fine 07/08/2009 (for Vitali Okorokov)
// 1. Points for linear function 
   const int nPoints1 = 3;
   Double_t ux1[]={1.0,2.0,3.0}; 
   Double_t uy1[]={1.15,1.84,3.02};
   Double_t exl1[] = {0.0,0.0,0.0};
   Double_t eyl1[] = {0.10,0.10,0.10};
   Double_t exh1[] = {0.0,0.0,0.0};
   Double_t eyh1[] = {0.10,0.10,0.10};
   Double_t exld1[] = {0.0,0.0,0.0};
   Double_t eyld1[] = {0.0,0.0,0.0};
   Double_t exhd1[] = {0.0,0.0,0.0};
   Double_t eyhd1[] = {0.0,0.0,0.0};
   
// 2. Points for quadratic function 
   const int nPoints2 = 5;
   Double_t x2[]={1.0,1.25,1.50,1.75,2.0}; 
   Double_t y2[]={0.95,1.56,2.22,3.07,4.02};
   Double_t exl2[] = {0.0,0.0,0.0,0.0,0.0};
   Double_t eyl2[] = {0.10,0.10,0.10,0.10,0.10};
   Double_t exh2[] = {0.0,0.0,0.0,0.0,0.0};
   Double_t eyh2[] = {0.10,0.10,0.10,0.10,0.10};
   Double_t exld2[] = {0.0,0.0,0.0,0.0,0.0};
   Double_t eyld2[] = {0.0,0.0,0.0,0.0,0.0};
   Double_t exhd2[] = {0.0,0.0,0.0,0.0,0.0};
   Double_t eyhd2[] = {0.0,0.0,0.0,0.0,0.0};
//------------------------------------------------------------------
//		Approximation functions
// 1. Linear function 
//------------------------------------------------------------------- 

Double_t myFunc_1(Double_t *xp, Double_t *par) 
{
// ------- approximation function -------
   Double_t x=*xp;
   return par[0]+par[1]*x;
}

//-------------------------------------------------------------------
// 2. Quadratic function 
//------------------------------------------------------------------- 
Double_t myFunc_2(Double_t *xp, Double_t *par) 
{
   // ------- approximation function -------
   Double_t x=*xp;
   return par[0]+par[1]*x+par[2]*x*x;
}

//______________________________________________________________________________
Double_t fcn1(Int_t &npar, Double_t *gin, Double_t *par, Int_t iflag)
{
   const Int_t nbins = nPoints1;
   Int_t i;

//calculate chisquare
   Double_t chisq = 0;
   for (i=0;i<  nbins; i++) {
     Double_t delta  = (uy1[i]-myFunc_1(&ux1[i],par))/eyl1[i];
     chisq += delta*delta;
   }
   return chisq;
}

//______________________________________________________________________________
Double_t fcn2(Int_t &npar, Double_t *gin, Double_t *par, Int_t iflag)
{ 
   const Int_t nbins = nPoints2;
   Int_t i;

//calculate chisquare
   Double_t chisq = 0;
   for (i=0;i<  nbins; i++) {
     Double_t delta  = (y2[i]-myFunc_2(&x2[i],par))/eyl2[i];
     chisq += delta*delta;
   }
   return chisq;
}
//______________________________________________________________________________
void fcn(Int_t &npar, Double_t *gin, Double_t &f, Double_t *par, Int_t iflag)
{
   f = fcn1(npar,gin,par, iflag) + fcn2(npar,gin, par, iflag);
}
//-------------------------------------------------------------------
//
void fit()
//void SimultFit_MyTest()
{
// Vitaly A. Okorokov, 07/08/2009
// Simultaneus fit of linear and quadratic function.
// My test code. 
//
    gROOT->Reset();
    gROOT->SetStyle("Plain");
// 0 -- no Grid for corresponding axis
    gStyle->SetPadGridX(0);
    gStyle->SetPadGridY(0);
// 0 -- no Date for canvas (1-default-date is printed in the bottom/left corner)
    gStyle->SetOptDate(0);
    gStyle->SetOptFit(1111);
//
  TCanvas *c1 = new TCanvas("c1","Simultaneous fit: test",800,800);
  c1->SetBorderSize(0);
  c1->SetFillStyle(4000);
//
// ------------>Primitives in pad: c1_1
    TPad *c1_1 = new TPad("c1_1", "c1_1",0.01,0.01,0.99,0.99);
    c1_1->Draw();
    c1_1->cd();
    c1_1->SetTopMargin(0.03);
    c1_1->SetBottomMargin(0.14);
    c1_1->SetRightMargin(0.06);
    c1_1->SetLeftMargin(0.11);
//
  TGraphBentErrors *gr1 = new TGraphBentErrors(3,ux1,uy1,exl1,exh1,eyl1,eyh1,exld1,exhd1,eyld1,eyhd1);
  gr1->SetMarkerStyle(20);
  gr1->SetMarkerSize(1.4);
  gr1->SetMarkerColor(kBlue);
  gr1->SetLineColor(kBlue);
  gr1->SetFillColor(0);
  gr1->SetLineWidth(2);
//
  TGraphBentErrors *gr2 = new TGraphBentErrors(5,x2,y2,exl2,exh2,eyl2,eyh2,exld2,exhd2,eyld2,eyhd2);
  gr2->SetMarkerStyle(21);
  gr2->SetMarkerSize(1.2);
  gr2->SetMarkerColor(kRed);
  gr2->SetLineColor(kRed);
  gr2->SetFillColor(0);
  gr2->SetLineWidth(2);
//
// create a MultiGraph
  TMultiGraph *mg1 = new TMultiGraph();
// absent of the MultiGraph title
   mg1->SetTitle("");   
   mg1->Add(gr1);
   mg1->Add(gr2);
// draw this MultiGraph
   mg1->Draw("AP");
// X-axis
   TAxis *axisx1=mg1->GetXaxis();
   axisx1->SetTickLength(0.025);
   axisx1->SetLabelSize(0.045);
   axisx1->SetTitle("x");
   axisx1->SetTitleSize(0.05);
   axisx1->SetTitleOffset(0.9);
   axisx1->Draw();
   gPad->SetTickx();
   axisx1->Draw();
   gPad->SetTickx();
// Y-axis
   TAxis *axisy1=mg1->GetYaxis();
   axisy1->SetRangeUser(0.,5.);
   axisy1->SetLabelSize(0.055);
   axisy1->SetTickLength(0.017);
// axisy1->SetNdivisions(804);
   axisy1->SetTitle("y");
   axisy1->SetTitleSize(0.065);
   axisy1->SetTitleOffset(0.65);
   axisy1->Draw();
   gPad->SetTicky();
//---------------------------------------------------------------------
//     Blocks for fitting procedure
//---------------------------------------------------------------------
   const Int_t npars = 3;
   TMinuit *gMinuit = new TMinuit(npars);  //initialize TMinuit with a maximum of 3 params
   gMinuit->SetFCN(fcn);

   Double_t arglist[10];
   Int_t ierflg = 0;

   arglist[0] = 1;
   gMinuit->mnexcm("SET ERR", arglist ,1,ierflg);

// Set starting values and step sizes for parameters
   static Double_t vstart[] = {0.7, 1.5 ,  0 };
   static Double_t step[]   = {0.1, 0.1 , 0.1};
   gMinuit->mnparm(0, "a0", vstart[0], step[0], 0,0,ierflg);
   gMinuit->mnparm(1, "a1", vstart[1], step[1], 0,0,ierflg);
   gMinuit->mnparm(2, "a2", vstart[2], step[2], 0,0,ierflg);
//   gMinuit->mnparm(3, "a4", vstart[3], step[3], 0,0,ierflg);

// Now ready for minimization step
   arglist[0] = 500;
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
   Int_t i = 0;
   Double_t currentPar[npars] = {0};
   for (i=0; i< npars;i++) {
      gMinuit->mnpout(i, chnam, currentPar[i], errl, xlolim, xuplim, iuint);
   }
// 1. Linear function
     TF1 *fun_1 = new TF1("fun_1",myFunc_1,1.0,3.0,2); 	

     fun_1->SetParameters(currentPar);
     fun_1->SetLineColor(kBlue);
     fun_1->SetLineStyle(1);
     fun_1->SetLineWidth(2);
//     mg1->Add(new TGraph(fun_1));
     fun_1->Draw("same");
     
     TF1 *fun_2=new TF1("fun_2",myFunc_2,1.0,2.0,3);
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
  }
















