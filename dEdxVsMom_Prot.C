Int_t ican;
TCanvas *cnv[100];
Int_t ii;
Float_t gPmin;
Float_t gPmax;
TF1* gdEdx_pion;
TF1* gdEdx_kaon;
TF1* gdEdx_prot;
TF1* gdEdx_muon;
TF1* gdEdx_elec;
TF1* gdEdx_deut;

void TPstyle();
Double_t FindMaximum(Double_t* array,Int_t N);
Double_t Univ(Double_t *x, Double_t *par);
void FindFitRange(TH1F* H,Float_t& xmin,Float_t& xmax,Float_t fraction);
TF1* GetGaussFitPol(TH1D* h, Float_t p, const Double_t* params);
TF1* GetQuintupleGaussFit(TH1D* h,Float_t p,const Double_t* params);
TF1* GetTripleGaussFit(TH1D* h,Float_t p,const Double_t* params);
TF1* GetDoubleGaussFit(TH1D* h,const Double_t* params);
void DrawSlice(TH1D* h,TCanvas* canv);
void DyVsY2(Int_t nSlicesInit, TH2F* h2);
// void DyVsY(Int_t nSlices, TH2F* h2);

//__________________________________________________________________________________
void FitSlices(Int_t runid=1234, Int_t nSlices=15, Int_t saveCalibParams=0)
{
  // VTPC1: range == 1   
  // VTPC2: range == 2   
  //Int_t run=168;
  //Int_t run=168;

  gStyle->SetOptStat(0);
  gStyle->SetOptFit(1);

  gPmin=5.5;
  gPmax=30.0;
  //gPmax=7.0;

  if(!gdEdx_pion){
    ///// pion dEdx
    Float_t pmin = 0.1;
    Float_t pmax = 50.0;
    // from InvMass selection 1350-1450:
    //Double_t A = 0.0776135;
    //Double_t B = 7.11011e-05;
    // from InvMass selection 652-812:
    Double_t A = 0.0778881;
    Double_t B = 6.71445e-05;

    ///// Proton dEdx
    TString beta2 = "(x*x/(x*x+0.880354496))";
    gdEdx_prot = new TF1("gdEdx_prot",Form("[0]/%s * (TMath::Log(%s/([1]*(1-%s))) - 2*%s)",beta2.Data(),
						beta2.Data(),beta2.Data(),beta2.Data()),pmin,pmax);
    gdEdx_prot->SetParNames("A","B");
    gdEdx_prot->SetParameters(A,B); 

    beta2 = "(x*x/(x*x+0.01947983))";
    gdEdx_pion = new TF1("gdEdx_pion",Form("[0]/%s * (TMath::Log(%s/([1]*(1-%s))) - 2*%s)",beta2.Data(),
					   beta2.Data(),beta2.Data(),beta2.Data()),pmin,pmax);
    gdEdx_pion->SetParNames("A","B"); 
    gdEdx_pion->SetParameters(gdEdx_prot->GetParameter(0),gdEdx_prot->GetParameter(1)); 

    ///// Kaon dEdx
    beta2 = "(x*x/(x*x+0.24371698))";
    gdEdx_kaon = new TF1("gdEdx_kaon",Form("[0]/%s * (TMath::Log(%s/([1]*(1-%s))) - 2*%s)",beta2.Data(),
						beta2.Data(),beta2.Data(),beta2.Data()),pmin,pmax);
    gdEdx_kaon->SetParNames("A","B");
    gdEdx_kaon->SetParameters(gdEdx_pion->GetParameter(0),gdEdx_pion->GetParameter(1)); 
    

    ///// muon dEdx, mass = 0.1056583745 (mass^2 = 0.011163692102)
    beta2 = "(x*x/(x*x+0.0111636921))";
    gdEdx_muon = new TF1("gdEdx_muon",Form("[0]/%s * (TMath::Log(%s/([1]*(1-%s))) - 2*%s)",beta2.Data(),
						beta2.Data(),beta2.Data(),beta2.Data()),pmin,pmax);
    gdEdx_muon->SetParNames("A","B");
    gdEdx_muon->SetParameters(gdEdx_pion->GetParameter(0),gdEdx_pion->GetParameter(1));

    ///// electron dEdx, mass = 0.0005109989461 (mass^2 = 2.611199269e-07)
    beta2 = "(x*x/(x*x+2.61119927e-07))";
    gdEdx_elec = new TF1("gdEdx_elec",Form("[0]/%s * (TMath::Log(%s/([1]*(1-%s))) - 2*%s)",beta2.Data(),
						beta2.Data(),beta2.Data(),beta2.Data()),pmin,pmax);
    gdEdx_elec->SetParNames("A","B");
    gdEdx_elec->SetParameters(gdEdx_pion->GetParameter(0),gdEdx_pion->GetParameter(1));

    ///// deuteron dEdx, mass = 1.875612928 (mass^2 = 3.51792385568)
    //// alpha dEdx, mass= 3.727379378 (mass^2 = 13.893357)
    beta2 = "(x*x/(x*x+3.517923856))";
    //beta2 = "(x*x/(x*x+13.893357))";
    gdEdx_deut = new TF1("gdEdx_deut",Form("[0]/%s * (TMath::Log(%s/([1]*(1-%s))) - 2*%s)",beta2.Data(),
						beta2.Data(),beta2.Data(),beta2.Data()),pmin,pmax);
    gdEdx_deut->SetParNames("A","B");
    gdEdx_deut->SetParameters(gdEdx_pion->GetParameter(0),gdEdx_pion->GetParameter(1));
  }


  TString dir = "LambAnaModule";  

  //TString fineName=Form("d0ana_run%d_p42_ap9.root",runid);
  //TString fineName=Form("d0ana_run%d_p42_ap1.root",runid);

  TString fineName=Form("d0ana_run1234_p5_ap6.root");
  // TString fineName2=Form("d0ana_run1395_p35_ap8.root",runid);

  TFile* InputFile = TFile::Open(Form("../root_files/%s",fineName.Data()));
  if(!InputFile)return;
  // TFile* InputFile2 = TFile::Open(Form("../root_files/xela_p35/%s",fineName2.Data()));
  // if(!InputFile2)return;

  
  //TH2F* hr = (TH2F*)InputFile->Get(Form("%s/hdEdx_mtpc",dir.Data()));
  TH2F* hr = (TH2F*)InputFile->Get(Form("%s/hProtons",dir.Data()));
  TH2F* hr2 = (TH2F*)InputFile->Get(Form("%s/hPbar",dir.Data()));
  // TH2F* Hr = (TH2F*)InputFile2->Get(Form("%s/hProtons",dir.Data()));
  // TH2F* Hr2 = (TH2F*)InputFile2->Get(Form("%s/hPbar",dir.Data()));

  hr->Add(hr2);
  // hr->Add(Hr);
  // hr->Add(Hr2);

  TH2F* hrplot = (TH2F*)hr->Clone("hrplot");
  //TH2F* hrplot = new TH2F(*hr);

  //DyVsY(nSlices,hr);
  DyVsY2(nSlices,hr);
  
  TCanvas* c1 = new TCanvas ("dEdx","dEdx",100,100,650,560);

  c1->Divide(1,1);
  c1->cd(1); 
 
  TList* list = hr->GetListOfFunctions();
  TMultiGraph* mg1 = new TMultiGraph();
  
  TGraphErrors* mean     = (TGraphErrors*)list->FindObject(Form("Mean"));
  TGraphErrors* Sigma    = (TGraphErrors*)list->FindObject(Form("Sigma"));
  TGraphErrors* SigmaVsMean    = (TGraphErrors*)list->FindObject(Form("SigmaVsMean"));
  TGraphErrors* SigmaVsBeta    = (TGraphErrors*)list->FindObject(Form("SigmaVsBeta"));
  TGraphErrors* SigmaVsBetaProt    = (TGraphErrors*)list->FindObject(Form("SigmaVsBetaProt"));
  //TGraphErrors* meanPsigma     = (TGraphErrors*)list->FindObject(Form("MeanPSigma"));
  //TGraphErrors* meanMsigma     = (TGraphErrors*)list->FindObject(Form("MeanMSigma"));

  //mg1->Add(mean);
  mg1->Add(SigmaVsBeta);
  //mg1->Add(SigmaVsBetaProt);
  //mg1->Add(SigmaVsMean);
  //mg1->Add(Sigma);


  hrplot->GetXaxis()->SetTitle("p [GeV/c]");
  hrplot->GetYaxis()->SetTitle("dEdx [arbitrary]");

  hrplot->GetListOfFunctions()->Add(mean);
  //hrplot->GetListOfFunctions()->Add(SigmaVsMean);
  //hrplot->GetListOfFunctions()->Add(meanMsigma);


  hrplot->Draw();

  
  TCanvas* c2 = new TCanvas ("sigmaVsmean","",200,200,650,560);
  
  c2->Divide(1,1);
  c2->cd(1); 
  mg1->Draw("a");
  
  
}

//_______________________________________________________________________
// void DyVsY(Int_t nSlices, TH2F* h2)
// {
//   // this is general method to perform pid on m2 versus momentum
//   // plot. The results are saved in file at path location.
//   // path should point to the particular (Tof/Rich) selector
//   // pid directory
  
//   Int_t countsInMostSideBin = 100; // this has to be adjusted
//   Int_t binsToMerge = 1;
  
//   /*  
//   Int_t ibin = h2->GetYaxis()->FindBin(0.5);
//   TH1F* pdist = (TH1F*)h2->ProjectionX("pdist",ibin,h2->GetYaxis()->GetNbins());
  
//   TList* list = h2->GetListOfFunctions();
  
//   Float_t pmin,pmax;
//   Float_t fraction = (Float_t)countsInMostSideBin/pdist->GetMaximum();
//   FindFitRange(pdist,pmin,pmax,fraction);
  
//   //Int_t nSlices = (pdist->FindBin(pmax)-pdist->FindBin(pmin))/binsToMerge;
//   */

//   TH1D* slice[100];
 
//   Float_t pmin = gPmin;
//   Float_t pmax = gPmax;
  

//   cout<<"DyVsY: number of slices is: "<<nSlices<<endl;
  
  
//   Float_t pstep = (pmax-pmin)/nSlices;
//   cout<<"initial hbw/2 is "<<pstep/2.<<endl;
    

//   for ( Int_t i = 0 ; i < nSlices ; i++ ) {
//     Float_t pL = pmin + pstep * i;
//     Float_t pU = pmin + pstep * (i+1);
//     Int_t binL = h2->GetXaxis()->FindBin(pL);
//     Int_t binU = h2->GetXaxis()->FindBin(pU); 
//     Float_t p = 0.5*(pL+pU);
    
//     slice[i] = (TH1D*)h2->ProjectionY(Form("slice%d",i),binL,binU);
//     slice[i] -> SetTitle(Form("p range: %f - %f GeV/c",pL,pU));
//     slice[i]->Sumw2();
//     slice[i]->Rebin(binsToMerge);
//     slice[i]->Scale(1./(Float_t)binsToMerge);
//     cout<<"i="<<i<<"  pL="<<pL<<"   pU="<<pU<<" max="<<slice[i]->GetMaximum()<<endl;
//     if(slice[i]->GetMaximum()<3){
//       slice[i]->Rebin(2);
//       slice[i]->Scale(1./(Float_t)2);
//     }
//   }
  

//   cout<<"Starting initial fit of slices"<<endl;

//   Float_t momPart[100];
//   Float_t momPartErr[100];
//   Float_t meanPart[100];
//   Float_t meanPartErr[100];
//   Float_t sigmaPart[100];
//   Float_t sigmaPartProt[100];
//   Float_t sigmaPartErr[100];
//   Float_t betaPart[100];
//   Float_t betaPartProt[100];
//   Float_t betaPartErr[100];
  
//   Float_t xmin;
//   Float_t xmax;
//   for ( Int_t i = 0 ; i < nSlices ; i++ ) {
//     Float_t pL = pmin + pstep * i;
//     Float_t pU = pmin + pstep * (i+1);
//     Int_t binL = h2->GetXaxis()->FindBin(pL);
//     Int_t binU = h2->GetXaxis()->FindBin(pU); 
//     Float_t p = 0.5*(pL+pU);
    
//     if(i==0)xmin = pL;
//     if(i==(nSlices-1))xmax = pU;

//     TH1D* h = slice[i];
    
//     Int_t maxBin = h->GetMaximumBin();
//     Float_t cons = h->GetMaximum();
//     Float_t mean = h->GetBinCenter(maxBin);
//     Float_t sigm = 1.3;
    
//     Float_t consB = h->GetBinContent(maxBin/2);    
//     Float_t meanB = h->GetBinCenter(maxBin);
//     Float_t sigmB = 100*sigm;
    
//     Double_t params[6];
    
//     params[0] = cons;
//     params[1] = mean;
//     params[2] = sigm;
//     params[3] = consB;
//     params[4] = meanB;
//     params[5] = sigmB;
    
//     cout<<"i="<<i<<"  pL="<<pL<<"   pU="<<pU<<endl;

//     //cout<<"params: "<<cons<<" "<< mean<<" "<<sigm<<" "<<consB<<" "<< meanB<<" "<<sigmB<<endl;

//     //TF1* func = GetDoubleGaussFit(slice[i],params);

//     //TF1* func = GetTripleGaussFit(slice[i],0.5*(pL+pU),params);
//     //TF1* func = GetQuintupleGaussFit(slice[i],0.5*(pL+pU),params);
//     TF1* func = GetGaussFitPol(slice[i],0.5*(pL+pU),params);
    
//     Bool_t drawSlices=kTRUE;
//     if(drawSlices){
//       cnv[i] = new TCanvas (Form("slice%d",i),Form("slice%d",i),5*i,5*i,400,560);
//       cnv[i]->SetTopMargin(0.02);
//       cnv[i]->SetRightMargin(0.05);
//       cnv[i]->SetLeftMargin(0.16);
//       cnv[i]->SetBottomMargin(0.15);
//       cnv[i]->Divide(1,1);    
      
//       DrawSlice(slice[i],cnv[i]);
//     }
    
//     //cout<<" func="<<func<<endl;

//     //cout<<"PidMass2: islice="<<i<<" p="<<p<<endl;
//     //cout<<"mass2 (p,k,prot): "<<meanPion<<" "<<meanKaon<<" "<<meanProt<<endl;
       
//     momPart[i]      = p;
//     momPartErr[i]   = p*0.005; // this should depend on field settings
//     meanPart[i]     = func->GetParameter(1);
//     meanPartErr[i]  = TMath::Abs(func->GetParError(1));
//     sigmaPart[i]    = TMath::Abs(func->GetParameter(2));
//     sigmaPartErr[i]    = TMath::Abs(func->GetParError(2));
//     betaPart[i] = p/TMath::Sqrt(p*p+0.01947983);
//     betaPartErr[i] = betaPart[i]*0.000001;

//     //if(i==0)sigmaPartErr[i]=sigmaPartErr[i] * 0.3;
//   }


//   cout<<"xmin="<<xmin<<"  xmax="<<xmax<<endl;

//   h2->GetListOfFunctions()->Clear();
    
//   TGraphErrors* Mean = new TGraphErrors(nSlices,momPart,meanPart,momPartErr,meanPartErr);        
//   TGraphErrors* Sigma = new TGraphErrors(nSlices,momPart,sigmaPart,momPartErr,sigmaPartErr);
//   TGraphErrors* SigmaVsMean = new TGraphErrors(nSlices,meanPart,sigmaPart,meanPartErr,sigmaPartErr);
//   TGraphErrors* SigmaVsBeta = new TGraphErrors(nSlices,betaPart,sigmaPart,betaPartErr,sigmaPartErr);
//   TGraphErrors* SigmaVsBetaProt = new TGraphErrors(nSlices,betaPartProt,sigmaPartProt,betaPartErr,sigmaPartErr);
  
//   Mean->SetName(Form("Mean"));
//   Sigma->SetName(Form("Sigma"));
//   SigmaVsMean->SetName("SigmaVsMean");
//   SigmaVsBeta->SetName("SigmaVsBeta");
//   SigmaVsBetaProt->SetName("SigmaVsBetaProt");
  
//   Mean->SetMarkerStyle(20);
//   Mean->SetMarkerSize(0.05);
//   Mean->SetMarkerColor(2);
//   Mean->SetLineColor(2);
//   Mean->SetLineWidth(3);

//   Sigma->SetMarkerStyle(20);
//   Sigma->SetMarkerSize(0.05);
//   Sigma->SetMarkerColor(4);
//   Sigma->SetLineColor(4);
//   Sigma->SetLineWidth(3);
  
//   SigmaVsMean->SetMarkerStyle(20);
//   SigmaVsMean->SetMarkerSize(0.05);
//   SigmaVsMean->SetMarkerColor(4);
//   SigmaVsMean->SetLineColor(4);
//   SigmaVsMean->SetLineWidth(3);

//   SigmaVsBeta->SetMarkerStyle(20);
//   SigmaVsBeta->SetMarkerSize(0.05);
//   SigmaVsBeta->SetMarkerColor(6);
//   SigmaVsBeta->SetLineColor(4);
//   SigmaVsBeta->SetLineWidth(3);

//   SigmaVsBetaProt->SetMarkerStyle(20);
//   SigmaVsBetaProt->SetMarkerSize(0.05);
//   SigmaVsBetaProt->SetMarkerColor(2);
//   SigmaVsBetaProt->SetLineColor(2);
//   SigmaVsBetaProt->SetLineWidth(3);
  
//   /*
//   TF1* meanFit = new TF1("meanFit","[0] + [1]*x + [2]*x*x + [3]*x*x*x",pmin,pmax);
//   meanFit->SetParNames("p0","p1","p2","p3");
//   meanFit->SetLineColor(2);
//   meanFit->SetLineWidth(4);
//   meanFit->SetParameters(-1.25,-0.016,-2.6e-5,0); 
//   Mean->Fit(meanFit,"w","",xmin-0.1,xmax+0.1);
//   */

//   TString beta2 = "(x*x/(x*x+0.01947983))";
//   TF1* meanFit_pion = new TF1("meanFit_pion",Form("[0]/%s * (TMath::Log(%s/([1]*(1-%s))) - %s)",beta2.Data(),
// 						  beta2.Data(),beta2.Data(),beta2.Data()),pmin,pmax);
//   meanFit_pion->SetParNames("A","B");
//   meanFit_pion->SetLineColor(3);
//   meanFit_pion->SetLineWidth(4);
//   meanFit_pion->SetParameters(0.07,2.8e-6); 
//   Mean->Fit(meanFit_pion,"w","",pmin-0.1,pmax+0.1);


//   beta2 = "(x*x/(x*x+0.880354496))";
//   TF1* meanFit_prot = new TF1("meanFit_prot",Form("[0]/%s * (TMath::Log(%s/([1]*(1-%s))) - %s)",beta2.Data(),
// 						  beta2.Data(),beta2.Data(),beta2.Data()),pmin,pmax);
//   meanFit_prot->SetParNames("A","B");
//   meanFit_prot->SetLineColor(4);
//   meanFit_prot->SetLineWidth(4);
//   meanFit_prot->SetParameters(meanFit_kaon->GetParameter(0),meanFit_kaon->GetParameter(1));

//   beta2 = "(x*x/(x*x+0.0111636921))";
//   TF1* meanFit_muon = new TF1("meanFit_muon",Form("[0]/%s * (TMath::Log(%s/([1]*(1-%s))) - %s)",beta2.Data(),
// 						  beta2.Data(),beta2.Data(),beta2.Data()),pmin,pmax);
//   meanFit_muon->SetParNames("A","B");
//   meanFit_muon->SetLineColor(8);
//   meanFit_muon->SetLineWidth(4);
//   meanFit_muon->SetParameters(meanFit_kaon->GetParameter(0),meanFit_kaon->GetParameter(1));
 
//   beta2 = "(x*x/(x*x+2.61119927e-07))";
//   TF1* meanFit_elec = new TF1("meanFit_elec",Form("[0]/%s * (TMath::Log(%s/([1]*(1-%s))) - %s)",beta2.Data(),
// 						  beta2.Data(),beta2.Data(),beta2.Data()),pmin,pmax);
//   meanFit_elec->SetParNames("A","B");
//   meanFit_elec->SetLineColor(kMagenta+2);
//   meanFit_elec->SetLineWidth(4);
//   meanFit_elec->SetParameters(meanFit_kaon->GetParameter(0),meanFit_kaon->GetParameter(1));

//   beta2 = "(x*x/(x*x+3.517923856))";
//   //beta2 = "(x*x/(x*x+13.893357))";
//   TF1* meanFit_deut = new TF1("meanFit_deut",Form("[0]/%s * (TMath::Log(%s/([1]*(1-%s))) - %s)",beta2.Data(),
// 						  beta2.Data(),beta2.Data(),beta2.Data()),pmin,pmax);
//   meanFit_deut->SetParNames("A","B");
//   meanFit_deut->SetLineColor(kRed+1);
//   meanFit_deut->SetLineWidth(4);
//   meanFit_deut->SetParameters(meanFit_kaon->GetParameter(0),meanFit_kaon->GetParameter(1));


  
//   // add graphs and fits to m2 histogram
  
//   //Mean->GetListOfFunctions()->Add(meanFit_pion);
//   Mean->GetListOfFunctions()->Add(meanFit_pion);
//   Mean->GetListOfFunctions()->Add(meanFit_prot);
//   //Mean->GetListOfFunctions()->Add(meanFit_muon);
//   Mean->GetListOfFunctions()->Add(meanFit_elec);
  
//   h2->GetListOfFunctions()->Add(Mean);
//   h2->GetListOfFunctions()->Add(Sigma);
//   h2->GetListOfFunctions()->Add(SigmaVsMean);
//   h2->GetListOfFunctions()->Add(SigmaVsBeta);
//   h2->GetListOfFunctions()->Add(SigmaVsBetaProt);

//   cout<<"functions added"<<endl;
  
  
// }

//_______________________________________________________________________
void DyVsY2(Int_t nSlicesInit, TH2F* h2)
{
  // this is general method to perform pid on m2 versus momentum
  // plot. The results are saved in file at path location.
  // path should point to the particular (Tof/Rich) selector
  // pid directory
  

  TH1D* slice[100];
 
  Float_t pmin = gPmin;
  Float_t pmax = gPmax;
  

  // cout<<"DyVsY2: number of slices is: "<<nSlicesInit<<endl;
  
  
  Float_t pstepInit = (pmax-pmin)/nSlicesInit;
  cout<<"initial hbw/2 is "<<pstepInit/2.<<endl;
    
  Float_t pstep = pstepInit;

  Float_t pL = pmin;
  Float_t pU;
  Int_t nSlices=0;

  // 
  Float_t pstep_array[] = {3.3, 3.3, 3.3, 3.3, 3.3, 3.3, 3.3, 3.3, 4, 5, 6};
  nSlicesInit = sizeof(pstep_array)/4;

  for ( Int_t i = 0 ; i < nSlicesInit ; i++ ) {
    if(i>0)pL = pU;
    // if(pL>7)pstep = 2*pstepInit;
    pU = pL + pstep_array[i];
    
    cout<<"DyVsY2: number of slices is: "<<nSlicesInit<<endl;
    
    //pU = pmin + pstep * (i+1);
    Int_t binL = h2->GetXaxis()->FindBin(pL);
    Int_t binU = h2->GetXaxis()->FindBin(pU); 
    Float_t p = 0.5*(pL+pU);
    // if(p>pmax)break;    
    nSlices++;
    slice[i] = (TH1D*)h2->ProjectionY(Form("slice%d",i),binL,binU);
    slice[i] -> SetTitle(Form("p range: %f - %f GeV/c",pL,pU));
    slice[i]->Sumw2();
    //slice[i]->Rebin(binsToMerge);
    //slice[i]->Scale(1./(Float_t)binsToMerge);
    cout<<"i="<<i<<"  pL="<<pL<<"   pU="<<pU<<" max="<<slice[i]->GetMaximum()<<endl;
    if(slice[i]->GetMaximum()<3){
      slice[i]->Rebin(2);
      //slice[i]->Scale(1./(Float_t)2);
    }
  }
  
  pmax = pU;

  cout<<"Starting initial fit of slices"<<endl;

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
  
  for ( Int_t i = 0 ; i < nSlicesInit ; i++ ) {
    //Float_t pL = pmin + pstep * i;
    //Float_t pU = pmin + pstep * (i+1);
    if(i>0)pL = pU;
    // if(pL>7)pstep = 2*pstepInit;
    pU = pL + pstep_array[i];
    Int_t binL = h2->GetXaxis()->FindBin(pL);
    Int_t binU = h2->GetXaxis()->FindBin(pU); 
    Float_t p = 0.5*(pL+pU);
    // if(p>pmax)break;    

    TH1D* h = slice[i];
    
    Int_t maxBin = h->GetMaximumBin();
    Float_t cons = h->GetMaximum();
    Float_t mean = h->GetBinCenter(maxBin);
    Float_t sigm = 1.3;
    
    Float_t consB = h->GetBinContent(maxBin/2);    
    Float_t meanB = h->GetBinCenter(maxBin);
    Float_t sigmB = 100*sigm;
    
    Double_t params[6];
    
    params[0] = cons;
    params[1] = mean;
    params[2] = sigm;
    params[3] = consB;
    params[4] = meanB;
    params[5] = sigmB;
    
    cout<<"i="<<i<<"  pL="<<pL<<"   pU="<<pU<<endl;

    TF1* func = GetGaussFitPol(slice[i],0.5*(pL+pU),params);
    
    Bool_t drawSlices=kTRUE;
    if(drawSlices){
      cnv[i] = new TCanvas (Form("slice%d",i),Form("slice%d",i),5*i,5*i,400,560);
      cnv[i]->SetTopMargin(0.02);
      cnv[i]->SetRightMargin(0.05);
      cnv[i]->SetLeftMargin(0.16);
      cnv[i]->SetBottomMargin(0.15);
      cnv[i]->Divide(1,1);    
      
      DrawSlice(slice[i],cnv[i]);
    }
    
       
    momPart[i]      = p;
    momPartErr[i]   = p*0.005; // this should depend on field settings
    meanPart[i]     = func->GetParameter(1);
    meanPartErr[i]  = TMath::Abs(func->GetParError(1));
    sigmaPart[i]    = TMath::Abs(func->GetParameter(2));
    sigmaPartErr[i]    = TMath::Abs(func->GetParError(2));
    betaPart[i] = p/TMath::Sqrt(p*p+0.880354496);
    betaPartErr[i] = betaPart[i]*0.000001;

  }


  h2->GetListOfFunctions()->Clear();
    
  TGraphErrors* Mean = new TGraphErrors(nSlices,momPart,meanPart,momPartErr,meanPartErr);        
  TGraphErrors* Sigma = new TGraphErrors(nSlices,momPart,sigmaPart,momPartErr,sigmaPartErr);
  TGraphErrors* SigmaVsMean = new TGraphErrors(nSlices,meanPart,sigmaPart,meanPartErr,sigmaPartErr);
  TGraphErrors* SigmaVsBeta = new TGraphErrors(nSlices,betaPart,sigmaPart,betaPartErr,sigmaPartErr);
  TGraphErrors* SigmaVsBetaProt = new TGraphErrors(nSlices,betaPartProt,sigmaPartProt,betaPartErr,sigmaPartErr);
  
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
  


  TString beta2 = "(x*x/(x*x+0.880354496))";

  // Tu parametryzujemy funkcję, ktora fitujemy Bethe-Blocha

  // Odtworzyć funkcję Betego Blocha
  // popatrzyć czy parametryzacja jest poprawna
  TF1* meanFit_prot = new TF1("meanFit_prot",Form("[0]/%s * (TMath::Log(%s/([1]*(1-%s))) - 2*%s)",beta2.Data(),
						  beta2.Data(),beta2.Data(),beta2.Data()),pmin,pmax);
  meanFit_prot->SetParNames("A","B");
  meanFit_prot->SetLineColor(4);
  meanFit_prot->SetLineWidth(4);
  meanFit_prot->SetParameters(gdEdx_prot->GetParameter(0),gdEdx_prot->GetParameter(1));
  Mean->Fit(meanFit_prot,"w","",pmin-0.1,pmax+0.1);
 
// Te piony są z protonu?
  beta2 = "(x*x/(x*x+0.01947983))";
  TF1* meanFit_pion = new TF1("meanFit_pion",Form("[0]/%s * (TMath::Log(%s/([1]*(1-%s))) - 2*%s)",beta2.Data(),
						  beta2.Data(),beta2.Data(),beta2.Data()),pmin,pmax);
  meanFit_pion->SetParNames("A","B");
  meanFit_pion->SetLineColor(3);
  meanFit_pion->SetLineWidth(4);
  meanFit_pion->SetParameters(meanFit_prot->GetParameter(0),meanFit_prot->GetParameter(1));

// te piony z oryginalnego fitu do pionow
  TF1* meanFit_pionPion = new TF1("meanFit_pionPion",Form("[0]/%s * (TMath::Log(%s/([1]*(1-%s))) - 2*%s)",beta2.Data(),
							  beta2.Data(),beta2.Data(),beta2.Data()),pmin,pmax);
  meanFit_pionPion->SetParNames("A","B");
  meanFit_pionPion->SetLineColor(6);
  meanFit_pionPion->SetLineWidth(4);
  meanFit_pionPion->SetParameters(0.0599219,5.62172e-07);
  //meanFit_pionPion->SetParameters(0.0779697,6.61051e-05);



  beta2 = "(x*x/(x*x+0.0111636921))";
  TF1* meanFit_muon = new TF1("meanFit_muon",Form("[0]/%s * (TMath::Log(%s/([1]*(1-%s))) - 2*%s)",beta2.Data(),
						  beta2.Data(),beta2.Data(),beta2.Data()),pmin,pmax);
  meanFit_muon->SetParNames("A","B");
  meanFit_muon->SetLineColor(8);
  meanFit_muon->SetLineWidth(4);
  meanFit_muon->SetParameters(meanFit_pion->GetParameter(0),meanFit_pion->GetParameter(1));
 
  beta2 = "(x*x/(x*x+2.61119927e-07))";
  TF1* meanFit_elec = new TF1("meanFit_elec",Form("[0]/%s * (TMath::Log(%s/([1]*(1-%s))) - 2*%s)",beta2.Data(),
						  beta2.Data(),beta2.Data(),beta2.Data()),pmin,pmax);
  meanFit_elec->SetParNames("A","B");
  meanFit_elec->SetLineColor(kMagenta+2);
  meanFit_elec->SetLineWidth(4);
  meanFit_elec->SetParameters(meanFit_pion->GetParameter(0),meanFit_pion->GetParameter(1));

  beta2 = "(x*x/(x*x+3.517923856))";
  TF1* meanFit_deut = new TF1("meanFit_deut",Form("[0]/%s * (TMath::Log(%s/([1]*(1-%s))) - 2*%s)",beta2.Data(),
						  beta2.Data(),beta2.Data(),beta2.Data()),pmin,pmax);
  meanFit_deut->SetParNames("A","B");
  meanFit_deut->SetLineColor(kRed+1);
  meanFit_deut->SetLineWidth(4);
  meanFit_deut->SetParameters(meanFit_pion->GetParameter(0),meanFit_pion->GetParameter(1));


  
  // add graphs and fits to m2 histogram
  
  Mean->GetListOfFunctions()->Add(meanFit_pion);
  Mean->GetListOfFunctions()->Add(meanFit_pionPion);
  Mean->GetListOfFunctions()->Add(meanFit_prot);
  //Mean->GetListOfFunctions()->Add(meanFit_elec);
  
  h2->GetListOfFunctions()->Add(Mean);
  h2->GetListOfFunctions()->Add(Sigma);
  h2->GetListOfFunctions()->Add(SigmaVsMean);
  h2->GetListOfFunctions()->Add(SigmaVsBeta);
  h2->GetListOfFunctions()->Add(SigmaVsBetaProt);

  cout<<"functions added"<<endl;
  
  
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
TF1* GetDoubleGaussFit(TH1D* h,const Double_t* params)
{
  // h - histograms to fit
  // params - pointer to list of parameters 
  // n - number of Gaussians to consider
  // The method create function fit it to histogram and return.

  TString g1 = "[0]*exp(-0.5*((x-[1])/[2])^2)";
  TString g2 = "[3]*exp(-0.5*((x-[4])/[5])^2)";
  TString name = h->GetName();
  name.Append("_Fit");
  cout<<"name: "<<name.Data()<<endl;
  TF1* func = new TF1(name.Data(),Form("%s + %s",g1.Data(),g2.Data()),-10.,10.);  
  //TF1* func = new TF1(name.Data(),Form("%s + %s",g1.Data(),g2.Data()),-10.,10.,6);  

  func->SetParameters(params);
  func->SetParLimits(0,0,1000);
  func->SetParLimits(2,0.5,2.0);

  func->SetParLimits(3,0,500);
  func->SetParLimits(4,-5,5);
  func->SetParLimits(5,10,50);

  func->SetParNames("const","mean","sigma","constB","meanB","sigmaB");
  
  //limits on mean, sigma and constant  
  Float_t rmin = -19.;
  Float_t rmax =  19.;
 
  //h->Rebin(4);
  //h->Scale(0.25);    
  
  //h->Fit(func,"Q0","0",rmin,rmax);
  h->Fit(func,"w0","",rmin,rmax);
  h->GetListOfFunctions()->Add(func);
  
  return func;
}

//______________________________________________________________
TF1* GetTripleGaussFit(TH1D* h,Float_t p,const Double_t* params)
{
  // h - histograms to fit
  // params - pointer to list of parameters 
  // n - number of Gaussians to consider
  // The method create function fit it to histogram and return.

  TString g1 = "[0]*exp(-0.5*((x-[1])/[2])^2)";
  TString g2 = "[3]*exp(-0.5*((x-[4])/[5])^2)";
  TString g3 = "[6]*exp(-0.5*((x-[7])/[8])^2)";
  TString name = h->GetName();
  name.Append("_Fit");
  cout<<"name: "<<name.Data()<<endl;
  TF1* func = new TF1(name.Data(),Form("%s + %s + %s",g1.Data(),g2.Data(),g3.Data()),0.7,1.8);  
  //TF1* func = new TF1(name.Data(),Form("%s + %s",g1.Data(),g2.Data()),-10.,10.,6);  

  Float_t dd = 0.01;
  Float_t sig = 0.06336;
  func->SetParameters(params);
  func->SetParLimits(0,0,100000);
  func->SetParLimits(1,gdEdx_pion->Eval(p)-dd,gdEdx_pion->Eval(p)+2*dd);
  //func->SetParLimits(2,sig-0.009,sig+0.002);
  func->SetParLimits(2,sig-0.009,sig+0.01);

  func->SetParLimits(3,0,50000);
  //func->FixParameter(3,0);
  func->SetParLimits(4,gdEdx_kaon->Eval(p)-dd,gdEdx_kaon->Eval(p)+dd);
  func->SetParLimits(5,sig-0.009,sig+0.002);

  func->SetParLimits(6,0,50000);
  func->SetParLimits(7,gdEdx_prot->Eval(p)-10*dd,gdEdx_prot->Eval(p)+dd);
  //func->SetParLimits(8,sig-0.009,sig+0.002);
  func->SetParLimits(8,sig-0.009,sig+0.01);

  func->SetParNames("constPion","meanPion","sigmaPion","constKaon","meanKaon","sigmaKaon",
		    "constP","meanP","sigmaP");
  
  //limits on mean, sigma and constant  
  Float_t rmin = 0.7;
  Float_t rmax =  1.8;
 
  //h->Rebin(4);
  //h->Scale(0.25);    
  
  //h->Fit(func,"Q0","0",rmin,rmax);
  h->Fit(func,"w0","",rmin,rmax);
  h->GetListOfFunctions()->Add(func);
  
  // add componet functions
  TF1* fpions = new TF1("pions",Form("%s",g1.Data()),0.7,1.8);  
  fpions->SetParameters(func->GetParameter(0),func->GetParameter(1),func->GetParameter(2));
  fpions->SetLineColor(3);
  
  TF1* fkaons = new TF1("kaons",Form("%s",g1.Data()),0.7,1.8);  
  fkaons->SetParameters(func->GetParameter(3),func->GetParameter(4),func->GetParameter(5));
  fkaons->SetLineColor(6);
  
  TF1* fprots = new TF1("prots",Form("%s",g1.Data()),0.7,1.8);  
  fprots->SetParameters(func->GetParameter(6),func->GetParameter(7),func->GetParameter(8));
  fprots->SetLineColor(4);
  
  h->GetListOfFunctions()->Add(fpions);
  h->GetListOfFunctions()->Add(fkaons);
  h->GetListOfFunctions()->Add(fprots);

 return func;
}

//______________________________________________________________
TF1* GetQuintupleGaussFit(TH1D* h,Float_t p,const Double_t* params)
{
  // h - histograms to fit
  // params - pointer to list of parameters 
  // n - number of Gaussians to consider
  // The method create function fit it to histogram and return.

  TString g1 = "[0]*exp(-0.5*((x-[1])/[2])^2)";
  TString g2 = "[3]*exp(-0.5*((x-[4])/[5])^2)";
  TString g3 = "[6]*exp(-0.5*((x-[7])/[8])^2)";
  TString g4 = "[9]*exp(-0.5*((x-[10])/[11])^2)";
  TString g5 = "[12]*exp(-0.5*((x-[13])/[14])^2)";
  TString name = h->GetName();
  name.Append("_Fit");
  cout<<"name: "<<name.Data()<<endl;
  TF1* func = new TF1(name.Data(),Form("%s + %s + %s + %s + %s",g1.Data(),g2.Data(),g3.Data(),g4.Data(),g5.Data()),0.7,1.8);  
 

  Float_t dd = 0.01;
  //Float_t sig = 0.06336;
  Float_t sig = 0.06;
  Float_t sigKaon = 0.054;
  func->SetParameters(params);
  // initial setup of params
  func->SetParLimits(0,0,900000);
  //func->SetParLimits(1,gdEdx_pion->Eval(p)-dd,gdEdx_pion->Eval(p)+2*dd);
  func->SetParLimits(1,gdEdx_pion->Eval(p)-3*dd,gdEdx_pion->Eval(p)+3*dd);
  func->SetParLimits(2,sig-0.005,sig+0.005);

  func->SetParLimits(3,0,90000);
  //func->FixParameter(3,0);
  //func->SetParLimits(4,gdEdx_kaon->Eval(p)-dd,gdEdx_kaon->Eval(p)+dd);
  func->SetParLimits(4,gdEdx_kaon->Eval(p)-3*dd,gdEdx_kaon->Eval(p)+3*dd);
  func->SetParLimits(5,sigKaon-0.001,sigKaon+0.001);

  func->SetParLimits(6,0,900000);
  //func->SetParLimits(7,gdEdx_prot->Eval(p)-10*dd,gdEdx_prot->Eval(p)+dd);
  func->SetParLimits(7,gdEdx_prot->Eval(p)-5*dd,gdEdx_prot->Eval(p)+3*dd);
  func->SetParLimits(8,sig-0.004,sig+0.004);

  func->SetParLimits(9,0,90000);
  func->SetParLimits(10,gdEdx_pion->Eval(p)+3*sig-3*dd,gdEdx_pion->Eval(p)+3*sig+10*dd);
  func->SetParLimits(11,sig-0.001,sig+0.001);

  func->SetParLimits(12,0,9000);
  func->SetParLimits(13,gdEdx_prot->Eval(p)-3*sig-10*dd,gdEdx_prot->Eval(p)-3*sig+3*dd);
  func->SetParLimits(14,sig-0.001,sig+0.001);

  func->SetParNames("constPion","meanPion","sigmaPion","constKaon","meanKaon","sigmaKaon",
		    "constP","meanP","sigmaP");
  func->SetParName(9,"constElect");  
  func->SetParName(10,"meanElect");  
  func->SetParName(11,"sigmaElect");  
  func->SetParName(12,"constDeut");  
  func->SetParName(13,"meanDeut");  
  func->SetParName(14,"sigmaDeut");  

  //limits on mean, sigma and constant  
  Float_t rmin = 0.7;
  Float_t rmax =  1.8;
 
  //h->Rebin(4);
  //h->Scale(0.25);    
  
  //h->Fit(func,"Q0","0",rmin,rmax);
  h->Fit(func,"w0","",rmin,rmax);
  
  ////////// final setup of params
  func->SetParLimits(1,func->GetParameter(1)-3*dd,func->GetParameter(1)+3*dd);
  func->SetParLimits(2,func->GetParameter(2)-0.005,func->GetParameter(2)+0.005);

  Float_t sigKf = (func->GetParameter(2) + func->GetParameter(8))/2.;
  func->SetParLimits(4,func->GetParameter(4)-3*dd,func->GetParameter(4)+3*dd);
  func->SetParLimits(5,sigKf-0.001,sigKf+0.001);

  func->SetParLimits(7,func->GetParameter(7)-3*dd,func->GetParameter(7)+3*dd);
  func->SetParLimits(8,func->GetParameter(8)-0.005,func->GetParameter(8)+0.005);

  func->SetParLimits(10,func->GetParameter(10)-3*dd,func->GetParameter(10)+10*dd);
  //func->SetParLimits(10,1.55-5*dd,1.55+2*dd);
  func->SetParLimits(11,func->GetParameter(2)-0.005,func->GetParameter(2)+0.0001);

  func->SetParLimits(13,func->GetParameter(13)-10*dd,func->GetParameter(13)+3*dd);
  func->SetParLimits(14,func->GetParameter(8)-0.0001,func->GetParameter(8)+0.005);
  
  h->Fit(func,"w0","",rmin,rmax);


  h->GetListOfFunctions()->Add(func);
  
  // add componet functions
  TF1* fpions = new TF1("pions",Form("%s",g1.Data()),0.7,1.8);  
  fpions->SetParameters(func->GetParameter(0),func->GetParameter(1),func->GetParameter(2));
  fpions->SetLineColor(3);
  
  TF1* fkaons = new TF1("kaons",Form("%s",g1.Data()),0.7,1.8);  
  fkaons->SetParameters(func->GetParameter(3),func->GetParameter(4),func->GetParameter(5));
  fkaons->SetLineColor(6);
  
  TF1* fprots = new TF1("prots",Form("%s",g1.Data()),0.7,1.8);  
  fprots->SetParameters(func->GetParameter(6),func->GetParameter(7),func->GetParameter(8));
  fprots->SetLineColor(4);

  TF1* fmuons = new TF1("muons",Form("%s",g1.Data()),0.7,1.8);  
  fmuons->SetParameters(func->GetParameter(9),func->GetParameter(10),func->GetParameter(11));
  fmuons->SetLineColor(8);

  TF1* fdeuts = new TF1("deuts",Form("%s",g1.Data()),0.7,1.8);  
  fdeuts->SetParameters(func->GetParameter(12),func->GetParameter(13),func->GetParameter(14));
  fdeuts->SetLineColor(kMagenta+2);
  
  h->GetListOfFunctions()->Add(fpions);
  h->GetListOfFunctions()->Add(fkaons);
  h->GetListOfFunctions()->Add(fprots);
  h->GetListOfFunctions()->Add(fmuons);
  h->GetListOfFunctions()->Add(fdeuts);

 return func;
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

  func->SetParLimits(0,1,2000);
  func->SetParLimits(1,gdEdx_prot->Eval(p)-3*dd, gdEdx_prot->Eval(p)+10*dd);
  func->SetParLimits(2,sig-0.005,sig+0.015);

  //func->SetParLimits(3,0,200);
  //func->SetParLimits(4,-0.5,0.5);
  //func->SetParLimits(5,-0.3,0);
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

  func->SetParLimits(1,func->GetParameter(1)-1*dd, func->GetParameter(1)+1*dd);
  func->SetParLimits(2,func->GetParameter(2)-0.005,func->GetParameter(2)+0.035);

  h->Fit(func,"w0","",rmin,rmax);
  
  h->GetListOfFunctions()->Add(func);
  
  return func;
}

//___________________________________________________________
void FindFitRange(TH1F* H,Float_t& xmin,Float_t& xmax,Float_t fraction)
{
  Float_t Max = H->GetMaximum();
  Int_t MaxBin = H->GetMaximumBin();
  Int_t StartBin,StopBin;
  for(Int_t i=1;i<H->GetNbinsX()+1;i++){
    Float_t content=H->GetBinContent(i);
    //cout<<i<<" "<<content/Max<<" "<<accepar<<" "<<MaxBin<<endl;
    if(i==1 && (content/Max)>fraction)StartBin=i; // to support cases like T5H2 vs statNo.
    if((content/Max)<fraction && i<MaxBin)StartBin=i+1;
    if((content/Max)>fraction && i>MaxBin)StopBin=i;
  }
  xmin  = H->GetXaxis() -> GetBinCenter(StartBin);
  xmax  = H->GetXaxis() -> GetBinCenter(StopBin);
  cout<<"pmin="<<xmin<<",   pmax="<<xmax<<endl;
} 

//____________________________________________________________________
Double_t Univ(Double_t *x, Double_t *par)
{
  Double_t p = x[0];
  Double_t A = par[0]; 
  Double_t B = par[1];     
  Double_t C = par[2];     
  Double_t D = par[3];
  Double_t E = par[4];
    
  Double_t f = A * TMath::Power(p-B,C) *
    TMath::Exp(D*TMath::Power(p-B,E));

  return f;
  
}

//__________________________________________________________________________
Double_t FindMaximum(Double_t* array,Int_t N)
{

  Double_t max=0;
  for(Int_t i=0;i<N;i++){
    if(array[i]>max){
      max=array[i];
    }
  }
  return max;
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

//____________________________________________________________________
//
// $Log: RichPid.C,v $
// Revision 1.1  2006/12/30 11:21:42  ufstasze
// initial release
//
// Revision 1.1  2006/12/08 14:09:19  ufstasze
// Initial release
//
