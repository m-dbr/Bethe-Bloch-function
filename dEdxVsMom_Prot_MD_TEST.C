Double_t BetheBloch_p(Double_t *x, Double_t *par);
Double_t BetheBloch_pi(Double_t *x, Double_t *par);
void TPstyle();

//__________________________________________________________________________________
void Test_BetheBloch(Double_t X0=1.44, Double_t X1=4, Double_t miu=3.77, Double_t s=1.53) 
{

  gStyle->SetOptStat(0);
  gStyle->SetOptFit(1);

    // W nowej parameteryzacji mamy 4 parametry: s, miu, X0 i X1
    // test_B = BetheBloch(Double_t *x, Double_t *par);
  Double_t par[] = {X0, X1, miu, s};

  // PROTONY
  TF1* Bethe_Bloch_p_func = new TF1("Bethe Bloch_p", BetheBloch_p, -0.5, 2, 4); 
    // set the parameters to the mean and RMS of the histogram
  Bethe_Bloch_p_func->SetParameters(par);
  Bethe_Bloch_p_func->SetParNames("X0", "X1", "miu", "s");

  // Eval wartosc dla minimum - uzyc zamiast rysowania

  // PIONY
  TF1* Bethe_Bloch_pi_func = new TF1("Bethe Bloch_pi", BetheBloch_pi, -0.5, 2, 4);
  // set the parameters to the mean and RMS of the histogram
  Bethe_Bloch_pi_func->SetParameters(par);
  Bethe_Bloch_pi_func->SetParNames("X0", "X1", "miu", "s");

  TH2F* frame = new TH2F("Frame", "", 1000, -0.5, 2, 1000, 0.8, 2);

  TCanvas* c1 = new TCanvas ("dEdx test","",200, 200, 650, 560);


  c1->Divide(1,1);
  c1->cd(1);

  Bethe_Bloch_p_func->SetLineWidth(2);
  Bethe_Bloch_p_func->SetLineColor(4);

  Bethe_Bloch_pi_func->SetLineWidth(2);
  Bethe_Bloch_pi_func->SetLineColor(2);
  
  frame->Draw();
  Bethe_Bloch_p_func->Draw("same");
  Bethe_Bloch_pi_func->Draw("same");

  
}



//__________________________________________________________________________
Double_t BetheBloch_p(Double_t *x, Double_t *par) 
{

  // x + 4 paramtery Bethego-Blocha:
  // Double_t p = x[0];

  Double_t log10p = x[0];
  Double_t p = TMath::Power(10, log10p);

  Double_t X0 = par[0]; // interval limits
  Double_t X1 = par[1];
  Double_t miu = par[2]; // betha*gamma value where Bethe-Bloch function gets minimum
  Double_t s = par[3]; // tam, gdzie wysyca sie Bethe Bloch dla duzych pedow
  
  std::cout << "------------------------" << std::endl;
  std::cout << "PARAMETRY: " << std::endl;
  std::cout << "p =: " << p << std::endl;
  std::cout << "X0 =: " << X0 << std::endl;
  std::cout << "X1 =: " << X1 << std::endl;
  std::cout << "miu =: " << miu << std::endl;
  std::cout << "s =: " << s << std::endl;

  // do funkcji wartości:
  Double_t mass = 0.9382720799426998;  // proton rest mass 0.9382720799426998
  Double_t pp = p*p;
  std::cout << "PP: " << pp << std::endl;
  Double_t mm = mass*mass;
    std::cout << "MM: " << mm << std::endl;

  Double_t ppmm = pp+mm;
      std::cout << "ppmm: " << ppmm << std::endl;

  Double_t sqppmm = TMath::Sqrt(ppmm);
    std::cout << "sqppmm: " << sqppmm << std::endl;

  Double_t ok = p/sqppmm;
      std::cout << "OK: " << ok << std::endl;


  Double_t betha = p/TMath::Sqrt(p*p + mass*mass); // betha = p/sqrt(p^2 + m^2)
  Double_t gamma = 1/TMath::Sqrt(1 - betha*betha); // gamma = 1/sqrt(1 - beta^2)
  Double_t alpha = -2;
  Double_t X = TMath::Log10(betha * gamma);   // betha*gamma=p/mc
  Double_t B = miu/(TMath::Sqrt(1 + miu));
  Double_t E0 = TMath::Power(TMath::Power(miu, 2), -1);

  // pojawia sie logarytm z ujemnej wartosci dodalam Abs
  Double_t b = miu*miu + TMath::Log(TMath::Abs(1 - B));

  Double_t M = (X1 - X0) * TMath::Power( (s/E0 - b + 1) / (2 * TMath::Log(10)), -1);
  Double_t XA = X0 - (X0 - X1)/M;
  Double_t a = 2 * TMath::Log(10) * TMath::Power(miu * TMath::Power((X1 - X0), (1 - M)), -1);
  Double_t delta = 0;

  std::cout << "ZALEZNOSCI: " << std::endl;
  std::cout << "betha =: " << betha << std::endl;
  std::cout << "gamma =: " << gamma << std::endl;
  std::cout << "X =: " << X << std::endl;
  std::cout << "B =: " << B << std::endl;
  std::cout << "E0 =: " << E0 << std::endl;
  std::cout << "b =: " << b << std::endl;
  std::cout << "M =: " << M << std::endl;
  std::cout << "XA =: " << XA << std::endl;
  std::cout << "a =: " << a << std::endl;
  std::cout << "poczatkowa delta =: " << delta << std::endl;

  std::cout << "/n WCHODZE DO IFA" << std::endl;
  std::cout << "X: " << X << " X0: " << X0 << " X1: " << X1 << std::endl;


  if (X < X0) {
    delta = 0;
    std::cout << " jestem w pierwszym " << std::endl;
    std::cout << " delta = " << delta << std::endl;}
  else if (X > X0 || X < X1){
    // musialam dac Abs na pierwszy argument Power 
    delta = 2 * TMath::Log(10) * (X - XA) + a * TMath::Power((X1 - X), M);
    std::cout << " jestem w drugim " << std::endl;
    std::cout << " delta = " << delta << std::endl;}
  else if ( X1 < X ){
    delta = 2 * TMath::Log(10) * (X - XA);
    std::cout << " jestem w trzecim " << std::endl;
    std::cout << " delta = " << delta << std::endl;}

  std::cout << "delta =: " << delta << std::endl;

  Double_t dEdx = E0 * TMath::Power(betha, alpha) * (b + 2 * TMath::Log(gamma) - TMath::Power(betha, 2) - delta);
  
  std::cout << "dEdx = " << dEdx << std::endl;

  return dEdx;
}

//__________________________________________________________________________
Double_t BetheBloch_pi(Double_t *x, Double_t *par) 
{

  // x + 4 paramtery Bethego-Blocha:
  // Double_t p = x[0];

  Double_t log10p = x[0];
  Double_t p = TMath::Power(10, log10p);

  Double_t X0 = par[0]; // interval limits
  Double_t X1 = par[1];
  Double_t miu = par[2]; // betha*gamma value where Bethe-Bloch function gets minimum
  Double_t s = par[3]; // tam, gdzie wysyca sie Bethe Bloch dla duzych pedow
  
  std::cout << "------------------------" << std::endl;
  std::cout << "PARAMETRY: " << std::endl;
  std::cout << "p =: " << p << std::endl;
  std::cout << "X0 =: " << X0 << std::endl;
  std::cout << "X1 =: " << X1 << std::endl;
  std::cout << "miu =: " << miu << std::endl;
  std::cout << "s =: " << s << std::endl;

  // do funkcji wartości:
  Double_t mass = 0.139;  // pion rest mass ????
  Double_t pp = p*p;
  std::cout << "PP: " << pp << std::endl;
  Double_t mm = mass*mass;
    std::cout << "MM: " << mm << std::endl;

  Double_t ppmm = pp+mm;
      std::cout << "ppmm: " << ppmm << std::endl;

  Double_t sqppmm = TMath::Sqrt(ppmm);
    std::cout << "sqppmm: " << sqppmm << std::endl;

  Double_t ok = p/sqppmm;
      std::cout << "OK: " << ok << std::endl;


  Double_t betha = p/TMath::Sqrt(p*p + mass*mass); // betha = p/sqrt(p^2 + m^2)
  Double_t gamma = 1/TMath::Sqrt(1 - betha*betha); // gamma = 1/sqrt(1 - beta^2)
  Double_t alpha = -2;
  Double_t X = TMath::Log10(betha * gamma);   // betha*gamma=p/mc
  Double_t B = miu/(TMath::Sqrt(1 + miu));
  Double_t E0 = TMath::Power(TMath::Power(miu, 2), -1);

  // pojawia sie logarytm z ujemnej wartosci dodalam Abs
  Double_t b = miu*miu + TMath::Log(TMath::Abs(1 - B));

  Double_t M = (X1 - X0) * TMath::Power( (s/E0 - b + 1) / (2 * TMath::Log(10)), -1);
  Double_t XA = X0 - (X0 - X1)/M;
  Double_t a = 2 * TMath::Log(10) * TMath::Power(miu * TMath::Power((X1 - X0), (1 - M)), -1);
  Double_t delta = 0;

  std::cout << "ZALEZNOSCI: " << std::endl;
  std::cout << "betha =: " << betha << std::endl;
  std::cout << "gamma =: " << gamma << std::endl;
  std::cout << "X =: " << X << std::endl;
  std::cout << "B =: " << B << std::endl;
  std::cout << "E0 =: " << E0 << std::endl;
  std::cout << "b =: " << b << std::endl;
  std::cout << "M =: " << M << std::endl;
  std::cout << "XA =: " << XA << std::endl;
  std::cout << "a =: " << a << std::endl;
  std::cout << "poczatkowa delta =: " << delta << std::endl;

  std::cout << "/n WCHODZE DO IFA" << std::endl;
  std::cout << "X: " << X << " X0: " << X0 << " X1: " << X1 << std::endl;


  if (X < X0) {
    delta = 0;
    std::cout << " jestem w pierwszym " << std::endl;
    std::cout << " delta = " << delta << std::endl;}
  else if (X > X0 || X < X1){
    // musialam dac Abs na pierwszy argument Power 
    delta = 2 * TMath::Log(10) * (X - XA) + a * TMath::Power(TMath::Abs((X1 - X)), M);
    std::cout << " jestem w drugim " << std::endl;
    std::cout << " delta = " << delta << std::endl;}
  else if ( X1 < X ){
    delta = 2 * TMath::Log(10) * (X - XA);
    std::cout << " jestem w trzecim " << std::endl;
    std::cout << " delta = " << delta << std::endl;}

  std::cout << "delta =: " << delta << std::endl;

  Double_t dEdx = E0 * TMath::Power(betha, alpha) * (b + 2 * TMath::Log(gamma) - TMath::Power(betha, 2) - delta);
  
  std::cout << "dEdx = " << dEdx << std::endl;

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