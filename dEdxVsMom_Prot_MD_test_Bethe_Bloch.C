Double_t BetheBloch_p(Double_t *x, Double_t *par);
Double_t BetheBloch_pi(Double_t *x, Double_t *par);
void Plot_BetheBloch(Double_t X0, Double_t X1, Double_t miu, Double_t s);
void TPstyle();

//__________________________________________________________________________________
void Test_BetheBloch() 
{

  Double_t log10p=1.2;

  gStyle->SetOptStat(0);
  gStyle->SetOptFit(1);

    // W nowej parameteryzacji mamy 4 parametry: s, miu, X0 i X1
    // test_B = BetheBloch(Double_t *x, Double_t *par);
  Double_t par[] = {1.44, 4, 3.77, 1.53, 0.9382720799426998};

  TF1* Bethe_Bloch_func = new TF1("Bethe Bloch", BetheBloch_p, -0.5, 3, 5); // <- ostatnie trzy cyfry to liczba parametrow
    // set the parameters to the mean and RMS of the histogram
  Bethe_Bloch_func->SetParameters(par);
  Bethe_Bloch_func->SetParNames("X0","X1","miu","s","mass");
  Bethe_Bloch_func->FixParameter(4, 0.9382720799426998);


  Double_t log10p_min = TMath::Log10(3.77);
  cout<<" log10p in min: "<<log10p_min<<endl;

  Double_t dEdx_min = Bethe_Bloch_func->Eval(log10p_min);
  Double_t dEdx = Bethe_Bloch_func->Eval(log10p);

  std::cout << std::endl;  


}

//__________________________________________________________________________________
void Plot_BetheBloch() 
{

  gStyle->SetOptStat(0);
  gStyle->SetOptFit(1);

    // W nowej parameteryzacji mamy 4 parametry: s, miu, X0 i X1
    // test_B = BetheBloch(Double_t *x, Double_t *par);
  Double_t par_p[] = {1.44, 4, 3.77, 1.53, 0.9382720799426998};

  TF1* Bethe_Bloch_func_p = new TF1("Bethe Bloch p", BetheBloch_p, -0.5, 3, 5); // <- ostatnie trzy cyfry to liczba parametrow
    // set the parameters to the mean and RMS of the histogram
  Bethe_Bloch_func_p->SetParameters(par_p);
  Bethe_Bloch_func_p->SetParNames("X0","X1","miu","s","mass");
  Bethe_Bloch_func_p->FixParameter(4, 0.9382720799426998);

  Double_t par_pi[] = {1.44, 4, 3.77, 1.53, 1.1349766};

  TF1* Bethe_Bloch_func_pi = new TF1("Bethe Bloch pi", BetheBloch_pi, -0.5, 3, 5); // <- ostatnie trzy cyfry to liczba parametrow
    // set the parameters to the mean and RMS of the histogram
  Bethe_Bloch_func_pi->SetParameters(par_pi);
  Bethe_Bloch_func_pi->SetParNames("X0","X1","miu","s","mass");
  Bethe_Bloch_func_pi->FixParameter(4, 1.1349766);


  TH2F* frame = new TH2F("Frame", "", 1000, -0.5, 2.3, 1000, 0.8, 2);

  TCanvas* c1 = new TCanvas ("dEdx test","",200,200,700,560);

  c1->Divide(1,1);
  TPad* pad = (TPad*)c1->cd(1);
  pad->SetGridx(1);
  pad->SetGridy(1);

  Bethe_Bloch_func_p->SetLineWidth(2);
  Bethe_Bloch_func_p->SetLineColor(4);

  Bethe_Bloch_func_pi->SetLineWidth(2);
  Bethe_Bloch_func_pi->SetLineColor(2);
  
  frame->Draw();
  Bethe_Bloch_func_p->Draw("same");
  Bethe_Bloch_func_pi->Draw("same");
  
  cout<<" momentum in minimum = "<<TMath::Power(10,0.585)<<endl;

}



//__________________________________________________________________________
Double_t BetheBloch_p(Double_t *x, Double_t *par) 
{

  // x + 4 paramtery Bethego-Blocha:
  //Double_t p = x[0];
  Double_t log10p = x[0];
  Double_t p = TMath::Power(10,log10p);

  Double_t X0 = par[0]; // interval limits
  Double_t X1 = par[1];
  Double_t miu = par[2]; // betha*gamma value where Bethe-Bloch function gets minimum
  Double_t s = par[3]; // tam, gdzie wysyca sie Bethe Bloch dla duzych pedow
  Double_t mass = par[4];

  std::cout << "------------------------" << std::endl;
  std::cout << "PARAMETRY: " << std::endl;
  std::cout << "p =: " << p << std::endl;
  std::cout << "X0 =: " << X0 << std::endl;
  std::cout << "X1 =: " << X1 << std::endl;
  std::cout << "miu =: " << miu << std::endl;
  std::cout << "s =: " << s << std::endl;
  
  // do funkcji wartości:
  // Double_t mass = 0.9382720799426998;  // proton rest mass 0.9382720799426998
  Double_t pp = p*p;
  //std::cout << "PP: " << pp << std::endl;
  Double_t mm = mass*mass;
  //std::cout << "MM: " << mm << std::endl;

  Double_t ppmm = pp+mm;
  //std::cout << "ppmm: " << ppmm << std::endl;

  Double_t sqppmm = TMath::Sqrt(ppmm);
  //std::cout << "sqppmm: " << sqppmm << std::endl;

  Double_t ok = p/sqppmm;
  // std::cout << "OK: " << ok << std::endl;

// BETHA WYCHOODZI 1 I PRZEZ TO GAMMA MA DZIELENIE PRZEZ 0 !!!!!!
  Double_t betha = p/TMath::Sqrt(p*p + mass*mass); // betha = p/sqrt(p^2 + m^2)
  Double_t gamma = 1/TMath::Sqrt(1 - betha*betha); // gamma = 1/sqrt(1 - beta^2)
  Double_t alpha = -2;

  std::cout << "\n ____ betha = " << betha << std::endl;
  std:: cout << "\n ____ gamma = " << gamma << std::endl;
  std:: cout << "\n ____ mass = " << mass << std::endl;
  Double_t X = TMath::Log10(betha * gamma);   // betha*gamma=p/mc
  Double_t B = miu/TMath::Sqrt(1 + miu*miu);                                    // PS: tu był błąd
  //cout<<" B="<<B<<" "<<miu<<" "<<TMath::Sqrt(1 + miu*miu)<<endl;  

  Double_t E0 = 1./(miu*miu);

  // pojawia sie logarytm z ujemnej wartosci dodalam Abs
  Double_t b = miu*miu + TMath::Log(TMath::Abs(1 - B*B));                      // PS: tu był bląd

  //Double_t M = (X1 - X0) * TMath::Power( (s/E0 - b + 1) / (2 * TMath::Log(10)), -1);
  //cout<<" M="<<M<<endl;
  Double_t M = (X1 - X0) / ((s/E0 - b + 1)/(2 * TMath::Log(10)) - X0);          // PS: tu był bląd 
  //cout<<" M="<<M<<endl;

  Double_t XA = X0 - (X0 - X1)/M;
  //Double_t a = 2 * TMath::Log(10.) * TMath::Power(miu * TMath::Power((X1 - X0), (1 - M)), -1);
  //cout<<" a="<<a<<endl;
  //Double_t a = 2. * TMath::Log(10.) * TMath::Power(X1-X0,1.-M) / miu;
  Double_t a = 2. * TMath::Log(10.) * TMath::Power(X1-X0,1.-M) / M;     // w pracy Gabora jest błąd w definicji "a": 
                                                                        // Aby dostać dobre dE/dx w formule 
                                                                        // na "a" należy miu zamienić na M
  cout<<" a="<<a<<endl;


  Double_t delta = 0;

  /*
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
  */
  std::cout << "/n WCHODZE DO IFA" << std::endl;
  std::cout << "X: " << X << " X0: " << X0 << " X1: " << X1 << std::endl;
  
 

  if (X < X0) {
    delta = 0;
    std::cout << " jestem w pierwszym " << std::endl;
    std::cout << " delta = " << delta << std::endl;
  }
  else if (X >= X0 || X < X1){
    // musialam dac Abs na pierwszy argument Power 
    delta = 2 * TMath::Log(10.) * (X - XA) + a * TMath::Power(X1 - X, M);
    std::cout << " jestem w drugim " << std::endl;
    std::cout << " delta = " << delta << std::endl;
  }
  else if ( X1 <= X ){
    delta = 2 * TMath::Log(10) * (X - XA);
    //delta = 2 * TMath::Log(10.) * (X - XA) + a * TMath::Power(X1 - X, M);
    std::cout << " jestem w trzecim " << std::endl;
    std::cout << " delta = " << delta << std::endl;
  }

    //std::cout << "delta =: " << delta << std::endl;

  Double_t dEdx = E0 * TMath::Power(betha, alpha) * (b + 2*TMath::Log(gamma) - betha*betha - delta);
  
  //std::cout << "dEdx = " << dEdx << std::endl;

  return dEdx;
}

//__________________________________________________________________________
Double_t BetheBloch_pi(Double_t *x, Double_t *par) 
{

  // x + 4 paramtery Bethego-Blocha:
  //Double_t p = x[0];
  Double_t log10p = x[0];
  Double_t p = TMath::Power(10,log10p);

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
  Double_t mass = 0.13957;  // pion mass
  Double_t pp = p*p;
  //std::cout << "PP: " << pp << std::endl;
  Double_t mm = mass*mass;
  //std::cout << "MM: " << mm << std::endl;

  Double_t ppmm = pp+mm;
  //std::cout << "ppmm: " << ppmm << std::endl;

  Double_t sqppmm = TMath::Sqrt(ppmm);
  //std::cout << "sqppmm: " << sqppmm << std::endl;

  Double_t ok = p/sqppmm;
  // std::cout << "OK: " << ok << std::endl;


  Double_t betha = p/TMath::Sqrt(p*p + mass*mass); // betha = p/sqrt(p^2 + m^2)
  Double_t gamma = 1/TMath::Sqrt(1 - betha*betha); // gamma = 1/sqrt(1 - beta^2)
  Double_t alpha = -2;
  Double_t X = TMath::Log10(betha * gamma);   // betha*gamma=p/mc
  Double_t B = miu/TMath::Sqrt(1 + miu*miu);                                    // PS: tu był błąd
  //cout<<" B="<<B<<" "<<miu<<" "<<TMath::Sqrt(1 + miu*miu)<<endl;  

  Double_t E0 = 1./(miu*miu);

  // pojawia sie logarytm z ujemnej wartosci dodalam Abs
  Double_t b = miu*miu + TMath::Log(TMath::Abs(1 - B*B));                      // PS: tu był bląd

  //Double_t M = (X1 - X0) * TMath::Power( (s/E0 - b + 1) / (2 * TMath::Log(10)), -1);
  //cout<<" M="<<M<<endl;
  Double_t M = (X1 - X0) / ((s/E0 - b + 1)/(2 * TMath::Log(10)) - X0);          // PS: tu był bląd 
  //cout<<" M="<<M<<endl;

  Double_t XA = X0 - (X0 - X1)/M;
  //Double_t a = 2 * TMath::Log(10.) * TMath::Power(miu * TMath::Power((X1 - X0), (1 - M)), -1);
  //cout<<" a="<<a<<endl;
  Double_t a = 2. * TMath::Log(10.) * TMath::Power(X1-X0,1.-M) / M;
  cout<<" a="<<a<<endl;


  Double_t delta = 0;

  /*
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
  */
  std::cout << "/n WCHODZE DO IFA" << std::endl;
  std::cout << "X: " << X << " X0: " << X0 << " X1: " << X1 << std::endl;
 

  if (X < X0) {
    delta = 0;
    std::cout << " jestem w pierwszym " << std::endl;
    std::cout << " delta = " << delta << std::endl;
  }
  else if (X >= X0 || X < X1){
    // musialam dac Abs na pierwszy argument Power 
    delta = 2 * TMath::Log(10.) * (X - XA) + a * TMath::Power(X1 - X, M);
    std::cout << " jestem w drugim " << std::endl;
    std::cout << " delta = " << delta << std::endl;
  }
  else if ( X1 <= X ){
    delta = 2 * TMath::Log(10) * (X - XA);
    //delta = 2 * TMath::Log(10.) * (X - XA) + a * TMath::Power(X1 - X, M);
    std::cout << " jestem w trzecim " << std::endl;
    std::cout << " delta = " << delta << std::endl;
  }

    //std::cout << "delta =: " << delta << std::endl;

  Double_t dEdx = E0 * TMath::Power(betha, alpha) * (b + 2*TMath::Log(gamma) - betha*betha - delta);
  
  //std::cout << "dEdx = " << dEdx << std::endl;

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
