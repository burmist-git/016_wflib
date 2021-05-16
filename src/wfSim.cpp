//my
#include "wfSim.hh"

//root
#include "TRandom3.h"
#include "TCanvas.h"
#include "TGraph.h"

//C, C++
#include <iostream>
#include <fstream>
#include <assert.h>
#include <iomanip>

using namespace std;

// 0.0<digitTime0<5.0
wfSim::wfSim(TRandom3 *rnd, int nDigitPoint, double dTimeDigit, double digitTime0){
  //cout<<"wfSim::wfSim()"<<endl;//
  //_maxDigitPoint = 100000;
  //if(_nDigitPoint>_maxDigitPoint){
  //cout<<"ERROR --> _nDigitPoint>_maxDigitPoint "<<endl
  //<<"          _nDigitPoint   = "<<_nDigitPoint<<endl 
  //<<"          _maxDigitPoint = "<<_maxDigitPoint<<endl;
  //assert(0);
  //}
  _nDigitPoint = nDigitPoint;
  _rnd = rnd;
  _digitTime0 = digitTime0;
  _dTimeDigit = dTimeDigit;
  _wfT = new double[_nDigitPoint];
  _wfA = new double[_nDigitPoint];
  for( Int_t i = 0; i<_nDigitPoint; i++){
    _wfT[i] = 0.0;
    _wfA[i] = 0.0;
  }
  _nPiontsMcpPM_SP = 0;
  _nPiontsMcpPM_CT = 0;
  _nPiontsMcpPM_CS = 0;
  _constDelay = 10.0;

}

wfSim::~wfSim(){
  delete _wfT;
  delete _wfA;
}

void wfSim::CleanWf(){
  for( Int_t i = 0; i<_nDigitPoint; i++){
    //_wfT[i] = 0.0;
    _wfA[i] = 0.0;
  }
}

int wfSim::getLinearPointID( double x, double xxMin, double xxMax, int nPoints){
  if(nPoints<2){
    cout<<" ERROR --> nPoints<2 "<<endl
	<<"           nPoints = "<<nPoints<<endl;
    assert(0);
  }
  if(xxMin>=xxMax){
    cout<<" ERROR --> xxMin>=xxMax "<<endl
	<<"           xxMin = "<<xxMin<<endl
	<<"           xxMax = "<<xxMax<<endl;
    assert(0);
  }
  int iID = (int)floor((x-xxMin)/((xxMax - xxMin)/(nPoints-1)));
  if(iID < 0 || iID >= nPoints)
    return -999;
  return iID;
}

void wfSim::SetIdealShape_SP(int nPiontsMcpPM_SP, double *xxMcpPM_SP, 
			     double *yyMcpPM_SP, double timeMax_SP){
  _nPiontsMcpPM_SP = nPiontsMcpPM_SP;
  _xxMcpPM_SP = xxMcpPM_SP;
  _yyMcpPM_SP = yyMcpPM_SP;
  _timeMax_SP = timeMax_SP;
  //for(int i = 0;i<nPiontsMcpPM_SP;i++){
  //cout<<setw(20)<<_xxMcpPM_SP[i]<<setw(20)<<_yyMcpPM_SP[i]<<endl;
  //}
}

void wfSim::SetIdealShape_CT(int nPiontsMcpPM_CT, double *xxMcpPM_CT,
			     double *yyMcpPM_CT, double timeMax_CT){
  _nPiontsMcpPM_CT = nPiontsMcpPM_CT;
  _xxMcpPM_CT = xxMcpPM_CT;
  _yyMcpPM_CT = yyMcpPM_CT;
  _timeMax_CT = timeMax_CT;
}

void wfSim::SetIdealShape_CS(int nPiontsMcpPM_CS, double *xxMcpPM_CS,
			     double *yyMcpPM_CS, double timeMax_CS){
  _nPiontsMcpPM_CS = nPiontsMcpPM_CS;
  _xxMcpPM_CS = xxMcpPM_CS;
  _yyMcpPM_CS = yyMcpPM_CS;
  _timeMax_CS = timeMax_CS;
}

void wfSim::genMCPPMT_SP_WF( double timeTrue){
  int i = 0;
  double x;
  double ampl = genAmplitudeDistMCPMPT_SP();
  timeTrue = timeTrue + _constDelay;
  for( i = 0; i<_nDigitPoint; i++){
    _wfT[i] = _dTimeDigit*i;
    x = _wfT[i] + _digitTime0;
    _wfA[i] = _wfA[i] + mcpPMShapeSP( x, timeTrue)*ampl;
    //_wfA[i] = _wfA[i] + 1.0;
    //if(_wfA[i]>1.0)
    //cout<<_wfA[i]<<endl;
  }
}

void wfSim::genMCPPMT_CT_WF( double timeTrue){
  int i = 0;
  double x;
  double ampl = genAmplitudeDistMCPMPT_CT();
  timeTrue = timeTrue + _constDelay;
  for( i = 0; i<_nDigitPoint; i++){
    _wfT[i] = _dTimeDigit*i;
    x = _wfT[i] + _digitTime0;
    _wfA[i] = _wfA[i] + mcpPMShapeCT( x, timeTrue)*ampl;
  }
}

void wfSim::genMCPPMT_CS_WF( double timeTrue){
  int i = 0;
  double x;
  double ampl = genAmplitudeDistMCPMPT_CS();
  timeTrue = timeTrue + _constDelay;
  //cout<<"_nDigitPoint "<<_nDigitPoint<<endl;
  for( i = 0; i<_nDigitPoint; i++){
    _wfT[i] = _dTimeDigit*i;
    x = _wfT[i] + _digitTime0;
    _wfA[i] = _wfA[i] + mcpPMShapeCS( x, timeTrue)*ampl;
    //cout<<"i "<<i<<endl;
  }
}

void wfSim::genTriangleWF( double timeTrue, double basis){
  int i = 0;
  double x;
  double ampl = genAmplitudeDistGauss(0.4, 0.05);
  timeTrue = timeTrue + 10.0;
  for( i = 0; i<_nDigitPoint; i++){
    _wfT[i] = _dTimeDigit*i;
    x = _wfT[i] + _digitTime0;
    _wfA[i] = _wfA[i] + triangleShape( x, timeTrue, basis)*ampl;
  }
}

void wfSim::genNinoWF( double timeTrue, double amplitude, double width){
  int i = 0;
  double x;
  double ampl = 1.0;
  double k = 4.0;
  double xL = timeTrue;
  double xR = timeTrue + width;
  timeTrue = timeTrue + 10.0;
  for( i = 0; i<_nDigitPoint; i++){
    _wfT[i] = _dTimeDigit*i;
    x = _wfT[i] + _digitTime0;
    _wfA[i] = _wfA[i] + amplitude/(TMath::Exp(k*(xL-x)) + 1.0)/(TMath::Exp(-k*(xR-x)) + 1.0);
    //_wfA[i] = _wfA[i] + amplitude/(TMath::Exp(k*(xL-x)) + 1.0);
    //cout<<_wfA[i]<<endl;
  }
}

void wfSim::genGaussWF( double timeTrue, double sigma){
  int i = 0;
  double x;
  //LB 09.02.2011
  //test was proposed by Achille (parametrisation of the global crosstalk)
  //double ampl = genAmplitudeDistGauss(0.2, 0.05);  
  //double ampl = _rnd->Uniform(-0.005,0.0);
  //double ampl = _rnd->Uniform(0.0,1.0);
  double ampl = 0.2;
  timeTrue = timeTrue + 10.0;
  for( i = 0; i<_nDigitPoint; i++){
    _wfT[i] = _dTimeDigit*i;
    x = _wfT[i] + _digitTime0;
    _wfA[i] = _wfA[i] + TMath::Gaus( x, timeTrue, sigma)*ampl;
    //cout<<_wfA[i]<<endl;
  }
}

void wfSim::Draw(){
  Double_t amplMax =  1.2;
  Double_t amplMin = -0.2;

  TCanvas *c1 = new TCanvas("c1","wavefor",10,10,1000,800);
  c1->SetFillColor(kWhite);

  TGraph *gr1 = new TGraph( _nDigitPoint, _wfT, _wfA); 

  gr1->SetMaximum(amplMax);
  gr1->SetMinimum(amplMin);

  gr1->SetMarkerStyle(20);
  gr1->SetMarkerSize(1);

  gr1->Draw("APL");

  //gr1->SaveAs("wfSim_wf.C");

}
//gaussian
//double wfSim::getRealAmplitude( double tt, double meanT){
//double aa = TMath::Gaus( tt, meanT, 0.5);
//if(usbConst::rmsOfTheNoise>0.0)
//aa = aa + _rnd->Gaus(0.0,usbConst::rmsOfTheNoise);
//if(aa< 1.0e-20)
//aa = 0.0;
//return aa;
//}

//triangle
double wfSim::triangleShape( double tt, double timeTrue, double basis){
  double meanT = timeTrue;

  double aa = 0.0;

  double kk = -999.0;
  double bb = -999.0;

  if( basis <= 0.0){
    cout<<"basis <= 0.0"<<endl
	<<"basis  = "<<basis<<endl;
    assert(0);
  }
  
  if( tt < (meanT - basis/2.0) || tt > (meanT + basis/2.0) ){
    kk = 0.0;
    bb = 0.0;
  }
  else if(tt >= (meanT - basis/2.0) && tt <= meanT){
    kk = 2.0/basis;
    bb = 1 - 2.0*meanT/basis;
  }
  else if(tt <= (meanT + basis/2.0) && tt >= meanT){
    kk = -2.0/basis;
    bb = 1 + 2.0*meanT/basis;  
  }
  aa = kk*tt + bb;

  return aa;
}

double wfSim::genAmplitudeDistGauss(double meanA, double rmsA){
  return TMath::Abs(_rnd->Gaus( meanA, rmsA));
}

void wfSim::generateNoiseGauss(double noiseRMS){
  for( int i = 0; i<_nDigitPoint; i++){
    _wfA[i] = _wfA[i] + _rnd->Gaus(0.0,noiseRMS);
  }
}

double wfSim::mcpPMShapeSP( double x, double timeTrue){
  if(_nPiontsMcpPM_SP < 2 ){
    cout<<" ERROR--> _nPiontsMcpPM_SP  < 2 "<<endl;
    assert(0);
  }
  x = (x-timeTrue) + _timeMax_SP;
  double xxMin = _xxMcpPM_SP[0];
  double xxMax = _xxMcpPM_SP[(_nPiontsMcpPM_SP-1)];
  if(x<=xxMin){
    return 0.0;
  }
  else if(x>=xxMax){
    return 0.0;
  }
  else{
    //for(int i = 0;i<(_nPiontsMcpPM_SP-1);i++){
    int i = getLinearPointID( x, xxMin, xxMax, _nPiontsMcpPM_SP);
    if(i!=-999){
      if(x>=_xxMcpPM_SP[i] && x<=_xxMcpPM_SP[i+1]){
	//cout<<linInterpol( x, _yyMcpPM_SP[i+1], _yyMcpPM_SP[i], _xxMcpPM_SP[i+1], _xxMcpPM_SP[i])<<endl;
	//return _yyMcpPM_SP[100];
	return linInterpol( x, _yyMcpPM_SP[i+1], _yyMcpPM_SP[i], _xxMcpPM_SP[i+1], _xxMcpPM_SP[i]);
	//}
      }
      else if((i-1)>=0 && (i+1)<_nPiontsMcpPM_SP){
	if(x>=_xxMcpPM_SP[i-1] && x<=_xxMcpPM_SP[i])
	  return linInterpol( x, _yyMcpPM_SP[i], _yyMcpPM_SP[i-1], _xxMcpPM_SP[i], _xxMcpPM_SP[i-1]);
	if(x>=_xxMcpPM_SP[i+1] && x<=_xxMcpPM_SP[i+2])
	  return linInterpol( x, _yyMcpPM_SP[i+2], _yyMcpPM_SP[i+1], _xxMcpPM_SP[i+2], _xxMcpPM_SP[i+1]);	
      }
      else if(i==0){
	if(TMath::Abs(x-_xxMcpPM_SP[1])<(_xxMcpPM_SP[1] - _xxMcpPM_SP[0])/100.0)
	  return _yyMcpPM_SP[1];
      }
      else{
	cout<<"  i                = "<<i<<endl
	    <<"  x                = "<<x<<endl
	    <<"  _xxMcpPM_SP[i]   = "<<_xxMcpPM_SP[i]<<endl
	    <<"  _xxMcpPM_SP[i+1] = "<<_xxMcpPM_SP[i+1]<<endl;
	assert(0);
      }
    }
  }
  return 0.0;
}

double wfSim::mcpPMShapeCT( double x, double timeTrue){
  if(_nPiontsMcpPM_CT == 0 ){
    cout<<" ERROR--> _nPiontsMcpPM_CT == 0 "<<endl;
    assert(0);
  }
  double xxMin = _xxMcpPM_CT[0];
  double xxMax = _xxMcpPM_CT[(_nPiontsMcpPM_CT-1)];
  x = (x-timeTrue) + _timeMax_CT;
  if(x<xxMin){
    return 0.0;
  }
  else if(x>xxMax){
    return 0.0;
  }
  else{
    //for(int i = 0;i<(_nPiontsMcpPM_SP-1);i++){
    int i = getLinearPointID( x, xxMin, xxMax, _nPiontsMcpPM_CT);
    if(i!=-999){
      if(x>=_xxMcpPM_CT[i] && x<=_xxMcpPM_CT[i+1]){
	return linInterpol( x, _yyMcpPM_CT[i+1], _yyMcpPM_CT[i], _xxMcpPM_CT[i+1], _xxMcpPM_CT[i]);
      }
      else if((i-1)>=0 && (i+1)<_nPiontsMcpPM_CT){
	if(x>=_xxMcpPM_CT[i-1] && x<=_xxMcpPM_CT[i])
	  return linInterpol( x, _yyMcpPM_CT[i], _yyMcpPM_CT[i-1], _xxMcpPM_CT[i], _xxMcpPM_CT[i-1]);
	if(x>=_xxMcpPM_CT[i+1] && x<=_xxMcpPM_CT[i+2])
	  return linInterpol( x, _yyMcpPM_CT[i+2], _yyMcpPM_CT[i+1], _xxMcpPM_CT[i+2], _xxMcpPM_CT[i+1]);	
      }
      else if(i==0){
	if(TMath::Abs(x-_xxMcpPM_CT[1])<(_xxMcpPM_CT[1] - _xxMcpPM_CT[0])/100.0)
	  return _yyMcpPM_CT[1];
      }
      else{
	cout<<"  i                = "<<i<<endl
	    <<"  x                = "<<x<<endl
	    <<"  _xxMcpPM_CT[i]   = "<<_xxMcpPM_CT[i]<<endl
	    <<"  _xxMcpPM_CT[i+1] = "<<_xxMcpPM_CT[i+1]<<endl;
	assert(0);
      }
    }
  }
  return 0.0;
}

double wfSim::mcpPMShapeCS( double x, double timeTrue){
  if(_nPiontsMcpPM_CS == 0 ){
    cout<<" ERROR--> _nPiontsMcpPM_CS == 0 "<<endl;
    assert(0);
  }
  double xxMin = _xxMcpPM_CS[0];
  double xxMax = _xxMcpPM_CS[(_nPiontsMcpPM_CS-1)];
  x = (x-timeTrue) + _timeMax_CS;
  if(x<xxMin){
    return 0.0;
  }
  else if(x>xxMax){
    return 0.0;
  }
  else{
    int i = getLinearPointID( x, xxMin, xxMax, _nPiontsMcpPM_CS);
    if(i!=-999){
      if(x>=_xxMcpPM_CS[i] && x<=_xxMcpPM_CS[i+1]){
	return linInterpol( x, _yyMcpPM_CS[i+1], _yyMcpPM_CS[i], _xxMcpPM_CS[i+1], _xxMcpPM_CS[i]);
      }
      else if((i-1)>=0 && (i+1)<_nPiontsMcpPM_CS){
	if(x>=_xxMcpPM_CS[i-1] && x<=_xxMcpPM_CS[i])
	  return linInterpol( x, _yyMcpPM_CS[i], _yyMcpPM_CS[i-1], _xxMcpPM_CS[i], _xxMcpPM_CS[i-1]);
	if(x>=_xxMcpPM_CS[i+1] && x<=_xxMcpPM_CS[i+2])
	  return linInterpol( x, _yyMcpPM_CS[i+2], _yyMcpPM_CS[i+1], _xxMcpPM_CS[i+2], _xxMcpPM_CS[i+1]);	
      }
      else if(i==0){
	if(TMath::Abs(x-_xxMcpPM_CS[1])<(_xxMcpPM_CS[1] - _xxMcpPM_CS[0])/100.0)
	  return _yyMcpPM_CS[1];
      }
      else{
	cout<<"  i                = "<<i<<endl
	    <<"  x                = "<<x<<endl
	    <<"  _xxMcpPM_CS[i]   = "<<_xxMcpPM_CS[i]<<endl
	    <<"  _xxMcpPM_CS[i+1] = "<<_xxMcpPM_CS[i+1]<<endl;
	assert(0);
      }
    }
  }
  return 0.0;
}

double wfSim::linInterpol( double x, double y2, double y1, double x2, double x1){
  double dx = (x1-x2);
  if(dx==0.0 || x2<x1){
    cout<<endl<<" ERROR--> x2<x1 || (x1 - x2)==0.0"<<endl
	<<" x1 = "<<x1<<endl
	<<" x2 = "<<x2<<endl
	<<" y1 = "<<y1<<endl
	<<" y2 = "<<y2<<endl;
    assert(0);
  }
  double k = (y1 - y2)/(x1-x2);
  double b = y1 - k*x1;
  if((y1 - y2)==0.0){
    //cout<<"  Warning --> (y1 - y2)==0.0 in wfSim::linInterpol"<<endl;
    return y1;
  }
  return (k*x + b);
}

double wfSim::genAmplitudeDistMCPMPT_SP(){
  const int _nP_SP = 22;
  double _ampl_SP[_nP_SP];
  double _prob_SP[_nP_SP];
  _ampl_SP[0]=0.04225;  _prob_SP[0]=0;
  _ampl_SP[1]=0.06775;  _prob_SP[1]=171;
  _ampl_SP[2]=0.09325;  _prob_SP[2]=731;
  _ampl_SP[3]=0.11875;  _prob_SP[3]=803;
  _ampl_SP[4]=0.14425;  _prob_SP[4]=884;//
  _ampl_SP[5]=0.16975;  _prob_SP[5]=814;
  _ampl_SP[6]=0.19525;  _prob_SP[6]=750;
  _ampl_SP[7]=0.22075;  _prob_SP[7]=641;
  _ampl_SP[8]=0.24625;  _prob_SP[8]=549;
  _ampl_SP[9]=0.27175;  _prob_SP[9]=412;
  _ampl_SP[10]=0.29725; _prob_SP[10]=321;
  _ampl_SP[11]=0.32275; _prob_SP[11]=215;
  _ampl_SP[12]=0.34825; _prob_SP[12]=139;
  _ampl_SP[13]=0.37375; _prob_SP[13]=92;
  _ampl_SP[14]=0.39925; _prob_SP[14]=54;
  _ampl_SP[15]=0.42475; _prob_SP[15]=34;
  _ampl_SP[16]=0.45025; _prob_SP[16]=19;
  _ampl_SP[17]=0.47575; _prob_SP[17]=8;
  _ampl_SP[18]=0.50125; _prob_SP[18]=6;
  _ampl_SP[19]=0.52675; _prob_SP[19]=4;
  _ampl_SP[20]=0.55225; _prob_SP[20]=2;
  _ampl_SP[21]=0.57775; _prob_SP[21]=0;
  double xx;
  double yy;
  double accP;
  int i = 0;
  bool Ok = false;
  if(_nP_SP<2){
    cout<<" ERROR --> _nP_SP<2"<<endl
	<<"           _nP_SP = "<<_nP_SP<<endl;
    assert(0);
  }
  while(!Ok){
    xx = _rnd->Uniform( _ampl_SP[0], _ampl_SP[_nP_SP-1]);
    yy = _rnd->Uniform(     0.0,    _prob_SP[4]);
    for(i=0;i<(_nP_SP-1);i++){
      if(_ampl_SP[i]<=xx && _ampl_SP[i+1]>=xx ){
	accP = linInterpol(xx,_prob_SP[i+1],_prob_SP[i],_ampl_SP[i+1],_ampl_SP[i]);
	i = _nP_SP-1;
      }
    }
    if(yy<accP){
      Ok = true;
      return xx;
    }
  }
  return 0.0;
}

double wfSim::genAmplitudeDistMCPMPT_CT(){
  const int nP = 21;
  double ampl[nP];
  double prob[nP];
  ampl[0]=0.00825; prob[0]=0;
  ampl[1]=0.00995; prob[1]=389;
  ampl[2]=0.01165; prob[2]=786;//
  ampl[3]=0.01335; prob[3]=559;
  ampl[4]=0.01505; prob[4]=431;
  ampl[5]=0.01675; prob[5]=286;
  ampl[6]=0.01845; prob[6]=177;
  ampl[7]=0.02015; prob[7]=127;
  ampl[8]=0.02185; prob[8]=93;
  ampl[9]=0.02355; prob[9]=65;
  ampl[10]=0.02525; prob[10]=57;
  ampl[11]=0.02695; prob[11]=56;
  ampl[12]=0.02865; prob[12]=29;
  ampl[13]=0.03035; prob[13]=31;
  ampl[14]=0.03205; prob[14]=21;
  ampl[15]=0.03375; prob[15]=15;
  ampl[16]=0.03545; prob[16]=19;
  ampl[17]=0.03715; prob[17]=13;
  ampl[18]=0.03885; prob[18]=11;
  ampl[19]=0.04055; prob[19]=1;
  ampl[20]=0.04225; prob[20]=0;
  double xx;
  double yy;
  double accP;
  int i = 0;
  bool Ok = false;
  if(nP<2){
    cout<<" ERROR --> nP<2"<<endl
	<<"           nP = "<<nP<<endl;
    assert(0);
  }
  while(!Ok){
    xx = _rnd->Uniform( ampl[0], ampl[nP-1]);
    yy = _rnd->Uniform(     0.0,    prob[2]);
    for(i=0;i<(nP-1);i++){
      if(ampl[i]<=xx && ampl[i+1]>=xx ){
	accP = linInterpol(xx,prob[i+1],prob[i],ampl[i+1],ampl[i]);
	i = nP-1;
      }
    }
    if(yy<accP){
      Ok = true;
      return xx;
    }
  }
  return 0.0;
}

double wfSim::genAmplitudeDistMCPMPT_CS(){
  const int nP = 33;
  double ampl[nP];
  double prob[nP];
  ampl[0]=0.00825; prob[0]=0;
  ampl[1]=0.00995; prob[1]=149;
  ampl[2]=0.01165; prob[2]=301;
  ampl[3]=0.01335; prob[3]=270;
  ampl[4]=0.01505; prob[4]=229;
  ampl[5]=0.01675; prob[5]=197;
  ampl[6]=0.01845; prob[6]=150;
  ampl[7]=0.02015; prob[7]=138;
  ampl[8]=0.02185; prob[8]=127;
  ampl[9]=0.02355; prob[9]=108;
  ampl[10]=0.02525; prob[10]=97;
  ampl[11]=0.02695; prob[11]=83;
  ampl[12]=0.02865; prob[12]=69;
  ampl[13]=0.03035; prob[13]=56;
  ampl[14]=0.03205; prob[14]=52;
  ampl[15]=0.03375; prob[15]=39;
  ampl[16]=0.03545; prob[16]=34;
  ampl[17]=0.03715; prob[17]=31;
  ampl[18]=0.03885; prob[18]=29;
  ampl[19]=0.04055; prob[19]=27;
  ampl[20]=0.04225; prob[20]=25;
  ampl[21]=0.04395; prob[21]=23;
  ampl[22]=0.04565; prob[22]=16;
  ampl[23]=0.04735; prob[23]=17;
  ampl[24]=0.04905; prob[24]=18;
  ampl[25]=0.05075; prob[25]=24;
  ampl[26]=0.05245; prob[26]=14;
  ampl[27]=0.05415; prob[27]=20;
  ampl[28]=0.05585; prob[28]=14;
  ampl[29]=0.05755; prob[29]=15;
  ampl[30]=0.05925; prob[30]=11;
  ampl[31]=0.06095; prob[31]=5;
  ampl[32]=0.06265; prob[32]=0;
  double xx;
  double yy;
  double accP;
  int i = 0;
  bool Ok = false;
  if(nP<2){
    cout<<" ERROR --> nP<2"<<endl
	<<"           nP = "<<nP<<endl;
    assert(0);
  }
  while(!Ok){
    xx = _rnd->Uniform( ampl[0], ampl[nP-1]);
    yy = _rnd->Uniform(     0.0,    prob[2]);
    for(i=0;i<(nP-1);i++){
      if(ampl[i]<=xx && ampl[i+1]>=xx ){
	accP = linInterpol(xx,prob[i+1],prob[i],ampl[i+1],ampl[i]);
	i = nP-1;
      }
    }
    if(yy<accP){
      Ok = true;
      return xx;
    }
  }
  return 0.0;
}
