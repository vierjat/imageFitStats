#include <iostream>
#include <sstream>
#include <string.h>
#include <stdio.h>
#include <sys/time.h>
#include <time.h>
#include <inttypes.h>
#include <fstream>
#include <unistd.h>
#include <getopt.h>    /* for getopt_long; standard getopt is in unistd.h */
#include <vector>
#include <algorithm>
#include <ctime>
#include <climits>
#include <cmath>
#include <iomanip>

#include "TF1.h"
#include "TFile.h"
#include "TCanvas.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TGraph.h"
#include "TNtuple.h"
#include "TStyle.h"
#include "TText.h"
#include "TApplication.h"
#include "TFitResultPtr.h"
#include "TFitResult.h"

#include "fitsio.h"

#include "globalConstants.h"

using namespace std;

bool fileExist(const char *fileName){
  ifstream in(fileName,ios::in);
  
  if(in.fail()){
    //cout <<"\nError reading file: " << fileName <<"\nThe file doesn't exist!\n\n";
    in.close();
    return false;
  }
  
  in.close();
  return true;
}

/*========================================================
  ASCII progress bar
==========================================================*/
void showProgress(unsigned int currEvent, unsigned int nEvent) {

  const int nProgWidth=50;

  if ( currEvent != 0 ) {
    for ( int i=0;i<nProgWidth+8;i++)
      cout << "\b";
  }

  double percent = (double) currEvent/ (double) nEvent;
  int nBars = (int) ( percent*nProgWidth );

  cout << " |";
  for ( int i=0;i<nBars-1;i++)
    cout << "=";
  if ( nBars>0 )
    cout << ">";
  for ( int i=nBars;i<nProgWidth;i++)
    cout << " ";
  cout << "| " << setw(3) << (int) (percent*100.) << "%";
  cout << flush;

}


void printHelp(const char *exeName, bool printFullHelp=true){
  
  if(printFullHelp){
    cout << bold;
    cout << endl;
    cout << "This program computes the mean and the sigma of the pedestal for a given rectangular region.\n";
    cout << normal;
  }
  cout << "==========================================================================\n";
  cout << yellow;
  cout << "\nUsage:\n";
  cout << "\t" << exeName << " <fits file name> Xmin Ymin  Xmax Ymax \n";
  cout << normal;
  cout << blue;
  cout << "For any problems or bugs contact Javier Tiffenberg <javiert@fnal.gov>\n\n";
  cout << normal;
  cout << "==========================================================================\n\n";
}


void computeParameters(const double* sArray, const long npix, double &mean, double &sigma){
  
  TNtuple pixNT("pixels","pixels","pix");
  
  for(long c=0;c<npix;++c) pixNT.Fill(sArray[c]);
  
  pixNT.Draw("pix","","goff");
  
  Double_t meanAll = pixNT.GetHistogram()->GetMean();
  Double_t sigmaAll = pixNT.GetHistogram()->GetRMS();
  
  ostringstream range1OSS;
  range1OSS << "TMath::Abs(pix-" << meanAll << ")<" << sigmaAll;
  
  pixNT.Draw("pix",range1OSS.str().c_str(),"goff");
  
  TFitResultPtr firstFit = pixNT.GetHistogram()->Fit("gaus","S Q");
  
  Double_t meanRange1 = firstFit->Parameter(1);
  Double_t sigmaRange1 = firstFit->Parameter(2);
  
//  TCanvas c;
  ostringstream range2OSS;
  range2OSS << "TMath::Abs(pix-" << meanRange1 << ")<" << sigmaRange1*3;
  pixNT.Draw("pix",range2OSS.str().c_str(),"");
  TFitResultPtr secondFit = pixNT.GetHistogram()->Fit("gaus","S Q");
  
  mean  = secondFit->Parameter(1);
  sigma = secondFit->Parameter(2);
//  c.WaitPrimitive();
  
  return;
}

struct ext_t{
  int e;
  double s1;
  double s2;
  ext_t(int ext, double sig1, double sig2){ ext = e; s1 = sig1; s2 = sig2;};
};


int computeStats(string inFile, int xMin, int yMin, int xMax, int yMax){
  
  cout << "computeStats\n";
  int status = 0;
  int single = 0;
  
  int singleHdu = -1;
  
  const char* inF = inFile.c_str();
  fitsfile *infPtr;   /* FITS file pointers defined in fitsio.h */
  /* Open the input file */
  fits_open_file(&infPtr, inF, READONLY, &status);
  if (status != 0){
    return(status);
  }
  
  int nhdu = 0;
  fits_get_num_hdus(infPtr, &nhdu, &status);
  
//   if (singleHdu>0){
//     single = 1; /* Copy only a single HDU if a specific extension was given */
//   }
  
  
  for (int n=1; n<=nhdu; ++n)  /* Main loop through each extension */
  {
    if (single){
      n = singleHdu;
    }
    
    /* get image dimensions and total number of pixels in image */
    int hdutype, bitpix, bytepix, naxis = 0, anynul;
    long naxes[9] = {1, 1, 1, 1, 1, 1, 1, 1, 1};
    double nulval = 0.;
    fits_movabs_hdu(infPtr, n, &hdutype, &status);
    for (int i = 0; i < 9; ++i) naxes[i] = 1;
    fits_get_img_param(infPtr, 9, &bitpix, &naxis, naxes, &status);
    long totpix = naxes[0] * naxes[1];
    
    /* Don't try to process data if the hdu is empty */    
    if (hdutype != IMAGE_HDU || naxis == 0 || totpix == 0){
      if(single) break;
      continue;
    }
    
    
    bytepix = abs(bitpix) / 8;
    if(bytepix!=1 && bytepix!=2 && bytepix!=4 && bytepix!=8){
      return -1000;
    }
    
//     if(gVerbosity){
//       cout << bold << "\rProcessing HDU: " << n << normal << flush;
//       showProgress(0,1);
//     }
    
    double mean1=0;
    double sigma1=0;
    {
      long inc[2]={1,1};
      long fpixel[2]={xMin,yMin};
      long lpixel[2]={xMax,yMax};
      const long npix = (lpixel[0]-fpixel[0]+1) * (lpixel[1]-fpixel[1]+1);
      double* sArray = new double[npix];
      fits_read_subset(infPtr, TDOUBLE, fpixel, lpixel, inc, &nulval, sArray, &anynul, &status);
      if (status != 0){
        return(status);
      }
      computeParameters(sArray,npix, mean1, sigma1);
      delete[] sArray;
    }
    
    
    cout << setw(3) << n << "\t" << mean1 << "\t" << sigma1 << endl;
    
    if (single) break;
  }
  
  /* Close the input file */
  fits_close_file(infPtr, &status);
  
  cout << green << "Stats computed.\n\n" << normal;
  return(status);
}

void checkArch(){
  if(sizeof(float)*CHAR_BIT!=32 || sizeof(double)*CHAR_BIT!=64){
    cout << red;
    cout << "\n ========================================================================================\n";
    cout << "   WARNING: the size of the float and double variables is non-standard in this computer.\n";
    cout << "   The program may malfunction or produce incorrect results\n";
    cout << " ========================================================================================\n";
    cout << normal;
  }
}

int processCommandLineArgs(const int argc, char *argv[], int &singleHdu, vector<string> &inFileList, string &outFile){
  
  if(argc == 1) return 1;
  
  bool outFileFlag = false;
  int opt=0;
  while ( (opt = getopt(argc, argv, "o:s:vVhH?")) != -1) {
    switch (opt) {
    case 'o':
      if(!outFileFlag){
        outFile = optarg;
        outFileFlag = true;
      }
      else{
        cerr << red << "\nError, can not set more than one output file!\n\n" << normal;
        return 2;
      }
      break;
    case 's':
      if(singleHdu<0){
        singleHdu = atoi(optarg);
      }
      else{
        cerr << red << "\nError, can not set more than one HDU!\n\n"  << normal;
        return 2;
      }
      break;
    case 'V':
    case 'v':
      gVerbosity = 1;
      break;
    case 'h':
    case 'H':
    default: /* '?' */
      return 1;
    }
  }
  
  if(!outFileFlag){
    cerr << red << "\nError: output filename missing.\n" << normal;
    return 2;
  }

  inFileList.clear();
  for(int i=optind; i<argc; ++i){
    inFileList.push_back(argv[i]);
    if(!fileExist(argv[i])){
      cout << red << "\nError reading input file: " << argv[i] <<"\nThe file doesn't exist!\n\n" << normal;
      return 1;
    }
  }
  
  return 0;
}

int main(int argc, char *argv[])
{
  
  TApplication myapp("myapp", 0, 0);
  
  checkArch(); //Check the size of the double and float variables.
  
  time_t start,end;
  double dif;
  time (&start);
  
  /* Do the actual processing */
  vector<ext_t> extData;
  
  int xMin, yMin, xMax, yMax;
  
  xMin = atoi(argv[2]);
  yMin = atoi(argv[3]);
  xMax = atoi(argv[4]);
  yMax = atoi(argv[5]);
  cout << "["<< xMin << ", "<< yMin << "]"<< " to ["<< xMax << ", "<< yMax << "]\n";
  
  int status = computeStats(argv[1], xMin, yMin, xMax, yMax);
  
  
  /* Report */
  time (&end);
  dif = difftime (end,start);
  if(gVerbosity) 
    cout << green << "All done!\n" << bold << "-> It took me " << dif << " seconds to do it!\n\n" << normal;

  return status;
}
