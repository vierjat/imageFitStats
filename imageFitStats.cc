#include <string.h>
#include <stdio.h>
#include "fitsio.h"

#include <iostream>
#include <sstream>

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

#include "globalConstants.h"

using namespace std;

int deleteFile(const char *fileName){
  cout << yellow;
  cout << "Will overwrite: " << fileName << endl << endl;
  cout << normal;
  return unlink(fileName);
}

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


void printCopyHelp(const char *exeName, bool printFullHelp=false){
  
  if(printFullHelp){
    cout << bold;
    cout << endl;
    cout << "This program computes the median images from N input fit files.\n";
    cout << "It handles all the available HDUs. The HDUs in the output fit file\n";
    cout << "will be:\n";
    cout << " * float (32bits) for:  int8, int16, int32 and float input images\n";
    cout << " * double (64bits) for: double input images\n";
    cout << normal;
  }
  cout << "==========================================================================\n";
  cout << yellow;
  cout << "\nUsage:\n";
  cout << "  "   << exeName << " <input file 1> .. <input file N> -o <output filename> [-v for verbosity] \n\n";
  cout << normal;
  cout << blue;
  cout << "For any problems or bugs contact Javier Tiffenberg <javiert@fnal.gov>\n\n";
  cout << normal;
  cout << "==========================================================================\n\n";
}

string bitpix2TypeName(int bitpix){
  
  string typeName;
  switch(bitpix) {
      case BYTE_IMG:
          typeName = "BYTE(8 bits)";
          break;
      case SHORT_IMG:
          typeName = "SHORT(16 bits)";
          break;
      case LONG_IMG:
          typeName = "INT(32 bits)";
          break;
      case FLOAT_IMG:
          typeName = "FLOAT(32 bits)";
          break;
      case DOUBLE_IMG:
          typeName = "DOUBLE(64 bits)";
          break;
      default:
          typeName = "UNKNOWN";
  }
  return typeName;
}


void sortAndComputeMedian( vector< vector <double> > &vLinePix, vector<double> &vMedian ){
  
  int npix = vLinePix.size();
  vMedian.clear();
  vMedian.resize( npix );
  
  for(int c=0;c<npix;++c) std::sort(vLinePix[c].begin(), vLinePix[c].end());
  int nElementsPerPix = vLinePix[0].size();
  bool isOdd = !!(nElementsPerPix & 1);
  
  int m = nElementsPerPix/2;
  if( isOdd )
    for(int c=0;c<npix;++c) vMedian[c] = vLinePix[c][m];
  else{
    for(int c=0;c<npix;++c){
      vMedian[c] = ( vLinePix[c][m-1] + vLinePix[c][m] )/2.;
    }
  }

}


double computeMean(const double* sArray, const long npix){
  
  double sum=0;
  for(long c=0;c<npix;++c) sum+=sArray[c];
  return sum/npix;
}

double computeSigma(const double* sArray, const long npix){
  
  double sum=0;
  for(long c=0;c<npix;++c) sum+=sArray[c];
  double mean = sum/npix;
  
  double var=0;
  for(long c=0;c<npix;++c){
    double diff=(sArray[c]-mean);
    var+=diff*diff;
  }
  var/=npix;
  return sqrt(var);
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
  
  
 
  
  const int nHDUsToProcess = (single>0)? 1 : nhdu;
  
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
      mean1 = computeMean(sArray,npix);
      sigma1 = computeSigma(sArray,npix);
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
