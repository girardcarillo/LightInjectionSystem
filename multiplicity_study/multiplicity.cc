// Author Cloé Girard-Carillo girardcarillo@lal.in2p3.fr
// Macro to compare amplitudes of secondary and primary LEDs, OM by OM
// root 'comparison.cc()'
#include <limits>
#include <string>
#include <iostream>
#include <fstream>

#include <TTree.h>
#include <TProfile.h>
#include <TFile.h>
#include <TObject.h>
#include <TH1.h>
#include <TH2.h>

#include "/home/girardcarillo/Workdir/SNPlot/RootDisplay.h"
#include "/home/girardcarillo/Workdir/SNPlot/EventTree.h"

using namespace std ;

const int column_tot_number=20 ;
const int row_tot_number=13 ;
const int area_number=5 ;

TH1F* GenerateMultiplicity(string filename) ;
void select_area(int area, int *xmin, int *xmax, int *ymin, int *ymax) ;
int transform_index(int x, int y) ;
void DrawAmplSpect(TH1F *histo[column_tot_number][row_tot_number],string fiber) ;

void multiplicity(string filename1,string run1, string filename2, string run2, string filename3, string run3){

  gStyle->SetOptStat(0) ;
  gStyle->SetOptTitle(0) ;
  TCanvas *c0 = new TCanvas("c0","c0",10,10,2000,1000) ;

  TH1F* multiplicity1 = GenerateMultiplicity(filename1) ;
  config_histo1D(multiplicity1,"","Event multiplicity","",1,1,MultiPlotColors(0)) ;
  multiplicity1->SetFillColorAlpha(MultiPlotColors(0),0.5) ;
  multiplicity1->SetTitle(run1.c_str()) ;

  TH1F* multiplicity2 = GenerateMultiplicity(filename2) ;
  config_histo1D(multiplicity2,"SAME","Event multiplicity","",1,1,MultiPlotColors(1)) ;
  multiplicity2->SetFillColorAlpha(MultiPlotColors(1),0.5) ;
  multiplicity2->SetTitle(run2.c_str()) ;

  // TH1F* multiplicity3 = GenerateMultiplicity(filename3) ;
  // config_histo1D(multiplicity3,"SAME","Event multiplicity","",1,1,MultiPlotColors(2)) ;
  // multiplicity3->SetFillColorAlpha(MultiPlotColors(2),0.5) ;
  // multiplicity3->SetTitle(run3.c_str()) ;

  multiplicity1->GetYaxis()->SetRangeUser(0,7000) ;

  c0->BuildLegend(0.77,0.79,0.87,0.88) ;
  c0->SaveAs("plots/mult.pdf") ;
}


TH1F* GenerateMultiplicity(string filename){


  //data file
  TFile *DataFile = new TFile(filename.c_str(),"READ") ;
  TTree *theTree = nullptr ;
  DataFile->GetObject("T",theTree) ;

  if (DataFile->IsOpen()) {
    cout << "File " << filename << " opened sucessfully" << endl ;
  }

  theTree = (TTree*)DataFile->Get("DataCut") ;

  theTree->SetBranchAddress("calo_peak",&calo_peak) ;
  theTree->SetBranchAddress("calo_row",&calo_row) ;
  theTree->SetBranchAddress("calo_column",&calo_column) ;
  theTree->SetBranchAddress("calo_multiplicity", &calo_multiplicity) ;

  TH1F* histo = new TH1F("mult","mult",60,0,60) ;

  int total_event = -1 ;
  for (Long64_t i=0 ;i<theTree->GetEntries() ;i++) {

    theTree->GetEntry(i) ;
    if (i%10000==0) cout << "event " << i << endl ;


    histo->Fill(calo_multiplicity) ;

    total_event = i ;
    if (i > 40000) {
      cout << "\033[1;31mwarning break at \033[0m" << i << endl ;
      break ;
    }


  }// end loop on event tree

  return histo ;

}

void select_area(int area, int *xmin, int *xmax, int *ymin, int *ymax){
  if (area == 1) {
    *xmin = 0 ; *xmax = 9 ;
    *ymin = 8 ; *ymax = 12 ;
  }
  if (area == 2) {
    *xmin = 10 ; *xmax = 19 ;
    *ymin = 8 ; *ymax = 12 ;
  }
  if (area == 3) {
    *xmin = 0 ; *xmax = 5 ;
    *ymin = 0 ; *ymax = 7 ;
  }
  if (area == 4) {
    *xmin = 6 ; *xmax = 13 ;
    *ymin = 0 ; *ymax = 7 ;
  }
  if (area == 5) {
    *xmin = 14 ; *xmax = 19 ;
    *ymin = 0 ; *ymax = 7 ;
  }
}

int transform_index(int x, int y){

  int index = -1 ;

  int matrice [column_tot_number][row_tot_number] ;
  for (int i = 0; i < column_tot_number; ++i) {
    for (int j = 0; j < row_tot_number; ++j) {
      matrice[i][j] = i*row_tot_number+j ;
    }
  }

  if (x<column_tot_number && y<row_tot_number) {
    index = matrice[x][y] ;
  }
  else {
    cout << "Function transform_index(int channel, int slot): Bad channel or slot number!!" << endl ;
  }

  return index ;
}
