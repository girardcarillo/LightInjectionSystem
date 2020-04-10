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

int GenerateAmplitudes(string filename, TH1F *histo[column_tot_number][row_tot_number]) ;
void select_area(int area, int *xmin, int *xmax, int *ymin, int *ymax) ;
int transform_index(int x, int y) ;
void DrawAmplSpect(TH1F *histo[column_tot_number][row_tot_number],string fiber) ;
// void FitAmplSpect(TH1F *histo[column_tot_number][row_tot_number]) ;
void FitAmplSpect(TH1F *histo, double fit_parameters[6]) ;

void comparison(string filename_prim, string filename_sec,bool enable_drawing = 0){

  TH1F *hamplSpect_prim[column_tot_number][row_tot_number] ;
  TH1F *hamplSpect_sec[column_tot_number][row_tot_number] ;
  int total_event_prim = GenerateAmplitudes(filename_prim,hamplSpect_prim) ;
  int total_event_sec = GenerateAmplitudes(filename_sec,hamplSpect_sec) ;

  TH1F *hmean_ampl_tot = new TH1F("total","",20, -2000, 0) ;
  TH1F *hsigma_ampl_tot = new TH1F("total","",20, 0, 45) ;

  TH2D *h2mean = new TH2D ("Mean of deltat t distribution","", column_tot_number+2, -1, column_tot_number+1, row_tot_number+2, -1, row_tot_number+1) ;

  TH1F *h2ampl_tot = new TH1F("Mean amplitude","", 100, 0, 2000) ;

  TProfile *hampl[area_number] ;
  TH1F *hratio[area_number] ;
  TH1F *h2ampl[area_number] ;
  TH1F *hmean_ampl[area_number] ;
  TH1F *hsigma_ampl[area_number] ;
  for (int i=0 ;i<area_number ;i++) {
    hampl[i] = new TProfile(Form("Area %d",i+1),"", 100, 0, 260, 0, 5000) ;
    hratio[i] = new TH1F(Form("Area %d",i+1),"",40, 0, 0) ;
    h2ampl[i] = new TH1F(Form("Area %d",i+1),"", 100, 0, 2000) ;
    hmean_ampl[i] = new TH1F(Form("Area %d",i+1),"",20, -2000, 0) ;
    hsigma_ampl[i] = new TH1F(Form("Area %d",i+1),"",20, 0, 45) ;
  }

  for (int area = 1; area < area_number+1; ++area) {

    int imin,imax,jmin,jmax ;
    select_area(area,&imin,&imax,&jmin,&jmax) ;


    for (int i = imin ; i < imax+1 ; ++i) {
      for (int j = jmin ; j< jmax+1 ; ++j) {
        if (j!=12) {// on retire la rangée d'OM du haut qui compte plus et qui n'a pas de fibre sec
          if (j!=0) {// on retire la rangée d'OM du bas qui n'a pas de fibre sec

            double fit_parameters_prim[6] ;
            double fit_parameters_sec[6] ;
            FitAmplSpect(hamplSpect_prim[i][j],fit_parameters_prim) ;
            FitAmplSpect(hamplSpect_sec[i][j],fit_parameters_sec) ;

            if (fit_parameters_prim[0]!=0) {

              hratio[area-1]->Fill(fit_parameters_prim[0]/fit_parameters_sec[0]);
              h2ampl[area-1]->Fill(-fit_parameters_prim[0],-fit_parameters_sec[0]) ;
              hmean_ampl[area-1]->Fill(fit_parameters_prim[0]) ;
              hsigma_ampl[area-1]->Fill(fit_parameters_prim[1]) ;
              hmean_ampl_tot->Fill(fit_parameters_prim[0]) ;
              hsigma_ampl_tot->Fill(fit_parameters_prim[1]) ;
              h2ampl_tot->Fill(-fit_parameters_prim[0],-fit_parameters_sec[0]) ;
              h2mean->SetBinContent(i+2,j+2,fit_parameters_prim[0]) ;
            }

          }
        }
      }
    }
  }

  // gStyle->SetOptFit(1) ;
  // gStyle->SetOptStat(0) ;

  // TCanvas *c4 = new TCanvas("c4","c4",10,10,2000,1000) ;
  // for (int i = 0; i < column_tot_number; ++i) {
  //   for (int j = 0; j < row_tot_number; ++j) {
  //     hamplSpect_prim[i][j]->Draw() ;
  //     c4->SaveAs(Form("tmp_%d_%d.png",i,j)) ;
  //   }
  // }

  // h2mean->Draw("COLORTEXT") ;
  // h2ampl_tot->Draw("E") ;


  //TCanvas *c11 = new TCanvas("c11","c11",10,10,2000,1000) ;

  // THStack *hs = new THStack("hs","");
  // hs->Add(hamplSpect_prim[0][1]) ;
  // hamplSpect_prim[0][1]->SetLineColor(kMagenta+2) ;
  // hs->Add(hamplSpect_prim[0][2]) ;
  // hamplSpect_prim[0][2]->SetLineColor(kOrange+7) ;
  // hs->Add(hamplSpect_prim[0][3]) ;
  // hamplSpect_prim[0][3]->SetLineColor(kCyan+2) ;
  // hs->Draw("nostack");
  // hs->GetXaxis()->SetTitle("Amplitude (mV)") ;
  // c11->BuildLegend(0.82,0.66,0.89,0.89) ;
  // c11->SaveAs("plots/ex_ampl_spect.pdf") ;



  // TCanvas *c0 = new TCanvas("c0","c0",10,10,2000,1000) ;

  // THStack *hs = new THStack("hs","");
  // hs->Add(hamplSpect_prim[0][1]) ;
  // hamplSpect_prim[0][1]->SetLineColor(kMagenta+2) ;
  // hs->Add(hamplSpect_prim[0][2]) ;
  // hamplSpect_prim[0][2]->SetLineColor(kOrange+7) ;
  // hs->Add(hamplSpect_prim[0][3]) ;
  // hamplSpect_prim[0][3]->SetLineColor(kCyan+2) ;
  // hs->Draw("nostack");
  // hs->GetXaxis()->SetTitle("Amplitude (mV)") ;
  // c0->BuildLegend(0.82,0.66,0.89,0.89) ;
  // c0->SaveAs("plots/ex_ampl_spect.pdf") ;





  if (enable_drawing) {


    // // generate a set of energy spectra, one for each group of OM (unlighted by the same LED)
    TCanvas *c10 = new TCanvas("c10","c10",10,10,2000,1000) ;

    // here for central region
    THStack *hs = new THStack("hs","");
    int counter_color = 0 ;
    for (int j=6 ;j<14 ;j++) {
      for (int i=0 ;i<2 ;i++) {
        hs->Add(hamplSpect_prim[j][i]) ;
        hamplSpect_prim[j][i]->SetLineColor(MultiPlotColors(counter_color)) ;
        counter_color++ ;
      }
    }

    hs->Draw("nostack");
    hs->GetXaxis()->SetTitle("Amplitude (mV)") ;
    hs->GetYaxis()->SetTitle("#Counts") ;
    c10->BuildLegend(0.82,0.28,0.89,0.89) ;
    c10->SaveAs("plots/ampl_area4.pdf") ;


    TCanvas *c6 = new TCanvas("c6","c6",10,10,2000,1000) ;
    hsigma_ampl[0]->Draw() ;
    config_histo1D(hsigma_ampl[0],"","Sigma (mV)","",1,1,MultiPlotColors(0)) ;
    hsigma_ampl[0]->SetLineColor(MultiPlotColors(0)) ;
    hsigma_ampl[0]->SetFillColorAlpha(MultiPlotColors(0),0.35) ;
    hsigma_ampl[0]->GetYaxis()->SetRangeUser(0,25) ;
    for (int i = 1; i < 5; ++i) {
      if (hsigma_ampl[i]->GetEntries()!=0) {
        config_histo1D(hsigma_ampl[i],"SAME","Sigma (mV)","",1,1,MultiPlotColors(i)) ;
        hsigma_ampl[i]->SetLineColor(MultiPlotColors(i)) ;
        hsigma_ampl[i]->SetFillColorAlpha(MultiPlotColors(i),0.35) ;
      }
    }
    hsigma_ampl_tot->Draw("SAME") ;
    hsigma_ampl_tot->SetLineColor(15) ;
    hsigma_ampl_tot->SetName("Total") ;
    c6->BuildLegend(0.897,0.703,0.999,0.998) ;
    c6->SaveAs("plots/1d_sigma.pdf") ;

    gStyle->SetOptStat(0) ;
    TCanvas *c5 = new TCanvas("c5","c5",10,10,2000,1000) ;
    config_histo1D(hmean_ampl[0],"","Mean amplitude (mV)","",1,1,MultiPlotColors(0)) ;
    hmean_ampl[0]->SetFillColorAlpha(MultiPlotColors(0),0.35) ;
    hmean_ampl[0]->GetYaxis()->SetRangeUser(0,16) ;
    for (int i = 1; i < 5; ++i) {
      if (hmean_ampl[i]->GetEntries()!=0) {
        config_histo1D(hmean_ampl[i],"SAME","Mean amplitude (mV)","",1,1,MultiPlotColors(i)) ;
        hmean_ampl[i]->SetFillColorAlpha(MultiPlotColors(i),0.35) ;
      }
    }
    hmean_ampl_tot->Draw("SAME") ;
    hmean_ampl_tot->SetLineColor(15) ;
    hmean_ampl_tot->SetName("Total") ;
    c5->BuildLegend(0.897,0.703,0.999,0.998) ;
    c5->SaveAs("plots/1d_mean.pdf") ;


    TCanvas *c4 = new TCanvas("c4","c4",10,10,2000,1000) ;
    gStyle->SetOptStat(0) ;

    config_histo1D(h2ampl[0],"","Mean amplitude (mV) primary","Mean amplitude (mV) secondary",1,1,MultiPlotColors(0)) ;

    for (int area = 2; area < area_number+1; ++area) {
      config_histo1D(h2ampl[area-1],"SAME","Mean amplitude (mV) primary","Mean amplitude (mV) secondary",1,1,MultiPlotColors(area-1)) ;
    }

    c4->BuildLegend(0.897,0.703,0.999,0.998) ;

    TF1 *f1 = new TF1("f1","0.85*x",0,2000) ;
    f1->SetLineWidth(1) ;
    f1->SetLineStyle(1) ;
    f1->SetLineColor(kGray+2) ;
    f1->Draw("SAME") ;

    c4->SaveAs("plots/amplitude_prim_sec.pdf") ;

    double maximum_y_ratio = -1 ;
    double maximum_x_ratio = -1 ;
    for (int i = 0; i < area_number; ++i) {

      if (hratio[i]->GetMaximum() > maximum_y_ratio) {
        maximum_y_ratio = hratio[i]->GetMaximum() ;
      }

      if (hratio[i]->GetXaxis()->GetBinCenter(hratio[i]->GetMaximumBin()) > maximum_x_ratio) {
        maximum_x_ratio = hratio[i]->GetXaxis()->GetBinCenter(hratio[i]->GetMaximumBin());
      }
    }


    TCanvas *c0 = new TCanvas("c0","c0",10,10,2000,1000) ;
    gStyle->SetOptStat(0) ;

    config_histo1D(hratio[0],"E","Ratio","",1,1,MultiPlotColors(0)) ;
    hratio[0]->GetXaxis()->SetRangeUser(0, maximum_x_ratio+0.9*maximum_x_ratio);
    hratio[0]->GetYaxis()->SetRangeUser(0,maximum_y_ratio+0.9*maximum_y_ratio) ;

    for (int area = 2; area < area_number+1; ++area) {
      config_histo1D(hratio[area-1],"ESAME","Ratio","",1,1,MultiPlotColors(area-1)) ;
    }

    // config_histo1D(hratio[1],"E","Ratio","Counts",1,1,MultiPlotColors(1)) ;
    // config_histo1D(hratio[3],"ESAME","Ratio","Counts",1,1,MultiPlotColors(0)) ;

    c0->BuildLegend(0.897,0.703,0.999,0.998) ;

    TLine *line = new TLine(1,0,1,maximum_y_ratio+0.9*maximum_y_ratio);
    line->SetLineStyle(2) ;
    line->SetLineColor(kGray+2) ;
    line->Draw("SAME") ;

    c0->SaveAs("plots/ratio.pdf") ;

    DrawAmplSpect(hamplSpect_prim,"prim") ;
    DrawAmplSpect(hamplSpect_sec,"sec") ;

  }


}// end macro

void FitAmplSpect(TH1F *histo, double fit_parameters[6]){

  double range_neg = histo->GetMean()-histo->GetStdDev()*10 ;
  double range_pos = histo->GetMean()+histo->GetStdDev()*10 ;
  TF1 *f1 = new TF1("f1","gaus",range_neg,range_pos) ;
  histo->Fit("f1","QR0") ;

  bool fitting = 0 ;
  if (f1->GetParameter(1)<-50-2*f1->GetParameter(2)) {
    if (f1->GetChisquare()/f1->GetNDF()<4&&f1->GetChisquare()/f1->GetNDF()>0.01) {
      fitting = 1 ;
      // cout << "bool1 = " << fitting << endl ;
    }
    else {
      fitting = 0 ;
      cout << "histogram " << histo->GetName() << ": fit parameters not passing chi2 condition" << endl ;
      // cout << "bool0 = " << fitting << endl ;
    }
  }
  else {
    fitting = 0 ;
    cout << "histogram " << histo->GetName() << ": fit parameters not passing trigger conditions" << endl ;
    // cout << "bool0 = " << fitting << endl ;
  }


  if (fitting) {
    fit_parameters[0] = f1->GetParameter(1) ;
    fit_parameters[1] = f1->GetParameter(2) ;
    fit_parameters[2] = f1->GetParError(1) ;
    fit_parameters[3] = f1->GetParError(2) ;
    fit_parameters[4] = f1->GetChisquare() ;
    fit_parameters[5] = f1->GetNDF() ;
  }
  else {
    for (int i = 0; i < 6; ++i) {
      fit_parameters[i] = 0 ;
    }
  }

  delete f1 ;

}

void DrawAmplSpect(TH1F *histo[column_tot_number][row_tot_number],string fiber){

  for (int area = 1; area < area_number+1; ++area) {
    TCanvas *c1 = new TCanvas("c1","c1",10,10,2000,1000) ;
    gStyle->SetOptStat(0) ;

    int imin,imax,jmin,jmax ;
    select_area(area,&imin,&imax,&jmin,&jmax) ;

    int min_x = 1 ;
    int max_y = -1 ;
    for (int i = imin ; i < imax+1 ; ++i) {
      for (int j = jmin ; j< jmax+1 ; ++j) {
        if (histo[i][j]->GetMaximum() > max_y) {
          max_y = histo[i][j]->GetMaximum() ;
        }

        if (histo[i][j]->GetXaxis()->GetBinCenter(histo[i][j]->GetMaximumBin()) < min_x) {
          min_x = histo[i][j]->GetXaxis()->GetBinCenter(histo[i][j]->GetMaximumBin());
        }
      }
    }

    config_histo1D(histo[imin][jmin],"","Amplitude (mV)","Counts",1,1,MultiPlotColors(0,histo[imin][jmin])) ;
    histo[imin][jmin]->GetXaxis()->SetLimits(min_x+0.5*min_x, 0);
    histo[imin][jmin]->GetYaxis()->SetRangeUser(0,max_y+0.2*max_y) ;

    for (int i = imin ; i < imax+1 ; ++i) {
      for (int j = jmin ; j< jmax+1 ; ++j) {
        if (j!=12) {// on retire la rangée d'OM du haut qui compte plus
          if (i!=imin||j!=jmin) {

            config_histo1D(histo[i][j],"SAME","Amplitude (mV)","Counts",1,1,MultiPlotColors(transform_index(i,j),histo[i][j])) ;

          }
        }
      }
    }

    gStyle->SetOptTitle(0);

    c1->BuildLegend(0.909,0.030,0.953,0.984,Form("Area %d",area)) ;
    c1->SaveAs(Form("plots/Amplitudes_%s_%d.png",fiber.c_str(),area)) ;
  }

}

int GenerateAmplitudes(string filename, TH1F *histo[column_tot_number][row_tot_number]){

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

  vector <vector <vector <double> > > amplitudes(column_tot_number) ;
  for (int j=0 ;j<column_tot_number ;j++) {
    amplitudes[j]=vector<vector <double> >(row_tot_number) ;
  }

  int total_event = -1 ;
  for (Long64_t i=0 ;i<theTree->GetEntries() ;i++) {

    theTree->GetEntry(i) ;
    if (i%10000==0) cout << "event " << i << endl ;


    for (int j = 0 ; j < calo_row->size() ; ++j) {

      if (!isnan(calo_peak->at(j))) {// guaranty no pb with eventual waveforms not analysed by SNFEE

        amplitudes[calo_column->at(j)][calo_row->at(j)].push_back(calo_peak->at(j)) ;

      }
    }

    total_event = i ;

    // if (i > 1e3) {
    //   cout << "\033[1;31mwarning break at \033[0m" << i << endl ;
    //   break ;
    // }


  }// end loop on event tree

  double mean = 0. ;
  double rms = 0. ;
  double inf = 0. ;
  double sup = 0. ;
  Long64_t size = 0. ;
  int bin  = 0 ;

  for (int j=0 ;j<column_tot_number ;j++) {
    for (int i=0 ;i<row_tot_number ;i++) {

      //      if (amplitudes[j][i].size()) {

      mean = TMath::Mean(amplitudes[j][i].begin(),amplitudes[j][i].end()) ;
      size = amplitudes[j][i].size() ;
      rms =  TMath::RMS(amplitudes[j][i].begin(), amplitudes[j][i].end()) ;
      inf = mean-3*rms ;
      sup = mean+3*rms ;

      if (amplitudes[j][i].size()>100) {
        bin = size/100 ;
      }
      else {
        bin = 10;
      }

      histo[j][i] = new TH1F(Form("%d_%d",j,i),Form("%d_%d",j,i),bin,inf,sup) ;
      for (Long64_t k = 0; k < size; ++k) {
        if (amplitudes[j][i].at(k) > inf && amplitudes[j][i].at(k) < sup) {
          histo[j][i]->Fill(amplitudes[j][i].at(k)) ;
        }
      }
      //      }
    }
  }

  return total_event ;

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
