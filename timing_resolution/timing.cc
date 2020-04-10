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


bool select_PM(vector<int> *vect_column,vector<int> *vect_row,int selected_column,int selected_row,int *hit) ;
void GenerateTimeSpect(string filename, TH1F *histo[column_tot_number][row_tot_number],double charges[column_tot_number][row_tot_number], int coinc_col, int coinc_row) ;
void select_area(int area, int *xmin, int *xmax, int *ymin, int *ymax) ;
int transform_index(int x, int y) ;
void DrawTimeSpect(TH1F *histo[column_tot_number][row_tot_number],string fiber) ;
void FitTimeSpect(TH1F *histo, double fit_parameters[6]) ;

void timing(string filename,bool enable_drawing = 0){

  TH1F *hsigma = new TH1F("sigma","",50, 0, 2) ;
  TH1F *henergy1 = new TH1F("sigma","",50, 0, 2) ;
  TH1F *henergy2 = new TH1F("sigma","",50, 0, 2) ;
  TH1F *hsigma_energy = new TH1F("sigma","",50, 0, 2) ;
  TProfile* psigma_energy = new TProfile("sigma_count","sigma_count", 50, 0.2, 1, 0, 1) ;

  double charges[column_tot_number][row_tot_number] ;

  int loop_counter = 0 ;

  // //area 5
  // for (int col = 14; col < 20; ++col) {
  //   for (int row = 0; row < 8; ++row) {

  // //area4
  for (int col = 6; col < 14; ++col) {
    for (int row = 0; row < 8; ++row) {

  // //tests
  // for (int col = 19; col < 20; ++col) {
  //   for (int row = 5; row < 6; ++row) {

      cout << "Loop #" << loop_counter << endl ;

      double charges[column_tot_number][row_tot_number];

      TH1F *htime[column_tot_number][row_tot_number] ;
      GenerateTimeSpect(filename,htime,charges,col,row) ;
      //cout << charges[col][row] << endl ;


      TH2D *h2mean = new TH2D ("Mean of deltat t distribution","", column_tot_number+2, -1, column_tot_number+1, row_tot_number+2, -1, row_tot_number+1) ;
      TH2D *h2chi2 = new TH2D ("chi2 deltat t distribution","", column_tot_number+2, -1, column_tot_number+1, row_tot_number+2, -1, row_tot_number+1) ;
      TH2D *h2sigma = new TH2D ("sigma deltat t distribution","", column_tot_number+2, -1, column_tot_number+1, row_tot_number+2, -1, row_tot_number+1) ;


      int area = 5 ;
      int imin,imax,jmin,jmax ;
      select_area(area,&imin,&imax,&jmin,&jmax) ;

      for (int i = imin ; i < imax+1 ; ++i) {
        for (int j = jmin ; j< jmax+1 ; ++j) {
          if (j!=12) {// on retire les rangées d'OM du haut et du bas car pas la même résolution en temps
            if (j!=0) {

              double fit_parameters_prim[6] ;
              FitTimeSpect(htime[i][j],fit_parameters_prim) ;

              if (fit_parameters_prim[0]!=0) {

                h2mean->SetBinContent(i+2,j+2,fit_parameters_prim[0]) ;
                h2chi2->SetBinContent(i+2,j+2,fit_parameters_prim[4]/fit_parameters_prim[5]) ;
                h2sigma->SetBinContent(i+2,j+2,fit_parameters_prim[1]) ;

                hsigma->Fill(fit_parameters_prim[1]) ;
                hsigma_energy->Fill(fit_parameters_prim[1]/sqrt(charges[i][j])) ;
                psigma_energy->Fill(charges[i][j],pow(fit_parameters_prim[1],2)) ;
                //cout << charges[i][j] << " " << pow(fit_parameters_prim[1],2) << endl ;

              }
            }
          }
        }
      }

      if (enable_drawing) {


        // gStyle->SetOptFit(1) ;
        // for (int i = 0; i < column_tot_number; ++i) {
        //   for (int j = 0; j < row_tot_number; ++j) {
        //     if (htime[i][j]->GetEntries()!=0) {
        //       TCanvas *c1 = new TCanvas("c1","c1",10,10,2000,1000) ;
        //       config_histo1D(htime[i][j],"","#Delta t (ns)","#Counts",1,1,1) ;
        //       c1->SaveAs(Form("plots/times_%d_%d.png",i,j)) ;
        //     }
        //   }
        // }


        TCanvas *c2 = new TCanvas("c2","c2",10,10,2000,1000) ;
        h2mean->Draw("COLZTEXT") ;
        c2->SaveAs("plots/map_mean.pdf") ;

        TCanvas *c3 = new TCanvas("c3","c3",10,10,2000,1000) ;
        h2chi2->Draw("COLZTEXT") ;
        c3->SaveAs("plots/map_chi2.pdf") ;

        TCanvas *c4 = new TCanvas("c4","c4",10,10,2000,1000) ;
        h2sigma->Draw("COLZTEXT") ;
        c4->SaveAs("plots/map_sigma.pdf") ;

      }

      loop_counter++ ;
    }
  }

  // TCanvas *c5 = new TCanvas("c5","c5",10,10,2000,1000) ;
  // config_histo1D(hsigma,"","#sigma (ns)","#Counts",1,1,1) ;
  // c5->SaveAs("plots/sigma.png") ;

  // TCanvas *c6 = new TCanvas("c6","c6",10,10,2000,1000) ;
  // config_histo1D(hsigma_energy,"","#sigma/(#sqrt{E_{1}}+#sqrt{E_{2}})","#Counts",1,1,1) ;
  // c6->SaveAs("plots/sigma_energy.png") ;

  gStyle->SetOptStat("e") ;
  TCanvas *c7 = new TCanvas("c7","c7",10,10,2000,1000) ;
  config_profile(psigma_energy, "", "1/(Q_{1}+Q_{2})","#sigma^{2}","",1,MultiPlotColors(0)) ;
  c7->SaveAs("plots/p_sigma_energy.png") ;



}// end macro

void FitTimeSpect(TH1F *histo, double fit_parameters[6]){

  double range_neg = histo->GetMean()-histo->GetStdDev()*100 ;
  double range_pos = histo->GetMean()+histo->GetStdDev()*100 ;
  TF1 *f1 = new TF1("f1","gaus",range_neg,range_pos) ;
  histo->Fit("f1","QR") ;

  bool fitting = 0 ;
  if (f1->GetChisquare()/f1->GetNDF()<4&&f1->GetChisquare()/f1->GetNDF()>0.01) {
    fitting = 1 ;
    // cout << "bool1 = " << fitting << endl ;
  }
  else {
    fitting = 0 ;
    //cout << "histogram " << histo->GetName() << ": fit parameters not passing chi2 condition" << endl ;
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

void DrawTimeSpect(TH1F *histo[column_tot_number][row_tot_number],string fiber){

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
    // histo[imin][jmin]->GetXaxis()->SetLimits(min_x+0.5*min_x, 0); a remplacer par hstack
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

void GenerateTimeSpect(string filename, TH1F *histo[column_tot_number][row_tot_number],double charges[column_tot_number][row_tot_number], int coinc_col, int coinc_row){

  TH2D *h2coincidence = new TH2D ("counts in coincidence","", column_tot_number+2, -1, column_tot_number+1, row_tot_number+2, -1, row_tot_number+1) ;

  //data file
  TFile *DataFile = new TFile(filename.c_str(),"READ") ;
  TTree *theTree = nullptr ;
  DataFile->GetObject("T",theTree) ;

  if (DataFile->IsOpen()) {
    cout << "File " << filename << " opened sucessfully" << endl ;
  }

  theTree = (TTree*)DataFile->Get("DataCut") ;

  theTree->SetBranchAddress("calo_multiplicity", &calo_multiplicity) ;
  theTree->SetBranchAddress("trigger_id",&trigger_id) ;
  theTree->SetBranchAddress("calo_peak",&calo_peak) ;
  theTree->SetBranchAddress("calo_charge",&calo_charge) ;
  theTree->SetBranchAddress("calo_row",&calo_row) ;
  theTree->SetBranchAddress("calo_column",&calo_column) ;
  theTree->SetBranchAddress("calo_time",&calo_time) ;
  theTree->SetBranchAddress("channel_raw_tdc", &channel_raw_tdc) ;

  vector <vector <double> > counts(column_tot_number) ;
  vector <vector <double> > charges1(column_tot_number) ;
  vector <vector <double> > charges2(column_tot_number) ;
  vector <vector <vector <double> > > times(column_tot_number) ;

  for (int j=0 ;j<column_tot_number ;j++) {
    times[j] = vector<vector <double> >(row_tot_number) ;
    counts[j]=vector<double>(row_tot_number);
    charges1[j]=vector<double>(row_tot_number);
    charges2[j]=vector<double>(row_tot_number);
  }


  cout << theTree->GetEntries() << " entries" << endl ;
  for (Long64_t i=0 ;i<theTree->GetEntries() ;i++) {

    theTree->GetEntry(i) ;
    if (i%100000==0) cout << "event " << i << endl ;

    int hit=-1 ;
    bool selected_PMT = select_PM(calo_column,calo_row,coinc_col,coinc_row,&hit) ;

    if (selected_PMT) {


      for (int j = 0 ; j < calo_row->size() ; ++j) {

        if (!isnan(calo_time->at(j))&&!isnan(calo_time->at(hit))) {// guaranty no pb with eventual waveforms not analysed by SNFEE

          if (calo_row->size() > 1) {// For coincidences

            if (calo_column->at(j)>13&&calo_row->at(j)<8) {// Same area

              double dtdc = channel_raw_tdc->at(hit) - channel_raw_tdc->at(j) ;
              double dt=calo_time->at(hit)-calo_time->at(j) ;

              if (trigger_id->at(hit)==trigger_id->at(j)) {

                if (fabs(dt)<100) {
                  // pb with several calo_times really small/big -> to be investigated
                  times[calo_column->at(j)][calo_row->at(j)].push_back(dt) ;
                  h2coincidence->Fill(calo_column->at(j),calo_row->at(j)) ;
                }

                if (fabs(calo_charge->at(j))<50&&fabs(calo_charge->at(hit))<50) {
                  charges1[calo_column->at(j)][calo_row->at(j)] += calo_charge->at(j) ;
                  charges2[calo_column->at(j)][calo_row->at(j)] += calo_charge->at(hit) ;
                  counts[calo_column->at(j)][calo_row->at(j)] ++ ;
                }

              }
            }
          }
        }
      }
      //cout << charge2 << " " << charge1/coinc_counter << endl ;
    }

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
  gStyle->SetOptFit(1) ;

  for (int j=0 ;j<column_tot_number ;j++) {
    for (int i=0 ;i<row_tot_number ;i++) {

      if (charges1[j][i] != 0 && charges2[j][i] != 0) {
        charges1[j][i] = charges1[j][i]/counts[j][i] ;
        charges2[j][i] = charges2[j][i]/counts[j][i] ;
        charges[j][i] =  1./charges1[j][i]+1./charges2[j][i] ;
      }

      mean = TMath::Mean(times[j][i].begin(),times[j][i].end()) ;
      size = times[j][i].size() ;
      rms =  TMath::RMS(times[j][i].begin(), times[j][i].end()) ;
      inf = mean-5*rms ;
      sup = mean+5*rms ;

      if (times[j][i].size()>100) {
        bin = size/100 ;
      }
      else {
        bin = 10;
      }

      histo[j][i] = new TH1F(Form("%d_%d",j,i),Form("%d_%d",j,i),bin,inf,sup) ;
      for (Long64_t k = 0; k < size; ++k) {

        if (times[j][i].at(k) > inf && times[j][i].at(k) < sup) {

          histo[j][i]->Fill(times[j][i].at(k)) ;

        }
      }
    }
  }


  // TCanvas *c2 = new TCanvas("c2","c2",10,10,2000,1000) ;
  // h2coincidence->Draw("COLZTEXT") ;
  // c2->SaveAs("plots/coincidences.pdf") ;

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

bool select_PM(vector<int> *vect_column,vector<int> *vect_row,int selected_column,int selected_row,int *hit){
  bool flag_test=0 ;
  for (int i = 0 ; i < vect_row->size() ; ++i) {
    if (vect_column->at(i)==selected_column) {
      if (vect_row->at(i)==selected_row) {
        *hit=i ;
        flag_test=1 ;
      }
    }
  }
  return flag_test ;
}
