// Author Cloé Girard-Carillo girardcarillo@lal.in2p3.fr

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

typedef numeric_limits<double> dbl ;

const int row_tot_number=13 ;
const int column_tot_number=20 ;
const int area_number=5 ;

bool select_PM(vector<int> *vect_column,vector<int> *vect_row,int selected_column,int selected_row,int *hit) ;
void fit_delta_t(TH1F* histo, double fit_parameters[4]) ;
int transform_index(int x, int y) ;
int index_area(int area, int x, int y) ;
void select_area(int area, int *xmin, int *xmax, int *ymin, int *ymax) ;

void RootAnalyzer(string filename, string correctedTimesFilename, string str_wall, int coinc_col, int coinc_row, bool enable_drawing = 0){

  int breaking_event = 1e9 ;

  // corrected times from reflecto files
  TTree *tree_correctedTimes = new TTree("ntuple_correctedTimes","corrected times from reflecto") ;

  tree_correctedTimes->ReadFile(correctedTimesFilename.c_str(),"wall/I:col_cable/I:row_cable/I:CT/F") ;

  Int_t wall,col_cable,row_cable ;
  Float_t CT ;
  tree_correctedTimes->SetBranchAddress("wall",&wall) ;
  tree_correctedTimes->SetBranchAddress("col_cable",&col_cable) ;
  tree_correctedTimes->SetBranchAddress("row_cable",&row_cable) ;
  tree_correctedTimes->SetBranchAddress("CT",&CT) ;


  //data file
  TFile *DataFile = new TFile(filename.c_str(),"READ") ;
  TTree *theTree = nullptr ;
  DataFile->GetObject("T",theTree) ;

  if (DataFile->IsOpen()) {
    cout << "File " << filename << " opened sucessfully" << endl ;
  }

  theTree = (TTree*)DataFile->Get("DataCut") ;

  theTree->SetBranchAddress("trigger_id",&trigger_id) ;
  theTree->SetBranchAddress("calo_row",&calo_row) ;
  theTree->SetBranchAddress("calo_column",&calo_column) ;
  theTree->SetBranchAddress("calo_id",&calo_id) ;
  theTree->SetBranchAddress("calo_module",&calo_module) ;
  theTree->SetBranchAddress("channel_raw_tdc", &channel_raw_tdc) ;
  theTree->SetBranchAddress("calo_time",&calo_time) ;
  theTree->SetBranchAddress("calo_energy",&calo_energy) ;
  theTree->SetBranchAddress("calo_charge",&calo_charge) ;
  theTree->SetBranchAddress("calo_peak",&calo_peak) ;
  theTree->SetBranchAddress("calo_number",&calo_number) ;
  theTree->SetBranchAddress("calo_multiplicity", &calo_multiplicity) ;

  TH1F *hamplSpect[column_tot_number][row_tot_number] ;
  TH1F *htime[column_tot_number][row_tot_number] ;
  for (int j=0 ;j<column_tot_number ;j++) {
    for (int i=0 ;i<row_tot_number ;i++) {
      htime[j][i] = new TH1F(Form("coicidences between OMs [%d:%d]&[%d:%d] ",coinc_col,coinc_row,j,i),Form("coicidences between OMs [%d:%d]&[%d:%d]",coinc_col,coinc_row,j,i),10000,-100,100) ;
      hamplSpect[j][i] = new TH1F(Form("%d_%d",j,i),Form("%d_%d",j,i),100,0,0) ;
    }
  }

  TH2D *h2counts = new TH2D ("counts","", column_tot_number+2, -1, column_tot_number+1, row_tot_number+2, -1, row_tot_number+1) ;
  TH2D *h2coincidence = new TH2D ("counts in coincidence","", column_tot_number+2, -1, column_tot_number+1, row_tot_number+2, -1, row_tot_number+1) ;
  TH2D *h2mean = new TH2D ("Mean of deltat t distribution","", column_tot_number+2, -1, column_tot_number+1, row_tot_number+2, -1, row_tot_number+1) ;

  TH1F *hmultiplicity = new TH1F("Multiplicity","",250, 0, 250) ;
  TH1F *hcharge = new TH1F("Charge","",1000, 0, 20) ;
  TH1F *henergy = new TH1F("Energy","",200, 0, 25) ;
  TH1F *hsigma = new TH1F("Sigma","",50, 0, 0) ;
  TH1F *hamplitude = new TH1F("Amplitude","",1000, 0, 0) ;
  TH1F *htemps = new TH1F("Time","",500, 50, 550) ;
  TH1F *hpulseTime = new TH1F("Time of pulse","",500, 0,0) ;

  TProfile *hampl_count[area_number] ;
  TProfile *hampl[area_number] ;
  TH1F *hcount[area_number] ;
  for (int i=0 ;i<area_number ;i++) {
    hcount[i] = new TH1F(Form("Area %d",i+1),"",260, 0,259) ;
    hampl_count[i] = new TProfile(Form("Area %d",i+1),"", 100, 0, 0.025*breaking_event, 0, 5000) ;
  }




  for (Long64_t i=0 ;i<theTree->GetEntries() ;i++) {
    theTree->GetEntry(i) ;
    if (i%10000==0) cout << "event " << i << endl ;

    int hit=-1 ;
    bool selected_PMT = select_PM(calo_column,calo_row,coinc_col,coinc_row,&hit) ;

    hmultiplicity->Fill(calo_multiplicity) ;

    for (int j = 0 ; j < calo_row->size() ; ++j) {

      if (!isnan(calo_time->at(j))) {

        // hpulseTime->Fill(calo_time->at(j)-channel_raw_tdc->at(j)*6.25e-9) ; // faut car déjà
        // corrigé du tdc dans le builder
        htemps->Fill(calo_time->at(j)) ;
        henergy->Fill(calo_energy->at(j)) ;
        hamplSpect[calo_column->at(j)][calo_row->at(j)]->Fill(calo_peak->at(j)) ;
        hcharge->Fill(calo_charge->at(j)) ;
        hamplitude->Fill(calo_peak->at(j)) ;

        //counting number of events for each PM
        if (calo_column->at(j) != coinc_col || calo_row->at(j) != coinc_row) {

          if (str_wall == "it") {
            h2counts->SetBinContent(21-calo_column->at(j),calo_row->at(j)+2,h2counts->GetBinContent(21-calo_column->at(j),calo_row->at(j)+2)+1) ;
          }

          else {
            h2counts->Fill(calo_column->at(j),calo_row->at(j)) ;
          }

        }

        // Coincidences
        if (calo_row->size() > 1) {

          if (selected_PMT) {

            double dt=calo_time->at(hit)-calo_time->at(j) ;
            htime[calo_column->at(j)][calo_row->at(j)]->Fill(dt) ;

            if (calo_column->at(j) != coinc_col || calo_row->at(j) != coinc_row) {

              if (str_wall == "it") {
                h2coincidence->SetBinContent(21-calo_column->at(j),calo_row->at(j)+2,h2coincidence->GetBinContent(21-calo_column->at(j),calo_row->at(j)+2)+1) ;
              }

              else if (str_wall == "fr") {
                h2coincidence->Fill(calo_column->at(j),calo_row->at(j)) ;
              }

              else {
                cout << "Wrong wall id" << endl ;
                break ;
              }

            }
          }
        }
      }//end for loop on trigger_id
    }//end for loop on events

    if (i > breaking_event) {
      cout << "\033[1;31mwarning break at \033[0m" << i << endl ;
      break ;
    }


  }


  int OM_counter = 0 ;
  for (int i = 0 ; i < column_tot_number ; ++i) {
    for (int j = 0 ; j< row_tot_number ; ++j) {
      OM_counter++ ;


      if (htime[i][j]->GetEntries() == 0) {
        continue ;
      }

      else {

        // Fit Delta t ditributions
        if (i != coinc_col || j != coinc_row) {

          double fit_parameters[4] ;
          fit_delta_t(htime[i][j],fit_parameters) ;

          if (str_wall == "it") {
            h2mean->SetBinContent(21-i,j+2,fit_parameters[0]) ;
          }
          else {
            h2mean->SetBinContent(i+2,j+2,fit_parameters[0]) ;
          }

          hsigma->Fill(fit_parameters[2]) ;

        }
      }
    }
  }


  double mean_count = 0. ;
  int counter = 0 ;

  for (int area = 1; area < area_number+1; ++area) {

    int imin,imax,jmin,jmax ;
    select_area(area,&imin,&imax,&jmin,&jmax) ;

    for (int i = imin ; i < imax+1 ; ++i) {
      for (int j = jmin ; j< jmax+1 ; ++j) {
        if (j!=12) {// on retire la rangée d'OM du haut qui compte plus

          mean_count += h2counts->GetBinContent(i+2,j+2) ;
          counter++ ;

          int index_wall = transform_index(i,j) ;
          hcount[area-1]->SetBinContent(index_wall,h2counts->GetBinContent(i+2,j+2));
          hampl_count[area-1]->Fill(hamplSpect[i][j]->GetEntries(),-hamplSpect[i][j]->GetMean());
        }
      }
    }
  }




  // ///Drawing

  if (enable_drawing) {


    double maximum_ampl = -1 ;
    for (int i = 0; i < area_number; ++i) {
      if (hampl_count[i]->GetMaximum() > maximum_ampl) {
        maximum_ampl = hampl_count[i]->GetMaximum() ;
      }
    }

    TCanvas *c13 = new TCanvas("c13","c13",10,10,2000,1000) ;
    gStyle->SetOptStat(0) ;
    config_profile(hampl_count[0], "", "Counts", "Mean amplitude (mV)","",1,MultiPlotColors(0)) ;
    hampl_count[0]->GetYaxis()->SetRangeUser(0,maximum_ampl+0.1*maximum_ampl) ;
    for (int area = 2; area < area_number+1; ++area) {
      config_profile(hampl_count[area-1], "", "Counts", "Mean amplitude (mV)","SAME",1,MultiPlotColors(area-1)) ;
    }
    c13->BuildLegend(0.897,0.703,0.999,0.998) ;
    c13->SaveAs("plots/1Dampl_count.pdf") ;


    double maximum_count = -1 ;
    for (int i = 0; i < area_number; ++i) {
      if (hcount[i]->GetMaximum() > maximum_count) {
        maximum_count = hcount[i]->GetMaximum() ;
      }
    }

    TCanvas *c12 = new TCanvas("c12","c12",10,10,2000,1000) ;
    gStyle->SetOptStat(0) ;
    config_histo1D(hcount[0],"E","OM index","Counts",2,1,MultiPlotColors(0)) ;
    hcount[0]->LabelsOption("v");
    hcount[0]->SetMarkerStyle(8) ;
    hcount[0]->SetMarkerSize(1.3) ;
    hcount[0]->SetMarkerColor(MultiPlotColors(0)) ;
    hcount[0]->GetYaxis()->SetRangeUser(0,maximum_count+0.1*maximum_count) ;
    for (int area = 2; area < area_number+1; ++area) {
      config_histo1D(hcount[area-1],"ESAME","OM index","Counts",2,1,MultiPlotColors(area-1)) ;
      hcount[area-1]->SetMarkerStyle(area+18) ;
      hcount[area-1]->SetMarkerSize(1.3) ;
      hcount[area-1]->SetMarkerColor(MultiPlotColors(area-1)) ;
    }
    mean_count /= counter ;
    TLine *line = new TLine(0,mean_count,259,mean_count);
    line->SetLineStyle(2) ;
    line->SetLineColor(kGray+2) ;
    c12->BuildLegend(0.897,0.703,0.999,0.998) ;
    line->Draw("SAME") ;
    c12->SaveAs("plots/1Dcounts.pdf") ;


    gStyle->SetPaintTextFormat("1.2f") ;
    gStyle->SetOptStat(0) ;
    TCanvas *c1 = new TCanvas("c1","c1",10,10,2000,1000) ;
    config_histo2D(h2counts, "Number of events in each PMT", "Column","Row","COLZTEXT") ;
    if (str_wall == "it") {
      ReverseXAxis(h2counts) ;
    }
    c1->SaveAs("plots/counts.pdf") ;


    // // // generate a set of energy spectra, one for each column
    // TCanvas *c2 = new TCanvas("c2","c2",10,10,2000,1000) ;

    // for (int j=0 ;j<column_tot_number ;j++) {
    //   config_histo1D(hamplSpect[j][0],"","Amplitude","",3,1,1) ;
    //   hamplSpect[j][0]->GetYaxis()->SetRangeUser(0,1000) ;

    //   gStyle->SetOptTitle(0);
    //   TPaveLabel *title = new TPaveLabel(0.46,0.93,0.53,0.99,"new title");
    //   title->Draw();

    //   c2->SetTitle(Form("Energy Spectra for column %d",j)) ;
    //   int counter_color = 1 ;
    //   for (int i=1 ;i<row_tot_number ;i++) {
    //     counter_color++ ;
    //     config_histo1D(hamplSpect[j][i],"SAME","Amplitude","",3,1,MultiPlotColors(counter_color)) ;
    //   }
    //   c2->BuildLegend(0.82,0.35,0.89,0.89) ;
    //   c2->SaveAs(Form("plots/energy_col_%d.png",j)) ;
    // }

    // c10->BuildLegend(0.82,0.35,0.89,0.89) ;
    // c10->SaveAs("plots/energy_LEDgroup_cent.png") ;
    // //


    // TCanvas *c4 = new TCanvas("c4","c4",10,10,2000,1000) ;

    // for (int i = 0 ; i < column_tot_number ; ++i) {
    //   for (int j = 0 ; j < row_tot_number ; ++j) {
    //     if (htime[i][j]) {

    //       gStyle->SetOptStat(1) ;
    //       gStyle->SetOptFit(1) ;
    //       TString n = TString::Format("coinc%d_%d", i,j) ;
    //       htime[i][j]->Draw() ;
    //       htime[i][j]->GetXaxis()->SetRangeUser(htime[i][j]->GetMean()-htime[i][j]->GetStdDev()*10,htime[i][j]->GetMean()+htime[i][j]->GetStdDev()*10) ;
    //       htime[i][j]->GetXaxis()->SetTitle("Delta t (ns)") ;
    //       c4->SaveAs("fit_plots/"+n+".jpg") ;

    //     }
    //   }
    // }

    gStyle->SetOptStat(1) ;
    TCanvas *c3 = new TCanvas("c3","c3",10,10,2000,1000) ;
    c3->Divide(1,2) ;
    c3->cd(1) ;
    config_histo1D(hcharge,"","Charge","",2,1,1) ;
    c3->cd(2) ;
    config_histo1D(henergy,"","Energy (MeV)","",2,1,1) ;
    c3->SaveAs("plots/energy_charge.pdf") ;

    TCanvas *c5 = new TCanvas("c5","c5",10,10,2000,1000) ;
    config_histo1D(hmultiplicity,"","Multiplicity","",2,1,1) ;
    c5->SaveAs("plots/multiplicity.pdf") ;

    gStyle->SetPaintTextFormat("1.f") ;
    gStyle->SetOptStat(0) ;
    TCanvas *c6 = new TCanvas("c6","c6",10,10,2000,1000) ;
    h2coincidence->Draw("COLZTEXT") ;
    if (str_wall == "it") {
      ReverseXAxis(h2coincidence) ;
    }
    c6->SaveAs("plots/coincidences.pdf") ;

    gStyle->SetPaintTextFormat("1.2f") ;
    TCanvas *c7 = new TCanvas("c7","c7",10,10,2000,1000) ;
    h2mean->Draw("COLZTEXT") ;
    if (str_wall == "it") {
      ReverseXAxis(h2mean) ;
    }

    c7->SaveAs("plots/mean.pdf") ;

    gStyle->SetOptStat(1) ;
    TCanvas *c8 = new TCanvas("c8","c8",10,10,2000,1000) ;
    config_histo1D(hsigma,"","Sigma (ns)","",2,1,1) ;
    c8->SaveAs("plots/sigma.pdf") ;

    gStyle->SetOptStat(1) ;
    TCanvas *c9 = new TCanvas("c9","c9",10,10,2000,1000) ;
    c9->Divide(2,2) ;
    c9->cd(1) ;
    config_histo1D(hamplitude,"","Amplitude (mV)","",1,1,1) ;
    c9->cd(2) ;
    //config_histo1D(htemps,"","time (s)","",2,1,1) ;
    config_histo1D(hcharge,"","Charge","",1,1,1) ;
    c9->cd(3) ;
    config_histo1D(hpulseTime,"","Pulse time (ns)","",1,1,1) ;
    c9->cd(4) ;
    config_histo1D(hmultiplicity,"","Multiplicity","",1,1,1) ;
    c9->SaveAs("plots/4_plots.pdf") ;

  }



  theTree->ResetBranchAddresses() ;

  // TCanvas *c11 = new TCanvas("c11","c11",10,10,2000,1000) ;
  // gStyle->SetOptStat(0) ;
  // gStyle->SetOptTitle(0) ;
  // config_histo1D(hamplSpect[8][1],"","Amplitude","",1,1,kCyan+2) ;
  // hamplSpect[8][1]->GetXaxis()->SetLimits(-500,0) ;
  // config_histo1D(hamplSpect[13][0],"SAME","Amplitude","",1,1,kMagenta+2) ;
  // config_histo1D(hamplSpect[9][6],"SAME","Amplitude","",1,1,kOrange+7) ;
  // c11->BuildLegend(0.82,0.66,0.89,0.89) ;

  // c11->SaveAs("plots/test_ampl.pdf") ;

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

void fit_delta_t(TH1F* histo, double fit_parameters[4]){

  double range_neg = histo->GetMean()-histo->GetStdDev()*10 ;
  double range_pos = histo->GetMean()+histo->GetStdDev()*10 ;
  TF1 *f1 = new TF1("f1","gaus",range_neg,range_pos) ;
  f1->SetParameter(1,histo->GetMean()) ;
  f1->SetParameter(2,histo->GetRMS()) ;
  histo->Fit("f1","QR") ;

  fit_parameters[0] = f1->GetParameter(1) ;
  fit_parameters[1] = f1->GetParError(1) ;
  fit_parameters[2] = f1->GetParameter(2) ;
  fit_parameters[3] = f1->GetParError(2) ;

  delete f1 ;
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

int index_area(int area, int x, int y){
  int index = -1 ;
  int imin,imax,jmin,jmax ;
  select_area(area,&imin,&imax,&jmin,&jmax) ;

  int irange = imax-imin ;
  int jrange = jmax-jmin ;

  int matrice_area[irange+1][jrange+1] ;

  for (int i = 0; i < irange+1; ++i) {
    for (int j = 0; j < jrange+1; ++j) {
      matrice_area[i][j] = i*(jmax+1)+j ;
    }
  }

  if (area > 0 && area < 6) {
    index = matrice_area[x][y] ;
  }
  else {
    cout << "Bad area number!!" << endl ;
  }

  return index ;
}
