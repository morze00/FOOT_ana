#include"libs.hh"
#include <TMathBase.h>
using namespace std;
const double NSIGMA = 3;
bool is_good_strip(UInt_t i)
{
    if(i==2) return false;
    if(i>609 && i<615) return false;
    if((i%64)>62 || (i%64)<3) return false;
    if(i==273 || i==276) return false;
    return true;
}
void FOOT_ana(int firstEvent, int max_events)
{
    TApplication* theApp = new TApplication("App", 0, 0);
    TH2D * h2d_cal = new TH2D("h2d_cal","After simple pedestal subtraction",640,1,641,1200,-50,50);
    TH1D * h1d_ped= new TH1D("h1d_ped","After simple pedestal subtraction",110,-10,100);
    TH1D * h1d_cal = new TH1D("h1d_cal","After baseline correction",10000,0,100);
    TH1I * h1d_mul = new TH1I("h1d_mul","Multiplicity of good hits(strip)",20,0,20);
    TH2D * h2d_e1_vs_e2 = new TH2D("h2d_e1_vs_e2","Strip 1 vs stip 2 energy",1000,-100,100,1000,-100,100);
    TH1D * h1d_strip_no = new TH1D("strip_no","Strip number",640,1,641);
    TH1D * h2d_ev[25];
    TString hname;
    for(auto i=0; i<25; i++){
        hname = Form("cosmic hit %d",2*i+firstEvent);
        h2d_ev[i] = new TH1D(hname, hname, 640,1,641);
    }
    //============== Drawing raw data ===============
    TChain * ch = new TChain("h101");
    //ch->Add("cosmicrun-19-01-holddelay-6p5-us.root");
    ch->Add("../cosmicrun-21-01-holddelay-5us.root");
    TCanvas * c1 = new TCanvas("c1","c1",1000,1000);
    c1->Divide(2,2);
    c1->cd(1);
    ch->Draw("FOOT1E:FOOT1I>>h(640,1,641,500,0,500)","Entry$%2==0","colz");
    TH2D * h2d = (TH2D*)gDirectory->Get("h")->Clone("h2d");
    h2d->GetXaxis()->SetTitle("Strip No.");
    h2d->GetYaxis()->SetTitle("ADC channel");
    h2d->Draw("colz");
    //============== Slicing, fitting, drawing params ===============
    TF1* foo = new TF1("foo","gaus",0,1000);
    h2d->FitSlicesY(foo,1,640,0,"QNR",0);
    TH1D * h2d_1 = (TH1D*)gDirectory->Get("h2d_1")->Clone("h2d_1");
    TH1D * h2d_2 = (TH1D*)gDirectory->Get("h2d_2")->Clone("h2d_2");
    h2d_1->SetMarkerStyle(kFullCircle);
    h2d_1->SetMarkerSize(0.3);
    h2d_1->SetMarkerColor(kRed);
    h2d_1->SetLineColor(kRed);
    h2d_1->Draw("same");
    c1->cd(2);
    h2d_1->Draw();
    c1->cd(3);
    h2d_2->SetMarkerStyle(kFullCircle);
    h2d_2->SetMarkerSize(0.3);
    h2d_2->SetMarkerColor(kRed);
    h2d_2->SetLineColor(kRed);
    h2d_2->GetYaxis()->SetRangeUser(0,10);
    h2d_2->Draw();
    //============== Writing pedestals ===============
    auto pedfilename = "pedestal.dat";
    std::ofstream fout;
    fout.open(pedfilename);
    for(Int_t i = 0; i<640; i++){
        fout << "1  " << i+1 << "  " <<  h2d_1->GetBinContent(i+1) << "  " <<  h2d_2->GetBinContent(i+1) << endl;
        cout << "1  " << i+1 << "  " <<  h2d_1->GetBinContent(i+1) << "  " <<  h2d_2->GetBinContent(i+1) << endl;
    }
    fout.close();
    //============== Reading pedestals ===============
    ifstream pedfile(pedfilename,ifstream::in);
    if ( !pedfile.is_open() ){
        cout << "Cannot open Ped file" << endl;
    }
    int Nlines = 0;
    std::string line;
    while( std::getline(pedfile, line) )
        Nlines++;
    pedfile.clear();
    pedfile.seekg( 0, std::ios::beg );
    Int_t DetId_[Nlines];
    Int_t StripId_[Nlines];
    Double_t Ped_[Nlines];
    Double_t Sig_[Nlines];
    for(Int_t i=0 ; i<Nlines ; i++){
        pedfile >> DetId_[i] >> StripId_[i] >> Ped_[i] >> Sig_[i];
    }
    pedfile.close();
    //============== Applying pedestals to root data ===============
    UInt_t  FOOT1I[640];
    UInt_t  FOOT1E[640];
    ch->SetBranchAddress("FOOT1I",FOOT1I);
    ch->SetBranchAddress("FOOT1E",FOOT1E);
    Int_t Nevents = ch->GetEntries();
    if(max_events>0) Nevents = max_events; 
    Int_t     mul_good_hits = 0;
    Int_t     mul_good_events = 0;
    Double_t  mean_ssd=0 ;
    Double_t  asic_offset[10];
    Double_t  asic_offset_fine[10];
    Double_t  high_mul_events =0 ;
    Double_t  signal = 0;
    Double_t  signal_sum = 0;
    int       counter_asic =0;
    Int_t     stat=0;
    Double_t  estrip[2];
    int i=0;
    for(int ev=firstEvent; ev<firstEvent+Nevents; ev++)
    {
        cout << "\r-- Event # : " << ev << flush;
        if(ev%2==1) continue;
        ch->GetEntry(ev);
        //======== Global base line correction for this event 
        mean_ssd=0;
        stat=0;
        for(i=0; i<640; i++)
        {
            if(!is_good_strip(FOOT1I[i])) continue; 
            signal = FOOT1E[i] - Ped_[i];
            stat++;
            mean_ssd += signal;
            h1d_ped->Fill(signal);
            //h2d_cal->Fill(FOOT1I[i],signal);
        }
        mean_ssd = mean_ssd/stat;
        if(mean_ssd>10){
            cout << "\nMean ssd = " << mean_ssd;
            continue;
        }
        //======== Coarse base line correction for individual asics
        stat=0;            
        counter_asic=0;            
        for(i=0; i<10; i++){  asic_offset[i]=0; }//reset asic baselines
        for(i=0; i<640; i++)
        {
            if(is_good_strip(FOOT1I[i]))
            {
                signal = FOOT1E[i] - Ped_[i] - mean_ssd;
                stat++;
                asic_offset[counter_asic] += signal;
            }
            if((FOOT1I[i]%64)==0) 
            {
                asic_offset[counter_asic] /= stat;
                //cout << "\n Calculated asic_offset = " << asic_offset[counter_asic] << " in " << counter_asic << " asic" <<  endl; 
                counter_asic++;
                stat=0;
            }
        }
        //======== Fine base line correction for individual asics
        stat=0;            
        counter_asic=0;            
        for(i=0; i<10; i++){  asic_offset_fine[i]=0; }//reset asic baselines
        for(i=0; i<640; i++)
        {
            if(is_good_strip(FOOT1I[i]))
            {
                signal = FOOT1E[i] - Ped_[i] - mean_ssd - asic_offset[counter_asic];
                if(fabs(signal) < (NSIGMA * Sig_[i]) )
                {
                    stat++;
                    asic_offset_fine[counter_asic] += signal;
                }
            }
            if((FOOT1I[i]%64)==0) 
            {
                asic_offset_fine[counter_asic] /= stat;
                //cout << "\n Calculated asic_offset = " << asic_offset[counter_asic] << " in " << counter_asic << " asic" <<  endl; 
                counter_asic++;
                stat=0;
            }
        }
        //======== Get number of "good" hits
        counter_asic=0;
        mul_good_hits=0;
        double _asic;
        for(i=0; i<640; i++)
        {
            if(FOOT1I[i]%64 == 1 && FOOT1I[i]>1) counter_asic++;
            if(!is_good_strip(FOOT1I[i])) continue;
            signal = FOOT1E[i] - Ped_[i] - mean_ssd - asic_offset[counter_asic] - asic_offset_fine[counter_asic];
            //if(mul_good_hits>10)
            //if(signal>1000)
            //if(fabs(signal) > (NSIGMA * Sig_[i]) )
            //{
            //    cout << "\n Signal " << signal << " in strip " << FOOT1I[i];
            //    cout << "\n Raw data = " <<  FOOT1E[i];
            //    cout << "\n Ped_[i] = " << Ped_[i];
            //    cout << "\n Mean ssd = " << mean_ssd;
            //    cout << "\n asic_offset = " << asic_offset[counter_asic] << " in " << counter_asic << " asic" <<  endl; 
            //}
            //if(fabs(signal) > (NSIGMA * 1.9) ) mul_good_hits++;
            if(signal > (NSIGMA * Sig_[i]) ) mul_good_hits++;
            //if(signal > (NSIGMA * 1.9) ) mul_good_hits++;
            //if(signal < (-1)*(NSIGMA * Sig_[i]) ) mul_good_hits++;
        }
        if(mul_good_hits>10)
        {
            cout << "\n--High multiplitcity!" << mul_good_hits << endl;
            high_mul_events++;
            //continue;
        }
        h1d_mul->Fill(mul_good_hits);
        //if(mul_good_hits!=2) continue;
        if(mul_good_hits<1 || mul_good_hits>2) continue;
        //if(mul_good_hits<1) continue;
        //if(mul_good_hits<10) continue;
        //========= Filling final histograms
        counter_asic=0;
        signal_sum=0;
        estrip[0]=-10000; 
        estrip[1]=-10000;
        for(i=0; i<640; i++)
        {
            if((FOOT1I[i]%64) == 1 && FOOT1I[i]>1) counter_asic++;
            if(!is_good_strip(FOOT1I[i])) continue;
            signal = FOOT1E[i] - Ped_[i] - mean_ssd - asic_offset[counter_asic] - asic_offset_fine[counter_asic];
            if(mul_good_events<25)
            {
                h2d_ev[mul_good_events]->Fill(FOOT1I[i],signal);
                //if(fabs(signal) > (NSIGMA * Sig_[i]))
                //    h2d_ev[mul_good_events]->GetXaxis()->SetRangeUser(FOOT1I[i]-40, FOOT1I[i]+40);
            }
            h2d_cal->Fill(FOOT1I[i], signal);
            //if( signal > (NSIGMA * 1.9) )
            if( fabs(signal) > (NSIGMA * Sig_[i]) )
            {
                signal_sum += signal;
                if(estrip[0]==(-10000)) estrip[0] = signal;
                else estrip[1] = signal;
                //h1d_cal->Fill(signal);
                h1d_strip_no->Fill(FOOT1I[i]);
            }
        }
        if(mul_good_hits==2) h2d_e1_vs_e2->Fill(estrip[0], estrip[1]);
        h1d_cal->Fill(signal_sum);
        mul_good_events++;
    }//end of eventloop
    cout << "\n--Identified " << mul_good_events << " good events" << " and " << high_mul_events<< "chrismas-tree events" <<  endl;
    c1->cd(4);
    h2d_cal->Draw("colz");
    TCanvas * cmul = new TCanvas("cmul","cmul",1000,1000);
    cmul->cd(1);
    h1d_mul->Draw();
    TCanvas * c2 = new TCanvas("c2","c2",1000,1000);
    h1d_cal->Draw();
    //h1d_cal->Fit("gaus","","",-10,10);
    //double bin_width = 500./110.;
    //h1d_ped->Sumw2();
    //h1d_ped->Scale(1.3*h1d_cal->Integral(1,20*bin_width)/h1d_ped->Integral(1,20*bin_width));
    //h1d_ped->Draw("same L");
    //h1d_ped->Fit("gaus","","",-10,10);
    //gPad->SetLogy();
    TCanvas * c_events =  new TCanvas("c_events","c_events",2000,1400);
    //c_events->Divide(5,5,0,0);
    c_events->Divide(5,5);
    for(auto i=0; i<25; i++)
    {
        c_events->cd(i+1);
        h2d_ev[i]->SetBarWidth(1);
        h2d_ev[i]->Draw("HISTO B");
        if(i%2==0)
        {
            h2d_ev[i]->SetFillColor(kBlue);
        }
        else
        {
            h2d_ev[i]->SetFillColor(kRed);
        }
        h2d_ev[i]->GetYaxis()->SetRangeUser(-20., 50.);
    }
    TCanvas * c3 = new TCanvas("c3","c3",1000,1000);
    h2d_e1_vs_e2->Draw("colz");
    TCanvas * c4 = new TCanvas("c4","c4",1000,1000);
    h1d_strip_no->Sumw2();
    h1d_strip_no->Scale(1./26./3600.);//26 hours run
    h1d_strip_no->Draw();
    theApp->Run();
    return;
}
int main(Int_t argc, Char_t* argv[])
{
    gRandom = new TRandom3();
    gRandom->SetSeed(0);
    gROOT->Macro("rootlogon.C");
    gStyle->SetPalette(kRainBow);
    FOOT_ana(0,-1);
    return 0;
}
