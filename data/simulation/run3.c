#include <TFile.h>
#include <TGraphErrors.h>
#include <TCanvas.h>
#include <TAxis.h>
#include <TSystem.h>
#include <iostream>

void run3() {
    
    TFile *_file0 = TFile::Open("simulation_BPix_BmI_SEC1_LYR1_phase1_newL1.root", "READ");

    if (!_file0 || _file0->IsZombie()) {
        std::cerr << "Error opening simulation_BPix_BmI_SEC1_LYR1_phase1_newL1.root" << std::endl;
        return;
    }

    TFile *_file1 = TFile::Open("/afs/cern.ch/user/s/singhsh/PixelMonitoring/plots/currents/i_leak_run3_vs_fill_number.root", "READ");

    if (!_file1 || _file1->IsZombie()) {
        std::cerr << "Error opening i_leak_run3_vs_fill_number.root!" << std::endl;
        return;
    }

    
    _file0->ls();
    _file1->ls();

    
    TGraphErrors *graph1 = (TGraphErrors*)_file0->Get("I_leak");
    TGraphErrors *graph2 = (TGraphErrors*)_file1->Get("i_leak_run3");

    if (graph1 && graph2) {
        
        TCanvas *canvas1 = new TCanvas("canvas", "Overlapping Graphs", 800, 600);
        
       
        graph1->SetTitle("Leakage Current");
        graph1->GetXaxis()->SetTitle("Fills");
        graph1->GetYaxis()->SetTitle("Leakage Current");

       
        graph1->SetMarkerStyle(21);  
        graph1->SetMarkerColor(kRed);  
        graph1->Draw("AP");  

        
        graph2->SetMarkerStyle(22);  
        graph2->SetMarkerColor(kBlue);  
        graph2->Draw("P SAME");  

       
        canvas1->Update();

       
        canvas1->SaveAs("leakageCurrents_graphs.png");
    } else {
        std::cerr << "Error: Could not retrieve the graphs!" << std::endl;
    }

    
    _file0->Close();
    _file1->Close();
}
 
