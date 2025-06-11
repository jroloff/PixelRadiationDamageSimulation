#include <TFile.h>
#include <TGraphErrors.h>
#include <TCanvas.h>
#include <TAxis.h>
#include <TSystem.h>
#include <iostream>
#include <TLegend.h>

//normalizeGraph Function
//write here if graph needs to be normalized
 
void run3() {
    
    TFile *_file0 = TFile::Open("/afs/cern.ch/user/s/singhsh/PixelRadiationDamageSimulation/testFile.root", "READ");

    if (!_file0 || _file0->IsZombie()) {
        std::cerr << "Error opening testFile.root" << std::endl;
        return;
    }

    TFile *_file1 = TFile::Open("/afs/cern.ch/user/s/singhsh/PixelMonitoring/plots/currents/i_leak_7920_9300_graphs.root", "READ");

    if (!_file1 || _file1->IsZombie()) {
        std::cerr << "Error opening i_leak_7920_9300_graphs.root" << std::endl;
        return;
    }

    
    _file0->ls();
    _file1->ls();

    
    TGraphErrors *graph1 = (TGraphErrors*)_file0->Get("I_leak_per_module_data");
    TGraphErrors *graph2 = (TGraphErrors*)_file0->Get("I_leak_volume");
  //TGraphErrors *graph3 = (TGraphErrors*)_file1->Get("i_leak_7920_9300_graph_2_fluence");
 // TGraphErrors *graph4 = (TGraphErrors*)_file1->Get("i_leak_7920_9300_graph_3_fluence");
 // TGraphErrors *graph5 = (TGraphErrors*)_file1->Get("i_leak_7920_9300_graph_4_fluence");

       if (graph1 && graph2 ){
       // Normalize all graphs 1-5
       // normalizeGraph(graph1);
       // normalizeGraph(graph2);
       // normalizeGraph(graph3);
       // normalizeGraph(graph4);
       // normalizeGraph(graph5);

    // double m = 0.1;  
   // Subtracting constant  from graph1
     int n = graph1->GetN();
     double *x1 = graph1->GetX();
     double *y1 = graph1->GetY();
     double c = y1[0] + 0.01; 

     TGraphErrors *graph1Sub = new TGraphErrors(n);
      for (int i = 0; i < n; ++i)
      {
      graph1Sub->SetPoint(i, x1[i], y1[i] -c); 
      }
 
   
      TCanvas *canvas1 = new TCanvas("canvas", "Overlapping Graphs", 800, 600);      
        graph1Sub->SetTitle("Leakage Current");
        graph1Sub->GetXaxis()->SetTitle("Date");
        graph1Sub->GetYaxis()->SetTitle("Leakage Current");
        graph1Sub->SetLineWidth(1);
	graph1Sub->SetMarkerStyle(22);
        graph1Sub->SetMarkerSize(0.8);	
        graph1Sub->SetMarkerColor(kBlack);
	// making sure x axis dates are displayed correctly
        graph1Sub->GetXaxis()->SetTimeDisplay(1);                   
        graph1Sub->GetXaxis()->SetTimeFormat("%d/%m/%Y");            
        graph1Sub->GetXaxis()->SetTimeOffset(0, "gmt");              
	
        graph1Sub->Draw("AP");  

        
        graph2->SetMarkerStyle(22);  
        graph2->SetMarkerColor(kBlue);  
        graph2->Draw("P SAME");

	//graph3->SetMarkerStyle(22);
        //graph3->SetMarkerColor(kRed);
       //graph3->Draw("P SAME");
        
        TLegend *legend = new TLegend(0.15, 0.75, 0.4, 0.88);
        legend->SetBorderSize(0);
        legend->SetFillStyle(0);
        legend->SetTextFont(42);
        legend->SetTextSize(0.035);
        legend->AddEntry(graph1Sub, "Leakage current (Data)", "lp");
        legend->AddEntry(graph2, "Leakage current (Simulation)", "lp");
       // legend->AddEntry(graph3, "Data layer 2", "lp");
       //legend->AddEntry(graph4, "Data layer 3", "lp");
       //legend->AddEntry(graph5, "Data layer 4", "lp");
        legend->Draw();
        canvas1->SetLeftMargin(0.15);       
        canvas1->Update();

       
        canvas1->SaveAs("leakageCurrents_graphs.png");
        }
      else {
        std::cerr << "Error: Could not retrieve the graphs!" << std::endl;
    }

    
    _file0->Close();
    _file1->Close();
}
 
