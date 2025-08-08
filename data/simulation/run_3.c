#include <TFile.h>
#include <TGraphErrors.h>
#include <TCanvas.h>
#include <TAxis.h>
#include <TSystem.h>
#include <iostream>
#include <TLegend.h>
#include <TLatex.h>
#include <TStyle.h>
#include <TGraphAsymmErrors.h>

void run_3() {
    gStyle->SetOptStat(0);
    gStyle->SetTextFont(42);
    gStyle->SetTitleFont(42, "XYZ");
    gStyle->SetLabelFont(42, "XYZ");

    TFile *_file0 = TFile::Open("/afs/cern.ch/user/s/singhsh/PixelRadiationDamageSimulation/TestFile.root", "READ");

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
    TGraphErrors *graph2 = (TGraphErrors*)_file0->Get("I_leak");
    TGraphErrors *tempGraph = (TGraphErrors*)_file0->Get("temperature");


     if (graph1 && graph2 && tempGraph ){



// Subtracting constant  from graph1
     int n = graph2->GetN();
     double *x1 = graph2->GetX();
     double *y1 = graph2->GetY();
     double c = y1[0] + 1.5;

     TGraphErrors *graph2Sub = new TGraphErrors(n);
      for (int i = 0; i < n; ++i)
      {
      graph2Sub->SetPoint(i, x1[i], y1[i] + c);
      }


      TCanvas *canvas1 = new TCanvas("canvas", "Overlapping Graphs", 1600, 1200);
      TPad *padT = new TPad("padT", "Temperature", 0, 0.8, 1, 1);
      TPad *padM = new TPad("padM", "Main plot", 0, 0.3, 1, 0.8);
      TPad *padR = new TPad("padR", "Ratio", 0, 0 , 1, 0.3);
      padT->SetBottomMargin(0.02);
      padM->SetTopMargin(0.02);
      padM->SetBottomMargin(0.02);
      padR->SetTopMargin(0.02);
      padR->SetBottomMargin(0.3);
      padT->SetLeftMargin(0.15);
      padM->SetLeftMargin(0.15);
      padR->SetLeftMargin(0.15);


      padT->Draw();
      padM->Draw();
      padR->Draw();


 // Top pad with termprature panel


      padT->cd();

      tempGraph->SetMarkerStyle(20);
      tempGraph->SetLineColor(kRed + 2);
      tempGraph->SetLineWidth(1);
      tempGraph->GetYaxis()->SetTitle("T [K]");
      tempGraph->GetYaxis()->SetTitleSize(0.09);
      tempGraph->GetYaxis()->SetTitleOffset(0.5);
      tempGraph->GetYaxis()->SetLabelSize(0.07);
      tempGraph->SetTitle("");
      tempGraph->GetXaxis()->SetLabelSize(0);  // hide x-axis labels
      tempGraph->GetXaxis()->SetTitle("");
      tempGraph->Draw("AL");
      TLatex *tempLabel = new TLatex();
      tempLabel->SetNDC();
      tempLabel->SetTextFont(42);
      tempLabel->SetTextSize(0.08);
      tempLabel->DrawLatex(0.35, 0.7, "Sensor Temperature");
// Middle pad with Leakage current Data and simulation
       padM->cd();
       graph1->SetTitle("");

       graph1->SetLineWidth(1);
       graph1->SetMarkerStyle(22);
       graph1->SetMarkerSize(0.8);
       graph1->SetMarkerColor(kBlack);
       graph1->GetYaxis()->SetTitle("Leakage Current [mA]");
       graph1->GetYaxis()->SetRangeUser(0, 3);
       graph1->GetYaxis()->SetTitleOffset(1.1);
       graph1->GetYaxis()->SetTitleSize(0.04);
       graph1->GetYaxis()->SetLabelSize(0.03);
       graph1->GetXaxis()->SetLabelSize(0);  // hide x-axis labels
       graph1->GetXaxis()->SetTitle("");
       graph1->Draw("AP");

       graph2->SetMarkerStyle(22);      
       graph2->SetMarkerColor(kBlue);
       graph2->Draw("P SAME");
       const double Eg = 1.21;             
       const double kB = 8.617e-5;         
       
 // Creating the error band around graph2

       TGraphAsymmErrors* band = new TGraphAsymmErrors(graph2->GetN());
       for (int k = 0; k < graph2->GetN(); ++k) {
            double x, y;
            graph2->GetPoint(k, x, y);

 // Get corresponding temperature
            double xT, T;
            tempGraph->GetPoint(k, xT, T);

 // Computing leakage current at T ± 2K using scaling factor
            double scale_plus = pow((T+2)/T, 2) * exp(-Eg/(2*kB) * (1/(T+2) - 1/T));
            double scale_minus = pow((T-2)/T, 2) * exp(-Eg/(2*kB) * (1/(T-2) - 1/T));

            double I_plus = y * scale_plus;
            double I_minus = y * scale_minus;

            double errYHigh = I_plus - y;
            double errYLow = y - I_minus;

            band->SetPoint(k, x, y);
            band->SetPointError(k, 0, 0, errYLow, errYHigh); 
}
        band->SetFillColor(kGreen - 10);
        band->SetLineColor(kGreen + 2);
	band->SetFillColorAlpha(kGreen - 10, 0.5);
        band->SetMarkerColor(kGreen + 2);
        band->SetFillStyle(1001);  

        band->Draw("3");        

       TLegend *legend = new TLegend(0.25, 0.75, 0.45, 0.88);

       legend->SetBorderSize(0);
       legend->SetFillStyle(0);
       legend->SetTextFont(42);
       legend->SetTextSize(0.035);
       legend->AddEntry(graph1, "Leakage current (Data) 2017-2018 Layer 1", "lp");
       legend->AddEntry(graph2, "Leakage current (Simulation)", "lp");
       legend->AddEntry(band, "Temperature Uncertainty \xB1 2 K", "f");
       legend->Draw();

// Down pad with ratio panel
     padR->cd();

     TGraphErrors* ratioGraph = new TGraphErrors();

    for (int j = 0; j < graph2->GetN(); ++j) {
        double x_data, y_data, x_sim, y_sim;
        graph1->GetPoint(j, x_data, y_data);
        graph2->GetPoint(j, x_sim, y_sim);
        if (y_data != 0) {
        ratioGraph->SetPoint(j, x_data, y_data / y_sim);
        }
    }
      ratioGraph->SetMarkerStyle(21);
      ratioGraph->SetMarkerColor(kBlue+2);
      ratioGraph->SetMarkerSize(0.6);
      ratioGraph->SetTitle("");
      ratioGraph->GetYaxis()->SetTitle("Data/Model");
      ratioGraph->GetYaxis()->SetTitleSize(0.07);
      ratioGraph->GetYaxis()->SetTitleOffset(0.6);
      ratioGraph->GetYaxis()->SetLabelSize(0.06);
      ratioGraph->GetXaxis()->SetLabelSize(0.07);
      ratioGraph->GetXaxis()->SetTitleSize(0.08);
      ratioGraph->GetXaxis()->SetTitleOffset(1.0);
      ratioGraph->GetXaxis()->SetTimeDisplay(1);
      ratioGraph->GetXaxis()->SetTimeFormat("%d/%m/%Y");
      ratioGraph->GetXaxis()->SetTimeOffset(0, "gmt");
      ratioGraph->Draw("AP");


       canvas1->Update();


       canvas1->SaveAs("leakageCurrents_graphs.pdf");
       }
      else {
        std::cerr << "Error: Could not retrieve the graphs!" << std::endl;
    }


    _file0->Close();
    _file1->Close();
}

