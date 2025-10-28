#include <TFile.h>
#include <TGraphErrors.h>
#include <TCanvas.h>
#include <TAxis.h>
#include <TStyle.h>
#include <TLegend.h>
#include <TLatex.h>
#include <iostream>

void data_vs_fill() {
    // Style settings
    gStyle->SetOptStat(0);
    gStyle->SetTitleFont(42, "XYZ");
    gStyle->SetLabelFont(42, "XYZ");
    gStyle->SetTextFont(42);

    // Open ROOT file
    TFile *file = TFile::Open("/afs/cern.ch/user/s/singhsh/PixelRadiationDamageSimulation/testFile.root", "READ");
    if (!file || file->IsZombie()) {
        std::cerr << "Error opening testFile.root" << std::endl;
        return;
    }

    // Get the graph
    TGraphErrors *graph = (TGraphErrors*)file->Get("Temperature_vs_fill_data");
    if (!graph) {
        std::cerr << "Graph I_leak_vs_fill_data not found!" << std::endl;
        return;
    }

    // Create canvas
    TCanvas *c = new TCanvas("c", "Temperature vs Fill", 900, 600);
    c->SetGrid();

    // Improve graph appearance
    graph->SetTitle("Tempwerature vs Fill Number");
    graph->GetYaxis()->SetTitle("Temperature [K]");
    graph->GetXaxis()->SetTitle("Fill number");
    graph->SetMarkerStyle(20);
    graph->SetMarkerSize(0.7);
    graph->SetMarkerColor(kBlue+2);
    graph->SetLineColor(kBlue+2);
    graph->SetLineWidth(2);

    // Draw graph
    graph->Draw("AP");

    // Axis tweaks
    graph->GetXaxis()->SetTitleSize(0.045);
    graph->GetYaxis()->SetTitleSize(0.045);
    graph->GetXaxis()->SetLabelSize(0.04);
    graph->GetYaxis()->SetLabelSize(0.04);

    // Save plot
    c->SaveAs("Ileak_vs_fill.png");
}

