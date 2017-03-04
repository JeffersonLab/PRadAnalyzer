//============================================================================//
// A class contains a few root canvas                                         //
//                                                                            //
// Chao Peng                                                                  //
// 02/27/2016                                                                 //
//============================================================================//

#include <QLayout>

#include "TSystem.h"
#include "TStyle.h"
#include "TColor.h"
#include "TF1.h"
#include "TH1.h"
#include "TH2.h"
#include "TAxis.h"

#include "HistCanvas.h"
#include "QRootCanvas.h"

#define HIST_FONT_SIZE 0.07
#define HIST_LABEL_SIZE 0.07

HistCanvas::HistCanvas(QWidget *parent) : QWidget(parent)
{
    layout = new QGridLayout(this);

    bkgColor = new TColor(200, 1, 1, 0.96);
    gStyle->SetTitleFontSize(HIST_FONT_SIZE);
    gStyle->SetStatFontSize(HIST_FONT_SIZE);
}

void HistCanvas::AddCanvas(int row, int column, int color)
{
    QRootCanvas *newCanvas = new QRootCanvas(this);
    canvases.push_back(newCanvas);
    fillColors.push_back(color);

    // add canvas in vertical layout
    layout->addWidget(newCanvas, row, column);
    newCanvas->SetFillColor(bkgColor->GetNumber());
    newCanvas->SetFrameFillColor(10); // white
}

void HistCanvas::UpdateHist(int index, TH1 *hist, bool auto_range)
{
    if(!hist || index < 0 || index >= canvases.size())
        return;

    canvases[index]->cd();
    canvases[index]->SetGrid();

    //gPad->SetLogy();

    if(auto_range) {
        int firstBin = hist->FindFirstBinAbove(0,1)*0.7;
        int lastBin = hist->FindLastBinAbove(0,1)*1.3;

        hist->GetXaxis()->SetRange(firstBin, lastBin);
    }

    hist->GetXaxis()->SetLabelSize(HIST_LABEL_SIZE);
    hist->GetYaxis()->SetLabelSize(HIST_LABEL_SIZE);

    hist->SetFillColor(fillColors[index]);
    hist->Draw();

    canvases[index]->Refresh();
}

// show the histogram in first slot, try a Gaussian fit with given parameters
void HistCanvas::UpdateHist(int index, TH1 *hist, int range_min, int range_max)
{
    if(!hist || index < 0 || index >= canvases.size())
        return;

    canvases[index]->cd();
    canvases[index]->SetGrid();

    //gPad->SetLogy();

    hist->GetXaxis()->SetRange(range_min, range_max);

    hist->GetXaxis()->SetLabelSize(HIST_LABEL_SIZE);
    hist->GetYaxis()->SetLabelSize(HIST_LABEL_SIZE);

    hist->SetFillColor(fillColors[index]);
    hist->Draw();

    canvases[index]->Refresh();
}

void HistCanvas::UpdateHist(int index, TH2 *hist)
{
    if(!hist || index < 0 || index >= canvases.size())
        return;

    canvases[index]->cd();

    hist->GetXaxis()->SetLabelSize(HIST_LABEL_SIZE);
    hist->GetYaxis()->SetLabelSize(HIST_LABEL_SIZE);

    hist->Draw("colz");

    canvases[index]->Refresh();
}
